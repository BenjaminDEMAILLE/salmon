//! `--gpu` has to reach the path a bare `salmon quant` actually takes.
//!
//! 2.6.0 made deterministic mode the default, and that path returns from `quant`
//! before the one-pass path builds its alignment backend. A `--gpu` run was
//! therefore accepted, printed nothing, and quantified on the CPU: the flag was
//! silently inert on every default run. These tests pin the flag at the binary
//! level, in both builds.
//!
//! With the `gpu` feature the run must announce the backend it acquired and
//! still produce a `quant.sf` byte-identical to the CPU `--fullLengthAlignment`
//! run (the backend changes where alignment happens, not the result). Without
//! the feature `--gpu` must be a hard error rather than a no-op.

use std::path::{Path, PathBuf};
use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

fn sequence(mut seed: u64, len: usize) -> String {
    (0..len)
        .map(|_| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            "ACGT".as_bytes()[(seed >> 33) as usize & 3] as char
        })
        .collect()
}

fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let a = sequence(3, 3000);
    let b = sequence(17, 3000);
    let fasta = dir.join("txp.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n")).unwrap();

    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(100);
    for i in 0..200 {
        let source = if i % 2 == 0 { &a } else { &b };
        let start = (i * 13) % (source.len() - 300);
        let mate1 = &source[start..start + 100];
        let mate2: String = source[start + 200..start + 300]
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect();
        r1.push_str(&format!("@f{i}\n{mate1}\n+\n{qual}\n"));
        r2.push_str(&format!("@f{i}\n{mate2}\n+\n{qual}\n"));
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    (fasta, p1, p2)
}

fn build_index(dir: &Path, fasta: &Path) -> PathBuf {
    let index = dir.join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());
    index
}

#[cfg(feature = "gpu")]
#[test]
fn gpu_flag_reaches_the_default_path_and_does_not_change_quant_sf() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = fixture(dir.path());
    let index = build_index(dir.path(), &fasta);

    let run = |out: &Path, extra: &[&str]| -> String {
        let o = Command::new(SALMON)
            .args(["quant", "-i"])
            .arg(&index)
            .args(["-l", "A", "-1"])
            .arg(&r1)
            .arg("-2")
            .arg(&r2)
            .args(["-p", "2", "-o"])
            .arg(out)
            .args(extra)
            .output()
            .unwrap();
        assert!(o.status.success(), "salmon quant {extra:?} failed");
        String::from_utf8_lossy(&o.stderr).into_owned()
    };

    let gpu_out = dir.path().join("gpu");
    let stderr = run(&gpu_out, &["--gpu"]);
    // The flag has to actually select a backend on the default path. Either
    // message proves it was consulted; which one depends on the machine.
    assert!(
        stderr.contains("GPU alignment backend ready") || stderr.contains("no usable GPU adapter"),
        "--gpu never reached the backend on the default path; stderr was:\n{stderr}"
    );

    // Same scores, wherever they were computed.
    let cpu_out = dir.path().join("cpu");
    run(&cpu_out, &["--fullLengthAlignment"]);
    let gpu_sf = std::fs::read_to_string(gpu_out.join("quant.sf")).unwrap();
    let cpu_sf = std::fs::read_to_string(cpu_out.join("quant.sf")).unwrap();
    assert_eq!(
        gpu_sf, cpu_sf,
        "--gpu produced a different quant.sf than the CPU full-length path"
    );
}

#[cfg(not(feature = "gpu"))]
#[test]
fn gpu_flag_is_a_hard_error_without_the_feature() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = fixture(dir.path());
    let index = build_index(dir.path(), &fasta);

    let out = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "A", "-1"])
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .args(["-p", "2", "-o"])
        .arg(dir.path().join("gpu"))
        .arg("--gpu")
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "--gpu succeeded in a build without GPU support; it must not be silently inert"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("without GPU support"),
        "expected a rebuild hint, got:\n{stderr}"
    );
}
