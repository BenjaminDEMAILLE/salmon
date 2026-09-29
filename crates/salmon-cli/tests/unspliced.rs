//! `salmon index --unspliced` / `salmon quant` on an unspliced index, end to
//! end through the binary: option wiring and validation.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

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

struct Fx {
    txome: PathBuf,
    genome: PathBuf,
    gtf: PathBuf,
}

/// One two-exon gene on a 4 kb chromosome, its transcript, and the GTF.
fn fixture(dir: &Path) -> Fx {
    let chr = sequence(7, 4000);
    let tx = format!("{}{}", &chr[500..1000], &chr[2000..2500]);
    let txome = dir.join("txome.fa");
    std::fs::write(&txome, format!(">T1\n{tx}\n")).unwrap();
    let genome = dir.join("genome.fa");
    std::fs::write(&genome, format!(">chr1\n{chr}\n")).unwrap();
    let gtf = dir.join("ann.gtf");
    std::fs::write(
        &gtf,
        "chr1\tt\texon\t501\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
         chr1\tt\texon\t2001\t2500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();
    Fx { txome, genome, gtf }
}

fn salmon(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(SALMON).args(args).output().unwrap()
}

fn os(s: &str) -> &std::ffi::OsStr {
    s.as_ref()
}

#[test]
fn read_length_sets_the_flank_and_the_table_is_written() {
    let tmp = tempfile::tempdir().unwrap();
    let fx = fixture(tmp.path());
    let idx = tmp.path().join("idx");
    let out = salmon(&[
        os("index"),
        os("-t"),
        fx.txome.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("1"),
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
        os("--readLength"),
        os("51"),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(idx.join("info.json")).unwrap()).unwrap();
    assert_eq!(info["unspliced"]["mode"], "intron");
    assert_eq!(info["unspliced"]["flank"], 50);
    assert_eq!(
        std::fs::read_to_string(idx.join("t2g_3col.tsv")).unwrap(),
        "T1\tG1\tS\nG1-I\tG1\tU\n"
    );
}

#[test]
fn unspliced_options_are_validated() {
    let tmp = tempfile::tempdir().unwrap();
    let fx = fixture(tmp.path());
    let idx = tmp.path().join("idx");
    let base = [
        os("index"),
        os("-t"),
        fx.txome.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("1"),
    ];
    let run = |extra: &[&std::ffi::OsStr]| {
        let mut a = base.to_vec();
        a.extend_from_slice(extra);
        let o = salmon(&a);
        assert!(!o.status.success(), "{extra:?} should fail");
        String::from_utf8_lossy(&o.stderr).into_owned()
    };
    // --genome and --gtf go together
    run(&[
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
    ]);
    // --flank and --readLength are two spellings of one setting
    run(&[os("--flank"), os("10"), os("--readLength"), os("51")]);
    // a mode clap does not know
    run(&[os("--unspliced"), os("splici")]);
    // intron mode needs a flank
    let e = run(&[
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
    ]);
    assert!(e.contains("flank"), "{e}");
    // annotation inputs without a mode
    let e = run(&[
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
    ]);
    assert!(e.contains("without --unspliced"), "{e}");
}
