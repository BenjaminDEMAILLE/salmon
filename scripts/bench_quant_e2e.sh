#!/usr/bin/env bash
# Reproducible A/B/n harness for full quantification, mapping included.
#
# scripts/bench_em_rad.sh replays a RAD and therefore excludes mapping, which is
# the phase most build-level changes (allocator, hasher, mapping kernels) move.
# This harness runs `salmon quant` from FASTQ so that phase is in the
# measurement, and sweeps the surfaces a claim has to survive before it counts
# as an end-to-end win: selective alignment and sketch, with and without decoys,
# with and without mapping output, with and without bootstraps, across threads.
#
# Each arm is a separate salmon binary. Arms are interleaved within a repeat
# (rotated, so no arm is always first), results are appended to results.tsv, and
# quant.sf is hashed per row so arms that were meant to compute the same thing
# can be proven to have done so.
#
# Usage:
#   scripts/bench_quant_e2e.sh \
#     --arm mimalloc=/path/to/salmon-mimalloc \
#     --arm snmalloc=/path/to/salmon-snmalloc \
#     --index /path/to/index -1 reads_1.fq.gz -2 reads_2.fq.gz \
#     --out bench-run [--decoy-index DIR] [--surfaces sa,sketch] \
#     [--threads 16,32] [--repeats 5] [--bootstraps 32] [--cpu-order 0,2,4,6] \
#     [--keep-outputs] [-- <extra salmon quant args>]
#
# Surfaces (--surfaces, comma separated, default sa,sketch):
#   sa            selective alignment, main index
#   sketch        --sketch, main index
#   sa-decoy      selective alignment, decoy-aware index (--decoy-index)
#   sketch-decoy  --sketch, decoy-aware index (--decoy-index)
#   sa-mappings   selective alignment + --writeMappings (SAM output path)
#   sa-bam        selective alignment + --writeBam
#   sa-boot       selective alignment + --numBootstraps (--bootstraps, default 32)

set -euo pipefail

usage() {
  sed -n '2,31p' "$0"
  exit "${1:-0}"
}

ARM_NAMES=()
ARM_BINS=()
INDEX=""
DECOY_INDEX=""
MATES1=""
MATES2=""
UNMATED=""
OUT=""
SURFACES="sa,sketch"
THREADS="16"
REPEATS=5
BOOTSTRAPS=32
CPU_ORDER=""
KEEP_OUTPUTS=0
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arm)
      [[ "$2" == *=* ]] || { echo "--arm expects name=path, got: $2" >&2; exit 2; }
      ARM_NAMES+=("${2%%=*}"); ARM_BINS+=("${2#*=}"); shift 2 ;;
    --index) INDEX="$2"; shift 2 ;;
    --decoy-index) DECOY_INDEX="$2"; shift 2 ;;
    -1|--mates1) MATES1="$2"; shift 2 ;;
    -2|--mates2) MATES2="$2"; shift 2 ;;
    -r|--unmated) UNMATED="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --surfaces) SURFACES="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --repeats) REPEATS="$2"; shift 2 ;;
    --bootstraps) BOOTSTRAPS="$2"; shift 2 ;;
    --cpu-order) CPU_ORDER="$2"; shift 2 ;;
    --keep-outputs) KEEP_OUTPUTS=1; shift ;;
    -h|--help) usage 0 ;;
    --) shift; EXTRA=("$@"); break ;;
    *) echo "unknown argument: $1" >&2; usage 2 ;;
  esac
done

(( ${#ARM_NAMES[@]} >= 2 )) || { echo "at least two --arm name=path are required" >&2; exit 2; }
for bin in "${ARM_BINS[@]}"; do
  [[ -x "$bin" ]] || { echo "arm binary is not executable: $bin" >&2; exit 2; }
done
[[ -d "$INDEX" ]] || { echo "--index does not exist: $INDEX" >&2; exit 2; }
[[ -n "$OUT" ]] || { echo "--out is required" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo "refusing to reuse output directory: $OUT" >&2; exit 2; }
if [[ -n "$UNMATED" ]]; then
  [[ -f "$UNMATED" ]] || { echo "--unmated does not exist: $UNMATED" >&2; exit 2; }
else
  [[ -f "$MATES1" && -f "$MATES2" ]] || { echo "-1/-2 are required (or -r)" >&2; exit 2; }
fi

IFS=',' read -r -a surface_list <<< "$SURFACES"
for surface in "${surface_list[@]}"; do
  case "$surface" in
    sa|sketch|sa-mappings|sa-bam|sa-boot) ;;
    sa-decoy|sketch-decoy)
      [[ -d "$DECOY_INDEX" ]] || { echo "surface $surface needs --decoy-index" >&2; exit 2; } ;;
    *) echo "unknown surface: $surface" >&2; usage 2 ;;
  esac
done

# Portability: this runs on both the Linux benchmark hosts and macOS laptops,
# and the two disagree about `time`, `sha256sum` and CPU inventory.
if /usr/bin/time -f '%e' true >/dev/null 2>&1; then
  TIME_FLAVOR="gnu"
else
  TIME_FLAVOR="bsd"
fi

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

cpu_inventory() {
  if command -v lscpu >/dev/null 2>&1; then
    lscpu
  else
    sysctl -n machdep.cpu.brand_string hw.physicalcpu hw.logicalcpu hw.memsize
  fi
}

mkdir -p "$OUT"
cpu_inventory > "$OUT/cpu.txt" 2>/dev/null || true
uname -a > "$OUT/uname.txt"

for i in "${!ARM_NAMES[@]}"; do
  name="${ARM_NAMES[$i]}"; bin="${ARM_BINS[$i]}"
  "$bin" --version > "$OUT/arm-$name-version.txt" 2>&1 || true
  # Records what the binary *is*, not what the flags claimed: builds that differ
  # only in a global allocator are otherwise indistinguishable in the results.
  RUST_LOG=debug "$bin" swim 2>&1 | grep -i "global allocator" \
    > "$OUT/arm-$name-allocator.txt" || true
done

{
  echo "index=$INDEX"
  echo "decoy_index=$DECOY_INDEX"
  echo "mates1=$MATES1"
  echo "mates2=$MATES2"
  echo "unmated=$UNMATED"
  echo "surfaces=$SURFACES"
  echo "threads=$THREADS"
  echo "repeats=$REPEATS"
  echo "bootstraps=$BOOTSTRAPS"
  echo "cpu_order=$CPU_ORDER"
  echo "time_flavor=$TIME_FLAVOR"
  for i in "${!ARM_NAMES[@]}"; do echo "arm=${ARM_NAMES[$i]}:${ARM_BINS[$i]}"; done
  printf 'extra='; printf ' %q' "${EXTRA[@]+"${EXTRA[@]}"}"; echo
} > "$OUT/config.txt"

printf 'arm\tsurface\tthreads\trepeat\torder\tindex_load_s\tmapping_s\tem_s\tposterior_s\toutput_s\twall_s\tuser_s\tsys_s\tmax_rss_kb\tnum_mapped\tpercent_mapped\tquant_sha256\n' > "$OUT/results.tsv"

affinity_for() {
  local n="$1"
  [[ -n "$CPU_ORDER" ]] || return 0
  local cpus selected
  IFS=',' read -r -a cpus <<< "$CPU_ORDER"
  (( ${#cpus[@]} >= n )) || { echo "--cpu-order has fewer than $n CPUs" >&2; exit 2; }
  selected="${cpus[0]}"
  for ((i=1; i<n; i++)); do selected+=",${cpus[i]}"; done
  printf '%s' "$selected"
}

phase_seconds() {
  local phase="$1" log="$2"
  sed 's/\x1b\[[0-9;]*m//g' "$log" \
    | grep -oE "phase=\"${phase}\" elapsed_s=[0-9.eE+-]+" \
    | tail -n 1 \
    | sed 's/.*elapsed_s=//'
}

json_scalar() {
  local key="$1" file="$2"
  [[ -f "$file" ]] || return 0
  sed -nE "s/.*\"${key}\"[[:space:]]*:[[:space:]]*([^,}]+).*/\1/p" "$file" | tail -n 1
}

# Normalizes both `time` dialects to one line: wall<TAB>user<TAB>sys<TAB>rss_kb.
# GNU %M is already KiB; BSD reports maximum resident set size in bytes.
parse_timing() {
  local file="$1"
  if [[ "$TIME_FLAVOR" == "gnu" ]]; then
    tail -n 1 "$file"
  else
    awk '
      $1 == "real" { wall = $2 }
      $1 == "user" { user = $2 }
      $1 == "sys"  { sys  = $2 }
      /maximum resident set size/ { rss = $1 / 1024 }
      END { printf "%s\t%s\t%s\t%d\n", wall, user, sys, rss }
    ' "$file"
  fi
}

surface_index() {
  case "$1" in
    sa-decoy|sketch-decoy) printf '%s' "$DECOY_INDEX" ;;
    *) printf '%s' "$INDEX" ;;
  esac
}

run_one() {
  local arm="$1" binary="$2" surface="$3" threads="$4" repeat="$5" order="$6"
  local tag="${surface}-t${threads}-r${repeat}-o${order}-${arm}"
  local run_out="$OUT/$tag" log="$OUT/$tag.log" timing="$OUT/$tag.time"

  local surface_args=()
  case "$surface" in
    sketch|sketch-decoy) surface_args=(--sketch) ;;
    sa-mappings) surface_args=(--writeMappings "$run_out.mappings.sam") ;;
    sa-bam) surface_args=(--writeBam "$run_out.mappings.bam") ;;
    sa-boot) surface_args=(--numBootstraps "$BOOTSTRAPS") ;;
  esac

  local read_args=()
  if [[ -n "$UNMATED" ]]; then
    read_args=(-r "$UNMATED")
  else
    read_args=(-1 "$MATES1" -2 "$MATES2")
  fi

  local affinity prefix=()
  affinity="$(affinity_for "$threads")"
  [[ -z "$affinity" ]] || prefix=(taskset -c "$affinity")

  local time_cmd=()
  if [[ "$TIME_FLAVOR" == "gnu" ]]; then
    time_cmd=(/usr/bin/time -f '%e\t%U\t%S\t%M' -o "$timing")
  else
    time_cmd=(/usr/bin/time -lp -o "$timing")
  fi

  # Both flavors write their report to $timing, so salmon keeps stderr to
  # itself and the two streams never have to be untangled.
  "${time_cmd[@]}" \
    env RUST_LOG=salmon::timing=info "${prefix[@]+"${prefix[@]}"}" \
    "$binary" quant -i "$(surface_index "$surface")" -l A -p "$threads" \
    "${read_args[@]}" --sigDigits 9 -o "$run_out" \
    "${surface_args[@]+"${surface_args[@]}"}" "${EXTRA[@]+"${EXTRA[@]}"}" 2> "$log"

  local wall user sys rss
  IFS=$'\t' read -r wall user sys rss < <(parse_timing "$timing")

  local meta="$run_out/aux_info/meta_info.json"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$arm" "$surface" "$threads" "$repeat" "$order" \
    "$(phase_seconds index_load "$log")" \
    "$(phase_seconds mapping "$log")" \
    "$(phase_seconds em_bias "$log")" \
    "$(phase_seconds posterior "$log")" \
    "$(phase_seconds output "$log")" \
    "$wall" "$user" "$sys" "$rss" \
    "$(json_scalar num_mapped "$meta")" \
    "$(json_scalar percent_mapped "$meta")" \
    "$(sha256_of "$run_out/quant.sf")" >> "$OUT/results.tsv"

  # Mapping output is the surface that fills a disk: one SAM per run, per arm,
  # per repeat. Hashes and timings are already recorded by this point.
  if (( ! KEEP_OUTPUTS )); then
    rm -f "$run_out.mappings.sam" "$run_out.mappings.bam"
    rm -rf "$run_out/bootstraps"
  fi
}

IFS=',' read -r -a thread_counts <<< "$THREADS"
arm_count="${#ARM_NAMES[@]}"

# Warm the page cache over index and reads on every surface before recording:
# a cold first arm otherwise reads as a regression. Warmups are deliberately
# outside results.tsv.
warm_threads="${thread_counts[0]}"
for surface in "${surface_list[@]}"; do
  run_one "${ARM_NAMES[0]}" "${ARM_BINS[0]}" "$surface" "$warm_threads" 0 0
done
sed -i.bak $'/\t0\t0\t/d' "$OUT/results.tsv" && rm -f "$OUT/results.tsv.bak"

for surface in "${surface_list[@]}"; do
  for threads in "${thread_counts[@]}"; do
    for ((repeat=1; repeat<=REPEATS; repeat++)); do
      # Rotate which arm goes first, so no arm systematically inherits the
      # machine state the previous one left behind.
      for ((slot=0; slot<arm_count; slot++)); do
        idx=$(( (slot + repeat - 1) % arm_count ))
        run_one "${ARM_NAMES[$idx]}" "${ARM_BINS[$idx]}" \
          "$surface" "$threads" "$repeat" "$((slot + 1))"
      done
    done
  done
done

echo "wrote $OUT/results.tsv"
