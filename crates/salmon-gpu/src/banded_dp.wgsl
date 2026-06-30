// Banded affine-gap DP, one alignment per invocation.
//
// A line-by-line mirror of `reference::run_dp`: same band, same recurrence, same
// `mqe`-or-max extraction, all in integer arithmetic. Because each invocation
// scores one independent task and there is no cross-task reduction, the result
// is deterministic and bit-identical to the CPU reference.
//
// Per-invocation it keeps three target-length rows (H_prev, H_cur, F_col) in
// function-private storage, capped at MAX_T; the host routes any task whose
// target exceeds MAX_T to the CPU reference instead.

const NEG_INF: i32 = -0x40000000;
const HALF_NEG: i32 = -0x20000000; // NEG_INF / 2, the "is -inf" guard
const WILDCARD: u32 = 4u;
const MAX_T: u32 = 384u;

struct Params {
    match_score: i32,
    mismatch_pen: i32,
    gap_open_pen: i32,
    gap_extend_pen: i32,
    bandwidth: i32,
    n_tasks: u32,
};

// q_off/q_len/t_off/t_len for each task, packed as a vec4<u32>.
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> tasks: array<vec4<u32>>;
@group(0) @binding(2) var<storage, read> qbuf: array<u32>; // one DNA5 base per element
@group(0) @binding(3) var<storage, read> tbuf: array<u32>;
@group(0) @binding(4) var<storage, read_write> scores: array<i32>;

fn sat_sub(a: i32, b: i32) -> i32 {
    if (a <= HALF_NEG) {
        return NEG_INF;
    }
    return a - b;
}

fn sub_score(qb: u32, tb: u32) -> i32 {
    if (qb == WILDCARD || tb == WILDCARD) {
        return -params.gap_extend_pen;
    }
    if (qb == tb) {
        return params.match_score;
    }
    return -params.mismatch_pen;
}

var<private> h_prev: array<i32, MAX_T>;
var<private> h_cur: array<i32, MAX_T>;
var<private> f_col: array<i32, MAX_T>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let task = gid.x;
    if (task >= params.n_tasks) {
        return;
    }
    let tmeta = tasks[task];
    let q_off = tmeta.x;
    let q_len = tmeta.y;
    let t_off = tmeta.z;
    let t_len = tmeta.w;

    if (q_len == 0u || t_len == 0u || t_len > MAX_T) {
        scores[task] = NEG_INF;
        return;
    }

    let qlen = i32(q_len);
    let tlen = i32(t_len);
    let go = params.gap_open_pen;
    let ge = params.gap_extend_pen;
    let w = params.bandwidth;

    // Seed H_prev with the leading-deletion edge: H(-1, it) = -(go + (it+1)*ge).
    for (var it: u32 = 0u; it < t_len; it = it + 1u) {
        h_prev[it] = -(go + (i32(it) + 1) * ge);
        f_col[it] = NEG_INF;
    }

    var global_max: i32 = NEG_INF;
    var mqe: i32 = NEG_INF;

    for (var iq: u32 = 0u; iq < q_len; iq = iq + 1u) {
        let iqi = i32(iq);
        var diag_src: i32 = select(-(go + iqi * ge), 0, iq == 0u);
        var h_left: i32 = -(go + (iqi + 1) * ge);
        var e_run: i32 = NEG_INF;
        let qb = qbuf[q_off + iq];

        for (var it: u32 = 0u; it < t_len; it = it + 1u) {
            let iti = i32(it);
            let r = iqi + iti;

            // band bounds (ksw2 anti-diagonal band)
            var st: i32 = 0;
            var en: i32 = tlen - 1;
            if (st < r - qlen + 1) { st = r - qlen + 1; }
            if (en > r) { en = r; }
            let lo = (r - w + 1) >> 1u;
            let hi = (r + w) >> 1u;
            if (st < lo) { st = lo; }
            if (en > hi) { en = hi; }
            let in_band = (iti >= st) && (iti <= en);

            let diag = diag_src;
            diag_src = h_prev[it];

            if (!in_band) {
                h_cur[it] = NEG_INF;
                f_col[it] = NEG_INF;
                e_run = NEG_INF;
                h_left = NEG_INF;
                continue;
            }

            let s = sub_score(qb, tbuf[t_off + it]);
            var m: i32 = NEG_INF;
            if (diag > HALF_NEG) {
                m = diag + s;
            }

            e_run = max(sat_sub(h_left, go + ge), sat_sub(e_run, ge));
            let f = max(sat_sub(h_prev[it], go + ge), sat_sub(f_col[it], ge));
            f_col[it] = f;

            var h = max(max(m, e_run), f);
            if (h < NEG_INF) { h = NEG_INF; }
            h_cur[it] = h;
            h_left = h;

            if (h > global_max) { global_max = h; }
            if (iq + 1u == q_len && h > mqe) { mqe = h; }
        }

        // swap h_prev <- h_cur for the next row
        for (var it: u32 = 0u; it < t_len; it = it + 1u) {
            h_prev[it] = h_cur[it];
        }
    }

    var result: i32;
    if (mqe > NEG_INF) {
        result = mqe;
    } else {
        result = max(global_max, 0);
    }
    scores[task] = result;
}
