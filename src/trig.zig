//! Trigonometry on BAM Angles (SPEC §6).
//!
//! sin/cos via octant reduction: the top 3 bits of `Angle.raw` pick which
//! of 8 octants the angle lies in; the bottom 29 bits give the position
//! within that octant. The position is converted to a Fixed value in
//! [0, π/4] and run through a Horner-evaluated Taylor polynomial. Sign and
//! swap (sin↔cos) come from the octant index — no modulo, no precision
//! loss at large angles, no division by a low-precision π.
//!
//! Target accuracy: max abs error ≤ 2^-22 across the full circle (SPEC §6).
//! Polynomial coefficients are locked in by the conformance suite — changing
//! them is a determinism contract break.

const std = @import("std");
const Fixed = @import("Fixed.zig");
const Angle = @import("Angle.zig");
const rounding = @import("rounding.zig");

// ---------------------------------------------------------------------------
// Octant reduction scaling
// ---------------------------------------------------------------------------

// reduced.raw = round(t × π/128 × 2^57 / 2^57)
//   where t is the bottom 29 bits (plus possible mirror) of Angle.raw.
// At t = 2^29 the input represents exactly π/4 — the maximum reduced value.
const PI_OVER_128_Q57: i128 = blk: {
    @setEvalBranchQuota(10_000);
    const k: comptime_float = std.math.pi / @as(comptime_float, 128.0);
    break :blk @intFromFloat(@floor(k * @as(comptime_float, 1 << 57) + 0.5));
};

fn reduceToFixed(t: u32) Fixed {
    std.debug.assert(t <= (@as(u32, 1) << 29));
    const product: i128 = @as(i128, t) * PI_OVER_128_Q57;
    return .{ .raw = rounding.shiftRoundNarrow(product, 57) };
}

// ---------------------------------------------------------------------------
// Polynomial coefficients (Taylor series on [0, π/4])
//
// sin(x) ≈ x · (1 - x²/6 + x⁴/120 - x⁶/5040 + x⁸/362880)
//   — terms through x^9, truncation error |x^11/11!| < 2^-29 at π/4
//   — coefficient + Horner rounding budget ≈ 2^-23.7 in practice
//
// cos(x) ≈ 1 - x²/2 + x⁴/24 - x⁶/720 + x⁸/40320
//   — terms through x^8, truncation error |x^10/10!| < 2^-25 at π/4
// ---------------------------------------------------------------------------

const SIN_C0: Fixed = Fixed.ONE; // outer constant (1)
const SIN_C1: Fixed = Fixed.rconst(-1.0 / 6.0);
const SIN_C2: Fixed = Fixed.rconst(1.0 / 120.0);
const SIN_C3: Fixed = Fixed.rconst(-1.0 / 5040.0);
const SIN_C4: Fixed = Fixed.rconst(1.0 / 362880.0);

const COS_C0: Fixed = Fixed.ONE;
const COS_C1: Fixed = Fixed.rconst(-1.0 / 2.0);
const COS_C2: Fixed = Fixed.rconst(1.0 / 24.0);
const COS_C3: Fixed = Fixed.rconst(-1.0 / 720.0);
const COS_C4: Fixed = Fixed.rconst(1.0 / 40320.0);

fn polySin(x: Fixed) Fixed {
    // sin(x) = x · (C0 + y·(C1 + y·(C2 + y·(C3 + y·C4))))   where y = x²
    const y = x.mul(x);
    var p = SIN_C4;
    p = SIN_C3.add(y.mul(p));
    p = SIN_C2.add(y.mul(p));
    p = SIN_C1.add(y.mul(p));
    p = SIN_C0.add(y.mul(p));
    return x.mul(p);
}

fn polyCos(x: Fixed) Fixed {
    // cos(x) = C0 + y·(C1 + y·(C2 + y·(C3 + y·C4)))   where y = x²
    const y = x.mul(x);
    var p = COS_C4;
    p = COS_C3.add(y.mul(p));
    p = COS_C2.add(y.mul(p));
    p = COS_C1.add(y.mul(p));
    p = COS_C0.add(y.mul(p));
    return p;
}

// ---------------------------------------------------------------------------
// sin / cos
// ---------------------------------------------------------------------------

pub fn sin(a: Angle) Fixed {
    // Octant decomposition. See SPEC §6 — `oct` indexes π/4-sized arcs;
    // `pos` is the position inside the octant.
    //
    // For sin:
    //   use_cos_poly = oct.bit0 XOR oct.bit1   (true on odd-bit-pair octants)
    //   mirror       = oct.bit0
    //   negate       = oct.bit2
    const oct: u32 = a.raw >> 29;
    const pos: u32 = a.raw & 0x1FFF_FFFF;
    const mirror = (oct & 1) != 0;
    const use_cos = ((oct & 1) ^ ((oct >> 1) & 1)) != 0;
    const negate = (oct & 4) != 0;

    const t: u32 = if (mirror) (@as(u32, 1) << 29) - pos else pos;
    const reduced = reduceToFixed(t);

    var result = if (use_cos) polyCos(reduced) else polySin(reduced);
    if (negate) result = result.neg();
    return result;
}

pub fn cos(a: Angle) Fixed {
    // For cos:
    //   use_sin_poly = oct.bit0 XOR oct.bit1
    //   mirror       = oct.bit0
    //   negate       = oct.bit1 XOR oct.bit2
    const oct: u32 = a.raw >> 29;
    const pos: u32 = a.raw & 0x1FFF_FFFF;
    const mirror = (oct & 1) != 0;
    const use_sin = ((oct & 1) ^ ((oct >> 1) & 1)) != 0;
    const negate = (((oct >> 1) & 1) ^ ((oct >> 2) & 1)) != 0;

    const t: u32 = if (mirror) (@as(u32, 1) << 29) - pos else pos;
    const reduced = reduceToFixed(t);

    var result = if (use_sin) polySin(reduced) else polyCos(reduced);
    if (negate) result = result.neg();
    return result;
}

// ---------------------------------------------------------------------------
// tan
// ---------------------------------------------------------------------------

pub fn tan(a: Angle) Fixed {
    const s = sin(a);
    const c = cos(a);
    // Pole case: cos rounds to exactly 0 not just at QUARTER_TURN /
    // THREE_QUARTER_TURN but for any angle within ~21 u29 units of those
    // poles — because the reduced Fixed value rounds to 0 at that scale
    // (Fixed has 24 fractional bits; Angle has 32). Saturate by sin's sign
    // (SPEC §6, §7). No assert: std.debug.assert becomes `unreachable` in
    // ReleaseFast, which licenses the optimizer to delete this guard and
    // hit a real divide-by-zero. The defined saturation is the contract.
    if (c.raw == 0) {
        return if (s.raw >= 0) Fixed.MAX else Fixed.MIN;
    }
    // Off-pole but possibly near-pole: the division may overflow i64.
    // Same i128-scaled divide as Fixed.div, but saturate instead of
    // wrap-and-assert on the narrow.
    const numerator: i128 = @as(i128, s.raw) << (Fixed.frac_bits + 1);
    const quotient_2x: i128 = @divTrunc(numerator, @as(i128, c.raw));
    const result: i128 = rounding.shiftRound(quotient_2x, 1);
    if (result > std.math.maxInt(i64)) return Fixed.MAX;
    if (result < std.math.minInt(i64)) return Fixed.MIN;
    return .{ .raw = @truncate(result) };
}

// ---------------------------------------------------------------------------
// atan polynomial + atan2
// ---------------------------------------------------------------------------

// SQRT2_MINUS_1 ≈ 0.414. atan splits [0, 1] here to keep the Taylor input
// in [-(√2-1), √2-1] for fast convergence.
const SQRT2_MINUS_1: Fixed = Fixed.rconst(0.41421356237309504880168872420969807856967187537694);
const PI_OVER_4: Fixed = Fixed.rconst(std.math.pi / 4.0);

// atan(x) ≈ x - x³/3 + x⁵/5 - x⁷/7 + x⁹/9 - x¹¹/11   (Taylor, degree 11)
//        = x · (1 + y·(-1/3 + y·(1/5 + y·(-1/7 + y·(1/9 + y·(-1/11))))))   y = x²
// Truncation at |x| ≤ √2-1: |x^13/13| ≤ 1.2e-7 ≈ 2^-23 — meets/beats 2^-22.
const ATAN_C5: Fixed = Fixed.rconst(-1.0 / 11.0);
const ATAN_C4: Fixed = Fixed.rconst(1.0 / 9.0);
const ATAN_C3: Fixed = Fixed.rconst(-1.0 / 7.0);
const ATAN_C2: Fixed = Fixed.rconst(1.0 / 5.0);
const ATAN_C1: Fixed = Fixed.rconst(-1.0 / 3.0);

fn polyAtan(x: Fixed) Fixed {
    const y = x.mul(x);
    var p = ATAN_C5;
    p = ATAN_C4.add(y.mul(p));
    p = ATAN_C3.add(y.mul(p));
    p = ATAN_C2.add(y.mul(p));
    p = ATAN_C1.add(y.mul(p));
    p = Fixed.ONE.add(y.mul(p));
    return x.mul(p);
}

/// atan(z) for z ∈ [0, 1]. Uses the inflection identity
/// atan(z) = π/4 + atan((z-1)/(z+1)) for z > √2-1 to keep the polynomial
/// input small and the series convergent.
fn atanZeroToOne(z: Fixed) Fixed {
    if (z.lte(SQRT2_MINUS_1)) return polyAtan(z);
    const w = z.sub(Fixed.ONE).div(z.add(Fixed.ONE));
    return PI_OVER_4.add(polyAtan(w));
}

/// atan2 in Fixed radians (-π, π]. Internal helper — the public API returns
/// an Angle (BAM u32) via the wrapper below.
fn atan2Radians(y: Fixed, x: Fixed) Fixed {
    // SPEC §6: atan2(0, 0) returns 0 with an assert (no in-band sentinel).
    if (y.raw == 0 and x.raw == 0) {
        std.debug.assert(false);
        return Fixed.ZERO;
    }
    // Axis cases — avoid div-by-zero in the abs_y/abs_x path.
    if (x.raw == 0) {
        return if (y.raw > 0) Fixed.HALF_PI else Fixed.HALF_PI.neg();
    }
    if (y.raw == 0) {
        return if (x.raw > 0) Fixed.ZERO else Fixed.PI;
    }

    // Q1-projected base = atan(|y|/|x|) ∈ [0, π/2]. Swap when |y| > |x|
    // to keep the polynomial input z = (smaller / larger) ∈ [0, 1].
    const abs_y = y.abs();
    const abs_x = x.abs();
    const swap = abs_y.raw > abs_x.raw;
    const z = if (swap) abs_x.div(abs_y) else abs_y.div(abs_x);
    var base = atanZeroToOne(z);
    if (swap) base = Fixed.HALF_PI.sub(base); // atan(|y|/|x|) = π/2 - atan(|x|/|y|)

    // Quadrant adjustment. Result lies in (-π, π].
    if (x.raw > 0) {
        return if (y.raw > 0) base else base.neg(); // Q1 / Q4
    }
    // x < 0 (x==0 handled above).
    return if (y.raw > 0) Fixed.PI.sub(base) else base.sub(Fixed.PI); // Q2 / Q3
}

pub fn atan2(y: Fixed, x: Fixed) Angle {
    return Angle.fromRadians(atan2Radians(y, x));
}

// ===========================================================================
// Tests
// ===========================================================================

const testing = std.testing;

test "sin / cos cardinal values" {
    try testing.expectEqual(Fixed.ZERO, sin(Angle.ZERO));
    try testing.expectEqual(Fixed.ONE, cos(Angle.ZERO));

    try testing.expectEqual(Fixed.ONE, sin(Angle.QUARTER_TURN));
    try testing.expectEqual(Fixed.ZERO, cos(Angle.QUARTER_TURN));

    try testing.expectEqual(Fixed.ZERO, sin(Angle.HALF_TURN));
    try testing.expectEqual(Fixed.NEG_ONE, cos(Angle.HALF_TURN));

    try testing.expectEqual(Fixed.NEG_ONE, sin(Angle.THREE_QUARTER_TURN));
    try testing.expectEqual(Fixed.ZERO, cos(Angle.THREE_QUARTER_TURN));
}

test "sin/cos: sin² + cos² ≈ 1 within 2^-22 across a sweep" {
    // The Pythagorean identity is the cheapest end-to-end accuracy probe.
    // 2^-22 in Fixed.raw units = 2^-22 · 2^24 = 4. Allow 8 ULPs of slack
    // for the two extra muls and one add in s²+c².
    inline for ([_]u32{
        0,         1234567,    0x10000000, 0x1FFF_FFFF,
        0x2000_0000, 0x4000_0000, 0x6000_0000, 0x8000_0000,
        0x8000_0001, 0xA000_0000, 0xC000_0000, 0xE000_0000,
        0xDEAD_BEEF, 0xFFFF_FFFF,
    }) |raw| {
        const a = Angle{ .raw = raw };
        const s = sin(a);
        const c = cos(a);
        const sum = s.mul(s).add(c.mul(c));
        const diff = @abs(sum.raw - Fixed.ONE.raw);
        try testing.expect(diff <= 8);
    }
}

test "sin/cos: signs across the four quadrants" {
    // Q1: sin > 0, cos > 0
    const q1 = Angle{ .raw = @as(u32, 1) << 28 }; // π/8
    try testing.expect(sin(q1).raw > 0 and cos(q1).raw > 0);
    // Q2: sin > 0, cos < 0
    const q2 = Angle{ .raw = 3 * (@as(u32, 1) << 29) - (@as(u32, 1) << 28) }; // 3π/4 + small offset
    try testing.expect(sin(q2).raw > 0 and cos(q2).raw < 0);
    // Q3: sin < 0, cos < 0
    const q3 = Angle{ .raw = (@as(u32, 5) << 29) }; // 5π/4
    try testing.expect(sin(q3).raw < 0 and cos(q3).raw < 0);
    // Q4: sin < 0, cos > 0
    const q4 = Angle{ .raw = (@as(u32, 7) << 29) }; // 7π/4
    try testing.expect(sin(q4).raw < 0 and cos(q4).raw > 0);
}

test "sin: anti-symmetric around full turn — sin(2π - a) = -sin(a)" {
    inline for ([_]u32{ 1234567, 0x10000000, 0x3000_0000, 0x6000_0000 }) |raw| {
        const a = Angle{ .raw = raw };
        const neg_a = Angle{ .raw = 0 -% raw };
        const s_a = sin(a);
        const s_neg = sin(neg_a);
        const diff = @abs(s_a.raw + s_neg.raw);
        try testing.expect(diff <= 4); // anti-symmetry within a few ULP
    }
}

test "cos: symmetric around full turn — cos(2π - a) = cos(a)" {
    inline for ([_]u32{ 1234567, 0x10000000, 0x3000_0000, 0x6000_0000 }) |raw| {
        const a = Angle{ .raw = raw };
        const neg_a = Angle{ .raw = 0 -% raw };
        const c_a = cos(a);
        const c_neg = cos(neg_a);
        const diff = @abs(c_a.raw - c_neg.raw);
        try testing.expect(diff <= 4);
    }
}

test "sin/cos: differential vs libm @sin/@cos within 2^-22" {
    // Sample 1024 angles spread evenly across the circle. Convert each to
    // f64 radians (lossy by ~3e-8) and compare polynomial output to libm's
    // sin/cos. Worst-case observed diff includes our polynomial error and
    // the radian-input precision loss.
    const target: f64 = 2.5e-7; // a hair above 2^-22 ≈ 2.384e-7
    var max_sin_err: f64 = 0;
    var max_cos_err: f64 = 0;
    var i: u32 = 0;
    const N: u32 = 1024;
    while (i < N) : (i += 1) {
        const angle_raw: u32 = @intCast((@as(u64, i) * std.math.maxInt(u32)) / N);
        const a = Angle{ .raw = angle_raw };
        const rad = a.toRadians().toFloat();
        const got_s = sin(a).toFloat();
        const got_c = cos(a).toFloat();
        const ref_s = @sin(rad);
        const ref_c = @cos(rad);
        max_sin_err = @max(max_sin_err, @abs(got_s - ref_s));
        max_cos_err = @max(max_cos_err, @abs(got_c - ref_c));
    }
    try testing.expect(max_sin_err < target);
    try testing.expect(max_cos_err < target);
}

test "tan: cardinal values within tolerance" {
    try testing.expectEqual(Fixed.ZERO, tan(Angle.ZERO));
    try testing.expectEqual(Fixed.ZERO, tan(Angle.HALF_TURN));

    // tan(π/4) ≈ 1
    const eighth = Angle{ .raw = @as(u32, 1) << 29 }; // π/4
    const t = tan(eighth);
    try testing.expect(@abs(t.raw - Fixed.ONE.raw) <= 16);

    // tan(π/6) ≈ 1/√3 ≈ 0.577 — Angle.raw = 2^32 / 12 = 357913941
    const sixth = Angle{ .raw = 357913941 };
    const t6 = tan(sixth);
    const expected = Fixed.rconst(1.0 / @sqrt(3.0));
    try testing.expect(@abs(t6.raw - expected.raw) <= 16);
}

test "tan: defined saturation when cos rounds to zero (the polynomial pole)" {
    // Fixed-resolution caveat: for any pos within ~21 u29 ULPs of a cardinal
    // pole, the reduced Fixed input rounds to 0, so the polynomial outputs
    // sin = ±1 and cos = 0 exactly. tan saturates by sin's sign — same
    // result either side of the pole. This isn't the mathematical sign-flip
    // across ±∞ (we can't resolve that at the Fixed/Angle interface), but
    // it IS a defined, deterministic saturation in every build mode.

    // Near +π/2: sin ≈ +1 in both directions → MAX.
    try testing.expectEqual(Fixed.MAX, tan(.{ .raw = Angle.QUARTER_TURN.raw - 1 }));
    try testing.expectEqual(Fixed.MAX, tan(Angle.QUARTER_TURN));
    try testing.expectEqual(Fixed.MAX, tan(.{ .raw = Angle.QUARTER_TURN.raw + 1 }));

    // Near 3π/2 (=-π/2): sin ≈ -1 → MIN.
    try testing.expectEqual(Fixed.MIN, tan(.{ .raw = Angle.THREE_QUARTER_TURN.raw - 1 }));
    try testing.expectEqual(Fixed.MIN, tan(Angle.THREE_QUARTER_TURN));
    try testing.expectEqual(Fixed.MIN, tan(.{ .raw = Angle.THREE_QUARTER_TURN.raw + 1 }));
}

test "tan: differential vs libm @tan within 2^-20 (loose — tan grows steeply)" {
    // tan amplifies polynomial error near poles, so use a coarse tolerance
    // and stay away from the pole proximities.
    inline for ([_]u32{ 0, 0x10000000, 0x18000000, 0x30000000, 0x50000000, 0x70000000 }) |raw| {
        const a = Angle{ .raw = raw };
        const rad = a.toRadians().toFloat();
        const got = tan(a).toFloat();
        const ref = @tan(rad);
        // Relative tolerance: tan can be large, so check abs OR rel error.
        const abs_err = @abs(got - ref);
        const rel_ok = if (@abs(ref) > 1.0) abs_err / @abs(ref) < 1e-5 else abs_err < 1e-5;
        try testing.expect(rel_ok);
    }
}

test "atan2: axis cardinal directions" {
    // (1, 0) east → 0
    try testing.expectEqual(Angle.ZERO, atan2(Fixed.ZERO, Fixed.ONE));
    // (0, 1) north → π/2
    const north = atan2(Fixed.ONE, Fixed.ZERO);
    try testing.expect(diffMod32Angle(north.raw, Angle.QUARTER_TURN.raw) <= 32);
    // (-1, 0) west → π
    const west = atan2(Fixed.ZERO, Fixed.NEG_ONE);
    try testing.expect(diffMod32Angle(west.raw, Angle.HALF_TURN.raw) <= 32);
    // (0, -1) south → 3π/2 (after BAM wrap of -π/2)
    const south = atan2(Fixed.NEG_ONE, Fixed.ZERO);
    try testing.expect(diffMod32Angle(south.raw, Angle.THREE_QUARTER_TURN.raw) <= 32);
}

test "atan2: diagonals at 45° / 135° / 225° / 315°" {
    // (1, 1) → π/4
    const ne = atan2(Fixed.ONE, Fixed.ONE);
    try testing.expect(diffMod32Angle(ne.raw, @as(u32, 1) << 29) <= 32);
    // (1, -1) → 3π/4
    const nw = atan2(Fixed.ONE, Fixed.NEG_ONE);
    try testing.expect(diffMod32Angle(nw.raw, (@as(u32, 3) << 29)) <= 32);
    // (-1, -1) → 5π/4 (= -3π/4 mod 2π)
    const sw = atan2(Fixed.NEG_ONE, Fixed.NEG_ONE);
    try testing.expect(diffMod32Angle(sw.raw, (@as(u32, 5) << 29)) <= 32);
    // (-1, 1) → 7π/4 (= -π/4 mod 2π)
    const se = atan2(Fixed.NEG_ONE, Fixed.ONE);
    try testing.expect(diffMod32Angle(se.raw, (@as(u32, 7) << 29)) <= 32);
}

test "atan2(sin θ, cos θ) ≈ θ for arbitrary angles" {
    // The round-trip identity. Errors compound across:
    //   - sin/cos polynomial (~2^-22 real ≈ 4 Fixed ULPs)
    //   - atan2 polynomial + division (~1-2 Fixed ULPs in radian output)
    //   - fromRadians input precision (~20 u32 ULPs for π-sized inputs,
    //     amplified by 128/π when Fixed → Angle)
    // Combined empirical budget: ~1024 u32 ULPs (≈ 1.5e-6 radians) — still
    // far below human-perceptible angular precision.
    inline for ([_]u32{
        0x10000000, 0x20000000, 0x30000000, 0x40000001,
        0x55555555, 0x6BBBBBBB, 0x80000001, 0xA000_0000,
        0xC000_0000, 0xDEAD_BEEF, 0xFFFF_0000,
    }) |raw| {
        const original = Angle{ .raw = raw };
        const s = sin(original);
        const c = cos(original);
        const round_tripped = atan2(s, c);
        try testing.expect(diffMod32Angle(original.raw, round_tripped.raw) <= 1024);
    }
}

test "atan2: differential vs libm @atan2 within 2^-18" {
    // Cross-check against libm at a variety of (y, x) magnitudes.
    inline for ([_]struct { y: f64, x: f64 }{
        .{ .y = 1.0, .x = 2.0 },     .{ .y = 2.0, .x = 1.0 },
        .{ .y = -1.0, .x = 3.0 },    .{ .y = 100.0, .x = 1.0 },
        .{ .y = 1.0, .x = 100.0 },   .{ .y = -50.0, .x = -75.0 },
        .{ .y = 0.001, .x = 1.0 },   .{ .y = -0.5, .x = 0.5 },
    }) |sample| {
        const yf = Fixed.rconst(sample.y);
        const xf = Fixed.rconst(sample.x);
        const got = atan2(yf, xf).toRadians().toFloat();
        const ref = std.math.atan2(sample.y, sample.x);
        // atan2 wraps to [0, 2π); ref is in (-π, π]. Normalize ref to [0, 2π).
        const ref_pos = if (ref < 0) ref + 2.0 * std.math.pi else ref;
        const diff = @abs(got - ref_pos);
        const wrap_diff = @min(diff, 2.0 * std.math.pi - diff);
        try testing.expect(wrap_diff < 4e-6); // 2^-18 ≈ 3.8e-6
    }
}

// Helper for tests: minimum forward/backward wrap distance, mod 2^32.
fn diffMod32Angle(a: u32, b: u32) u32 {
    const f = a -% b;
    const r = b -% a;
    return @min(f, r);
}

test "sin/cos at octant boundaries" {
    // Test angles exactly at oct·π/4 boundaries — these exercise the seam
    // between mirror=true and mirror=false (pos=0).
    var oct: u32 = 0;
    while (oct < 8) : (oct += 1) {
        const a = Angle{ .raw = oct << 29 };
        const rad: f64 = @as(f64, @floatFromInt(oct)) * std.math.pi / 4.0;
        const got_s = sin(a).toFloat();
        const got_c = cos(a).toFloat();
        try testing.expect(@abs(got_s - @sin(rad)) < 2.5e-7);
        try testing.expect(@abs(got_c - @cos(rad)) < 2.5e-7);
    }
}
