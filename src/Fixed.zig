//! The fpz Fixed-point type (SPEC §1, §3, §4).
//!
//! File-is-type: `@import("Fixed.zig")` returns the struct itself. The format
//! is determined by the comptime constant `frac_bits` (default Q40.24).
//!
//! Distinct struct (not an i64 alias), method-only API (Zig has no operator
//! overloading): the type system blocks mixing raw integers with Fixed values
//! and blocks accidental integer `+`/`*` on the raw bits.

const std = @import("std");
const rounding = @import("rounding.zig");

const Fixed = @This();

// ---------------------------------------------------------------------------
// Format
// ---------------------------------------------------------------------------

/// Bits below the binary point. Q40.24 by default — see SPEC §1 for the
/// range/precision tradeoff. Changing this is a breaking change that requires
/// re-baselining the conformance vectors.
pub const frac_bits: comptime_int = 24;
pub const whole_bits: comptime_int = 64 - frac_bits;

/// Range of integers exactly representable by `fromInt`: [min_int, max_int].
pub const max_int: i64 = (@as(i64, 1) << (whole_bits - 1)) - 1;
pub const min_int: i64 = -(@as(i64, 1) << (whole_bits - 1));

// ---------------------------------------------------------------------------
// The single field
// ---------------------------------------------------------------------------

raw: i64,

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const ZERO: Fixed = .{ .raw = 0 };
pub const ONE: Fixed = .{ .raw = @as(i64, 1) << frac_bits };
pub const HALF: Fixed = .{ .raw = @as(i64, 1) << (frac_bits - 1) };
pub const MIN: Fixed = .{ .raw = std.math.minInt(i64) };
pub const MAX: Fixed = .{ .raw = std.math.maxInt(i64) };
pub const NEG_ONE: Fixed = .{ .raw = -(@as(i64, 1) << frac_bits) };

pub const PI: Fixed = rconst(std.math.pi);
pub const TWO_PI: Fixed = rconst(std.math.tau);
pub const HALF_PI: Fixed = rconst(std.math.pi / 2.0);
pub const E: Fixed = rconst(std.math.e);
pub const LN2: Fixed = rconst(std.math.ln2);

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

/// Construct a Fixed from an integer in [min_int, max_int].
pub fn fromInt(i: i64) Fixed {
    std.debug.assert(i >= min_int and i <= max_int);
    return .{ .raw = i << frac_bits };
}

/// Construct a Fixed from a raw i64 representation. Use sparingly — most
/// callers should prefer `fromInt` or `rconst`. Useful for deserializing
/// already-validated values.
pub fn fromRaw(r: i64) Fixed {
    return .{ .raw = r };
}

/// Comptime float → Fixed. Runs the float math in the compiler; bakes the
/// integer result into the binary. SPEC §4: comptime-only by signature
/// (`comptime r: f64`), so calling with a runtime float fails to compile.
/// Round-to-nearest, ties away from zero.
pub fn rconst(comptime r: f64) Fixed {
    @setEvalBranchQuota(10_000);
    const scale: f64 = @floatFromInt(@as(i64, 1) << frac_bits);
    const scaled: f64 = r * scale;
    const rounded_f: f64 = if (scaled >= 0)
        @floor(scaled + 0.5)
    else
        -@floor(-scaled + 0.5);
    return .{ .raw = @intFromFloat(rounded_f) };
}

// ---------------------------------------------------------------------------
// Conversion / inspection
// ---------------------------------------------------------------------------

/// Truncate toward zero, returning the integer part.
pub fn toInt(self: Fixed) i64 {
    // Arithmetic right shift floors toward -∞; explicit trunc-toward-zero
    // requires sign-aware handling.
    if (self.raw >= 0) return self.raw >> frac_bits;
    // For negatives, divTrunc by 2^frac_bits gives trunc-toward-zero.
    return @divTrunc(self.raw, @as(i64, 1) << frac_bits);
}

/// Round to the nearest integer, ties away from zero.
pub fn roundToInt(self: Fixed) i64 {
    // Input |self.raw| <= 2^63, shift right by 24 → magnitude ≤ 2^39: fits i64.
    return @intCast(rounding.shiftRound(@as(i128, self.raw), frac_bits));
}

/// Lossy convert to f64. **Debug / display only.** Never call from a sim path.
pub fn toFloat(self: Fixed) f64 {
    const scale: f64 = @floatFromInt(@as(i64, 1) << frac_bits);
    return @as(f64, @floatFromInt(self.raw)) / scale;
}

// ---------------------------------------------------------------------------
// Arithmetic (wrap-and-assert — SPEC §3)
// ---------------------------------------------------------------------------

pub fn add(a: Fixed, b: Fixed) Fixed {
    const r = @addWithOverflow(a.raw, b.raw);
    std.debug.assert(r[1] == 0);
    return .{ .raw = r[0] };
}

pub fn sub(a: Fixed, b: Fixed) Fixed {
    const r = @subWithOverflow(a.raw, b.raw);
    std.debug.assert(r[1] == 0);
    return .{ .raw = r[0] };
}

pub fn mul(a: Fixed, b: Fixed) Fixed {
    // i64 × i64 → i128 is exact (max |product| = 2^126). Shift-right by
    // frac_bits with RNE, then narrow to i64 with overflow assert.
    const product: i128 = @as(i128, a.raw) * @as(i128, b.raw);
    return .{ .raw = rounding.shiftRoundNarrow(product, frac_bits) };
}

pub fn div(a: Fixed, b: Fixed) Fixed {
    std.debug.assert(b.raw != 0);
    // We want round((a/2^F) / (b/2^F)) * 2^F = round(a/b * 2^F).
    // Trick: scale `a` by 2^(F+1), divTrunc by `b`, then shiftRound by 1
    // to apply RNE-ties-away on the final bit. Single rounding step keeps
    // the determinism contract intact.
    const numerator: i128 = @as(i128, a.raw) << (frac_bits + 1);
    const quotient_2x: i128 = @divTrunc(numerator, @as(i128, b.raw));
    return .{ .raw = rounding.shiftRoundNarrow(quotient_2x, 1) };
}

pub fn neg(a: Fixed) Fixed {
    // For a.raw == i64.min the negation overflows; defined wrap is i64.min
    // itself. Assert in safe builds (loud bug), wrap in fast (deterministic).
    std.debug.assert(a.raw != std.math.minInt(i64));
    return .{ .raw = -%a.raw };
}

pub fn abs(a: Fixed) Fixed {
    std.debug.assert(a.raw != std.math.minInt(i64));
    return .{ .raw = if (a.raw < 0) -%a.raw else a.raw };
}

// ---------------------------------------------------------------------------
// Saturating arithmetic (SPEC §3) — for genuine singularities
// ---------------------------------------------------------------------------

pub fn addSat(a: Fixed, b: Fixed) Fixed {
    return .{ .raw = a.raw +| b.raw };
}

pub fn subSat(a: Fixed, b: Fixed) Fixed {
    return .{ .raw = a.raw -| b.raw };
}

pub fn mulSat(a: Fixed, b: Fixed) Fixed {
    const product: i128 = @as(i128, a.raw) * @as(i128, b.raw);
    const shifted: i128 = rounding.shiftRound(product, frac_bits);
    if (shifted > std.math.maxInt(i64)) return MAX;
    if (shifted < std.math.minInt(i64)) return MIN;
    return .{ .raw = @truncate(shifted) };
}

// ---------------------------------------------------------------------------
// Comparators
// ---------------------------------------------------------------------------

pub fn eql(a: Fixed, b: Fixed) bool {
    return a.raw == b.raw;
}

pub fn lt(a: Fixed, b: Fixed) bool {
    return a.raw < b.raw;
}

pub fn lte(a: Fixed, b: Fixed) bool {
    return a.raw <= b.raw;
}

pub fn cmp(a: Fixed, b: Fixed) std.math.Order {
    return std.math.order(a.raw, b.raw);
}

pub fn min(a: Fixed, b: Fixed) Fixed {
    return if (a.raw < b.raw) a else b;
}

pub fn max(a: Fixed, b: Fixed) Fixed {
    return if (a.raw > b.raw) a else b;
}

pub fn clamp(self: Fixed, lo: Fixed, hi: Fixed) Fixed {
    std.debug.assert(lo.raw <= hi.raw);
    return .{ .raw = std.math.clamp(self.raw, lo.raw, hi.raw) };
}

// ===========================================================================
// Tests
// ===========================================================================

const testing = std.testing;

test "Fixed has i64 layout (size + alignment)" {
    try testing.expectEqual(@sizeOf(i64), @sizeOf(Fixed));
    try testing.expectEqual(@alignOf(i64), @alignOf(Fixed));
}

test "constants: ZERO / ONE / HALF / NEG_ONE raw values" {
    try testing.expectEqual(@as(i64, 0), ZERO.raw);
    try testing.expectEqual(@as(i64, 1) << frac_bits, ONE.raw);
    try testing.expectEqual(@as(i64, 1) << (frac_bits - 1), HALF.raw);
    try testing.expectEqual(-(@as(i64, 1) << frac_bits), NEG_ONE.raw);
}

test "constants: MIN / MAX match i64 extremes" {
    try testing.expectEqual(std.math.minInt(i64), MIN.raw);
    try testing.expectEqual(std.math.maxInt(i64), MAX.raw);
}

// GOLDEN: each constant's raw value is part of the determinism contract
// (SPEC §4). Any drift means rconst rounded differently or comptime float
// math changed — a contract break that requires re-baselining the
// conformance vectors. Split into separate tests so a divergence on one
// constant reports the actual value of all the others.

test "PI bit-stable" {
    try testing.expectEqual(@as(i64, 52707179), PI.raw);
}
test "TWO_PI bit-stable" {
    try testing.expectEqual(@as(i64, 105414357), TWO_PI.raw);
}
test "HALF_PI bit-stable" {
    try testing.expectEqual(@as(i64, 26353589), HALF_PI.raw);
}
test "E bit-stable" {
    try testing.expectEqual(@as(i64, 45605201), E.raw);
}
test "LN2 bit-stable" {
    try testing.expectEqual(@as(i64, 11629080), LN2.raw);
}

test "constants: 2 × HALF_PI relates to PI within 1 ULP" {
    // PI's rconst rounds; HALF_PI doubled may differ by 0 or 1 in the raw.
    const doubled = HALF_PI.raw * 2;
    try testing.expect(@abs(doubled - PI.raw) <= 1);
}

test "fromInt / toInt round-trip across the representable range" {
    const samples = [_]i64{ 0, 1, -1, 2, -2, 100, -100, max_int, min_int, max_int - 1, min_int + 1 };
    for (samples) |i| {
        try testing.expectEqual(i, fromInt(i).toInt());
    }
}

test "fromInt scales by 2^frac_bits" {
    try testing.expectEqual(@as(i64, 0), fromInt(0).raw);
    try testing.expectEqual(@as(i64, 1) << frac_bits, fromInt(1).raw);
    try testing.expectEqual(@as(i64, 5) << frac_bits, fromInt(5).raw);
    try testing.expectEqual(-(@as(i64, 7) << frac_bits), fromInt(-7).raw);
}

test "toInt truncates toward zero (not RNE)" {
    // 1.5 → 1, -1.5 → -1
    try testing.expectEqual(@as(i64, 1), (Fixed{ .raw = (@as(i64, 1) << frac_bits) + HALF.raw }).toInt());
    try testing.expectEqual(@as(i64, -1), (Fixed{ .raw = -(@as(i64, 1) << frac_bits) - HALF.raw }).toInt());
    // 0.9 → 0, -0.9 → 0
    try testing.expectEqual(@as(i64, 0), (Fixed{ .raw = HALF.raw + HALF.raw / 2 }).toInt());
    try testing.expectEqual(@as(i64, 0), (Fixed{ .raw = -(HALF.raw + HALF.raw / 2) }).toInt());
}

test "roundToInt rounds half away from zero" {
    // 0.5 → 1, -0.5 → -1
    try testing.expectEqual(@as(i64, 1), HALF.roundToInt());
    try testing.expectEqual(@as(i64, -1), (Fixed{ .raw = -HALF.raw }).roundToInt());
    // 1.5 → 2, -1.5 → -2
    try testing.expectEqual(@as(i64, 2), (Fixed{ .raw = ONE.raw + HALF.raw }).roundToInt());
    try testing.expectEqual(@as(i64, -2), (Fixed{ .raw = -(ONE.raw + HALF.raw) }).roundToInt());
    // 0.49ish → 0
    try testing.expectEqual(@as(i64, 0), (Fixed{ .raw = HALF.raw - 1 }).roundToInt());
    try testing.expectEqual(@as(i64, 0), (Fixed{ .raw = -(HALF.raw - 1) }).roundToInt());
}

test "rconst matches expected rounding behavior" {
    // 0.5 (exact in binary)
    try testing.expectEqual(HALF, rconst(0.5));
    try testing.expectEqual(NEG_ONE, rconst(-1.0));
    try testing.expectEqual(ONE, rconst(1.0));
    try testing.expectEqual(ZERO, rconst(0.0));
    // Negation symmetry for exactly representable values
    try testing.expectEqual(-HALF.raw, rconst(-0.5).raw);
}

test "toFloat is approximate inverse of rconst" {
    // `rconst` requires a comptime float (deliberate — keeps runtime float
    // off sim paths). Use `inline for` so each `f` is comptime-known.
    inline for (.{ 0.0, 0.5, -0.5, 1.0, -1.0, 3.14159, -2.71828, 100.25 }) |f| {
        const x = rconst(f);
        const back = x.toFloat();
        try testing.expect(@abs(back - f) < 1e-6);
    }
}

test "add: identity, commutativity, known values" {
    try testing.expectEqual(ONE, ONE.add(ZERO));
    try testing.expectEqual(ONE, ZERO.add(ONE));
    try testing.expectEqual(ONE.add(ONE), ONE.add(ONE)); // determinism
    try testing.expectEqual(Fixed{ .raw = 5 }, (Fixed{ .raw = 2 }).add(.{ .raw = 3 }));
    try testing.expectEqual(Fixed{ .raw = -1 }, (Fixed{ .raw = 2 }).add(.{ .raw = -3 }));
    // Commutativity
    const a = rconst(1.25);
    const b = rconst(-3.75);
    try testing.expectEqual(a.add(b), b.add(a));
}

test "sub: a - a == ZERO, identity, known values" {
    try testing.expectEqual(ZERO, ONE.sub(ONE));
    try testing.expectEqual(ZERO, MAX.sub(MAX));
    try testing.expectEqual(ONE, ONE.sub(ZERO));
    try testing.expectEqual(Fixed{ .raw = -1 }, (Fixed{ .raw = 2 }).sub(.{ .raw = 3 }));
}

test "mul: identity, ONE × x == x, RNE at the narrowing shift" {
    try testing.expectEqual(ZERO, ZERO.mul(ZERO));
    try testing.expectEqual(ZERO, ONE.mul(ZERO));
    try testing.expectEqual(ONE, ONE.mul(ONE));
    try testing.expectEqual(HALF, HALF.mul(ONE));
    try testing.expectEqual(HALF, ONE.mul(HALF));
    // HALF × HALF == 0.25 exactly
    try testing.expectEqual(rconst(0.25), HALF.mul(HALF));
    // Sign flips
    try testing.expectEqual(NEG_ONE, ONE.mul(NEG_ONE));
    try testing.expectEqual(ONE, NEG_ONE.mul(NEG_ONE));
}

test "mul: exactly-half rounds away from zero (RNE-ties-away)" {
    // a = raw 1 (smallest positive ULP), b = HALF.
    // product = 2^23; shiftRound(2^23, 24) = exactly 0.5 → rounds away to 1.
    const ulp = Fixed{ .raw = 1 };
    try testing.expectEqual(@as(i64, 1), ulp.mul(HALF).raw);
    try testing.expectEqual(@as(i64, -1), ulp.mul(.{ .raw = -HALF.raw }).raw);
    // Just below half: rounds to 0.
    try testing.expectEqual(@as(i64, 0), ulp.mul(.{ .raw = HALF.raw - 1 }).raw);
    try testing.expectEqual(@as(i64, 0), ulp.mul(.{ .raw = -(HALF.raw - 1) }).raw);
    // Just above half: rounds away.
    try testing.expectEqual(@as(i64, 1), ulp.mul(.{ .raw = HALF.raw + 1 }).raw);
    try testing.expectEqual(@as(i64, -1), ulp.mul(.{ .raw = -(HALF.raw + 1) }).raw);
}

test "div: identity, divide-by-self, divide-by-ONE" {
    try testing.expectEqual(ONE, ONE.div(ONE));
    try testing.expectEqual(HALF, HALF.div(ONE));
    try testing.expectEqual(ONE, HALF.div(HALF));
    try testing.expectEqual(rconst(0.5), ONE.div(rconst(2.0)));
    try testing.expectEqual(rconst(-0.5), ONE.div(rconst(-2.0)));
    try testing.expectEqual(rconst(2.0), ONE.div(HALF));
}

test "div: 1/3 rounds correctly" {
    const third = ONE.div(fromInt(3));
    // 1/3 ≈ 0.333..., × 2^24 ≈ 5592405.33 → trunc would give 5592405,
    // RNE-away with fractional 0.33 also gives 5592405.
    try testing.expectEqual(@as(i64, 5592405), third.raw);
}

test "div: sign combinations" {
    const a = fromInt(7);
    const b = fromInt(2);
    try testing.expectEqual(rconst(3.5), a.div(b));
    try testing.expectEqual(rconst(-3.5), a.div(b.neg()));
    try testing.expectEqual(rconst(-3.5), a.neg().div(b));
    try testing.expectEqual(rconst(3.5), a.neg().div(b.neg()));
}

test "neg / abs" {
    try testing.expectEqual(ZERO, ZERO.neg());
    try testing.expectEqual(NEG_ONE, ONE.neg());
    try testing.expectEqual(ONE, NEG_ONE.neg());
    try testing.expectEqual(ONE, ONE.abs());
    try testing.expectEqual(ONE, NEG_ONE.abs());
    try testing.expectEqual(ZERO, ZERO.abs());
}

test "addSat / subSat clamp at i64 extremes" {
    try testing.expectEqual(MAX, MAX.addSat(ONE));
    try testing.expectEqual(MAX, ONE.addSat(MAX));
    try testing.expectEqual(MIN, MIN.subSat(ONE));
    try testing.expectEqual(MIN, MIN.addSat(NEG_ONE));
    // In-range stays exact
    try testing.expectEqual(Fixed{ .raw = 5 }, (Fixed{ .raw = 2 }).addSat(.{ .raw = 3 }));
}

test "mulSat clamps on overflow, exact otherwise" {
    // Squaring a value near the max overflows → saturates to MAX
    const big = fromInt(max_int);
    try testing.expectEqual(MAX, big.mulSat(big));
    // Negative × positive overflow → saturates to MIN
    try testing.expectEqual(MIN, big.mulSat(fromInt(min_int)));
    // In-range exact
    try testing.expectEqual(rconst(6.0), rconst(2.0).mulSat(rconst(3.0)));
    try testing.expectEqual(rconst(0.25), HALF.mulSat(HALF));
}

test "comparators: eql / lt / lte / cmp" {
    try testing.expect(ONE.eql(ONE));
    try testing.expect(!ONE.eql(HALF));
    try testing.expect(HALF.lt(ONE));
    try testing.expect(!ONE.lt(HALF));
    try testing.expect(!ONE.lt(ONE));
    try testing.expect(ONE.lte(ONE));
    try testing.expect(HALF.lte(ONE));
    try testing.expect(!ONE.lte(HALF));

    try testing.expectEqual(std.math.Order.eq, ONE.cmp(ONE));
    try testing.expectEqual(std.math.Order.lt, HALF.cmp(ONE));
    try testing.expectEqual(std.math.Order.gt, ONE.cmp(HALF));
}

test "min / max / clamp" {
    try testing.expectEqual(HALF, HALF.min(ONE));
    try testing.expectEqual(HALF, ONE.min(HALF));
    try testing.expectEqual(ONE, HALF.max(ONE));
    try testing.expectEqual(ONE, ONE.max(HALF));

    try testing.expectEqual(HALF, HALF.clamp(ZERO, ONE));
    try testing.expectEqual(ZERO, NEG_ONE.clamp(ZERO, ONE));
    try testing.expectEqual(ONE, rconst(2.0).clamp(ZERO, ONE));
}

test "add/sub round-trip property" {
    const samples = [_]Fixed{
        ZERO, ONE, HALF, NEG_ONE, rconst(3.14), rconst(-2.71), rconst(0.001),
    };
    for (samples) |a| {
        for (samples) |b| {
            try testing.expectEqual(a, a.add(b).sub(b));
            try testing.expectEqual(a, a.sub(b).add(b));
        }
    }
}

test "rconst differential vs runtime computation within tolerance" {
    // Cross-check rconst against (toFloat ∘ rconst): bit-exact round trip
    // for values that don't exceed ULP at this format.
    const f = 1.5; // exactly representable in binary
    try testing.expectEqual(f, rconst(f).toFloat());
    // Values not exactly representable should round-trip within 1 ULP.
    const approx = rconst(0.1).toFloat();
    try testing.expect(@abs(approx - 0.1) < 1.0 / @as(f64, @floatFromInt(@as(i64, 1) << frac_bits)));
}
