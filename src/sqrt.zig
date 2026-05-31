//! Integer-only square root for Fixed (SPEC §6).
//!
//! Digit-by-digit (restoring) integer sqrt over a u128 promotion of
//! `raw << frac_bits`. Exact — no Newton convergence-count guessing —
//! correctly rounded to ±0.5 ULP (beating the spec's ±1 ULP target,
//! because for integer inputs there are no exact half-way cases).
//!
//! Native integer math. No libm, no runtime float.

const std = @import("std");
const Fixed = @import("Fixed.zig");

/// Integer square root with round-to-nearest.
///
/// Returns `round(sqrt(n))`. For an integer `n`, the value `sqrt(n)` is
/// never exactly halfway between two integers, so RNE is unambiguous.
pub fn isqrtRne(n: u128) u64 {
    if (n == 0) return 0;

    // Restoring digit-by-digit isqrt. Process two bits at a time from the
    // top: `bit` is the squared-bit position currently being decided.
    var rem: u128 = n;
    var root: u128 = 0;
    var bit: u128 = @as(u128, 1) << 126;
    while (bit > n) bit >>= 2;

    while (bit != 0) {
        if (rem >= root + bit) {
            rem -= root + bit;
            root = (root >> 1) + bit;
        } else {
            root >>= 1;
        }
        bit >>= 2;
    }

    // After the loop: root = floor(sqrt(n)), rem = n - root^2.
    // We need round(sqrt(n)). Round up iff sqrt(n) >= root + 0.5
    //   iff n >= (root + 0.5)^2 = root^2 + root + 0.25
    //   iff rem >= root + 0.25  (since rem is integer)
    //   iff rem >  root          (integer comparison)
    if (rem > root) root += 1;
    return @intCast(root);
}

/// sqrt(x) for Fixed.
///
/// Precondition: `x.raw >= 0`. In safe builds the assert fires loudly.
/// In ReleaseFast the assert is compiled out; negative input returns
/// `Fixed.ZERO` deterministically (SPEC §7).
pub fn sqrt(x: Fixed) Fixed {
    std.debug.assert(x.raw >= 0);
    if (x.raw < 0) return Fixed.ZERO;
    const n: u128 = @as(u128, @intCast(x.raw)) << Fixed.frac_bits;
    return .{ .raw = @intCast(isqrtRne(n)) };
}

/// sqrt(x), checked variant returning an error union for callers that want
/// to handle the domain violation explicitly.
pub fn sqrtChecked(x: Fixed) error{Negative}!Fixed {
    if (x.raw < 0) return error.Negative;
    const n: u128 = @as(u128, @intCast(x.raw)) << Fixed.frac_bits;
    return .{ .raw = @intCast(isqrtRne(n)) };
}

// ===========================================================================
// Tests
// ===========================================================================

const testing = std.testing;

test "isqrtRne: zero" {
    try testing.expectEqual(@as(u64, 0), isqrtRne(0));
}

test "isqrtRne: exact on perfect squares" {
    try testing.expectEqual(@as(u64, 1), isqrtRne(1));
    try testing.expectEqual(@as(u64, 2), isqrtRne(4));
    try testing.expectEqual(@as(u64, 3), isqrtRne(9));
    try testing.expectEqual(@as(u64, 10), isqrtRne(100));
    try testing.expectEqual(@as(u64, 100), isqrtRne(10_000));
    try testing.expectEqual(@as(u64, 1_000_000), isqrtRne(1_000_000_000_000));
    // Powers of 2 with even exponent.
    try testing.expectEqual(@as(u64, 1) << 20, isqrtRne(@as(u128, 1) << 40));
    try testing.expectEqual(@as(u64, 1) << 50, isqrtRne(@as(u128, 1) << 100));
}

test "isqrtRne: rounds to nearest (small non-squares)" {
    // sqrt(2) ≈ 1.414 → 1
    try testing.expectEqual(@as(u64, 1), isqrtRne(2));
    // sqrt(3) ≈ 1.732 → 2
    try testing.expectEqual(@as(u64, 2), isqrtRne(3));
    // sqrt(5) ≈ 2.236 → 2
    try testing.expectEqual(@as(u64, 2), isqrtRne(5));
    // sqrt(6) ≈ 2.449 → 2
    try testing.expectEqual(@as(u64, 2), isqrtRne(6));
    // sqrt(7) ≈ 2.646 → 3
    try testing.expectEqual(@as(u64, 3), isqrtRne(7));
    // sqrt(8) ≈ 2.828 → 3
    try testing.expectEqual(@as(u64, 3), isqrtRne(8));
    // sqrt(12) ≈ 3.464 → 3
    try testing.expectEqual(@as(u64, 3), isqrtRne(12));
    // sqrt(13) ≈ 3.606 → 4
    try testing.expectEqual(@as(u64, 4), isqrtRne(13));
}

test "isqrtRne: monotone non-decreasing over [0, 10000]" {
    var n: u128 = 0;
    var prev: u64 = 0;
    while (n <= 10_000) : (n += 1) {
        const got = isqrtRne(n);
        try testing.expect(got >= prev);
        // |got^2 - n| should be at most got (within ±0.5 ULP at integer scale).
        const got2: u128 = @as(u128, got) * @as(u128, got);
        const diff = if (got2 > n) got2 - n else n - got2;
        try testing.expect(diff <= got);
        prev = got;
    }
}

test "isqrtRne: handles large u128 close to fpz upper bound (2^88)" {
    // The biggest input fpz sqrt will see: i64.max << 24 ≈ 2^87.
    const big: u128 = (@as(u128, 1) << 87) - 1;
    const r = isqrtRne(big);
    const r_sq: u128 = @as(u128, r) * @as(u128, r);
    const diff = if (r_sq > big) r_sq - big else big - r_sq;
    // Within ±r of big — the ±0.5 ULP guarantee at integer scale.
    try testing.expect(diff <= r);
}

test "Fixed.sqrt: cardinal values are exact" {
    try testing.expectEqual(Fixed.ZERO, sqrt(Fixed.ZERO));
    try testing.expectEqual(Fixed.ONE, sqrt(Fixed.ONE));
    try testing.expectEqual(Fixed.fromInt(2), sqrt(Fixed.fromInt(4)));
    try testing.expectEqual(Fixed.fromInt(3), sqrt(Fixed.fromInt(9)));
    try testing.expectEqual(Fixed.fromInt(10), sqrt(Fixed.fromInt(100)));
    // 0.25 = (0.5)^2
    try testing.expectEqual(Fixed.HALF, sqrt(Fixed.rconst(0.25)));
}

test "Fixed.sqrt: sqrt(x)^2 ≈ x within an error bound that scales with sqrt(x)" {
    // Error analysis: sqrt has at most 0.5 ULP error in sqrt-units. Squaring
    // propagates as d(y²) = 2y·dy, so the error in sx² is bounded by
    // ~2·sqrt(x) ULPs in Fixed.raw units. Add 4 ULPs slack for the mul's
    // own RNE rounding.
    inline for (.{ 0.25, 0.5, 1.0, 2.0, 3.0, 7.5, 100.0, 1_000.0, 1_000_000.0 }) |f| {
        const x = Fixed.rconst(f);
        const sx = sqrt(x);
        const sx_sq = sx.mul(sx);
        const diff: u64 = @abs(x.raw - sx_sq.raw);
        // 2 * sx_real = sx.raw >> (frac_bits - 1) — bound in Fixed.raw units.
        const bound: u64 = (@as(u64, @intCast(sx.raw)) >> (Fixed.frac_bits - 1)) + 4;
        try testing.expect(diff <= bound);
    }
}

test "Fixed.sqrt: monotone non-decreasing over small raw values" {
    var raw: i64 = 0;
    var prev: i64 = 0;
    while (raw <= 1024) : (raw += 1) {
        const r = sqrt(Fixed{ .raw = raw }).raw;
        try testing.expect(r >= prev);
        prev = r;
    }
}

test "Fixed.sqrt: differential vs @sqrt(f64) within 1e-6 relative error" {
    inline for (.{ 1.0, 2.0, 3.0, 7.5, 100.0, 1e6, 1e9 }) |f| {
        const x = Fixed.rconst(f);
        const my_sqrt = sqrt(x).toFloat();
        const ref = @sqrt(f);
        const rel_err = @abs(my_sqrt - ref) / ref;
        try testing.expect(rel_err < 1e-6);
    }
}

test "sqrtChecked: returns Negative for negative input, value otherwise" {
    try testing.expectError(error.Negative, sqrtChecked(Fixed.fromInt(-1)));
    try testing.expectError(error.Negative, sqrtChecked(Fixed{ .raw = -1 }));
    try testing.expectError(error.Negative, sqrtChecked(Fixed.MIN));

    try testing.expectEqual(Fixed.ZERO, try sqrtChecked(Fixed.ZERO));
    try testing.expectEqual(Fixed.ONE, try sqrtChecked(Fixed.ONE));
    try testing.expectEqual(Fixed.fromInt(7), try sqrtChecked(Fixed.fromInt(49)));
}
