//! Narrowing-shift rounding for the determinism contract (SPEC §3).
//!
//! Every narrowing right-shift in fpz goes through `shiftRound`. It implements
//! round-to-nearest, ties-away-from-zero, on an `i128` intermediate. Callers
//! narrow the result to `i64` themselves and apply their own overflow check —
//! this helper has one job (the rounding rule) and never silently truncates.

const std = @import("std");

/// Shift `num` right by `k` bits with round-to-nearest, ties-away-from-zero.
///
/// Semantics:
///   k == 0: identity (returns `num`).
///   k  > 0: returns round_ties_away(num / 2^k) as i128.
///
/// Result is `i128` so the caller decides how to narrow and how to detect
/// overflow on the narrow. The helper itself cannot overflow for any input
/// satisfying `|num| + 2^(k-1) <= i128.max` — trivially true at fpz call sites
/// (largest `num` ≈ 2^126 from i64×i64; largest bias 2^(k-1) ≪ 2^126).
pub inline fn shiftRound(num: i128, comptime k: u7) i128 {
    if (k == 0) return num;
    const half: i128 = @as(i128, 1) << (k - 1);
    const biased: i128 = if (num >= 0) num + half else num - half;
    return @divTrunc(biased, @as(i128, 1) << k);
}

/// Narrow `x` from i128 to i64 with wrap-and-assert semantics (SPEC §3).
///
/// In Debug/ReleaseSafe: asserts the value fits — overflow is a loud, locatable
/// bug. In ReleaseFast: the assert is compiled out and `@truncate` returns the
/// defined 2's-complement wrap. Every build mode is total and deterministic.
pub inline fn narrow(x: i128) i64 {
    std.debug.assert(x >= std.math.minInt(i64) and x <= std.math.maxInt(i64));
    return @truncate(x);
}

/// Combined `shiftRound` + `narrow` — the common idiom for any op whose i128
/// intermediate must collapse back to i64 (mul, div, etc.).
pub inline fn shiftRoundNarrow(num: i128, comptime k: u7) i64 {
    return narrow(shiftRound(num, k));
}

test "narrow accepts in-range values" {
    try std.testing.expectEqual(@as(i64, 0), narrow(0));
    try std.testing.expectEqual(@as(i64, 42), narrow(42));
    try std.testing.expectEqual(@as(i64, -42), narrow(-42));
    try std.testing.expectEqual(std.math.maxInt(i64), narrow(std.math.maxInt(i64)));
    try std.testing.expectEqual(std.math.minInt(i64), narrow(std.math.minInt(i64)));
}

test "shiftRoundNarrow composes the two helpers" {
    try std.testing.expectEqual(@as(i64, 2), shiftRoundNarrow(3, 1)); // 1.5 → 2
    try std.testing.expectEqual(@as(i64, -2), shiftRoundNarrow(-3, 1));
    try std.testing.expectEqual(@as(i64, 0), shiftRoundNarrow(1, 2)); // 0.25 → 0
}

test "shiftRound k=0 is identity" {
    try std.testing.expectEqual(@as(i128, 0), shiftRound(0, 0));
    try std.testing.expectEqual(@as(i128, 42), shiftRound(42, 0));
    try std.testing.expectEqual(@as(i128, -42), shiftRound(-42, 0));
    try std.testing.expectEqual(
        @as(i128, std.math.maxInt(i128)),
        shiftRound(std.math.maxInt(i128), 0),
    );
    try std.testing.expectEqual(
        @as(i128, std.math.minInt(i128)),
        shiftRound(std.math.minInt(i128), 0),
    );
}

test "shiftRound exact divisions stay exact" {
    try std.testing.expectEqual(@as(i128, 1), shiftRound(2, 1));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-2, 1));
    try std.testing.expectEqual(@as(i128, 1), shiftRound(4, 2));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-4, 2));
    try std.testing.expectEqual(@as(i128, 0), shiftRound(0, 8));
    try std.testing.expectEqual(@as(i128, 5), shiftRound(20, 2));
    try std.testing.expectEqual(@as(i128, -5), shiftRound(-20, 2));
    // Larger powers
    try std.testing.expectEqual(@as(i128, 1), shiftRound(@as(i128, 1) << 24, 24));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-(@as(i128, 1) << 24), 24));
}

test "shiftRound halves go away from zero" {
    // 0.5 → 1, -0.5 → -1
    try std.testing.expectEqual(@as(i128, 1), shiftRound(1, 1));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-1, 1));
    // 1.5 → 2, -1.5 → -2
    try std.testing.expectEqual(@as(i128, 2), shiftRound(3, 1));
    try std.testing.expectEqual(@as(i128, -2), shiftRound(-3, 1));
    // 2.5 → 3, -2.5 → -3
    try std.testing.expectEqual(@as(i128, 3), shiftRound(5, 1));
    try std.testing.expectEqual(@as(i128, -3), shiftRound(-5, 1));
    // Half at k=2 (i.e., 0.5 represented as 2/4)
    try std.testing.expectEqual(@as(i128, 1), shiftRound(2, 2));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-2, 2));
}

test "shiftRound rounds toward zero below the half" {
    // 0.25 → 0, -0.25 → 0
    try std.testing.expectEqual(@as(i128, 0), shiftRound(1, 2));
    try std.testing.expectEqual(@as(i128, 0), shiftRound(-1, 2));
    // 1.25 → 1, -1.25 → -1
    try std.testing.expectEqual(@as(i128, 1), shiftRound(5, 2));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-5, 2));
}

test "shiftRound rounds away from zero above the half" {
    // 0.75 → 1, -0.75 → -1
    try std.testing.expectEqual(@as(i128, 1), shiftRound(3, 2));
    try std.testing.expectEqual(@as(i128, -1), shiftRound(-3, 2));
    // 1.75 → 2, -1.75 → -2
    try std.testing.expectEqual(@as(i128, 2), shiftRound(7, 2));
    try std.testing.expectEqual(@as(i128, -2), shiftRound(-7, 2));
}

test "shiftRound matches a half-away oracle exhaustively for small inputs" {
    // Oracle: scale to f64 (exact for these small magnitudes), apply
    // round-half-away-from-zero by branching on sign and using floor(|x|+0.5).
    inline for ([_]u7{ 1, 2, 3, 4, 5, 8, 16, 24 }) |k| {
        const div: i128 = @as(i128, 1) << k;
        const div_f: f64 = @floatFromInt(div);
        var n: i128 = -2048;
        while (n <= 2048) : (n += 1) {
            const got = shiftRound(n, k);
            const f: f64 = @as(f64, @floatFromInt(n)) / div_f;
            const oracle: i128 = if (f >= 0)
                @as(i128, @intFromFloat(@floor(f + 0.5)))
            else
                -@as(i128, @intFromFloat(@floor(-f + 0.5)));
            std.testing.expectEqual(oracle, got) catch |err| {
                std.debug.print("mismatch: n={d} k={d} got={d} oracle={d}\n", .{ n, k, got, oracle });
                return err;
            };
        }
    }
}

test "shiftRound is monotonic non-decreasing in num" {
    inline for ([_]u7{ 1, 4, 12, 24 }) |k| {
        var prev: i128 = shiftRound(-10_000, k);
        var n: i128 = -9_999;
        while (n <= 10_000) : (n += 1) {
            const cur = shiftRound(n, k);
            try std.testing.expect(cur >= prev);
            prev = cur;
        }
    }
}

test "shiftRound handles Q40.24 mul-intermediate magnitudes" {
    // i64.min * i64.min == 2^126 — the worst-case product. Fits in i128.
    const max_product: i128 = @as(i128, std.math.minInt(i64)) * @as(i128, std.math.minInt(i64));
    try std.testing.expectEqual(@as(i128, 1) << 126, max_product);
    // After shift-right by 24, exactly 2^102.
    try std.testing.expectEqual(@as(i128, 1) << 102, shiftRound(max_product, 24));

    // Mixed-sign extremes: i64.min * i64.max = -(2^126 - 2^63).
    const mixed: i128 = @as(i128, std.math.minInt(i64)) * @as(i128, std.math.maxInt(i64));
    const mixed_expected: i128 = -((@as(i128, 1) << 126) - (@as(i128, 1) << 63));
    try std.testing.expectEqual(mixed_expected, mixed);
    // Shift by 24 — the low 24 bits of |mixed| are zero, so the result is exact.
    const shifted_expected: i128 = -((@as(i128, 1) << 102) - (@as(i128, 1) << 39));
    try std.testing.expectEqual(shifted_expected, shiftRound(mixed, 24));
}
