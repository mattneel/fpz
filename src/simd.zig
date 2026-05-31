//! SIMD batch ops on Fixed (SPEC §9).
//!
//! Scalar `Fixed` is the canonical reference. Every SIMD op here must
//! produce bit-identical results to its scalar counterpart for any input
//! within the no-overflow precondition — verified by random differential
//! tests below.
//!
//! What vectorizes cleanly:
//!   - add / sub via `+%` / `-%` (wrapping arithmetic — matches scalar's
//!     defined wrap on overflow).
//!   - compare (eql, lt) via Zig's native vector comparison operators.
//!
//! What doesn't vectorize:
//!   - mul. The scalar uses a 64×64→128 widening multiply, which is largely
//!     absent on common ISAs (AVX-512-DQ's vpmullq is low-64 only). The
//!     compiler will scalarize the i128 product per lane. The point of
//!     having `mulVec` here is bit-identity guarantee, not throughput — if
//!     SIMD-mul throughput ever becomes the bottleneck, the escape hatch
//!     (per SPEC §9) is a 32-bit companion format for hot arrays.

const std = @import("std");
const Fixed = @import("Fixed.zig");
const rounding = @import("rounding.zig");

/// Vector type for N Fixed values (stored as their raw i64 representations).
/// `@Vector(N, i64)` so the LLVM backend can map to native SIMD where the
/// op vectorizes.
pub fn Vec(comptime N: comptime_int) type {
    return @Vector(N, i64);
}

pub fn fromArray(comptime N: comptime_int, arr: [N]Fixed) Vec(N) {
    var v: Vec(N) = undefined;
    inline for (0..N) |i| v[i] = arr[i].raw;
    return v;
}

pub fn toArray(comptime N: comptime_int, v: Vec(N)) [N]Fixed {
    var arr: [N]Fixed = undefined;
    inline for (0..N) |i| arr[i] = .{ .raw = v[i] };
    return arr;
}

/// Splat a single Fixed across all lanes.
pub fn splat(comptime N: comptime_int, x: Fixed) Vec(N) {
    return @splat(x.raw);
}

// ---------------------------------------------------------------------------
// Element-wise arithmetic — wrapping semantics match scalar's defined wrap
// (SPEC §3). The scalar's assert-on-overflow doesn't apply here: SIMD ops
// can't reasonably surface a per-lane assertion, so callers stay inside the
// no-overflow precondition. Tested bit-identity below covers in-range inputs.
// ---------------------------------------------------------------------------

pub fn addVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) Vec(N) {
    return a +% b;
}

pub fn subVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) Vec(N) {
    return a -% b;
}

pub fn negVec(comptime N: comptime_int, a: Vec(N)) Vec(N) {
    const zero: Vec(N) = @splat(0);
    return zero -% a;
}

/// Element-wise multiply with RNE narrowing — bit-identical to scalar.
/// Compiler will scalarize the i128 multiply per lane on most ISAs.
pub fn mulVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) Vec(N) {
    var result: Vec(N) = undefined;
    inline for (0..N) |i| {
        const product: i128 = @as(i128, a[i]) * @as(i128, b[i]);
        result[i] = rounding.shiftRoundNarrow(product, Fixed.frac_bits);
    }
    return result;
}

// ---------------------------------------------------------------------------
// Element-wise comparisons (return per-lane bool vector)
// ---------------------------------------------------------------------------

pub fn eqlVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) @Vector(N, bool) {
    return a == b;
}

pub fn ltVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) @Vector(N, bool) {
    return a < b;
}

pub fn lteVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) @Vector(N, bool) {
    return a <= b;
}

pub fn minVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) Vec(N) {
    return @min(a, b);
}

pub fn maxVec(comptime N: comptime_int, a: Vec(N), b: Vec(N)) Vec(N) {
    return @max(a, b);
}

// ===========================================================================
// Tests — bit-identity to scalar for every op
// ===========================================================================

const testing = std.testing;

// Random samples kept small enough that scalar add/sub/mul won't overflow,
// so the equality check is meaningful (in safe mode the scalar would panic
// on overflow; we don't want to confound the bit-identity test with that).
const MAGNITUDE_BITS: u6 = 28; // |raw| < 2^28 keeps add/sub/mul exact

fn fillRandom(comptime N: comptime_int, random: std.Random) [N]Fixed {
    var arr: [N]Fixed = undefined;
    inline for (0..N) |i| {
        const r = random.intRangeAtMost(i64, -(@as(i64, 1) << MAGNITUDE_BITS), @as(i64, 1) << MAGNITUDE_BITS);
        arr[i] = .{ .raw = r };
    }
    return arr;
}

test "simd addVec: bit-identical to scalar add" {
    var prng = std.Random.DefaultPrng.init(0xABCDEF);
    const random = prng.random();
    inline for (.{ 2, 4, 8, 16 }) |N| {
        var iter: u32 = 0;
        while (iter < 64) : (iter += 1) {
            const a_arr = fillRandom(N, random);
            const b_arr = fillRandom(N, random);
            const a_vec = fromArray(N, a_arr);
            const b_vec = fromArray(N, b_arr);
            const c_arr = toArray(N, addVec(N, a_vec, b_vec));
            inline for (0..N) |i| {
                try testing.expectEqual(a_arr[i].add(b_arr[i]), c_arr[i]);
            }
        }
    }
}

test "simd subVec: bit-identical to scalar sub" {
    var prng = std.Random.DefaultPrng.init(0x123456);
    const random = prng.random();
    inline for (.{ 2, 4, 8, 16 }) |N| {
        var iter: u32 = 0;
        while (iter < 64) : (iter += 1) {
            const a_arr = fillRandom(N, random);
            const b_arr = fillRandom(N, random);
            const a_vec = fromArray(N, a_arr);
            const b_vec = fromArray(N, b_arr);
            const c_arr = toArray(N, subVec(N, a_vec, b_vec));
            inline for (0..N) |i| {
                try testing.expectEqual(a_arr[i].sub(b_arr[i]), c_arr[i]);
            }
        }
    }
}

test "simd mulVec: bit-identical to scalar mul" {
    var prng = std.Random.DefaultPrng.init(0xDEADBEEF);
    const random = prng.random();
    inline for (.{ 2, 4, 8 }) |N| {
        var iter: u32 = 0;
        while (iter < 64) : (iter += 1) {
            const a_arr = fillRandom(N, random);
            const b_arr = fillRandom(N, random);
            const a_vec = fromArray(N, a_arr);
            const b_vec = fromArray(N, b_arr);
            const c_arr = toArray(N, mulVec(N, a_vec, b_vec));
            inline for (0..N) |i| {
                try testing.expectEqual(a_arr[i].mul(b_arr[i]), c_arr[i]);
            }
        }
    }
}

test "simd negVec: bit-identical to scalar neg" {
    var prng = std.Random.DefaultPrng.init(0xFEED);
    const random = prng.random();
    inline for (.{ 4, 8 }) |N| {
        const a_arr = fillRandom(N, random);
        const a_vec = fromArray(N, a_arr);
        const c_arr = toArray(N, negVec(N, a_vec));
        inline for (0..N) |i| {
            try testing.expectEqual(a_arr[i].neg(), c_arr[i]);
        }
    }
}

test "simd eqlVec / ltVec: bit-identical to scalar" {
    var prng = std.Random.DefaultPrng.init(42);
    const random = prng.random();
    inline for (.{ 4, 8 }) |N| {
        var iter: u32 = 0;
        while (iter < 32) : (iter += 1) {
            const a_arr = fillRandom(N, random);
            const b_arr = fillRandom(N, random);
            const a_vec = fromArray(N, a_arr);
            const b_vec = fromArray(N, b_arr);
            const eq_vec = eqlVec(N, a_vec, b_vec);
            const lt_vec = ltVec(N, a_vec, b_vec);
            inline for (0..N) |i| {
                try testing.expectEqual(a_arr[i].eql(b_arr[i]), eq_vec[i]);
                try testing.expectEqual(a_arr[i].lt(b_arr[i]), lt_vec[i]);
            }
        }
    }
}

test "simd minVec / maxVec: bit-identical to scalar" {
    var prng = std.Random.DefaultPrng.init(0x55);
    const random = prng.random();
    inline for (.{ 4, 8 }) |N| {
        const a_arr = fillRandom(N, random);
        const b_arr = fillRandom(N, random);
        const a_vec = fromArray(N, a_arr);
        const b_vec = fromArray(N, b_arr);
        const min_arr = toArray(N, minVec(N, a_vec, b_vec));
        const max_arr = toArray(N, maxVec(N, a_vec, b_vec));
        inline for (0..N) |i| {
            try testing.expectEqual(a_arr[i].min(b_arr[i]), min_arr[i]);
            try testing.expectEqual(a_arr[i].max(b_arr[i]), max_arr[i]);
        }
    }
}

test "simd splat / fromArray / toArray round-trips" {
    const x = Fixed.fromInt(42);
    const v: Vec(4) = splat(4, x);
    const arr = toArray(4, v);
    inline for (0..4) |i| try testing.expectEqual(x, arr[i]);

    const samples = [_]Fixed{ Fixed.ONE, Fixed.HALF, Fixed.NEG_ONE, Fixed.ZERO };
    const v2 = fromArray(4, samples);
    const back = toArray(4, v2);
    inline for (0..4) |i| try testing.expectEqual(samples[i], back[i]);
}
