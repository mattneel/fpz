//! Angle — Binary Angle Measure (SPEC §5).
//!
//! File-is-type. A full turn maps onto the full u32 range, so angle arithmetic
//! wraps correctly by construction (u32 +%/-% is mod 2^32, which IS mod 2π).
//! This eliminates an entire class of angle-normalization bugs and removes any
//! need to range-reduce by a low-precision π.
//!
//! Resolution: 360° / 2^32 ≈ 8.4e-8 degrees — far finer than any trig
//! approximation's accuracy, so conversion ULP errors are below the noise.
//!
//! Trig (SPEC §6) consumes Angle directly: range reduction becomes bit
//! extraction (the high 2–3 bits select the quadrant/octant).

const std = @import("std");
const Fixed = @import("Fixed.zig");
const rounding = @import("rounding.zig");

const Angle = @This();

raw: u32,

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const ZERO: Angle = .{ .raw = 0 };
pub const QUARTER_TURN: Angle = .{ .raw = @as(u32, 1) << 30 }; // π/2
pub const HALF_TURN: Angle = .{ .raw = @as(u32, 1) << 31 }; // π
pub const THREE_QUARTER_TURN: Angle = .{ .raw = @as(u32, 3) << 30 }; // 3π/2

// ---------------------------------------------------------------------------
// Comptime scaling constants
//
// All conversions use a precomputed integer multiplier at high fixed-point
// precision (Q57–Q68). Comptime_float internally is f128 (112-bit mantissa)
// so the rounded multipliers carry far more precision than the u32 angle's
// 32 bits — conversion error is below the angle's own resolution.
// ---------------------------------------------------------------------------

/// `angle.raw = round(r.raw × (128/π) / 2^57) mod 2^32`.
/// Derivation: angle = r × (2^32 / 2π); substitute r.raw = r × 2^24 → angle
/// = r.raw × 2^8 / (2π) = r.raw × 128/π. We carry 128/π at Q57 so the i128
/// product `r.raw × K` has plenty of headroom (max ≈ 2^126).
const K_FROM_RAD_Q57: i128 = blk: {
    @setEvalBranchQuota(10_000);
    const k: comptime_float = @as(comptime_float, 128.0) / std.math.pi;
    break :blk @intFromFloat(@floor(k * @as(comptime_float, 1 << 57) + 0.5));
};

/// `r.raw = round(a.raw × (π/128) × 2^68) / 2^68`.
const K_TO_RAD_Q68: i128 = blk: {
    @setEvalBranchQuota(10_000);
    const k: comptime_float = std.math.pi / @as(comptime_float, 128.0);
    break :blk @intFromFloat(@floor(k * @as(comptime_float, 1 << 68) + 0.5));
};

/// `angle.raw = round(d.raw × (32/45) × 2^63 / 2^63) mod 2^32`.
/// (32/45 = (2^32/360)/2^24 — same derivation as radians but with 360 in
/// place of 2π.)
const K_FROM_DEG_Q63: i128 = blk: {
    const k: comptime_float = @as(comptime_float, 32.0) / @as(comptime_float, 45.0);
    break :blk @intFromFloat(@floor(k * @as(comptime_float, 1 << 63) + 0.5));
};

/// `d.raw = round(a.raw × (45/32) × 2^62) / 2^62`.
const K_TO_DEG_Q62: i128 = blk: {
    const k: comptime_float = @as(comptime_float, 45.0) / @as(comptime_float, 32.0);
    break :blk @intFromFloat(@floor(k * @as(comptime_float, 1 << 62) + 0.5));
};

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

pub fn fromRaw(r: u32) Angle {
    return .{ .raw = r };
}

pub fn fromRadians(r: Fixed) Angle {
    const product: i128 = @as(i128, r.raw) * K_FROM_RAD_Q57;
    const scaled: i128 = rounding.shiftRound(product, 57);
    // Take low 32 bits as u32; the AND-mask gives 2's-complement wrap for
    // negative `scaled` (e.g. scaled = -1 → 0xFFFF_FFFF, the "just before zero"
    // angle), matching the angle's natural modular semantics.
    return .{ .raw = @intCast(scaled & ((@as(i128, 1) << 32) - 1)) };
}

pub fn toRadians(a: Angle) Fixed {
    const product: i128 = @as(i128, a.raw) * K_TO_RAD_Q68;
    const scaled: i128 = rounding.shiftRound(product, 68);
    return .{ .raw = rounding.narrow(scaled) };
}

pub fn fromDegrees(d: Fixed) Angle {
    const product: i128 = @as(i128, d.raw) * K_FROM_DEG_Q63;
    const scaled: i128 = rounding.shiftRound(product, 63);
    return .{ .raw = @intCast(scaled & ((@as(i128, 1) << 32) - 1)) };
}

pub fn toDegrees(a: Angle) Fixed {
    const product: i128 = @as(i128, a.raw) * K_TO_DEG_Q62;
    const scaled: i128 = rounding.shiftRound(product, 62);
    return .{ .raw = rounding.narrow(scaled) };
}

// ---------------------------------------------------------------------------
// Arithmetic — wrap IS the semantic. No overflow path exists.
// ---------------------------------------------------------------------------

pub fn addAngle(a: Angle, b: Angle) Angle {
    return .{ .raw = a.raw +% b.raw };
}

pub fn subAngle(a: Angle, b: Angle) Angle {
    return .{ .raw = a.raw -% b.raw };
}

pub fn eql(a: Angle, b: Angle) bool {
    return a.raw == b.raw;
}

// ===========================================================================
// Tests
// ===========================================================================

const testing = std.testing;

test "Angle has u32 layout (size + alignment)" {
    try testing.expectEqual(@sizeOf(u32), @sizeOf(Angle));
    try testing.expectEqual(@alignOf(u32), @alignOf(Angle));
}

test "constants: ZERO / QUARTER / HALF / THREE_QUARTER raw values" {
    try testing.expectEqual(@as(u32, 0), ZERO.raw);
    try testing.expectEqual(@as(u32, 1) << 30, QUARTER_TURN.raw);
    try testing.expectEqual(@as(u32, 1) << 31, HALF_TURN.raw);
    try testing.expectEqual(@as(u32, 3) << 30, THREE_QUARTER_TURN.raw);
    // The four quadrants partition the u32 range exactly.
    try testing.expectEqual(@as(u32, 0), QUARTER_TURN.raw +% QUARTER_TURN.raw +% QUARTER_TURN.raw +% QUARTER_TURN.raw);
}

test "addAngle wraps a full turn back to zero" {
    try testing.expectEqual(ZERO, QUARTER_TURN.addAngle(THREE_QUARTER_TURN));
    try testing.expectEqual(ZERO, HALF_TURN.addAngle(HALF_TURN));
    try testing.expectEqual(QUARTER_TURN, HALF_TURN.addAngle(THREE_QUARTER_TURN));
    // Max u32 + 1 = 0 (the wrap is the point of BAM).
    try testing.expectEqual(ZERO, (Angle{ .raw = std.math.maxInt(u32) }).addAngle(.{ .raw = 1 }));
}

test "subAngle wraps zero minus something to its complement" {
    try testing.expectEqual(THREE_QUARTER_TURN, ZERO.subAngle(QUARTER_TURN));
    try testing.expectEqual(HALF_TURN, ZERO.subAngle(HALF_TURN));
    try testing.expectEqual(Angle{ .raw = std.math.maxInt(u32) }, ZERO.subAngle(.{ .raw = 1 }));
}

test "addAngle and subAngle are inverses for any pair" {
    const samples = [_]Angle{
        ZERO,                       QUARTER_TURN,                       HALF_TURN, THREE_QUARTER_TURN,
        .{ .raw = 0xDEADBEEF }, .{ .raw = 0xCAFEBABE }, .{ .raw = 1 }, .{ .raw = std.math.maxInt(u32) },
    };
    for (samples) |a| {
        for (samples) |b| {
            try testing.expectEqual(a, a.addAngle(b).subAngle(b));
            try testing.expectEqual(a, a.subAngle(b).addAngle(b));
        }
    }
}

test "toRadians: cardinal angles match expected Fixed values within 1 ULP" {
    try testing.expectEqual(@as(i64, 0), ZERO.toRadians().raw);

    const half_rad = HALF_TURN.toRadians();
    try testing.expect(@abs(half_rad.raw - Fixed.PI.raw) <= 1);

    const q_rad = QUARTER_TURN.toRadians();
    try testing.expect(@abs(q_rad.raw - Fixed.HALF_PI.raw) <= 1);

    const tq_rad = THREE_QUARTER_TURN.toRadians();
    const three_half_pi = Fixed.HALF_PI.raw + Fixed.PI.raw;
    try testing.expect(@abs(tq_rad.raw - three_half_pi) <= 1);
}

test "fromRadians: cardinal radians map to expected angles within input-bound tolerance" {
    // Tolerance budget: Fixed.PI is itself rounded to nearest representable
    // Fixed value, so it has up to 0.5 Fixed-ULP error vs real π. That error
    // scales by 128/π ≈ 40.7 when projected into u32 angle ULPs — so even a
    // perfectly-implemented fromRadians can be off from HALF_TURN by ~20 ULPs
    // just because Fixed.PI ≠ π. Tolerance of 32 covers the input-error
    // budget plus a few ULPs of conversion rounding.
    const cardinal_tolerance: u32 = 32;

    try testing.expectEqual(@as(u32, 0), fromRadians(Fixed.ZERO).raw);
    try testing.expect(diffMod32(fromRadians(Fixed.PI).raw, HALF_TURN.raw) <= cardinal_tolerance);
    try testing.expect(diffMod32(fromRadians(Fixed.HALF_PI).raw, QUARTER_TURN.raw) <= cardinal_tolerance);
    // 2π → ZERO (wraps via mod-2^32). Fixed.TWO_PI happens to round very
    // close to real 2π, so the diff here is tiny.
    try testing.expect(diffMod32(fromRadians(Fixed.TWO_PI).raw, ZERO.raw) <= cardinal_tolerance);
}

test "fromRadians / toRadians round-trip within 4 ULP (positive inputs in [0, 2π))" {
    // Restricted to non-negative radians: the BAM Angle has no sign, so a
    // negative input round-trips to its positive equivalent (separate test
    // below covers the modular-wrap behavior explicitly).
    inline for (.{ 0.0, 0.5, 1.0, 1.5, 3.0, 3.14, 4.5, 6.0 }) |f| {
        const original = Fixed.rconst(f);
        const back = fromRadians(original).toRadians();
        try testing.expect(@abs(original.raw - back.raw) <= 4);
    }
}

test "fromRadians wraps mod 2π: x and x + Fixed.TWO_PI give the same angle" {
    // Fixed.TWO_PI happens to round within ~0.03 Fixed-ULPs of real 2π (the
    // RNE error there is tiny — see the bit-stability tests). After scaling
    // by 128/π that adds maybe 1-2 u32 ULPs vs ang(x). Tolerance: 8 is
    // generous.
    const ang_one = fromRadians(Fixed.ONE);
    const ang_one_plus_2pi = fromRadians(Fixed.ONE.add(Fixed.TWO_PI));
    try testing.expect(diffMod32(ang_one.raw, ang_one_plus_2pi.raw) <= 8);

    const ang_two = fromRadians(Fixed.rconst(2.0));
    const ang_two_plus_2pi = fromRadians(Fixed.rconst(2.0).add(Fixed.TWO_PI));
    try testing.expect(diffMod32(ang_two.raw, ang_two_plus_2pi.raw) <= 8);
}

test "fromRadians: negative-input wrap is symmetric around HALF_TURN" {
    // ±Fixed.PI inputs land symmetrically around HALF_TURN because
    // -Fixed.PI.raw represents the real value -π-ε (where ε is Fixed.PI's
    // ~0.5-ULP error), which mod 2π is π-ε — distinct from +π+ε by 2ε on
    // the circle. So this is not a wrap-equivalence (the two inputs encode
    // genuinely different real angles) but a sanity check on the mod-2^32
    // wrap math: pos_raw + neg_raw should equal 2^32 modulo 2π.
    const pos = fromRadians(Fixed.PI).raw;
    const neg = fromRadians(.{ .raw = -Fixed.PI.raw }).raw;
    // pos + neg should wrap to 0 (since -π + π = 0).
    try testing.expect(@as(u32, pos +% neg) <= 4);
}

test "toDegrees: cardinal angles match 0/90/180/270 within 1 ULP" {
    try testing.expectEqual(@as(i64, 0), ZERO.toDegrees().raw);

    try testing.expect(@abs(HALF_TURN.toDegrees().raw - Fixed.fromInt(180).raw) <= 1);
    try testing.expect(@abs(QUARTER_TURN.toDegrees().raw - Fixed.fromInt(90).raw) <= 1);
    try testing.expect(@abs(THREE_QUARTER_TURN.toDegrees().raw - Fixed.fromInt(270).raw) <= 1);
}

test "fromDegrees: cardinal degrees map to expected angles within a few ULP" {
    try testing.expectEqual(@as(u32, 0), fromDegrees(Fixed.ZERO).raw);
    try testing.expect(diffMod32(fromDegrees(Fixed.fromInt(180)).raw, HALF_TURN.raw) <= 2);
    try testing.expect(diffMod32(fromDegrees(Fixed.fromInt(90)).raw, QUARTER_TURN.raw) <= 2);
    try testing.expect(diffMod32(fromDegrees(Fixed.fromInt(360)).raw, ZERO.raw) <= 2);
}

test "fromDegrees / toDegrees round-trip within 2 ULP" {
    inline for (.{ 0, 30, 45, 90, 180, 270, 359 }) |d| {
        const original = Fixed.fromInt(d);
        const back = fromDegrees(original).toDegrees();
        try testing.expect(@abs(original.raw - back.raw) <= 2);
    }
}

test "deg / rad consistency: 180° via fromDegrees ≈ π via fromRadians" {
    // fromDegrees(180) lands within ~1 ULP of HALF_TURN (180 and 32/45 are
    // both exact at our precision). fromRadians(Fixed.PI) lands within ~12
    // ULPs of HALF_TURN due to Fixed.PI's own rounding error vs real π.
    const from_180 = fromDegrees(Fixed.fromInt(180));
    const from_pi = fromRadians(Fixed.PI);
    try testing.expect(diffMod32(from_180.raw, from_pi.raw) <= 32);
}

test "addAngle is the modular sum (sanity at u32 corners)" {
    const a: Angle = .{ .raw = std.math.maxInt(u32) - 10 };
    const b: Angle = .{ .raw = 100 };
    try testing.expectEqual(@as(u32, 89), a.addAngle(b).raw);
}

// Helper for tests: minimum of forward / backward wrap distance, mod 2^32.
fn diffMod32(a: u32, b: u32) u32 {
    const forward = a -% b;
    const backward = b -% a;
    return @min(forward, backward);
}
