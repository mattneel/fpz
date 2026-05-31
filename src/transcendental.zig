//! Transcendentals on Fixed: exp, ln (Phase 6); log, pow (Phase 7).
//!
//! exp(x): range-reduce x = k·LN2 + r with |r| ≤ LN2/2, polynomial for
//!         exp(r), scale by 2^k via shift. Saturates on overflow.
//! ln(x):  extract integer power-of-two k from raw, polynomial for ln(m)
//!         on m ∈ [1, 2) via the atanh substitution, recombine k·LN2 + ln(m).
//!         Precondition x > 0; lnChecked returns an error union; the
//!         unchecked variant returns Fixed.MIN as a defined fallback.

const std = @import("std");
const Fixed = @import("Fixed.zig");

// ===========================================================================
// exp
// ===========================================================================

// exp(r) on |r| ≤ LN2/2 ≈ 0.347. Truncation through r^7/5040:
//   |r^8/40320| ≤ 4.4e-9 ≈ 2^-27.7 — far inside the 2^-22 target.
const EXP_C7: Fixed = Fixed.rconst(1.0 / 5040.0);
const EXP_C6: Fixed = Fixed.rconst(1.0 / 720.0);
const EXP_C5: Fixed = Fixed.rconst(1.0 / 120.0);
const EXP_C4: Fixed = Fixed.rconst(1.0 / 24.0);
const EXP_C3: Fixed = Fixed.rconst(1.0 / 6.0);
const EXP_C2: Fixed = Fixed.rconst(1.0 / 2.0);
const EXP_C1: Fixed = Fixed.ONE;
const EXP_C0: Fixed = Fixed.ONE;

fn polyExp(r: Fixed) Fixed {
    // 1 + r·(1 + r·(1/2 + r·(1/6 + r·(1/24 + r·(1/120 + r·(1/720 + r·1/5040))))))
    var p = EXP_C7;
    p = EXP_C6.add(r.mul(p));
    p = EXP_C5.add(r.mul(p));
    p = EXP_C4.add(r.mul(p));
    p = EXP_C3.add(r.mul(p));
    p = EXP_C2.add(r.mul(p));
    p = EXP_C1.add(r.mul(p));
    p = EXP_C0.add(r.mul(p));
    return p;
}

// Outer bounds chosen so x/LN2 stays well inside Fixed range. The true
// representability bound is x ≈ 27.03 (exp ≈ 5.5e11 ≈ Fixed.MAX); we cap
// slightly above so the saturate path is reachable.
const EXP_BOUND_HI: Fixed = Fixed.rconst(28.0);
const EXP_BOUND_LO: Fixed = Fixed.rconst(-28.0);

pub fn exp(x: Fixed) Fixed {
    if (x.raw > EXP_BOUND_HI.raw) return Fixed.MAX;
    if (x.raw < EXP_BOUND_LO.raw) return Fixed.ZERO;

    // Decompose x = k·LN2 + r, |r| ≤ LN2/2.
    const k_int: i64 = x.div(Fixed.LN2).roundToInt();
    const k_fixed = Fixed.fromInt(k_int);
    const r = x.sub(k_fixed.mul(Fixed.LN2));

    const expr = polyExp(r);

    // Scale by 2^k. Use an i128 intermediate so the saturate path is total.
    var result: i128 = expr.raw;
    if (k_int >= 0) {
        if (k_int > 63) return Fixed.MAX;
        result <<= @as(u7, @intCast(k_int));
    } else {
        const neg_k = -k_int;
        if (neg_k > 63) return Fixed.ZERO;
        result >>= @as(u7, @intCast(neg_k));
    }

    if (result > std.math.maxInt(i64)) return Fixed.MAX;
    if (result < std.math.minInt(i64)) return Fixed.MIN;
    return .{ .raw = @truncate(result) };
}

// ===========================================================================
// ln
// ===========================================================================

// ln(m) for m ∈ [1, 2):
//   ln(m) = 2·atanh(t) where t = (m-1)/(m+1) ∈ [0, 1/3]
//   = 2·(t + t³/3 + t⁵/5 + … + t¹³/13 + …)
// Truncation through t^13: 2·|t^15/15| at t=1/3 ≈ 9.3e-9 ≈ 2^-27 — safe.
const LN_C6: Fixed = Fixed.rconst(1.0 / 13.0);
const LN_C5: Fixed = Fixed.rconst(1.0 / 11.0);
const LN_C4: Fixed = Fixed.rconst(1.0 / 9.0);
const LN_C3: Fixed = Fixed.rconst(1.0 / 7.0);
const LN_C2: Fixed = Fixed.rconst(1.0 / 5.0);
const LN_C1: Fixed = Fixed.rconst(1.0 / 3.0);

fn polyLnAtanh(t: Fixed) Fixed {
    // 2·t·(1 + y/3 + y²/5 + y³/7 + y⁴/9 + y⁵/11 + y⁶/13)   where y = t²
    const y = t.mul(t);
    var p = LN_C6;
    p = LN_C5.add(y.mul(p));
    p = LN_C4.add(y.mul(p));
    p = LN_C3.add(y.mul(p));
    p = LN_C2.add(y.mul(p));
    p = LN_C1.add(y.mul(p));
    p = Fixed.ONE.add(y.mul(p));
    const tp = t.mul(p);
    return tp.add(tp); // ×2, exact (no rounding)
}

fn lnPositive(x: Fixed) Fixed {
    // Decompose x = 2^k · m with m ∈ [1, 2).
    //   raw_log2 = floor(log2(x.raw))
    //   k = raw_log2 - frac_bits
    //   m.raw = x.raw shifted so its leading 1-bit sits at frac_bits position
    const raw_unsigned: u64 = @intCast(x.raw);
    const raw_log2: i64 = 63 - @as(i64, @intCast(@clz(raw_unsigned)));
    const k_int: i64 = raw_log2 - Fixed.frac_bits;

    var m_raw: i64 = x.raw;
    if (k_int > 0) {
        m_raw = x.raw >> @as(u6, @intCast(k_int));
    } else if (k_int < 0) {
        m_raw = x.raw << @as(u6, @intCast(-k_int));
    }
    const m = Fixed{ .raw = m_raw };

    // t = (m - 1) / (m + 1) ∈ [0, 1/3]
    const t = m.sub(Fixed.ONE).div(m.add(Fixed.ONE));
    const ln_m = polyLnAtanh(t);

    // ln(x) = k·LN2 + ln(m)
    const k_ln2 = Fixed.LN2.mul(Fixed.fromInt(k_int));
    return k_ln2.add(ln_m);
}

pub fn ln(x: Fixed) Fixed {
    // SPEC §7: x ≤ 0 is a domain violation; defined ReleaseFast fallback is
    // Fixed.MIN. No std.debug.assert here — it would license the optimizer
    // to assume x > 0 and delete the guard in ReleaseFast.
    if (x.raw <= 0) return Fixed.MIN;
    return lnPositive(x);
}

pub fn lnChecked(x: Fixed) error{NonPositive}!Fixed {
    if (x.raw <= 0) return error.NonPositive;
    return lnPositive(x);
}

// ===========================================================================
// log, pow
// ===========================================================================

/// log_b(x) = ln(x) / ln(b). Returns Fixed.MIN as a defined fallback for
/// invalid domains (x ≤ 0, base ≤ 0, base == 1).
pub fn log(x: Fixed, base: Fixed) Fixed {
    if (x.raw <= 0) return Fixed.MIN;
    if (base.raw <= 0 or base.eql(Fixed.ONE)) return Fixed.MIN;
    const ln_x = lnPositive(x);
    const ln_base = lnPositive(base);
    return ln_x.div(ln_base);
}

/// Integer-exponent fast path (exact, by squaring).
fn powInt(base: Fixed, n: i64) Fixed {
    if (n == 0) return Fixed.ONE;
    if (n < 0) {
        if (base.raw == 0) return Fixed.MAX; // 0^(neg) → ∞, saturate
        const inv = Fixed.ONE.div(base);
        return powInt(inv, -n);
    }
    var result = Fixed.ONE;
    var current = base;
    var remaining: u64 = @intCast(n);
    while (remaining > 0) {
        if (remaining & 1 == 1) {
            result = result.mul(current);
        }
        remaining >>= 1;
        if (remaining > 0) {
            current = current.mul(current);
        }
    }
    return result;
}

/// pow(base, exponent). Exact integer-exponent path; otherwise exp(exp·ln(base))
/// with base > 0 required. Domain fallback: returns ZERO for non-integer
/// exponent with non-positive base (mathematically complex / undefined).
pub fn pow(base: Fixed, exponent: Fixed) Fixed {
    if (exponent.raw == 0) return Fixed.ONE;

    // Whole-number exponent detection — toInt truncates toward zero, so this
    // matches iff the exponent has zero fractional part.
    const exp_int_part = exponent.toInt();
    const exp_int_as_fixed = Fixed.fromInt(exp_int_part);
    if (exp_int_as_fixed.eql(exponent)) {
        return powInt(base, exp_int_part);
    }

    // Non-integer exponent: requires positive base (SPEC §6).
    if (base.raw <= 0) return Fixed.ZERO;
    return exp(exponent.mul(lnPositive(base)));
}

// ===========================================================================
// Tests
// ===========================================================================

const testing = std.testing;

test "exp: cardinal values" {
    try testing.expectEqual(Fixed.ONE, exp(Fixed.ZERO));

    // exp(1) ≈ e
    const e_got = exp(Fixed.ONE);
    try testing.expect(@abs(e_got.raw - Fixed.E.raw) <= 8);

    // exp(LN2) ≈ 2
    const two_got = exp(Fixed.LN2);
    try testing.expect(@abs(two_got.raw - Fixed.fromInt(2).raw) <= 8);

    // exp(-1) ≈ 1/e ≈ 0.36788
    const inv_e = exp(Fixed.ONE.neg());
    const expected_inv = Fixed.rconst(1.0 / std.math.e);
    try testing.expect(@abs(inv_e.raw - expected_inv.raw) <= 8);
}

test "exp: saturates on overflow / underflow inputs" {
    // x well above the representable range → MAX
    try testing.expectEqual(Fixed.MAX, exp(Fixed.rconst(100.0)));
    try testing.expectEqual(Fixed.MAX, exp(Fixed.fromInt(30)));
    // x well below → ZERO
    try testing.expectEqual(Fixed.ZERO, exp(Fixed.rconst(-100.0)));
    try testing.expectEqual(Fixed.ZERO, exp(Fixed.fromInt(-30)));
}

test "exp: differential vs @exp(f64) within Fixed-precision tolerance" {
    // Combined absolute + relative tolerance — Fixed.raw has only ~6e-8
    // absolute resolution, so small exp values (e.g. exp(-10) ≈ 4.5e-5)
    // can't beat ~1e-3 relative even with a perfect algorithm. The
    // absolute_floor term encodes that representation limit.
    const absolute_floor: f64 = 2.0 / @as(f64, @floatFromInt(@as(i64, 1) << Fixed.frac_bits));
    inline for (.{ 0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, -1.0, -5.0, -10.0 }) |f| {
        const x = Fixed.rconst(f);
        const got = exp(x).toFloat();
        const ref = @exp(f);
        const diff = @abs(got - ref);
        const tolerance = absolute_floor + 1e-5 * @abs(ref);
        try testing.expect(diff < tolerance);
    }
}

test "ln: cardinal values" {
    try testing.expectEqual(Fixed.ZERO, ln(Fixed.ONE));

    // ln(2) ≈ LN2
    const ln_2 = ln(Fixed.fromInt(2));
    try testing.expect(@abs(ln_2.raw - Fixed.LN2.raw) <= 8);

    // ln(e) ≈ 1
    const ln_e = ln(Fixed.E);
    try testing.expect(@abs(ln_e.raw - Fixed.ONE.raw) <= 8);

    // ln(0.5) ≈ -LN2
    const ln_half = ln(Fixed.HALF);
    try testing.expect(@abs(ln_half.raw + Fixed.LN2.raw) <= 8);

    // ln(4) ≈ 2·LN2
    const ln_4 = ln(Fixed.fromInt(4));
    try testing.expect(@abs(ln_4.raw - Fixed.LN2.raw * 2) <= 8);
}

test "ln: defined fallback for non-positive inputs" {
    try testing.expectEqual(Fixed.MIN, ln(Fixed.ZERO));
    try testing.expectEqual(Fixed.MIN, ln(Fixed.NEG_ONE));
    try testing.expectEqual(Fixed.MIN, ln(Fixed.fromInt(-100)));
}

test "lnChecked: error union for non-positive, value otherwise" {
    try testing.expectError(error.NonPositive, lnChecked(Fixed.ZERO));
    try testing.expectError(error.NonPositive, lnChecked(Fixed.NEG_ONE));
    try testing.expectEqual(Fixed.ZERO, try lnChecked(Fixed.ONE));
    // ln(2) within tolerance
    const ln_2 = try lnChecked(Fixed.fromInt(2));
    try testing.expect(@abs(ln_2.raw - Fixed.LN2.raw) <= 8);
}

test "ln: differential vs @log(f64) within 1e-5 relative error" {
    inline for (.{ 0.001, 0.1, 0.5, 1.0, 2.0, 10.0, 1000.0, 1e9 }) |f| {
        const x = Fixed.rconst(f);
        const got = ln(x).toFloat();
        const ref = @log(f);
        const rel = if (@abs(ref) > 1e-10) @abs(got - ref) / @abs(ref) else @abs(got - ref);
        try testing.expect(rel < 1e-5);
    }
}

test "exp(ln(x)) ≈ x within 1e-4 relative error" {
    // Round-trip accumulates two polynomial errors + the k·LN2 mismatch.
    inline for (.{ 0.5, 1.0, 2.0, 7.5, 100.0, 1e6 }) |f| {
        const x = Fixed.rconst(f);
        const round_tripped = exp(ln(x)).toFloat();
        const rel = @abs(round_tripped - f) / f;
        try testing.expect(rel < 1e-4);
    }
}

test "ln(exp(x)) ≈ x within tolerance" {
    inline for (.{ -5.0, -1.0, 0.0, 0.5, 1.0, 2.0, 5.0, 10.0 }) |f| {
        const x = Fixed.rconst(f);
        const round_tripped = ln(exp(x)).toFloat();
        try testing.expect(@abs(round_tripped - f) < 1e-4);
    }
}

test "log: identity log_b(b) = 1, log_b(1) = 0" {
    inline for (.{ 2.0, 3.0, 10.0, std.math.e }) |b| {
        const base = Fixed.rconst(b);
        const one_result = log(base, base);
        try testing.expect(@abs(one_result.raw - Fixed.ONE.raw) <= 16);
        const zero_result = log(Fixed.ONE, base);
        try testing.expect(@abs(zero_result.raw) <= 4);
    }
}

test "log: known integer values" {
    // log_2(8) = 3
    const log_8_2 = log(Fixed.fromInt(8), Fixed.fromInt(2));
    try testing.expect(@abs(log_8_2.raw - Fixed.fromInt(3).raw) <= 32);
    // log_10(100) = 2
    const log_100_10 = log(Fixed.fromInt(100), Fixed.fromInt(10));
    try testing.expect(@abs(log_100_10.raw - Fixed.fromInt(2).raw) <= 32);
    // log_2(1024) = 10
    const log_1024_2 = log(Fixed.fromInt(1024), Fixed.fromInt(2));
    try testing.expect(@abs(log_1024_2.raw - Fixed.fromInt(10).raw) <= 64);
}

test "log: defined fallback for invalid domains" {
    try testing.expectEqual(Fixed.MIN, log(Fixed.ZERO, Fixed.fromInt(2)));
    try testing.expectEqual(Fixed.MIN, log(Fixed.NEG_ONE, Fixed.fromInt(2)));
    try testing.expectEqual(Fixed.MIN, log(Fixed.fromInt(10), Fixed.ZERO));
    try testing.expectEqual(Fixed.MIN, log(Fixed.fromInt(10), Fixed.NEG_ONE));
    try testing.expectEqual(Fixed.MIN, log(Fixed.fromInt(10), Fixed.ONE)); // log base 1
}

test "pow: integer exponent fast path is exact" {
    try testing.expectEqual(Fixed.ONE, pow(Fixed.fromInt(7), Fixed.ZERO));
    try testing.expectEqual(Fixed.fromInt(7), pow(Fixed.fromInt(7), Fixed.ONE));
    try testing.expectEqual(Fixed.fromInt(49), pow(Fixed.fromInt(7), Fixed.fromInt(2)));
    try testing.expectEqual(Fixed.fromInt(1024), pow(Fixed.fromInt(2), Fixed.fromInt(10)));
    try testing.expectEqual(Fixed.fromInt(125), pow(Fixed.fromInt(5), Fixed.fromInt(3)));
    // Negative base, integer exp — sign by parity.
    try testing.expectEqual(Fixed.fromInt(-8), pow(Fixed.fromInt(-2), Fixed.fromInt(3)));
    try testing.expectEqual(Fixed.fromInt(16), pow(Fixed.fromInt(-2), Fixed.fromInt(4)));
}

test "pow: pow(x, 0) == ONE for any x including 0" {
    try testing.expectEqual(Fixed.ONE, pow(Fixed.ZERO, Fixed.ZERO));
    try testing.expectEqual(Fixed.ONE, pow(Fixed.ONE, Fixed.ZERO));
    try testing.expectEqual(Fixed.ONE, pow(Fixed.fromInt(-5), Fixed.ZERO));
    try testing.expectEqual(Fixed.ONE, pow(Fixed.MAX, Fixed.ZERO));
}

test "pow: negative integer exponent gives the reciprocal" {
    // pow(2, -1) = 0.5
    try testing.expectEqual(Fixed.HALF, pow(Fixed.fromInt(2), Fixed.fromInt(-1)));
    // pow(4, -2) = 1/16 = 0.0625
    const expected = Fixed.rconst(0.0625);
    try testing.expectEqual(expected, pow(Fixed.fromInt(4), Fixed.fromInt(-2)));
    // pow(0, -n) saturates
    try testing.expectEqual(Fixed.MAX, pow(Fixed.ZERO, Fixed.fromInt(-1)));
}

test "pow: fractional exponent via exp(exp · ln(base))" {
    // sqrt(2) = pow(2, 0.5)
    const sqrt2 = pow(Fixed.fromInt(2), Fixed.HALF);
    try testing.expect(@abs(sqrt2.toFloat() - @sqrt(2.0)) < 1e-4);
    // cbrt(8) = pow(8, 1/3) = 2
    const cbrt8 = pow(Fixed.fromInt(8), Fixed.rconst(1.0 / 3.0));
    try testing.expect(@abs(cbrt8.toFloat() - 2.0) < 1e-3);
    // 2^1.5 = 2 · sqrt(2) ≈ 2.828
    const two_pow_1_5 = pow(Fixed.fromInt(2), Fixed.rconst(1.5));
    try testing.expect(@abs(two_pow_1_5.toFloat() - 2.828427) < 1e-3);
}

test "pow: defined fallback for negative base with fractional exponent" {
    try testing.expectEqual(Fixed.ZERO, pow(Fixed.fromInt(-2), Fixed.HALF));
    try testing.expectEqual(Fixed.ZERO, pow(Fixed.fromInt(-1), Fixed.rconst(1.5)));
}
