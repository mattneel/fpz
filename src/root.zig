//! fpz — a deterministic fixed-point math library for Zig.
//!
//! See SPEC.md for the design contract. The single inviolable rule:
//! identical inputs produce bit-identical outputs on every target arch,
//! build mode, and toolchain. Precision, speed, and ergonomics are all
//! subordinate to that.

const std = @import("std");

pub const rounding = @import("rounding.zig");
pub const Fixed = @import("Fixed.zig");
pub const Angle = @import("Angle.zig");

const sqrt_mod = @import("sqrt.zig");
pub const sqrt = sqrt_mod.sqrt;
pub const sqrtChecked = sqrt_mod.sqrtChecked;

const trig_mod = @import("trig.zig");
pub const sin = trig_mod.sin;
pub const cos = trig_mod.cos;
pub const tan = trig_mod.tan;
pub const atan2 = trig_mod.atan2;

const transcendental_mod = @import("transcendental.zig");
pub const exp = transcendental_mod.exp;
pub const ln = transcendental_mod.ln;
pub const lnChecked = transcendental_mod.lnChecked;
pub const log = transcendental_mod.log;
pub const pow = transcendental_mod.pow;

pub const simd = @import("simd.zig");

test {
    std.testing.refAllDecls(@This());
    // refAllDecls follows pub consts, so sqrt_mod gets pulled in via the
    // re-exports above. Keep an explicit import here as a safety net in case
    // the public API ever stops referencing a module file directly.
    _ = sqrt_mod;
    _ = trig_mod;
    _ = transcendental_mod;
}
