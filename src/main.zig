//! Demo runner. The library lives in src/root.zig; this exists so
//! `zig build run` is wired up. Not on any sim path.

const std = @import("std");
const fpz = @import("fpz");

pub fn main() !void {
    _ = fpz;
    std.debug.print("fpz: deterministic fixed-point math library\n", .{});
}
