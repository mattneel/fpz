//! Conformance-vector generator.
//!
//! Runs the current fpz implementation at a fixed set of inputs and writes
//! the resulting raw output bits to `src/conformance.zig`. The committed
//! conformance.zig becomes the determinism contract (SPEC §8): any diff
//! between the file and the live implementation output fails CI.
//!
//! To re-baseline (e.g., after a polynomial change):
//!     zig build gen-conformance
//!     # review the diff carefully — it IS the contract change
//!     git add src/conformance.zig
//!     git commit -m "fpz: rebaseline conformance vectors (reason: ...)"

const std = @import("std");
const Io = std.Io;
const fpz = @import("fpz");

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    var file = try Io.Dir.cwd().createFile(io, "src/conformance.zig", .{});
    defer file.close(io);

    var buf: [16384]u8 = undefined;
    var fw: Io.File.Writer = .init(file, io, &buf);
    const w = &fw.interface;

    try w.writeAll(
        \\//! AUTO-GENERATED conformance golden vectors (SPEC §8).
        \\//!
        \\//! Run `zig build gen-conformance` to regenerate. Any diff between
        \\//! these expected outputs and the current implementation fails the
        \\//! conformance test — that's the determinism guarantee, made
        \\//! empirical and CI-enforced across every target arch, toolchain,
        \\//! and build mode (Debug / ReleaseSafe / ReleaseFast).
        \\
        \\const std = @import("std");
        \\const fpz = @import("root.zig");
        \\
        \\const AngleGolden = struct { raw: u32, expected: i64 };
        \\const FixedGolden = struct { input_raw: i64, expected_raw: i64 };
        \\const PairGolden = struct { y_raw: i64, x_raw: i64, expected_angle_raw: u32 };
        \\
        \\
    );

    // ---- sin / cos / tan goldens ----
    try emitAngleHeader(w, "sin_goldens", "sin");
    try emitAngleSamples(w, fpz.sin);
    try w.writeAll("};\n\n");

    try emitAngleHeader(w, "cos_goldens", "cos");
    try emitAngleSamples(w, fpz.cos);
    try w.writeAll("};\n\n");

    // Skip tan near poles (saturation regions) — only sample stable angles.
    try w.writeAll(
        \\pub const tan_goldens = [_]AngleGolden{
        \\
    );
    inline for ([_]u32{
        0,         0x05000000, 0x0A000000, 0x10000000,
        0x15000000, 0x18000000, 0x35000000, 0x55000000,
        0x70000000, 0x90000000, 0xA0000000, 0xB0000000,
        0xC8000000, 0xD0000000, 0xE0000000, 0xF8000000,
    }) |raw| {
        const got = fpz.tan(fpz.Angle{ .raw = raw }).raw;
        try w.print("    .{{ .raw = 0x{X:0>8}, .expected = {d} }},\n", .{ raw, got });
    }
    try w.writeAll("};\n\n");

    // ---- sqrt goldens ----
    try w.writeAll(
        \\pub const sqrt_goldens = [_]FixedGolden{
        \\
    );
    inline for ([_]i64{
        0, 1, 16, 256, 16777216, 33554432,
        62914560, 100000000, 1677721600, 16777216000,
        167772160000, 9223372036854775000,
    }) |raw| {
        const got = fpz.sqrt(.{ .raw = raw }).raw;
        try w.print("    .{{ .input_raw = {d}, .expected_raw = {d} }},\n", .{ raw, got });
    }
    try w.writeAll("};\n\n");

    // ---- exp goldens ----
    try w.writeAll(
        \\pub const exp_goldens = [_]FixedGolden{
        \\
    );
    inline for ([_]i64{
        // -10, -5, -1, -0.5, 0, 0.5, 1, 2, 5, 10, 20 — as Fixed.raw
        -167772160, -83886080, -16777216, -8388608, 0,
        8388608,    16777216,  33554432,  83886080, 167772160, 335544320,
    }) |raw| {
        const got = fpz.exp(.{ .raw = raw }).raw;
        try w.print("    .{{ .input_raw = {d}, .expected_raw = {d} }},\n", .{ raw, got });
    }
    try w.writeAll("};\n\n");

    // ---- ln goldens ----
    try w.writeAll(
        \\pub const ln_goldens = [_]FixedGolden{
        \\
    );
    inline for ([_]i64{
        // 0.001, 0.5, 1, 2, e, 10, 100, 1000, large
        16777,    8388608,  16777216, 33554432, 45605201,
        167772160, 1677721600, 16777216000, 1000000000000,
    }) |raw| {
        const got = fpz.ln(.{ .raw = raw }).raw;
        try w.print("    .{{ .input_raw = {d}, .expected_raw = {d} }},\n", .{ raw, got });
    }
    try w.writeAll("};\n\n");

    // ---- atan2 goldens ----
    try w.writeAll(
        \\pub const atan2_goldens = [_]PairGolden{
        \\
    );
    inline for ([_]struct { y: i64, x: i64 }{
        .{ .y = 0, .x = 16777216 },             // (0, 1)
        .{ .y = 16777216, .x = 0 },             // (1, 0)
        .{ .y = 0, .x = -16777216 },            // (0, -1)
        .{ .y = -16777216, .x = 0 },            // (-1, 0)
        .{ .y = 16777216, .x = 16777216 },      // (1, 1) → π/4
        .{ .y = 16777216, .x = -16777216 },     // (1, -1) → 3π/4
        .{ .y = -16777216, .x = -16777216 },    // (-1, -1) → 5π/4
        .{ .y = -16777216, .x = 16777216 },     // (-1, 1) → 7π/4
        .{ .y = 33554432, .x = 16777216 },      // (2, 1)
        .{ .y = 16777216, .x = 33554432 },      // (1, 2)
        .{ .y = -16777216, .x = 50331648 },     // (-1, 3)
        .{ .y = 8388608, .x = 100000000 },      // (0.5, ~6)
    }) |pair| {
        const got = fpz.atan2(.{ .raw = pair.y }, .{ .raw = pair.x }).raw;
        try w.print("    .{{ .y_raw = {d}, .x_raw = {d}, .expected_angle_raw = 0x{X:0>8} }},\n", .{ pair.y, pair.x, got });
    }
    try w.writeAll("};\n");

    try w.flush();
}

fn emitAngleHeader(w: anytype, name: []const u8, op: []const u8) !void {
    try w.print("/// Goldens for {s}. 64 angles spread across the u32 range.\n", .{op});
    try w.print("pub const {s} = [_]AngleGolden{{\n", .{name});
}

fn emitAngleSamples(w: anytype, op: fn (fpz.Angle) fpz.Fixed) !void {
    var i: u32 = 0;
    while (i < 64) : (i += 1) {
        const raw: u32 = @intCast((@as(u64, i) * std.math.maxInt(u32)) / 63);
        const got = op(.{ .raw = raw }).raw;
        try w.print("    .{{ .raw = 0x{X:0>8}, .expected = {d} }},\n", .{ raw, got });
    }
}
