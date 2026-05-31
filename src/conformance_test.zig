//! Conformance tests — verify the live implementation output matches the
//! committed goldens in `conformance.zig`. Any diff fails CI (SPEC §8).
//!
//! Tests run under the test runner's chosen optimize mode; the
//! `test-all-modes` build step runs them under Debug, ReleaseSafe, AND
//! ReleaseFast, locking in the SPEC §2 corollary that output is bit-
//! identical across modes.

const std = @import("std");
const fpz = @import("root.zig");
const conformance = @import("conformance.zig");

test "conformance: sin goldens" {
    for (conformance.sin_goldens) |g| {
        const got = fpz.sin(fpz.Angle{ .raw = g.raw });
        try std.testing.expectEqual(g.expected, got.raw);
    }
}

test "conformance: cos goldens" {
    for (conformance.cos_goldens) |g| {
        const got = fpz.cos(fpz.Angle{ .raw = g.raw });
        try std.testing.expectEqual(g.expected, got.raw);
    }
}

test "conformance: tan goldens" {
    for (conformance.tan_goldens) |g| {
        const got = fpz.tan(fpz.Angle{ .raw = g.raw });
        try std.testing.expectEqual(g.expected, got.raw);
    }
}

test "conformance: sqrt goldens" {
    for (conformance.sqrt_goldens) |g| {
        const got = fpz.sqrt(fpz.Fixed{ .raw = g.input_raw });
        try std.testing.expectEqual(g.expected_raw, got.raw);
    }
}

test "conformance: exp goldens" {
    for (conformance.exp_goldens) |g| {
        const got = fpz.exp(fpz.Fixed{ .raw = g.input_raw });
        try std.testing.expectEqual(g.expected_raw, got.raw);
    }
}

test "conformance: ln goldens" {
    for (conformance.ln_goldens) |g| {
        const got = fpz.ln(fpz.Fixed{ .raw = g.input_raw });
        try std.testing.expectEqual(g.expected_raw, got.raw);
    }
}

test "conformance: atan2 goldens" {
    for (conformance.atan2_goldens) |g| {
        const got = fpz.atan2(
            fpz.Fixed{ .raw = g.y_raw },
            fpz.Fixed{ .raw = g.x_raw },
        );
        try std.testing.expectEqual(g.expected_angle_raw, got.raw);
    }
}
