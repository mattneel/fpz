# fpz

> A deterministic fixed-point math library for Zig.

Identical inputs produce **bit-identical outputs** on every target architecture, build mode, and toolchain. Precision, speed, and ergonomics are all subordinate to that. See [SPEC.md](SPEC.md) for the full design contract.

## Why

Game simulations, lockstep multiplayer, replay systems, and reproducible scientific computing all need the same guarantee: *if A == B on machine X at frame N, then A == B on machine Y at frame N*. IEEE-754 floats can't give that — rounding modes, FMA fusion, libm vendor variance, and the compiler's freedom to reassociate all conspire against you.

`fpz` is the numeric substrate for that guarantee:

- Single Q40.24 fixed-point type. ±5.5×10¹¹ range, ~6×10⁻⁸ precision.
- Integer-only runtime. No libm, no runtime float, ever, on a sim path.
- Overflow is defined as wrap — asserted in safe builds, deterministic in `ReleaseFast`.
- Binary Angle Measure (BAM) for trig — modular by construction, no range reduction.
- Golden-vector conformance suite locks the contract empirically.

## Install

Add `fpz` as a dependency in your `build.zig.zon`:

```zig
.dependencies = .{
    .fpz = .{
        .url = "https://example.com/fpz/v0.1.0.tar.gz",
        .hash = "...",
    },
},
```

Then wire it in `build.zig`:

```zig
const fpz = b.dependency("fpz", .{ .target = target, .optimize = optimize });
exe.root_module.addImport("fpz", fpz.module("fpz"));
```

## Quick start

```zig
const fpz = @import("fpz");
const Fixed = fpz.Fixed;
const Angle = fpz.Angle;

// Construction
const x = Fixed.fromInt(42);
const y = Fixed.rconst(3.14);          // comptime-only — baked at compile time
const z = x.add(y).mul(Fixed.HALF);    // method-call chaining

// Trig via BAM angles (modular by construction)
const a = Angle.fromDegrees(Fixed.fromInt(45));
const s = fpz.sin(a);                  // ≈ √2/2
const c = fpz.cos(a);

// Transcendentals
const r = fpz.sqrt(Fixed.fromInt(2));  // ≈ 1.414
const e_value = fpz.exp(Fixed.ONE);    // = e
const two = fpz.pow(Fixed.fromInt(8), Fixed.rconst(1.0 / 3.0));  // ∛8 = 2

// SIMD batch (bit-identical to scalar)
const va = fpz.simd.fromArray(4, .{ Fixed.ONE, Fixed.HALF, Fixed.fromInt(2), Fixed.fromInt(3) });
const vb = fpz.simd.splat(4, Fixed.HALF);
const vc = fpz.simd.addVec(4, va, vb);
const arr = fpz.simd.toArray(4, vc);
```

## Public API

| Namespace | Surface |
| --- | --- |
| `fpz.Fixed` | `frac_bits`, `min_int`, `max_int`; `ZERO`, `ONE`, `HALF`, `NEG_ONE`, `MIN`, `MAX`, `PI`, `TWO_PI`, `HALF_PI`, `E`, `LN2`; `fromInt`, `fromRaw`, `rconst`, `toInt`, `roundToInt`, `toFloat`; `add`, `sub`, `mul`, `div`, `neg`, `abs`; `addSat`, `subSat`, `mulSat`; `eql`, `lt`, `lte`, `cmp`, `min`, `max`, `clamp` |
| `fpz.Angle` | `ZERO`, `QUARTER_TURN`, `HALF_TURN`, `THREE_QUARTER_TURN`; `fromRaw`, `fromRadians`, `toRadians`, `fromDegrees`, `toDegrees`; `addAngle`, `subAngle`, `eql` |
| `fpz` (top-level) | `sqrt`, `sqrtChecked`; `sin`, `cos`, `tan`, `atan2`; `exp`, `ln`, `lnChecked`, `log`, `pow` |
| `fpz.simd` | `Vec(N)`, `fromArray`, `toArray`, `splat`; `addVec`, `subVec`, `mulVec`, `negVec`; `eqlVec`, `ltVec`, `lteVec`, `minVec`, `maxVec` |
| `fpz.rounding` | `shiftRound`, `narrow`, `shiftRoundNarrow` — the single narrowing-shift RNE-ties-away helper used by every op |

`*Checked` variants return Zig error unions for domain violations instead of the defined fallback (`sqrtChecked → error.Negative`, `lnChecked → error.NonPositive`). The non-`Checked` form returns a deterministic value (e.g., `ln(0) → MIN`).

## Determinism contract

Per SPEC §2:

1. **Runtime is integer-only.** No libm, no runtime float, ever, on a sim path.
2. **Float is permitted only at comptime** (`rconst`, constant folding) and in test/debug display helpers (`toFloat`).
3. **Every op is total** — a defined result for all inputs, including overflow. Bare signed `+`/`-`/`*` are forbidden; ops use `@addWithOverflow` / i128 intermediates and either assert + return the defined wrap (`add`, `sub`, `mul`, `div`) or clamp (`addSat`, `subSat`, `mulSat`).
4. **Rounding is round-to-nearest, ties-away-from-zero**, applied identically at every narrowing shift via a single shared helper.
5. **Approximations are part of the contract.** The polynomial coefficients are frozen — changing them re-baselines the conformance vectors.
6. **Bit-exactness is tested, not assumed.**

**Corollary:** `Debug`, `ReleaseSafe`, and `ReleaseFast` produce **bit-identical** results. A divergence found by `ReleaseFast` fuzzing reproduces under a debugger at `-O0`, bit-for-bit. The `test-all-modes` build step enforces this on every run.

## Build steps

| Command | What it does |
| --- | --- |
| `zig build test` | Run all tests in the current optimize mode |
| `zig build test-all-modes` | Run tests under Debug, ReleaseSafe, AND ReleaseFast — fails if any disagree |
| `zig build test-cross-build` | Cross-compile tests for x86_64-linux, aarch64-linux, x86_64-windows, aarch64-macos, wasm32-wasi (build only — CI runs them) |
| `zig build gen-conformance` | Regenerate `src/conformance.zig` from current implementation output |
| `zig build run` | Run the demo executable |

## Conformance workflow

`src/conformance.zig` is the determinism artifact — a committed table of `(input, expected_raw)` pairs for every op. `src/conformance_test.zig` verifies the live implementation matches it bit-for-bit; any drift fails CI.

To re-baseline (e.g., after adopting a new polynomial):

```bash
zig build gen-conformance         # regenerates src/conformance.zig
zig build test-all-modes          # confirm the new baseline holds across modes
git diff src/conformance.zig      # review carefully — this IS the contract change
git commit -m "rebaseline conformance vectors: ..."
```

## Accuracy

| Function | Bound | Method |
| --- | --- | --- |
| `add`, `sub`, `neg`, `abs` | Exact | Integer arithmetic |
| `mul`, `div` | ±0.5 ULP | i128 intermediate + RNE narrow |
| `sqrt` | ±0.5 ULP | Restoring digit-by-digit u128 isqrt + RNE |
| `sin`, `cos` | ≤ 2⁻²² absolute | Octant reduction (bit-extract) + Taylor on [0, π/4] |
| `tan` | Bounded by `sin/cos`; saturates at poles | `sin(a) / cos(a)` with saturating narrow |
| `atan2` | ≤ ~2⁻¹⁸ absolute | Quadrant + swap + √2−1 inflection split + Taylor |
| `exp` | ≤ 2⁻²² relative (above ULP floor) | `k·LN2 + r` range reduction + Taylor, scale by 2^k |
| `ln` | ≤ 2⁻²² relative | Power-of-2 extraction + `atanh` series on mantissa ∈ [1, 2) |
| `log`, `pow` (fractional exp) | Derived from `ln` / `exp` | `ln(x)/ln(b)` / `exp(exp·ln(base))` |
| `pow` (integer exp) | Exact | Exponentiation by squaring |

The polynomial choices are frozen by the conformance vectors; they're documented in each module's header.

## Deviations from SPEC

- **`tan` at the polynomial pole — no assert, only saturate.** SPEC §6 says "saturate ... and assert". `std.debug.assert` becomes `unreachable` in `ReleaseFast`, which licenses the optimizer to delete the saturation branch as dead code and hit a real integer divide-by-zero. The defined saturation is the contract that survives in every build mode; an assert that violates the cross-mode bit-identity corollary (§2) is worse than no assert. `tan` near a pole returns `Fixed.MAX` or `Fixed.MIN` based on `sin`'s sign — same answer in Debug, ReleaseSafe, and ReleaseFast.

- **`ln` of non-positive input — no assert, only saturate.** Same reasoning as `tan`. `ln(x ≤ 0)` returns `Fixed.MIN`. Use `lnChecked` if you want an error union instead.

## Status

Phase 0–9 complete per SPEC.md. 110 tests pass under Debug, ReleaseSafe, and ReleaseFast (330/330 across modes). Library cross-compiles cleanly for x86_64-linux, aarch64-linux, x86_64-windows, aarch64-macos, wasm32-wasi.

Not yet covered: a CI workflow that actually runs the test binaries cross-arch under QEMU (the SPEC §8 "every target arch" guarantee is currently a local build-only check).
