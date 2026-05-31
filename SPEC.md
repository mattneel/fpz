# Fixed-Point Math Library — Specification (Zig)

> Numeric substrate for a deterministic simulation. One inviolable contract:
> **identical inputs produce bit-identical outputs on every target arch, build mode, and toolchain.**
> Precision, speed, and ergonomics are all subordinate to that.

---

## 0. Scope

A single-format, 64-bit fixed-point type with total arithmetic, native integer-only
transcendentals, and a binary-angle companion type. No runtime radix, no FFI, no runtime float.
The library is the reference oracle for the engine's determinism story; the cross-platform
conformance suite (§8) *is* the guarantee, not an assumption.

---

## 1. Representation

```zig
pub const frac_bits: comptime_int = 24;          // Q40.24 default
pub const whole_bits: comptime_int = 64 - frac_bits;

pub const Fixed = struct { raw: i64 };           // distinct type, NOT an i64 alias
```

- **Distinct struct, not `typedef i64`.** The type system then blocks mixing raw integers with
  fixed-point values and blocks accidental integer `+`/`*` on the bits. This is where the type
  safety lives — methods only, no operators (Zig has none to overload).
- **One format.** The codebase uses exactly one `Fixed`. The binary point is the single comptime
  constant `frac_bits` and never appears as a magic literal anywhere else. No per-quantity zoo,
  no runtime radix argument.
- **Default Q40.24:** integer range ±5.49×10¹¹, resolution 2⁻²⁴ ≈ 5.96×10⁻⁸. At 1 unit = 1 m that
  is ±5.5×10⁸ km of world at ~60 nm precision — chosen because for a game sim positional *range*
  matters more than sub-nanometer precision. Changing the split is a one-line edit + re-baselining
  the golden vectors (§8).

---

## 2. Determinism contract (the spine)

1. **Runtime is integer-only.** No libm, no runtime float, ever, on a sim path.
2. **Float is permitted only at comptime** (constant folding, §4) and in test/debug display helpers.
3. **Every operation is total** — a defined result for all inputs, including overflow. Bare signed
   `+`/`-`/`*` are forbidden: they are UB in unsafe builds, which licenses the optimizer to break
   reproducibility. Use the overflow builtins and define the result (§3).
4. **Rounding is exactly specified and uniform:** round-to-nearest, ties away from zero, applied
   identically at every narrowing shift. A rounding divergence is a determinism divergence.
5. **Approximations are part of the contract.** The chosen polynomial / iteration count for each
   transcendental is frozen; changing it is a breaking change that requires re-baselining golden
   vectors.
6. **Bit-exactness is tested, not assumed** (§8).

A corollary worth stating: because every op is defined regardless of whether asserts are compiled
in, `Debug`, `ReleaseSafe`, and `ReleaseFast` produce **bit-identical** results. A divergence found
in a fast fuzzing run reproduces under a debugger at `-O0`, bit-for-bit.

---

## 3. Core arithmetic

```zig
pub fn fromInt(i: i64) Fixed;
pub fn toInt(self: Fixed) i64;          // truncates toward zero
pub fn roundToInt(self: Fixed) i64;     // round-to-nearest, ties away from zero

pub fn add(a: Fixed, b: Fixed) Fixed;   // @addWithOverflow; assert no overflow; return defined wrap
pub fn sub(a: Fixed, b: Fixed) Fixed;
pub fn mul(a: Fixed, b: Fixed) Fixed;   // i128 intermediate, >>frac_bits w/ RNE, narrow + overflow assert
pub fn div(a: Fixed, b: Fixed) Fixed;   // assert b.raw != 0; (i128<<frac_bits)/b, RNE, narrow
pub fn neg(a: Fixed) Fixed;             // assert a.raw != INT64_MIN
pub fn abs(a: Fixed) Fixed;             // assert a.raw != INT64_MIN

pub fn addSat(a: Fixed, b: Fixed) Fixed; // +| ; clamp to range
pub fn subSat(a: Fixed, b: Fixed) Fixed;
pub fn mulSat(a: Fixed, b: Fixed) Fixed;

pub fn eql/lt/lte/cmp(...) ...;
pub fn min/max(a, b) Fixed;
pub fn clamp(self, lo, hi) Fixed;
```

**Overflow policy — wrap-and-assert.** `add`/`sub`/`mul`/`div` compute via `@addWithOverflow` /
`@mulWithOverflow` / i128 intermediate, **assert** the result didn't overflow (fires in
`Debug`/`ReleaseSafe`, compiled out in `ReleaseFast`), and **return the defined wrapped value** so
the op is total in every build mode. Overflow is therefore a loud, locatable bug in safe builds and
a deterministic wrap in fast builds — never UB.

**Saturating variants** (`*Sat`) exist for the specific call sites where clamping is the correct
semantics — i.e. genuine singularities, not general arithmetic. Default to wrap-and-assert
everywhere else so overflow surfaces rather than being silently absorbed.

**Rounding detail:** arithmetic `>>` floors (toward −∞). To get round-to-nearest-ties-away on a
narrowing shift, add the rounding bias with the operand's sign before the shift. Specify and reuse
one helper; do not re-derive it per call site.

**`mul`/`div` intermediates:** `i64 → i128`, multiply/shift, narrow to `i64`. The i128 product
cannot overflow for i64 operands; the only overflow check is on the narrow back to `i64`.

---

## 4. Constants

```zig
pub fn rconst(comptime r: f64) Fixed;   // comptime float literal -> baked integer; comptime-only

pub const ONE:     Fixed = ...;         // 1 << frac_bits
pub const HALF:    Fixed = ...;
pub const PI:      Fixed = rconst(3.14159265358979323846);
pub const TWO_PI:  Fixed = ...;
pub const HALF_PI: Fixed = ...;
pub const E:       Fixed = rconst(2.71828182845904523536);
pub const LN2:     Fixed = rconst(0.69314718055994530942);
```

- `rconst` is **comptime-only** — the float math runs in the compiler, the result is an integer
  literal in the binary. No runtime float.
- The baked constants are part of the determinism contract → the conformance suite (§8) verifies
  they are bit-identical across toolchains (compile-time float rounding could otherwise drift by a
  ULP between compilers and silently desync the sim).
- A runtime `fromFloat` is **not** provided on sim paths. A `toFloat` exists for debug/display only.

---

## 5. Angles — Binary Angle Measure (BAM)

```zig
pub const Angle = struct { raw: u32 };  // full turn = 2^32

pub fn fromRadians(r: Fixed) Angle;
pub fn toRadians(a: Angle) Fixed;
pub fn fromDegrees(d: Fixed) Angle;
pub fn toDegrees(a: Angle) Fixed;
pub fn addAngle(a: Angle, b: Angle) Angle;  // u32 +% ; always defined
pub fn subAngle(a: Angle, b: Angle) Angle;
```

- **A full turn maps onto the full `u32` range**, so angle arithmetic wraps correctly by
  construction (`u32` overflow *is* mod 2π). This eliminates an entire class of
  angle-normalization bugs and removes any need to range-reduce by a low-precision π.
- Resolution: 360° / 2³² ≈ 8.4×10⁻⁸ degrees — far finer than any transcendental's accuracy.
- **Trig consumes `Angle`** (§6). Range reduction is then *bit extraction*: the top 2–3 bits select
  the quadrant/octant directly, with no modulo and no precision loss at large angles. Radian
  convenience wrappers may exist but the canonical entry point is the BAM angle.

---

## 6. Transcendentals

For each: signature, method, domain, error budget. Approach is specified; bodies and coefficients
are implementation.

```zig
pub fn sin(a: Angle) Fixed;
pub fn cos(a: Angle) Fixed;
pub fn tan(a: Angle) Fixed;             // saturates at poles
pub fn atan2(y: Fixed, x: Fixed) Angle;

pub fn sqrt(x: Fixed) Fixed;            // precondition x.raw >= 0
pub fn exp(x: Fixed) Fixed;             // saturates on overflow
pub fn ln(x: Fixed) Fixed;              // precondition x.raw > 0
pub fn log(x: Fixed, base: Fixed) Fixed;
pub fn pow(base: Fixed, exp: Fixed) Fixed;

pub fn sqrtChecked(x: Fixed) error{Negative}!Fixed;
pub fn lnChecked(x: Fixed) error{NonPositive}!Fixed;
```

- **`sqrt`** — digit-by-digit (restoring) integer square root over a `u128` of `raw << frac_bits`;
  exact, no Newton convergence-count guessing. Target: correctly rounded to ±1 ULP. Precondition
  `x ≥ 0`; `sqrtChecked` returns an error union for callers that handle it.
- **`sin`/`cos`** — octant reduction via the high bits of the `Angle` → minimax (or Taylor on the
  reduced interval) polynomial on `[0, π/4]`, sign/swap from the octant index. Target max abs error
  ≤ 2⁻²² (state the achieved bound in the accuracy table, §8). This is intentionally far better than
  the crude two-coefficient approximations common in embedded libs.
- **`tan`** — `sin/cos`; at the poles, saturate to `MAX`/`MIN` (defined value) and assert. No
  in-band sentinel.
- **`atan2`** — quadrant from the signs of `(y, x)`, argument reduction (swap when `|y| > |x|`),
  core `atan` on `[0, 1]` via polynomial, result returned as a BAM `Angle` (so the caller's angle
  arithmetic stays wrap-correct). `atan2(0, 0)` returns a documented value (0) plus an assert.
- **`exp`** — `k = round(x / LN2)`, `r = x − k·LN2`, polynomial for `exp(r)` on the small residual,
  scale by `2^k` via shift. Saturate (don't wrap) on overflow.
- **`ln`** — normalize to a mantissa in `[1, 2)` by extracting the integer power-of-two `k`,
  `atanh`-style series / polynomial on the mantissa, recombine `ln = k·LN2 + poly`. Precondition
  `x > 0`; `lnChecked` for graceful handling.
- **`log(x, base)`** — `ln(x) / ln(base)`, with the same domain guards.
- **`pow(base, exp)`** — `exp == 0 → ONE`; integer-exponent fast path by squaring (exact) when
  `exp` is whole; otherwise `exp(exp · ln(base))` with `base > 0` guarded.

---

## 7. Error policy

- **Preconditions are asserts** (Tiger Style): active in `Debug`/`ReleaseSafe`, compiled out in
  `ReleaseFast`.
- **No in-band sentinels — ever.** A domain violation is a bug, surfaced loudly (assert) or as a
  Zig error union in the `*Checked` variant. It is never a representable value masquerading as a
  result. (Explicit anti-pattern: returning a raw `-1` from `sqrt`, or a width-dependent bit
  pattern from `ln(0)`.)
- **The `ReleaseFast` fallback for each domain-fallible op is specified and deterministic**, so even
  with asserts off the behavior is reproducible: `sqrt(neg) → 0`, `ln(≤0) → MIN`, `div by 0` is an
  assert-only invariant (callers must not pass zero; document it).

---

## 8. Testing & conformance (the determinism gate)

- **Golden-vector suite.** Committed input→output tables, run on every target arch + toolchain in
  CI. *Any* diff fails the build. This is the determinism guarantee made empirical — and it doubles
  as the check for constant-folding drift (§4) and any C-style IDB corners (there are none, since
  this is native Zig — that's the point of writing it ourselves).
- **Differential tests.** Compare against `f64`/libm at test time (float allowed in tests) within a
  documented tolerance, per function. Pins accuracy.
- **Property tests.** Round-trips (`fromInt`/`toInt`, radians↔BAM), monotonicity where applicable,
  identities within tolerance: `sin² + cos² ≈ 1`, `exp(ln x) ≈ x`, `sqrt(x)² ≈ x`,
  `atan2(sin θ, cos θ) ≈ θ`.
- **Accuracy table.** Maintain a max/avg-error table per function (documented + asserted bounds).
- **Edge/overflow tests.** Range boundaries, the `INT64_MIN` neg/abs corner, saturation clamps,
  pole behavior of `tan`, domain fallbacks of §7.

---

## 9. SIMD batch path

- **Scalar `Fixed` is canonical and the reference.** Batch ops operate on `@Vector(N, i64)` for hot
  arrays (positions, velocities) and **must produce bit-identical results to the scalar path**
  (tested).
- `add`/`sub`/`compare` vectorize trivially. **`mul` does not** — it is a widening 64×64→128 then a
  shift, and autovectorizers choke on the lane-width change. Write the explicit widening sequence;
  do not rely on the optimizer to discover it.
- **Honest constraint:** 64×64→128 SIMD multiply is largely absent on common ISAs (AVX-512-DQ's
  `vpmullq` is low-64 only), so vectorized fix-mul at Q40.24 is limited. If SIMD-mul throughput ever
  becomes the bottleneck, the lever is a 32-bit companion format for those specific arrays — out of
  scope under "one format," noted here as the escape hatch rather than a present feature.

---

## 10. Config & build

- Comptime config: `frac_bits` (default 24), an assert-on-overflow toggle, and the rounding mode if
  ever exposed (default: round-to-nearest, ties away from zero).
- One exported canonical `Fixed` and `Angle`. No instantiation ceremony at call sites.
- Determinism holds across all three build modes (see the corollary in §2).

---

## 11. Non-goals (scope fence)

- **No runtime or per-quantity radix.** One compile-time format.
- **No GPU path.** Scalar + SIMD on CPU only.
- **No vector/matrix/geometry types here.** Those are a separate module built atop this one.
- **No audio/wave/ADSR/printf-style formatting.** Not this library's job.
- **No runtime float on any sim path.** Comptime constants and debug display only.
