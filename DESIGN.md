# DESIGN: Wide-angle beam patterns on the sphere

## 1. Purpose and requirements

TART antennas are wide field-of-view elements: a single antenna sees essentially
the entire sky above its local horizon. We image the **full sphere** (HEALPix),
so the beam model must live natively on the sphere — not on a tangent/`(l,m)`
plane. We need a compact, smooth, differentiable representation of an antenna
**beam pattern** (its directional power or voltage response) with the following
properties:

1. **Full hemisphere of support.** A beam covers up to a 180° field of view —
   from horizon, through boresight, to the opposite horizon — and is embedded in
   a full-sphere map.
2. **Polynomial.** The angular response is represented by polynomials. This gives
   a small coefficient vector, cheap vectorised evaluation, analytic derivatives,
   and a well-posed least-squares fit.
3. **Pointable.** A beam is defined in a *local* frame attached to its boresight
   and is pointed anywhere on the sphere by a rotation, with no re-fitting.
4. **Zero behind the antenna.** The response is identically zero in the rear
   hemisphere (below the antenna ground plane / horizon), going smoothly to zero
   at the horizon.
5. **Composable.** Many overlapping beams, pointed in different directions, are
   combined into a single full-sphere map.

The key idea: a polynomial in the Cartesian direction cosines of a local frame is
intrinsically a function on the whole sphere (it is a truncated spherical-harmonic
expansion). We use that directly, with a horizon taper and a hard rear-hemisphere
clip to satisfy requirement 4. No map projection is ever introduced.

---

## 2. Coordinate conventions

### 2.1 The global sphere

Sky directions are unit vectors `ŝ ∈ S²`. The combined all-sky product is stored
on a **HEALPix** grid (consistent with the rest of the TART tooling:
`radio-imaging`, `tart2ms`) — equal-area pixels, fast spherical operations,
genuinely full-sphere.

### 2.2 The local beam frame

Each beam carries an orthonormal frame `(ê₁, ê₂, ê₃)`:

- `ê₃ = p̂` — the **boresight** (pointing direction).
- `ê₁, ê₂` — span the plane orthogonal to `p̂`; their roll about `p̂` fixes the
  azimuth origin and the polarisation/orientation of the element.

For a sky direction `ŝ` we form the three **direction cosines** in this frame:

```
u = ŝ · ê₁ = sinθ cosφ
v = ŝ · ê₂ = sinθ sinφ
w = ŝ · ê₃ = cosθ
```

`θ` is the angle from boresight, `φ` the azimuth about it. These are dot products,
**defined over the entire sphere** — `w` ranges over `[−1, +1]`, not just a disk.
The front hemisphere is `w ≥ 0` and the rear (behind the antenna) is `w < 0`. The
horizon is the great circle `w = 0`. `(u, v, w)` is just `ŝ` re-expressed in the
beam frame: `u² + v² + w² = 1` everywhere.

---

## 3. Representing a single beam

### 3.1 Polynomials on the sphere = spherical harmonics

A polynomial in `(u, v, w)` restricted to the unit sphere *is* a (finite)
spherical-harmonic expansion: homogeneous polynomials of degree `≤ N` in three
variables, restricted to `S²`, span exactly the spherical harmonics up to degree
`N`. So "use polynomials" and "use a truncated spherical-harmonic series on the
sphere" are the same statement. We adopt this as the model:

```
B(ŝ) = Σ_{ℓ=0}^{N} Σ_{m=−ℓ}^{ℓ}  a_ℓ^m  Y_ℓ^m(θ, φ)
       ≡  Σ_{i+j+k ≤ N}  c_{ijk}  u^i v^j w^k          (same object, two bases)
```

- The **spherical-harmonic form** (`a_ℓ^m`) is the natural *fitting/interchange*
  representation: `Y_ℓ^m` is orthonormal on the sphere, so least-squares fits are
  well-conditioned and coefficients are individually meaningful. The truncation
  degree `N` is the single accuracy/compactness knob.
- The **Cartesian-polynomial form** (`c_{ijk}`) is the natural *runtime*
  representation: `u, v, w` are just dot products, so evaluation and pointing are
  bare arithmetic, fully vectorisable over a HEALPix map, with no special
  functions. The two forms are related by a fixed linear change of basis.

The `m` index carries azimuthal structure: `m = 0` terms are azimuthally
symmetric; `m ≠ 0` terms capture E-plane/H-plane asymmetry, dipole lobes, etc.

This is intrinsically a full-sphere basis — there is no disk, no tangent plane,
and no horizon coordinate singularity. A beam smooth to a few degrees needs only
`N ≲ 8–12`.

### 3.2 Special case: azimuthally symmetric beams

If the element is rotationally symmetric about boresight, only `m = 0` survives
and the beam is a 1-D polynomial in `w = cosθ` — i.e. a **Legendre** series:

```
B(θ) = Σ_ℓ a_ℓ · P_ℓ(cosθ)            (w = cosθ ∈ [−1, 1])
```

the cheap, common case (a single radial profile). It is just the `m = 0` subset
of §3.1, handled by the same code.

### 3.3 Enforcing "zero behind the antenna"

A polynomial on the sphere cannot be *exactly* zero over an entire hemisphere
(a nonzero polynomial vanishes only on a measure-zero set), so we impose the rear
zero explicitly, with a smooth horizon crossing:

1. **Hard support.** Evaluation returns `0` whenever `w < 0`. This guarantees
   exact zero in the rear hemisphere, independent of the coefficients.

2. **Smooth horizon roll-off.** A raw series does not generally vanish at `w = 0`,
   so a bare clip would leave a step at the horizon. Factor the model:

   ```
   B(ŝ) = w^q · Q(ŝ)        with integer q ≥ 1,   for w ≥ 0;   else 0
   ```

   where `Q` is the spherical-harmonic / polynomial series of §3.1. Since
   `w = cosθ → 0` at the horizon, the `w^q` factor drives the beam (and, for
   `q ≥ 2`, its derivatives) smoothly to zero there. `q` sets the order of contact
   at the horizon (`q = 1` continuous, `q = 2` C¹, …). Note `w^q · Q` is itself a
   polynomial in `(u, v, w)`, so the model stays polynomial.

The combination — spherical-harmonic core for shape, `w^q` factor for the horizon
taper, hard clip for the rear hemisphere — is smooth on the front hemisphere,
exactly zero behind, and continuous across the horizon. All three pieces are
pointing-invariant (they depend on `ŝ` only through the frame dot products).

---

## 4. Pointing a beam

A beam fitted in a canonical frame (boresight at `+z`) is pointed by the rotation
`R ∈ SO(3)` mapping the canonical frame onto `(ê₁, ê₂, ê₃)`:

```
[u, v, w]ᵀ = Rᵀ ŝ                  (project sky direction into the beam frame)
B_pointed(ŝ) = B_canonical(u, v, w)
```

In the **Cartesian-polynomial** form, only the three dot products change; the
coefficient vector is invariant, so pointing is O(1) in storage and trivial to
apply to a whole HEALPix map. (Equivalently, in the spherical-harmonic form the
`a_ℓ^m` rotate by a Wigner-D matrix — useful analytically, but the Cartesian form
needs no such machinery.) The boresight `p̂` plus a roll angle fully specify `R`.

---

## 5. Evaluation algorithm

```python
def evaluate(beam, s_hat):           # s_hat: (..., 3) unit vectors (HEALPix pixels)
    # 1. project into the local frame (three dot products over the full sphere)
    u = s_hat @ beam.e1
    v = s_hat @ beam.e2
    w = s_hat @ beam.e3              # = cos(theta from boresight), in [-1, 1]

    # 2. rear hemisphere -> exactly zero
    front = w >= 0.0

    # 3. evaluate the polynomial core (Cartesian poly, or Y_lm sum)
    core = poly_eval(beam.coeffs, u, v, w)

    # 4. horizon taper + hard clip
    B = (w ** beam.q) * core
    return np.where(front, B, 0.0)
```

Vectorised over the pixel axis; `poly_eval` is a Horner/precomputed-monomial
evaluation of `Σ c_{ijk} u^i v^j w^k`.

---

## 6. Fitting beams to data

Input: samples `B_i` at directions `ŝ_i` (EM simulation — `.ffe`/NEC patterns — or
holographic / drone measurement), naturally given on the sphere.

1. Project each sample into the canonical frame → `(u_i, v_i, w_i)`; weight down or
   discard rear-hemisphere samples (`w_i < 0`).
2. Divide out the horizon factor: fit `core_i = B_i / w_i^q` (guard near
   `w = 0`), or fold `w^q` into the design matrix.
3. Build the design matrix from the basis functions up to degree `N`
   (`Y_ℓ^m(θ_i, φ_i)`, or the monomials `u^i v^j w^k`).
4. Solve the regularised least-squares problem
   `a = argmin ‖M a − core‖² + λ‖a‖²`. Spherical-harmonic orthonormality keeps
   `MᵀM` well-conditioned; small `λ` (or a degree cap) controls overfitting/ringing.
5. Validate: residual RMS vs `N` and on a held-out direction set; raise `N` until
   the residual plateaus. Convert the fitted `a_ℓ^m` to the Cartesian `c_{ijk}`
   for runtime.

Store the coefficients, degree `N`, taper power `q`, and the local frame.

---

## 7. Combining beams into a full sphere

Each beam `k` covers a hemisphere about its boresight `p̂_k` and is zero behind.
Neighbouring beams overlap. We combine them into one full-sky product on the
shared HEALPix grid.

### 7.1 Coverage / sensitivity map

```
A(ŝ) = Σ_k |B_k(ŝ)|²
```

Every `B_k` is zero behind its own antenna, so the sum is naturally restricted and
overlaps add. `A(ŝ)` shows which parts of the sky are seen and how well.

### 7.2 Optimal mosaic (overlap weighting)

To merge per-beam sky estimates `I_k(ŝ)` (or to form a combined effective beam),
use **inverse-variance / primary-beam weighting** — the standard radio
mosaicking estimator, optimal in the overlap regions:

```
            Σ_k  B_k(ŝ) I_k(ŝ) / σ_k²
I(ŝ)  =  ───────────────────────────────
            Σ_k  B_k(ŝ)²       / σ_k²
```

with `σ_k²` the noise variance of beam `k`. Low-response directions contribute
little; the denominator normalises the overlap.

### 7.3 Partition of unity (smooth blending)

For a seamless stitched all-sky weighting that sums to one wherever there is
coverage:

```
w_k(ŝ) = B_k(ŝ)² / Σ_j B_j(ŝ)²        (where Σ_j B_j² > 0)
```

The `{w_k}` form a partition of unity over the covered sphere; rear-hemisphere
zeros mean each beam participates only where it physically can. Directions with
`Σ_j B_j² = 0` are flagged as gaps.

---

## 8. Data structures and interchange format

A beam serialises to JSON (one file per fitted element, or a list):

```json
{
  "name": "tart_dipole_v3",
  "frequency_hz": 1.57542e9,
  "basis": "real_spherical_harmonic",
  "degree_N": 10,
  "horizon_power_q": 2,
  "boresight": [0.0, 0.0, 1.0],
  "roll_axis_e1": [1.0, 0.0, 0.0],
  "coeffs": [
    {"l": 0, "m": 0,  "a": 1.0000},
    {"l": 2, "m": 0,  "a": -0.4120},
    {"l": 2, "m": 2,  "a": 0.0310}
  ]
}
```

Runtime object:

```
Beam:
    e1, e2, e3        # orthonormal frame (e3 = boresight)
    coeffs            # spherical-harmonic a_lm and/or cached Cartesian c_ijk
    q                 # horizon taper power
    evaluate(s_hat)   # -> response on the sphere, zero behind antenna
```

---

## 9. Extensions

- **Polarisation.** Replace scalar `B` by a 2×2 complex **Jones matrix** beam
  (four spherical-harmonic expansions). §3–§7 carry over component-wise;
  combination generalises to the coherency/Mueller form.
- **Frequency dependence.** Fit coefficients on a frequency grid and interpolate,
  or make the leading `a_ℓ^m` low-order polynomials in frequency. The TART band is
  narrow, so a few samples suffice.
- **Beam-width scaling.** Width scales with wavelength; scale the angle (e.g. via a
  `cosθ` remap) by `f/f₀` before evaluation rather than re-fitting per frequency.

---

## 10. Summary

- Work in the **direction cosines** `(u, v, w)` of a local beam frame — three dot
  products defined over the **whole sphere** (`w ∈ [−1, 1]`), no projection.
- Represent each beam as a **truncated spherical-harmonic series** (equivalently a
  **polynomial in `(u, v, w)`**), multiplied by a **`w^q` horizon taper** and
  **hard-clipped to zero for `w < 0`**: smooth on the front hemisphere, exactly
  zero behind the antenna, continuous at the horizon.
- **Point** beams by rotating the three dot products; coefficients are
  pointing-invariant.
- **Fit** by orthogonal least squares; **combine** overlapping beams directly on
  the HEALPix sphere via inverse-variance / partition-of-unity weighting to build
  the full-sky product.
