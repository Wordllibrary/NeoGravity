# NeoGravity Verification Suite

Companion computational materials for **The Theory of NeoGravity & Ether Dynamics**

**John Salvatore Guagliardo**, Independent Researcher
ORCID: [0000-0003-0756-6886](https://orcid.org/0000-0003-0756-6886)
Archived DOI: [10.5281/zenodo.22080221](https://doi.org/10.5281/zenodo.22080221)

---

## Purpose

Every central mathematical claim in the series is reduced here to executable
code. A reviewer can verify the entire derivation chain in under a minute
rather than checking algebra by hand.

The suite is deliberately adversarial toward its own theory. It marks
`[FAIL]` on phrasings that overclaim, states plainly which constants are
fitted rather than derived, and identifies the results that are *not*
independent confirmation.

## Running

```bash
pip install sympy numpy
python3 verify_neogravity.py
```

Output is reproduced in `sample_output.txt`.

---

## Conventions (v2.0, canonical)

| Term | Definition |
|---|---|
| **Ether Plenum** (The Modern Ether) | The continuous, Lorentz-invariant medium proposed to fill space: quantum-vacuum zero-point energy, the Higgs field vacuum expectation value, and the universal electromagnetic background. |
| **Radiative Undertow** | The **inward**-directed reaction load returned by the Ether Plenum on matter, under Newton's Third Law, in response to the outward transport of mass-energy. Identified with gravitational attraction. |
| **Ether Density** (ρ_ether) | Local mass-energy density of the Plenum. **Increases** toward a gravitating mass. |
| **Refractive index** | n = √(ρ/ρ₀). A denser medium slows light, exactly as glass does. |

> **Note on v2.0.** Earlier drafts had the Plenum *rarefied* near mass with
> n = √(ρ₀/ρ). That convention produced identical results but inverted the
> physical picture. The compressive convention above is canonical: it yields
> the same w = 1/4, the same n(r), and the same falsifiable prediction, while
> matching both ordinary optical intuition and the inward direction of the
> Radiative Undertow.

---

## What is checked

| # | Claim | Method |
|---|---|---|
| 1 | Hydrostatic equilibrium of the Plenum | SymPy `dsolve` on a first-order ODE |
| 2 | The exponential is **not** a 1/r power law | Series expansion at large r |
| 3 | Calibration of w = 1/4 | First-order coefficient solved against measurement |
| 4 | Second-order divergence, 2.00 vs 1.75 | Taylor expansion of both refractive indices |
| 5 | Light deflection, three routes | Closed-form evaluation |
| 6 | Mercury perihelion advance | Secular advance formula |
| 7 | Gravitational redshift and Shapiro delay | Closed-form evaluation |
| 8 | Frame dragging: topology vs fitted prefactor | Explicit statement of what is matched |
| 9 | Two constraints on the Undertow mechanism | Radiation pressure and luminosity-to-mass ratios |
| 10 | PTA correlation functions | Hellings-Downs and scalar-transverse evaluated |
| 11 | SMBHB orbital frequencies | Kepler's third law across parameter space |
| — | Dimensional consistency | Unit audit of every composite quantity |

---

## Epistemic status, stated plainly

**DERIVED** — follows from the postulates with no free parameters:

- Inverse-square form of the Undertow (spherical shell geometry in ℝ³)
- Equilibrium density profile and refractive index form
- Mercury's perihelion advance
- Gravitational redshift
- Shapiro time delay

**CALIBRATED** — fixed by matching one observation:

- `w = 1/4` — from the measured solar light deflection
- `R_eff = 4920 m` — from the measured frame-dragging amplitude

**TOPOLOGY ONLY** — shared functional form, physically distinct mechanism:

- The frame-dragging dipole `sin(θ)/r²`

**CONSTRAINT** — required by observation, not supplied by the theory:

- The coupling tracks total mass-energy, not thermal luminosity

**OPEN** — not resolved by the framework:

- Derivation of `w = 1/4` from first principles
- Consistency with 10⁻¹⁷-level Lorentz-invariance bounds
- The entire strong-field regime
- Polarization content of Plenum oscillations (bears on PTA interpretation)

> ### The framework has two fitted constants, not one.
> Any statement to the contrary is incorrect and is flagged by this suite.

---

## Corrections enforced here

1. **`exp(A/r)` is not a 1/r power law.** The 1/r dependence sits in the
   exponent; a 1/r perturbation appears only at first order. Manuscripts
   claiming the integration "produces a 1/r profile" are flagged `[FAIL]`.

2. **Frame dragging is not "algebraically identical" to Stokes flow.**
   Stokes creeping flow requires viscosity and a no-slip boundary; Kerr is a
   vacuum solution. They share dipole topology only. The amplitude is matched
   by defining `R_eff`, so computing the ratio afterward returns 1 *by
   construction* — that is calibration, not derivation.

3. **Earth's moment of inertia is `I = 0.3307 M R²`**, not `(2/5)M R²`.
   Earth is centrally condensed. The uniform-sphere value overstates J by 21%
   and yields an incorrect `R_eff ≈ 5.24 km`. Canonical: **4920 m**.

4. **Light deflection is not independent confirmation.** `w` was calibrated
   on that measurement. The genuinely independent successes are Mercury's
   perihelion, gravitational redshift, and the Shapiro delay.

5. **Radiation pressure is not the Undertow.** Measured solar radiation
   pressure is ~6 × 10¹³ times too weak, and points outward rather than
   inward.

---

## Constants

CODATA 2018 and IAU nominal values, declared at the top of the script.

## Citation

> Guagliardo, J. S. (2026). *NeoGravity verification suite: Computational
> checks for the Theory of NeoGravity & Ether Dynamics* (Version 2.0)
> [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22080221

## License

Released for verification and reuse. Attribution appreciated.

