# NeoGravity Series: Independent Verification Suite

Companion computational materials for **The Theory of NeoGravity & Ether Dynamics**
(J. S. Guagliardo, 2026).

## Purpose

Every central mathematical claim in the series is reduced here to executable
code. A reviewer can verify the entire derivation chain in under a minute
rather than checking algebra by hand.

## Running

```bash
pip install sympy numpy
python3 verify_neogravity.py
```

## What is checked

| # | Claim | Method |
|---|---|---|
| 1 | Hydrostatic equilibrium solution | SymPy `dsolve` on the first-order ODE |
| 2 | Exponential is **not** a 1/r power law | Series expansion at large r |
| 3 | Calibration w = 1/4 | Solve first-order coefficient against measurement |
| 4 | Second-order divergence 2.00 vs 1.75 | Taylor expansion of both refractive indices |
| 5 | Light deflection, three routes | Closed-form evaluation |
| 6 | Mercury perihelion | Standard secular-advance formula |
| 7 | Redshift and Shapiro delay | Closed-form evaluation |
| 8 | Frame dragging: topology vs prefactor | Explicit statement of what is fitted |
| 9 | PTA correlation functions | Hellings-Downs and scalar-transverse evaluated |
| 10 | SMBHB orbital frequencies | Kepler's third law across parameter space |
| — | Dimensional consistency | Unit audit of every composite quantity |

## Epistemic status, stated plainly

The suite distinguishes three categories and prints the classification for
each result:

- **DERIVED** — follows from the postulates with no free parameters:
  equilibrium profile, refractive index form, perihelion advance,
  gravitational redshift, Shapiro delay.
- **CALIBRATED** — fixed by matching one observation:
  - `w = 1/4` from the measured solar light deflection
  - `R_eff` from the measured frame-dragging amplitude
- **TOPOLOGY ONLY** — shared functional form, physically distinct mechanism:
  the frame-dragging dipole `sin(theta)/r^2`.

**The series has two fitted constants, not one.** Any statement to the
contrary is incorrect and is flagged by this suite.

## Known corrections enforced here

1. `exp(-A/r)` is an exponential whose *exponent* depends on 1/r. It is not
   a 1/r power law. The 1/r perturbation appears only at first order.
2. Stokes creeping flow (viscous, no-slip boundary) and the Kerr vacuum
   solution share dipole topology but not physical mechanism. The amplitude
   is matched, not derived.
3. Earth's moment of inertia is `I = 0.3307 M R^2` (centrally condensed),
   not `(2/5) M R^2`. Canonical `R_eff = 4920 m`.

## Constants

CODATA 2018 and IAU nominal values, listed at the top of the script.

## License

Released for verification and reuse. Attribution appreciated.
