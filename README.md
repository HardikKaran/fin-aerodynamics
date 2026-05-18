# Hand Calculations & Methodology

**Author:** Hardik
**Purpose:** to set boundary conditions and sanity-check ranges for the CFD analysis.

---

## 1. Flight condition and effective angle of attack

The loads are evaluated at q-max:

- $M_\infty = 2.8$, $V_\infty = 924\ \text{m/s}$, $\rho_\infty = 0.3652\ \text{kg/m}^3$, $T_\infty = 226\ \text{K}$, $p_\infty = 33{,}980\ \text{Pa}$.
- $q_\infty = \tfrac{1}{2}\rho V^2 = 1.559 \times 10^5\ \text{Pa}$.

A lateral gust of 80 m/s (≈150 kt, conservative envelope) is superposed on the trajectory AoA of 10°. The gust induces an additional incidence:

$$\alpha_{\text{ind}} = \arctan\!\left(\frac{V_{\text{gust}}}{V_\infty}\right) \approx 4.95^\circ,
\qquad \alpha_{\text{eff}} = \alpha + \alpha_{\text{ind}} \approx 14.95^\circ.$$

This $\alpha_{\text{eff}}$ is the value used for all aero coefficients — it is the design-envelope angle, not the nominal trajectory angle.

## 2. Geometry and the LE-attachment criterion

The fin is a trapezoidal planform with a double-wedge (or wedge-leading-edge) section:

- Root chord $c_r = 0.6$ m, tip chord $c_t = 0.35$ m, span $b = 0.36$ m, thickness $t = 0.015$ m.
- LE sweep $\Lambda = 37.5^\circ$, wedge half-angle $\delta_{LE} = 10^\circ$.
- Planform area $S = \tfrac{1}{2}(c_r + c_t)b = 0.171\ \text{m}^2$.
- Mean aerodynamic chord $\bar{c} = \tfrac{2}{3}\dfrac{c_r^2 + c_r c_t + c_t^2}{c_r + c_t} = 0.486$ m.
- Aspect ratio $A_R = b^2/S = 0.758$.

The Mach number component normal to the LE is

$$M_{\perp,\,LE} = M_\infty \cos\Lambda = 2.221.$$

Because $M_{\perp,\,LE} > 1$, the LE shock attaches cleanly to the leading edge and 2D strip theory is admissible across the fin. This is the justification for treating the fin sectionally with oblique-shock / Prandtl–Meyer theory rather than needing a fully 3D conical-flow solution at the LE.

## 3. Sectional pressure coefficients via shock-expansion theory

The fin section is treated as a wedge at incidence $\alpha_{\text{eff}}$. Because $\alpha_{\text{eff}} > \delta_{LE}$ (15° > 10°), the lower surface still turns *into* the flow (compression → oblique shock) while the upper surface now turns *away* from it (expansion → Prandtl–Meyer fan). The deflection angles are:

$$\theta_{\text{lower}} = \alpha_{\text{eff}} + \delta_{LE} = 24.95^\circ,
\qquad \theta_{\text{upper}} = \alpha_{\text{eff}} - \delta_{LE} = 4.95^\circ.$$

### 3.1 Lower surface (oblique shock)

The $\theta$–$\beta$–$M$ relation (Ribera Vicent Eq. 3.16, Anderson Eq. 4.17):

$$\tan\theta = 2\cot\beta \, \frac{M_1^2 \sin^2\beta - 1}{M_1^2(\gamma + \cos 2\beta) + 2}.$$

This is implicit in $\beta$ for given $(M_1, \theta)$ and cannot be inverted in closed form. Hence the use of `scipy.optimize.brentq` (Brent's method) to bracket and solve for the weak-shock root numerically — it is a guaranteed-convergence, derivative-free root finder which is ideal here because the LHS minus RHS is continuous and changes sign across the root.

**Subtlety in the bracket.** $\theta(\beta)$ is not monotonic on $[\mu_1,\ \pi/2]$: it rises from zero at the Mach angle $\mu_1 = \arcsin(1/M_1)$, peaks at $\beta_{\max}$ (giving $\theta_{\max}$, the detachment limit), then falls back to zero at $\beta = \pi/2$ (normal shock). Naively bracketing on $[\mu_1, \pi/2]$ fails because both endpoints give $\theta - \theta_{\text{target}} < 0$, so `brentq` has no sign change to latch onto. The fix is to first locate $\beta_{\max}$ using `scipy.optimize.minimize_scalar` (bounded, Brent's algorithm on the negative of $\theta(\beta)$), then bracket the *weak* solution on $[\mu_1,\ \beta_{\max}]$. The strong solution would be on $[\beta_{\max},\ \pi/2]$ but is not physically realised for an attached wedge LE.

For $M_1 = 2.8$, $\theta = 24.95^\circ$: $\beta = 46.03^\circ$, $\beta_{\max} = 65.05^\circ$ (so the shock is comfortably attached — no detachment concern at this AoA).

The pressure jump across the oblique shock follows from the normal-shock relations applied to the normal Mach component (Ribera Vicent §3.3, Anderson Eq. 4.9):

$$\frac{p_2}{p_1} = 1 + \frac{2\gamma}{\gamma+1}\left(M_1^2 \sin^2\beta - 1\right) = 4.571.$$

### 3.2 Upper surface (Prandtl–Meyer expansion)

For an expansion turn $\theta_{\text{upper}}$ from $M_1$ to $M_2$ (Ribera Vicent Eq. 3.22–3.23):

$$\nu(M_2) = \nu(M_1) + \theta_{\text{upper}},
\qquad
\nu(M) = \sqrt{\tfrac{\gamma+1}{\gamma-1}}\,\arctan\sqrt{\tfrac{\gamma-1}{\gamma+1}(M^2-1)} - \arctan\sqrt{M^2-1}.$$

$\nu(M)$ is a monotonic function but again not analytically invertible, so `brentq` is used a second time to solve $\nu(M_2) = \nu(M_1) + \theta_{\text{upper}}$ for $M_2$. With $\nu(M_1) = 45.75^\circ$ and the turn of 4.95°, $\nu(M_2) = 50.70^\circ \Rightarrow M_2 = 3.049$. The pressure ratio across the expansion is isentropic:

$$\frac{p_2}{p_1} = \left(\frac{1 + \tfrac{\gamma-1}{2}M_1^2}{1 + \tfrac{\gamma-1}{2}M_2^2}\right)^{\gamma/(\gamma-1)} = 0.687.$$

### 3.3 Pressure coefficients and 2D sectional Cl / Cd

By the standard definition $C_p = (p - p_\infty)/q_\infty = \tfrac{2}{\gamma M_\infty^2}(p/p_\infty - 1)$ (Ribera Vicent §1, Anderson §11):

$$C_{p,\text{lower}} = 0.651, \qquad C_{p,\text{upper}} = -0.057.$$

Resolving the pressure-difference force into lift and lift-induced wave drag (Ribera Vicent Eq. 3.25 for the flat-plate analogue):

$$C_{\ell,\,2D} = (C_{p,L} - C_{p,U})\cos\alpha_{\text{eff}} = 0.684,$$
$$C_{d,\,\text{wave-lift}} = (C_{p,L} - C_{p,U})\sin\alpha_{\text{eff}} = 0.183.$$

For reference, Ackeret's linear theory (Ribera Vicent Eq. 3.29) on a flat plate gives $C_\ell = 4\alpha_{\text{eff}}/\sqrt{M^2-1} = 0.399$. The shock-expansion value sits ~70 % higher because the lower-face turn of nearly 25° is well outside the small-angle regime where Ackeret is valid (the notes warn the linear estimate is already 10 % low at 4°, Ribera Vicent §3.11).

## 4. Three-dimensional lift coefficient

The 2D value above is essentially a *root-strip* estimate (no tip relief, no spanwise pressure equalisation). A 3D estimate uses a finite-span supersonic lift-curve slope:

$$C_{L,\,3D} = \left(\frac{dC_L}{d\alpha}\right)_{3D}\!\alpha_{\text{eff}} = 1.479 \times 14.95^\circ \cdot \tfrac{\pi}{180} = 0.386.$$

Interpretation of the two $C_L$ values: they are two different methods applied to the same fin. The 2D shock-expansion value (0.684) captures nonlinear compression but ignores tip relief; the 3D linear-slope value (0.386) captures *finite-span 3D effects* but uses a linear lift curve (assumed). The truth lies between them. For structural sizing the 2D number is the conservative upper bound; for trim and performance the 3D number is more representative. **CFD will be the arbiter here**, and bracketing it between these two limits is a useful validation check.

## 5. Total sectional drag coefficient

Three contributions are summed:

- **Wave drag from lift** (already computed): $C_{d,\text{wave-lift}} = 0.183$.
- **Wave drag from thickness** (thin-aerofoil supersonic theory, Anderson §12.3): $\tau = t/\bar{c} = 0.0309$,
$$C_{d,\text{thick}} = \frac{4\tau^2}{\sqrt{M^2-1}} = 1.5 \times 10^{-3}.$$
- **Skin friction**, Prandtl–Schlichting 1/7-power law (Schlichting §21) corrected for compressibility by Frankl–Voishel (Gudmundsson Ch. 15):
$$C_f = \frac{0.0592}{Re^{0.2}}\bigl(1 + 0.144\,M^2\bigr)^{-0.65},
\qquad C_{d,\text{friction}} = 2 C_f = 3.3 \times 10^{-3},$$
with $Re_{\bar{c}} = \rho V \bar{c} / \mu = 5.16 \times 10^6$ (turbulent assumed throughout; transition near the LE at this Re is reasonable).

Total: $C_d \approx 0.187$, dominated (~98%) by lift-induced wave drag. Thickness and friction are essentially negligible at this AoA — a useful insight: at q-max the fin's drag bill is set almost entirely by the angle at which it slices the flow, not by its profile or finish.

## 6. Steady-state roll rate from canted fins

The 2° cant angle drives the rocket into a steady-state roll. The derivation follows OpenRocket §3.3 (Niskanen 2013), equating roll-forcing and roll-damping moments.

For one fin at cant $\delta$ and roll rate $\omega$, each spanwise strip sees a local incidence $\delta - \omega(r+y)/V$ (cant minus the roll-induced AoA). Integrating the strip lift × moment arm and setting net moment to zero gives:

$$\omega_{ss} = \frac{V \tan\delta}{y_{\text{eff}}},
\qquad y_{\text{eff}} = \frac{\displaystyle\int_0^b c(y)\,(r+y)^2\,dy}{\displaystyle\int_0^b c(y)\,(r+y)\,dy} = \frac{I_2}{I_1}.$$

The numerator and denominator are first and second moments of the chord distribution weighted by the (radial) moment arm — they map directly onto the roll-forcing and roll-damping sums in OpenRocket Eqs. 3.66 and 3.69. For a trapezoidal $c(y) = c_r + (c_t - c_r)y/b$, both integrals have closed form (OpenRocket Eq. 3.70), but `scipy.integrate.quad` is used instead because it generalises trivially to any planform (the same code works if the fin shape changes during design iteration). `quad` uses adaptive Gauss–Kronrod quadrature and is essentially exact for these smooth polynomial integrands.

Result: $y_{\text{eff}} = 0.326$ m, $\omega_{ss} = 99.1\ \text{rad/s} \approx 15.8\ \text{rev/s}$.

## 7. Centripetal load on a fin

With the fin spinning at $\omega_{ss}$ about the rocket centreline:

$$F_{\text{cent}} = m_f \, \omega_{ss}^2 \, r_{\text{CG}}.$$

$r_{\text{CG}}$ is the mass centroid of the fin from the rotation axis. For a trapezoid of uniform areal density,

$$\bar{y} = \frac{b}{3}\,\frac{c_r + 2c_t}{c_r + c_t} = 0.164\ \text{m},
\qquad r_{\text{CG}} = r + \bar{y} = 0.289\ \text{m},$$

giving $F_{\text{cent}} \approx 4256\ \text{N}$. 

## 8. Leading-edge stagnation temperature

From the isentropic stagnation relation:

$$T_0 = T_\infty\!\left(1 + \tfrac{\gamma-1}{2}M^2\right) = 580\ \text{K} \approx 307^\circ\text{C}.$$

The LE is hot but well below the softening point of common aerospace aluminium alloys (~573 K for 7075-T6 onset of property degradation). Carbon-fibre composites would warrant more careful checking.

---

## Areas of confusion and how they were resolved

1. **Non-monotonic $\theta$–$\beta$ relation.** The first attempt at `brentq` on $[\mu_1, \pi/2]$ failed because both endpoints have the same sign. Resolved by locating $\beta_{\max}$ first via `minimize_scalar`, then bracketing the weak-shock branch on $[\mu_1, \beta_{\max}]$. Still open: it would be cleaner to also bracket and report the strong-shock solution as a diagnostic, even though it is not physically realised at q-max.

2. **Two competing $C_L$ values.** The 2D shock-expansion result (0.684) and the 3D linear-slope result (0.386) differ by a factor of 1.8. Resolved by reframing them as *bracketing bounds* rather than competing point estimates — the 2D method handles nonlinearity, the 3D method handles finite span, and CFD must land between them. Still open: the source of $dC_L/d\alpha = 1.479$ is currently hard-coded; where did this come from?.

---

## References

1. Ribera Vicent, M. *Compressible Aerodynamics — AERO50001*, Imperial College London, 2025/26. (θ–β–M relation, oblique-shock pressure ratio, Prandtl–Meyer function, shock-expansion theory, Ackeret theory, pressure coefficient definition — §§3.3, 3.4, 3.9, 3.10, 3.11.)
2. Anderson, J. D. *Fundamentals of Aerodynamics*, McGraw-Hill, 6th ed., 2017. (Chs. 9–12 for the same material from the standard textbook; Ch. 12 for supersonic thin-aerofoil wave drag.)
3. Niskanen, S. *OpenRocket Technical Documentation*, v13.05, 2013. (§3.3 — roll-forcing and roll-damping coefficients, equilibrium roll frequency, Eqs. 3.66–3.74.)
4. Gudmundsson, S. *General Aviation Aircraft Design: Applied Methods and Procedures*, Butterworth-Heinemann, 2nd ed., 2022, Ch. 15 "Aircraft Drag Analysis". (Frankl–Voishel compressibility correction for turbulent skin friction.)
5. Schlichting, H. *Boundary-Layer Theory*, Springer, 8th ed., 2000. (Prandtl–Schlichting 1/7-power-law flat-plate skin friction; recovery temperature for turbulent boundary layers.)
6. NACA Report 1135, *Equations, Tables, and Charts for Compressible Flow*, Ames Research Staff, 1953. (Tabulated oblique-shock and Prandtl–Meyer values used as a cross-check of the numerical solves.)
7. `scipy.optimize.brentq` — Brent (1973), *Algorithms for Minimization Without Derivatives*, Prentice-Hall.
8. `scipy.integrate.quad` — Piessens et al., *QUADPACK: A Subroutine Package for Automatic Integration*, Springer, 1983.