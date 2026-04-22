---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Review

## 1. Checklist of Key Concepts

### Section 6.1: Phase Diagrams

1. **One-Component $P$–$T$ Phase Diagram**
   - **Single-phase regions** (solid, liquid, gas) are separated by **coexistence curves** on which two phases are simultaneously in equilibrium.
   - **Triple point:** the unique $(T_{tp}, P_{tp})$ at which solid, liquid, and gas coexist.
   - **Critical point:** the point at which the liquid–vapor coexistence curve *terminates*; beyond it, liquid and gas become indistinguishable (a **supercritical fluid**).
   - The solid–liquid line does not terminate in a critical point because solid and liquid differ qualitatively (long-range translational order), while liquid and gas differ only quantitatively.

2. **Open Systems and the Chemical Potential**
   - A closed system cannot describe two coexisting phases by itself; each phase must be treated as an **open system** whose amount of matter can change.
   - The Gibbs differential generalizes to

   ```{math}
   dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i,
   ```

   with the **chemical potential** defined by

   ```{math}
   \mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,\,n_{j\neq i}}.
   ```

   - $\mu_i$ is **intensive** — the "price per mole" of adding substance to a phase — which makes it the natural variable for comparing phases of different size.

3. **Pure Substance: $\mu = \bar g$**
   - For a one-component system the Euler relation gives $G = \mu n$, so the chemical potential equals the molar Gibbs free energy:

   ```{math}
   \mu = \frac{G}{n} \equiv \bar g.
   ```

   - This resolves the apparent paradox in Section 5.1 that $G$ alone appeared to vanish for a one-component system at fixed $N$.

4. **Phase-Equilibrium Condition**
   - Transferring $dn_\alpha$ moles from phase $\alpha$ to phase $\beta$ at fixed $T$ and $P$ gives $dG = (\mu_\alpha - \mu_\beta)\,dn_\alpha$. Requiring $dG = 0$ for arbitrary transfers yields

   ```{math}
   \mu_\alpha(T,P) = \mu_\beta(T,P).
   ```

   - Physically: matter flows from higher $\mu$ to lower $\mu$ and stops when they are equal. Because this is a single scalar constraint on two variables, it picks out a *curve* in the $(T,P)$ plane — the coexistence curve.

5. **Graphical View: $\bar g$ vs. $T$ at Fixed $P$**
   - From $dG = -S\,dT + V\,dP$: $(\partial \bar g/\partial T)_P = -\bar s$.
   - Every phase has *negative* slope; the larger $\bar s$, the steeper (more negative) the slope.
   - Since $\bar s_g > \bar s_l > \bar s_s$, the gas line is steepest, then liquid, then solid. Each intersection marks a first-order transition at which the stable (lowest-$\bar g$) phase changes.

6. **Graphical View: $\bar g$ vs. $P$ at Fixed $T$**
   - From $dG = -S\,dT + V\,dP$: $(\partial \bar g/\partial P)_T = \bar v$.
   - Every phase has *positive* slope; the larger $\bar v$, the steeper the slope.
   - Since $\bar v_g \gg \bar v_l \gtrsim \bar v_s$, compressing the system favors condensed phases over gas.
   - **Water anomaly:** $\bar v_s > \bar v_l$, so the ice curve rises *faster* than the liquid curve and increasing $P$ near $T_{fus}$ can melt ice.

7. **Latent Heat and Entropy Jump at a First-Order Transition**
   - A **first-order transition** is one at which the first derivatives of $G$ ($S$ and $V$) are discontinuous across the coexistence curve.
   - At coexistence, $\mu_\alpha = \mu_\beta$ gives $\Delta \bar g_{tr} = 0$; combined with $\Delta \bar g = \Delta \bar h - T\,\Delta \bar s$:

   ```{math}
   \Delta \bar h_{tr} = T_{tr}\,\Delta \bar s_{tr}.
   ```

   - One calorimetric measurement ($\Delta \bar h_{tr}$) therefore supplies both quantities.

8. **Trouton's Rule and Water as an Exception**
   - For many ordinary liquids at the normal boiling point,

   ```{math}
   \Delta \bar s_{vap} \approx 85.9\ \mathrm{J\,mol^{-1}\,K^{-1}} \approx 10.3\,R.
   ```

   - Approximate constancy reflects that the entropy cost of vaporization is dominated by the phase change itself (volume increase, loss of short-range order) rather than the liquid's chemical identity.
   - Hydrogen-bonded liquids (water, lower alcohols) exceed the Trouton value because the liquid is more ordered than "average," so vaporization releases more entropy. For water, $\Delta \bar s_{vap} \approx 109\ \mathrm{J\,mol^{-1}\,K^{-1}}$, roughly 25% above Trouton.

---

### Section 6.2: Clapeyron and Clausius–Clapeyron

1. **Derivation of the Clapeyron Equation**
   - Displacing along a coexistence curve preserves $\mu_\alpha = \mu_\beta$, so $d\mu_\alpha = d\mu_\beta$.
   - Using the molar Gibbs differential $d\mu = -\bar s\,dT + \bar v\,dP$ for each phase and collecting terms gives the **Clapeyron equation**:

   ```{math}
   \frac{dP}{dT} = \frac{\Delta\bar s_{tr}}{\Delta\bar v_{tr}} = \frac{\Delta\bar h_{tr}}{T\,\Delta\bar v_{tr}}.
   ```

   - Exact for any first-order transition between pure phases.

2. **Sign and Magnitude of $dP/dT$**
   - The **sign** is set by $\Delta\bar v$, since $\Delta\bar h$ and $T$ are positive at any "normal" transition. Most substances have $\bar v_l > \bar v_s$ (positive melting slope) and $\bar v_g \gg \bar v_l$ (positive boiling slope).
   - The **magnitude** is set by $\Delta\bar v$: solid–liquid transitions have very small $\Delta\bar v$ and therefore nearly vertical coexistence curves; liquid–vapor transitions have very large $\Delta\bar v$ and shallow coexistence curves.
   - This is why melting points are nearly pressure-independent at ordinary pressures while boiling points shift strongly.

3. **Water: Negative-Slope Ice–Water Line**
   - Because $\bar v_s > \bar v_l$ for water, $\Delta\bar v_{fus} < 0$ and therefore $dP/dT < 0$ along the ice–water line. Numerically, $dP/dT \approx -135\ \mathrm{bar\,K^{-1}}$, equivalently $dT/dP \approx -7.4\ \mathrm{mK\,bar^{-1}}$.
   - The effect is far too small to explain ice skating by "pressure melting"; surface friction and a pre-existing liquid layer dominate.

4. **Liquid–Vapor Approximation**
   - Two approximations collapse the Clapeyron equation into a far more useful form for the liquid–vapor line:
     - $\bar v_g \gg \bar v_l$, so $\Delta\bar v_{vap} \approx \bar v_g$ (for water at 1 atm the ratio exceeds $10^3$).
     - The vapor is treated as an ideal gas, $\bar v_g = RT/P_{sat}$.

5. **Clausius–Clapeyron Equation (Differential Form)**

   ```{math}
   \frac{d\ln P_{sat}}{dT} = \frac{\Delta\bar h_{vap}}{RT^2}.
   ```

   - The $T^{-2}$ dependence (rather than $T^{-1}$) traces to the two factors of $T$ introduced by the approximations: one from the explicit $T$ in Clapeyron, one from the $RT$ in $\bar v_g$.

6. **Clausius–Clapeyron Equation (Integrated Form, Constant $\Delta\bar h_{vap}$)**

   ```{math}
   \ln\!\left(\frac{P_2}{P_1}\right) = -\frac{\Delta\bar h_{vap}}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right).
   ```

   - **Two-point form:** two $(T,P)$ points give $\Delta\bar h_{vap}$ directly.
   - **Linearized form:** a plot of $\ln P_{sat}$ vs. $1/T$ should be a straight line with slope $-\Delta\bar h_{vap}/R$, *if* the constant-$\Delta\bar h_{vap}$ approximation holds.

7. **Diagnosing Breakdown and Recovering $\Delta\bar h_{vap}(T)$**
   - A $\ln P$-vs-$1/T$ plot can look straight even when $\Delta\bar h_{vap}$ is varying. A **residual plot** is the reliable diagnostic: systematic curvature (not random scatter) indicates breakdown of the constant-$\Delta\bar h_{vap}$ assumption.
   - A **local-slope estimate** from finite differences gives

   ```{math}
   \Delta\bar h_{vap}(T) = R\,T^2\,\frac{d\ln P_{sat}}{dT}.
   ```

   - $\Delta\bar h_{vap}(T)$ decreases with rising $T$ and vanishes at the critical point, where the liquid and vapor phases merge.

8. **Molecular Reading of $\Delta\bar h_{vap}$**
   - For water, $\Delta\bar h_{vap} \approx 40\ \mathrm{kJ\,mol^{-1}}$ is consistent with breaking roughly two hydrogen bonds per molecule (each worth $\sim 20\ \mathrm{kJ\,mol^{-1}}$) on going from liquid to vapor — the macroscopic shadow of the intermolecular interactions behind the partition functions of Chapter 2.

---

## 2. Checklist of Most Important Equations

Below is a unified list of the major equations from Sections 6.1–6.2.

A. **Open-System Gibbs Differential**

```{math}
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i.
```

- **Applicability**: any multicomponent system whose composition may change. Reduces to the closed-system form when all $dn_i = 0$.

---

B. **Chemical Potential (Definition)**

```{math}
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,\,n_{j\neq i}}.
```

- **Applicability**: intensive quantity conjugate to the amount of component $i$.

---

C. **Pure-Substance Identity**

```{math}
G = \mu n,\qquad \mu = \bar g.
```

- **Applicability**: one-component systems. Chemical potential and molar Gibbs free energy are the same object.

---

D. **Phase-Equilibrium Condition**

```{math}
\mu_\alpha(T,P) = \mu_\beta(T,P).
```

- **Applicability**: two phases of a pure substance in equilibrium at fixed $T$ and $P$. A single scalar constraint on two variables defines the coexistence *curve* in the $(T,P)$ plane.

---

E. **Slopes of $\bar g$ Along Natural Axes**

```{math}
\left(\frac{\partial \bar g}{\partial T}\right)_P = -\bar s,
\qquad
\left(\frac{\partial \bar g}{\partial P}\right)_T = \bar v.
```

- **Applicability**: pure substance, single phase. Used graphically to identify the stable phase (lowest $\bar g$) and to read off entropy and molar volume from slopes.

---

F. **Latent Heat / Entropy Jump at a First-Order Transition**

```{math}
\Delta \bar h_{tr} = T_{tr}\,\Delta \bar s_{tr}.
```

- **Applicability**: any first-order transition $\alpha \to \beta$ at coexistence $(T_{tr}, P_{tr})$. Follows directly from $\Delta \bar g_{tr} = 0$ on the coexistence curve.

---

G. **Trouton's Rule**

```{math}
\Delta \bar s_{vap} \approx 85.9\ \mathrm{J\,mol^{-1}\,K^{-1}} \approx 10.3\,R.
```

- **Applicability**: empirical benchmark at the *normal* boiling point for ordinary (non-hydrogen-bonded, non-associated) liquids. Substantial upward deviations indicate unusual order in the liquid.

---

H. **Clapeyron Equation**

```{math}
\frac{dP}{dT} = \frac{\Delta\bar s_{tr}}{\Delta\bar v_{tr}} = \frac{\Delta\bar h_{tr}}{T\,\Delta\bar v_{tr}}.
```

- **Applicability**: exact for any first-order transition between pure phases. The sign of $dP/dT$ is set by $\Delta\bar v$; the magnitude is dominated by $|\Delta\bar v|$ (small for solid–liquid, large for liquid–vapor).

---

I. **Clausius–Clapeyron Equation (Differential Form)**

```{math}
\frac{d\ln P_{sat}}{dT} = \frac{\Delta\bar h_{vap}}{RT^2}.
```

- **Applicability**: liquid–vapor coexistence under two approximations — $\bar v_g \gg \bar v_l$ and an ideal vapor.

---

J. **Clausius–Clapeyron Equation (Integrated Form)**

```{math}
\ln\!\left(\frac{P_2}{P_1}\right) = -\frac{\Delta\bar h_{vap}}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right).
```

- **Applicability**: liquid–vapor coexistence with $\Delta\bar h_{vap}$ approximately constant over the interval $[T_1, T_2]$. A plot of $\ln P_{sat}$ vs. $1/T$ is a straight line with slope $-\Delta\bar h_{vap}/R$.

---

K. **Local-Slope Estimate of $\Delta\bar h_{vap}(T)$**

```{math}
\Delta\bar h_{vap}(T) = R\,T^2\,\frac{d\ln P_{sat}}{dT}.
```

- **Applicability**: recovers the temperature-dependent enthalpy of vaporization from vapor-pressure data via finite differences. $\Delta\bar h_{vap}(T) \to 0$ as $T \to T_c$.
