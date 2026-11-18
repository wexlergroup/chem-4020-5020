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

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=6&amp;here=review"
  title="Course map, opened at Module 6"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-6.png
:class: hidden
:alt: Course map with Module 6 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6, which has no microscopic step. Below the map, one card for each of Sections 6.1 and 6.2 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## 1. Checklist of Key Concepts

### Section 6.1: Phase Diagrams

See [Section 6.1](section-01.md).

1. **One-Component $P$–$T$ Phase Diagram**
   - **Single-phase regions** (solid, liquid, gas) are separated by **coexistence curves** on which two phases are simultaneously in equilibrium.
   - **Triple point:** the unique $(T_{\mathrm{tp}}, P_{\mathrm{tp}})$ at which solid, liquid, and gas coexist.
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

   - $\mu_i$ is **intensive**: it is the "price per mole" of adding substance to a phase.

3. **Pure Substance: $\mu = \bar g$**
   - For a one-component system the Euler relation gives $G = \mu n$, so the chemical potential equals the molar Gibbs free energy:

   ```{math}
   \mu = \frac{G}{n} \equiv \bar g.
   ```

   - This is the Euler relation of [Section 5.1](../chapter-05/section-01.md), $G = \mu N$, written per mole.

4. **Phase-Equilibrium Condition**
   - Letting the amount of phase $\alpha$ change by $dn_\alpha$ at fixed $T$ and $P$, with $dn_\beta = -dn_\alpha$, gives $dG = (\mu_\alpha - \mu_\beta)\,dn_\alpha$. Requiring $dG = 0$ for an arbitrary small $dn_\alpha$ yields

   ```{math}
   \mu_\alpha(T,P) = \mu_\beta(T,P).
   ```

   - Physically: matter flows from higher $\mu$ to lower $\mu$ and stops when they are equal. Because this is a single scalar constraint on two variables, it picks out a *curve* in the $(T,P)$ plane — the coexistence curve.

5. **Graphical View: $\bar g$ vs. $T$ at Fixed $P$**
   - From $dG = -S\,dT + V\,dP$: $(\partial \bar g/\partial T)_P = -\bar s$.
   - Every phase has *negative* slope; the larger $\bar s$, the steeper (more negative) the slope.
   - Since $\bar s_{\mathrm{g}} > \bar s_{\mathrm{l}} > \bar s_{\mathrm{s}}$, the gas line is steepest, then liquid, then solid. Each intersection marks a first-order transition at which the stable (lowest-$\bar g$) phase changes.

6. **Graphical View: $\bar g$ vs. $P$ at Fixed $T$**
   - From $dG = -S\,dT + V\,dP$: $(\partial \bar g/\partial P)_T = \bar v$.
   - Every phase has *positive* slope; the larger $\bar v$, the steeper the slope.
   - Since $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}} \gtrsim \bar v_{\mathrm{s}}$, compressing the system favors condensed phases over gas.
   - **Water anomaly:** $\bar v_{\mathrm{s}} > \bar v_{\mathrm{l}}$ (ice is less dense than liquid water). Concept Check 5 of Section 6.1 asks what follows for the two curves.

7. **Latent Heat and Entropy Jump at a First-Order Transition**
   - A **first-order transition** is one at which the first derivatives of $G$ ($S$ and $V$) are discontinuous across the coexistence curve.
   - At coexistence, $\mu_\alpha = \mu_\beta$ gives $\Delta \bar g_{\mathrm{trs}} = 0$; combined with $\Delta \bar g = \Delta \bar h - T\,\Delta \bar s$:

   ```{math}
   \Delta \bar h_{\mathrm{trs}} = T_{\mathrm{trs}}\,\Delta \bar s_{\mathrm{trs}}.
   ```

   - One calorimetric measurement ($\Delta \bar h_{\mathrm{trs}}$) therefore supplies both quantities.

8. **Trouton's Rule and Water as an Exception**
   - For many ordinary liquids at the normal boiling point,

   ```{math}
   \Delta \bar s_{\mathrm{vap}} \approx 85.9\ \mathrm{J\,mol^{-1}\,K^{-1}} \approx 10.3\,R.
   ```

   - Approximate constancy reflects that the entropy gained on vaporization is dominated by the phase change itself (volume increase, loss of short-range order) rather than the liquid's chemical identity.
   - Hydrogen-bonded liquids (water, the lower alcohols, ammonia) lie well above the Trouton value. For water, $\Delta \bar s_{\mathrm{vap}} \approx 109\ \mathrm{J\,mol^{-1}\,K^{-1}}$, roughly 25% above Trouton. Liquids that boil at very low temperatures (neon, hydrogen, helium) lie well below it.

---

### Section 6.2: Clapeyron and Clausius–Clapeyron

See [Section 6.2](section-02.md).

1. **Derivation of the Clapeyron Equation**
   - Displacing along a coexistence curve preserves $\mu_\alpha = \mu_\beta$, so $d\mu_\alpha = d\mu_\beta$.
   - Using the molar Gibbs differential $d\mu = -\bar s\,dT + \bar v\,dP$ for each phase and collecting terms gives the **Clapeyron equation**:

   ```{math}
   \frac{dP}{dT} = \frac{\Delta\bar s_{\mathrm{trs}}}{\Delta\bar v_{\mathrm{trs}}} = \frac{\Delta\bar h_{\mathrm{trs}}}{T\,\Delta\bar v_{\mathrm{trs}}}.
   ```

   - Exact for any first-order transition between pure phases.

2. **Sign and Magnitude of $dP/dT$**
   - The **sign** is set by $\Delta\bar v$, since $\Delta\bar h$ and $T$ are positive at any "normal" transition. Most substances have $\bar v_{\mathrm{l}} > \bar v_{\mathrm{s}}$ (positive melting slope) and $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$ (positive boiling slope).
   - The **magnitude** is set by $\Delta\bar v$: solid–liquid transitions have very small $\Delta\bar v$ and therefore nearly vertical coexistence curves; liquid–vapor transitions have very large $\Delta\bar v$ and shallow coexistence curves.
   - This is why melting points are nearly pressure-independent at ordinary pressures while boiling points shift strongly.

3. **Water: Negative-Slope Ice–Water Line**
   - Because $\bar v_{\mathrm{s}} > \bar v_{\mathrm{l}}$ for water, $\Delta\bar v_{\mathrm{fus}} < 0$ and therefore $dP/dT < 0$ along the ice–water line. Numerically, $dP/dT \approx -135\ \mathrm{bar\,K^{-1}}$, equivalently $dT/dP \approx -7.4\ \mathrm{mK\,bar^{-1}}$.
   - The effect is far too small to explain ice skating by "pressure melting"; surface friction and a pre-existing liquid layer dominate.

4. **Liquid–Vapor Approximation**
   - Two approximations collapse the Clapeyron equation into a far more useful form for the liquid–vapor line:
     - $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$, so $\Delta\bar v_{\mathrm{vap}} \approx \bar v_{\mathrm{g}}$ (for water at 1 atm the ratio exceeds $10^3$).
     - The vapor is treated as an ideal gas, $\bar v_{\mathrm{g}} = RT/P_{\mathrm{sat}}$.

5. **Clausius–Clapeyron Equation (Differential Form)**

   ```{math}
   \frac{d\ln P_{\mathrm{sat}}}{dT} = \frac{\Delta\bar h_{\mathrm{vap}}}{RT^2}.
   ```

   - The $T^{-2}$ dependence (rather than $T^{-1}$) traces to two steps of the derivation, each of which brings in a factor of $T$.

6. **Clausius–Clapeyron Equation (Integrated Form, Constant $\Delta\bar h_{\mathrm{vap}}$)**

   ```{math}
   \ln\!\left(\frac{P_2}{P_1}\right) = -\frac{\Delta\bar h_{\mathrm{vap}}}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right).
   ```

   - **Two-point form:** two $(T,P)$ points give $\Delta\bar h_{\mathrm{vap}}$ directly.
   - **Linearized form:** a plot of $\ln P_{\mathrm{sat}}$ vs. $1/T$ should be a straight line with slope $-\Delta\bar h_{\mathrm{vap}}/R$, *if* the constant-$\Delta\bar h_{\mathrm{vap}}$ approximation holds.

7. **Diagnosing Breakdown and the Local Slope**
   - A $\ln P$-vs-$1/T$ plot can look straight even when $\Delta\bar h_{\mathrm{vap}}$ is varying. A **residual plot** is the reliable diagnostic: systematic curvature (not random scatter) indicates breakdown of the constant-$\Delta\bar h_{\mathrm{vap}}$ assumption.
   - The **local slope** from finite differences obeys the exact relation

   ```{math}
   R\,T^2\,\frac{d\ln P_{\mathrm{sat}}}{dT} = \frac{\Delta\bar h_{\mathrm{vap}}}{\Delta Z},
   \qquad
   \Delta Z \equiv \frac{P_{\mathrm{sat}}\,\Delta\bar v_{\mathrm{vap}}}{RT}.
   ```

   - Where the vapor is nearly ideal ($\Delta Z \approx 1$), the local slope estimates $\Delta\bar h_{\mathrm{vap}}(T)$; it exceeds $\Delta\bar h_{\mathrm{vap}}$ as the vapor becomes denser (by 8 % for water at 180 °C).
   - $\Delta\bar h_{\mathrm{vap}}(T)$ decreases with rising $T$. At the critical point, where the liquid and vapor phases merge, $\Delta\bar h_{\mathrm{vap}}$ and $\Delta\bar v_{\mathrm{vap}}$ vanish together.

8. **Molecular Reading of $\Delta\bar h_{\mathrm{vap}}$**
   - For water, $\Delta\bar h_{\mathrm{vap}} \approx 40\ \mathrm{kJ\,mol^{-1}}$ is consistent with breaking roughly two hydrogen bonds per molecule (each worth $\sim 20\ \mathrm{kJ\,mol^{-1}}$) on going from liquid to vapor — the macroscopic shadow of the intermolecular forces of [Section 1.4](../chapter-01/section-04.md).

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
\Delta \bar h_{\mathrm{trs}} = T_{\mathrm{trs}}\,\Delta \bar s_{\mathrm{trs}}.
```

- **Applicability**: any first-order transition $\alpha \to \beta$ at coexistence $(T_{\mathrm{trs}}, P_{\mathrm{trs}})$. Follows directly from $\Delta \bar g_{\mathrm{trs}} = 0$ on the coexistence curve.

---

G. **Trouton's Rule**

```{math}
\Delta \bar s_{\mathrm{vap}} \approx 85.9\ \mathrm{J\,mol^{-1}\,K^{-1}} \approx 10.3\,R.
```

- **Applicability**: empirical benchmark at the *normal* boiling point for ordinary (non-hydrogen-bonded, non-associated) liquids. Hydrogen-bonded liquids lie well above it, and liquids that boil at very low temperatures lie well below it.

---

H. **Clapeyron Equation**

```{math}
\frac{dP}{dT} = \frac{\Delta\bar s_{\mathrm{trs}}}{\Delta\bar v_{\mathrm{trs}}} = \frac{\Delta\bar h_{\mathrm{trs}}}{T\,\Delta\bar v_{\mathrm{trs}}}.
```

- **Applicability**: exact for any first-order transition between pure phases. The sign of $dP/dT$ is set by $\Delta\bar v$; the magnitude is dominated by $|\Delta\bar v|$ (small for solid–liquid, large for liquid–vapor).

---

I. **Clausius–Clapeyron Equation (Differential Form)**

```{math}
\frac{d\ln P_{\mathrm{sat}}}{dT} = \frac{\Delta\bar h_{\mathrm{vap}}}{RT^2}.
```

- **Applicability**: liquid–vapor coexistence under two approximations — $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$ and an ideal vapor.

---

J. **Clausius–Clapeyron Equation (Integrated Form)**

```{math}
\ln\!\left(\frac{P_2}{P_1}\right) = -\frac{\Delta\bar h_{\mathrm{vap}}}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right).
```

- **Applicability**: liquid–vapor coexistence with $\Delta\bar h_{\mathrm{vap}}$ approximately constant over the interval $[T_1, T_2]$. A plot of $\ln P_{\mathrm{sat}}$ vs. $1/T$ is a straight line with slope $-\Delta\bar h_{\mathrm{vap}}/R$.

---

K. **Local Slope of the Vapor-Pressure Curve**

```{math}
R\,T^2\,\frac{d\ln P_{\mathrm{sat}}}{dT} = \frac{\Delta\bar h_{\mathrm{vap}}}{\Delta Z},
\qquad
\Delta Z \equiv \frac{P_{\mathrm{sat}}\,\Delta\bar v_{\mathrm{vap}}}{RT}.
```

- **Applicability**: exact along the liquid–vapor coexistence curve. Where the vapor is nearly ideal ($\Delta Z \approx 1$), the local slope from finite differences of vapor-pressure data estimates the temperature-dependent $\Delta\bar h_{\mathrm{vap}}(T)$. As $T \to T_c$, $\Delta\bar h_{\mathrm{vap}} \to 0$ and $\Delta\bar v_{\mathrm{vap}} \to 0$ together.

## 3. Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. Liquid or vapor? Toluene near its boiling point

Draws on [Section 6.1](section-01.md).

Toluene boils at $T_{\mathrm b}=383.75\ \mathrm{K}$ under $1\ \mathrm{atm}$, where $\Delta\bar h_{\mathrm{vap}}=33.23\ \mathrm{kJ\,mol^{-1}}$. At that temperature the molar heat capacities are $\bar c_P=184.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for the liquid and $137.8\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for the vapor ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/fluid/)). The pressure is $1\ \mathrm{atm}$ throughout, and $\Delta\bar g_{\mathrm{vap}}=\bar g_{\mathrm g}-\bar g_{\mathrm l}=\mu_{\mathrm g}-\mu_{\mathrm l}$.

**(a)** Find $\Delta\bar s_{\mathrm{vap}}$ at $T_{\mathrm b}$, and compare it with the value of Trouton's rule.

**(b)** Treat $\Delta\bar h_{\mathrm{vap}}$ and $\Delta\bar s_{\mathrm{vap}}$ as independent of temperature. Find $\Delta\bar g_{\mathrm{vap}}$ at $363.75\ \mathrm{K}$ and at $403.75\ \mathrm{K}$, and name the stable phase at each temperature.

**(c)** Treat the two heat capacities as constant instead. Find $\Delta\bar h_{\mathrm{vap}}$, $\Delta\bar s_{\mathrm{vap}}$, and $\Delta\bar g_{\mathrm{vap}}$ at the same two temperatures. By how much does each $\Delta\bar g_{\mathrm{vap}}$ differ from the result of part (b)?

````{admonition} Solution
:class: dropdown

1. **Entropy of vaporization**

   At the boiling point the two phases coexist, so $\Delta\bar g_{\mathrm{vap}}=\Delta\bar h_{\mathrm{vap}}-T_{\mathrm b}\,\Delta\bar s_{\mathrm{vap}}=0$ and

   ```{math}
   \Delta\bar s_{\mathrm{vap}}=\frac{\Delta\bar h_{\mathrm{vap}}}{T_{\mathrm b}}=\frac{33\,230\ \mathrm{J\,mol^{-1}}}{383.75\ \mathrm{K}}=86.59\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   ```

   Trouton's rule gives $85.9\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Toluene lies $0.8\ \%$ above it.

2. **Constant $\Delta\bar h_{\mathrm{vap}}$ and $\Delta\bar s_{\mathrm{vap}}$**

   Because $\Delta\bar h_{\mathrm{vap}}=T_{\mathrm b}\,\Delta\bar s_{\mathrm{vap}}$,

   ```{math}
   \Delta\bar g_{\mathrm{vap}}(T)=\Delta\bar h_{\mathrm{vap}}-T\,\Delta\bar s_{\mathrm{vap}}=\Delta\bar s_{\mathrm{vap}}\,(T_{\mathrm b}-T).
   ```

   The two temperatures lie $20\ \mathrm{K}$ below and $20\ \mathrm{K}$ above $T_{\mathrm b}$, so

   ```{math}
   \begin{aligned}
   \Delta\bar g_{\mathrm{vap}}(363.75\ \mathrm{K})&=(86.593)(+20)=+1732\ \mathrm{J\,mol^{-1}},\\
   \Delta\bar g_{\mathrm{vap}}(403.75\ \mathrm{K})&=(86.593)(-20)=-1732\ \mathrm{J\,mol^{-1}}.
   \end{aligned}
   ```

   At $363.75\ \mathrm{K}$, $\mu_{\mathrm g}>\mu_{\mathrm l}$: vapor at $1\ \mathrm{atm}$ would condense spontaneously, and the liquid is the stable phase. At $403.75\ \mathrm{K}$, $\mu_{\mathrm g}<\mu_{\mathrm l}$, and the vapor is the stable phase.

3. **Constant heat capacities**

   Take the liquid from $T$ to $T_{\mathrm b}$, vaporize it there, and bring the vapor back to $T$. With $\Delta\bar c_P=\bar c_{P,\mathrm g}-\bar c_{P,\mathrm l}=137.8-184.6=-46.8\ \mathrm{J\,mol^{-1}\,K^{-1}}$, the heating integrals of [Section 6.1](section-01.md) give

   ```{math}
   \begin{aligned}
   \Delta\bar h_{\mathrm{vap}}(T)&=\Delta\bar h_{\mathrm{vap}}(T_{\mathrm b})+\Delta\bar c_P\,(T-T_{\mathrm b}),\\
   \Delta\bar s_{\mathrm{vap}}(T)&=\Delta\bar s_{\mathrm{vap}}(T_{\mathrm b})+\Delta\bar c_P\ln\frac{T}{T_{\mathrm b}}.
   \end{aligned}
   ```

   At $363.75\ \mathrm{K}$,

   ```{math}
   \begin{aligned}
   \Delta\bar h_{\mathrm{vap}}&=33\,230+(-46.8)(-20)=34\,166\ \mathrm{J\,mol^{-1}},\\
   \Delta\bar s_{\mathrm{vap}}&=86.593+(-46.8)\ln\frac{363.75}{383.75}\\
   &=86.593+2.505=89.098\ \mathrm{J\,mol^{-1}\,K^{-1}},\\
   \Delta\bar g_{\mathrm{vap}}&=34\,166-(363.75)(89.098)=+1757\ \mathrm{J\,mol^{-1}}.
   \end{aligned}
   ```

   The same steps at $403.75\ \mathrm{K}$ complete the table.

   | $T$ (K) | $\Delta\bar h_{\mathrm{vap}}$ ($\mathrm{kJ\,mol^{-1}}$) | $\Delta\bar s_{\mathrm{vap}}$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) | $\Delta\bar g_{\mathrm{vap}}$ ($\mathrm{J\,mol^{-1}}$) | Part (b) ($\mathrm{J\,mol^{-1}}$) |
   |---|---:|---:|---:|---:|
   | 363.75 | 34.17 | 89.10 | +1757 | +1732 |
   | 403.75 | 32.29 | 84.22 | −1708 | −1732 |

   Each $\Delta\bar g_{\mathrm{vap}}$ moves by about $+25\ \mathrm{J\,mol^{-1}}$, which is $1.4\ \%$ of its size. The enthalpy term and the entropy term each change by about $0.9\ \mathrm{kJ\,mol^{-1}}$, and the two changes nearly cancel in $\Delta\bar g_{\mathrm{vap}}$. The signs, and with them the stable phases, are the same as in part (b).

**Result.** $\Delta\bar s_{\mathrm{vap}}=86.59\ \mathrm{J\,mol^{-1}\,K^{-1}}$ at the boiling point. $\Delta\bar g_{\mathrm{vap}}$ is $+1.73\ \mathrm{kJ\,mol^{-1}}$ at $363.75\ \mathrm{K}$ (liquid stable) and $-1.73\ \mathrm{kJ\,mol^{-1}}$ at $403.75\ \mathrm{K}$ (vapor stable). The heat capacities change these values to $+1.76\ \mathrm{kJ\,mol^{-1}}$ and $-1.71\ \mathrm{kJ\,mol^{-1}}$.
````

### Problem 2. Latent heats of ammonia from its phase diagram

Draws on [Section 6.1](section-01.md) and [Section 6.2](section-02.md).

The table lists two points on the liquid–vapor coexistence curve of ammonia, with the tabulated enthalpy of vaporization at each ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/fluid/)).

| Point | $T$ (K) | $P_{\mathrm{sat}}$ (kPa) | $\Delta\bar h_{\mathrm{vap}}$ ($\mathrm{kJ\,mol^{-1}}$) |
|---|---:|---:|---:|
| Triple point | 195.50 | 6.056 | 25.36 |
| Normal boiling point | 239.83 | 101.325 | 23.33 |

Just below the triple point, the enthalpy of sublimation is $\Delta\bar h_{\mathrm{sub}}=31.2\ \mathrm{kJ\,mol^{-1}}$. Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Estimate $\Delta\bar h_{\mathrm{vap}}$ from the two pressures with the integrated Clausius–Clapeyron equation, and compare the estimate with the tabulated values.

**(b)** Find $\Delta\bar h_{\mathrm{fus}}$ at the triple point from $\Delta\bar h_{\mathrm{sub}}$ and the tabulated $\Delta\bar h_{\mathrm{vap}}$.

**(c)** Find the slopes $dP/dT$ of the sublimation curve and of the vaporization curve at the triple point. Treat the vapor as an ideal gas and neglect the molar volumes of the solid and the liquid.

**(d)** Near the triple point, the melting curve of ammonia follows $dT/dP=aT$ with $a=4\times10^{-5}\ \mathrm{atm^{-1}}$ [@haar1977]. Find $dP/dT$ of the melting curve at the triple point in $\mathrm{Pa\,K^{-1}}$, and compare it with the slopes of part (c). Then find $\Delta\bar v_{\mathrm{fus}}$ from the Clapeyron equation and the result of part (b).

````{admonition} Solution
:class: dropdown

1. **Enthalpy of vaporization from two points**

   Let point 1 be the triple point and point 2 the normal boiling point. The integrated Clausius–Clapeyron equation, solved for the enthalpy, is

   ```{math}
   \Delta\bar h_{\mathrm{vap}}=\frac{R\,\ln(P_2/P_1)}{1/T_1-1/T_2}.
   ```

   With the data of the table,

   ```{math}
   \begin{aligned}
   \ln\frac{P_2}{P_1}&=\ln\frac{101.325}{6.056}=2.8173,\\
   \frac{1}{T_1}-\frac{1}{T_2}&=\frac{1}{195.50\ \mathrm{K}}-\frac{1}{239.83\ \mathrm{K}}=9.4547\times10^{-4}\ \mathrm{K^{-1}},\\
   \Delta\bar h_{\mathrm{vap}}&=\frac{(8.3145)(2.8173)}{9.4547\times10^{-4}}=24.78\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

   The estimate lies between the tabulated values at the two ends of the interval, $25.36$ and $23.33\ \mathrm{kJ\,mol^{-1}}$.

2. **Enthalpy of fusion**

   At the triple point the solid, the liquid, and the vapor coexist at one temperature and one pressure. Enthalpy is a state function, so subliming the solid changes the enthalpy by the same amount as melting it and then vaporizing the liquid:

   ```{math}
   \Delta\bar h_{\mathrm{fus}}=\Delta\bar h_{\mathrm{sub}}-\Delta\bar h_{\mathrm{vap}}=31.2-25.36=5.8\ \mathrm{kJ\,mol^{-1}}.
   ```

3. **Slopes of the sublimation and vaporization curves**

   With $\Delta\bar v\approx\bar v_{\mathrm g}=RT/P$, the Clapeyron equation becomes $dP/dT=P\,\Delta\bar h/(RT^2)$ for either curve. At the triple point, $RT^2=(8.3145)(195.50)^2=3.1778\times10^{5}\ \mathrm{J\,K\,mol^{-1}}$ and $P=6056\ \mathrm{Pa}$, so

   ```{math}
   \begin{aligned}
   \left(\frac{dP}{dT}\right)_{\mathrm{sub}}&=\frac{(6056)(31\,200)}{3.1778\times10^{5}}=595\ \mathrm{Pa\,K^{-1}},\\
   \left(\frac{dP}{dT}\right)_{\mathrm{vap}}&=\frac{(6056)(25\,360)}{3.1778\times10^{5}}=483\ \mathrm{Pa\,K^{-1}}.
   \end{aligned}
   ```

   The sublimation curve is steeper than the vaporization curve by the factor $\Delta\bar h_{\mathrm{sub}}/\Delta\bar h_{\mathrm{vap}}=1.23$.

4. **Melting curve**

   The slope is the reciprocal of $dT/dP$:

   ```{math}
   \begin{aligned}
   \frac{dP}{dT}=\frac{1}{aT}&=\frac{1}{(4\times10^{-5}\ \mathrm{atm^{-1}})(195.50\ \mathrm{K})}\\
   &=128\ \mathrm{atm\,K^{-1}}=1.3\times10^{7}\ \mathrm{Pa\,K^{-1}}.
   \end{aligned}
   ```

   The melting curve is about $2\times10^{4}$ times steeper than the two curves of step 3. The Clapeyron equation, $dP/dT=\Delta\bar h_{\mathrm{fus}}/(T\,\Delta\bar v_{\mathrm{fus}})$, then gives

   ```{math}
   \begin{aligned}
   \Delta\bar v_{\mathrm{fus}}&=\frac{\Delta\bar h_{\mathrm{fus}}}{T\,(dP/dT)}=a\,\Delta\bar h_{\mathrm{fus}}\\
   &=\left(\frac{4\times10^{-5}}{101\,325\ \mathrm{Pa}}\right)(5800\ \mathrm{J\,mol^{-1}})\\
   &=2.3\times10^{-6}\ \mathrm{m^3\,mol^{-1}}=2.3\ \mathrm{cm^3\,mol^{-1}}.
   \end{aligned}
   ```

**Result.** $\Delta\bar h_{\mathrm{vap}}\approx24.78\ \mathrm{kJ\,mol^{-1}}$ from the two pressures, and $\Delta\bar h_{\mathrm{fus}}=5.8\ \mathrm{kJ\,mol^{-1}}$. At the triple point the slopes are $595\ \mathrm{Pa\,K^{-1}}$ (sublimation), $483\ \mathrm{Pa\,K^{-1}}$ (vaporization), and $1.3\times10^{7}\ \mathrm{Pa\,K^{-1}}$ (melting), with $\Delta\bar v_{\mathrm{fus}}=2.3\ \mathrm{cm^3\,mol^{-1}}$.
````

### Problem 3. A sublimation curve and the triple-point pressure

Draws on [Section 6.2](section-02.md).

Sulfur hexafluoride, the insulating gas of high-voltage switchgear, sublimes at atmospheric pressure: the solid and the vapor coexist at $209.3\ \mathrm{K}$ and $1.01325\ \mathrm{bar}$. Its enthalpy of sublimation is $\Delta\bar h_{\mathrm{sub}}=23.3\ \mathrm{kJ\,mol^{-1}}$, and its triple point lies at $T_{\mathrm{tp}}=223.56\ \mathrm{K}$ ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/)). Treat the vapor as an ideal gas, neglect the molar volume of the solid, and take $\Delta\bar h_{\mathrm{sub}}$ as constant. Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Start from the Clapeyron equation and show that the sublimation curve obeys the Clausius–Clapeyron equation with $\Delta\bar h_{\mathrm{sub}}$ in place of $\Delta\bar h_{\mathrm{vap}}$. Find the slope $dP/dT$ of the sublimation curve at $1.01325\ \mathrm{bar}$.

**(b)** Estimate the triple-point pressure $P_{\mathrm{tp}}$, and compare it with the measured $2.314\ \mathrm{bar}$.

**(c)** Can liquid sulfur hexafluoride be a stable phase at $1.01325\ \mathrm{bar}$? Give the reason.

````{admonition} Solution
:class: dropdown

1. **Clausius–Clapeyron equation for sublimation**

   The Clapeyron equation holds for any first-order transition, so for solid to vapor

   ```{math}
   \frac{dP}{dT}=\frac{\Delta\bar h_{\mathrm{sub}}}{T\,\Delta\bar v_{\mathrm{sub}}}.
   ```

   The two approximations of the problem give $\Delta\bar v_{\mathrm{sub}}=\bar v_{\mathrm g}-\bar v_{\mathrm s}\approx\bar v_{\mathrm g}=RT/P$. Then

   ```{math}
   \frac{dP}{dT}=\frac{P\,\Delta\bar h_{\mathrm{sub}}}{RT^2}
   \qquad\Longrightarrow\qquad
   \frac{d\ln P}{dT}=\frac{\Delta\bar h_{\mathrm{sub}}}{RT^2}.
   ```

2. **Slope at atmospheric pressure**

   ```{math}
   \frac{dP}{dT}=\frac{(101\,325\ \mathrm{Pa})(23\,300\ \mathrm{J\,mol^{-1}})}{(8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}})(209.3\ \mathrm{K})^2}=6.48\times10^{3}\ \mathrm{Pa\,K^{-1}}.
   ```

   The sublimation pressure rises by about $0.065\ \mathrm{bar}$ for each kelvin.

3. **Triple-point pressure**

   The triple point lies on the sublimation curve. Integrate from $T_1=209.3\ \mathrm{K}$, $P_1=1.01325\ \mathrm{bar}$ to $T_{\mathrm{tp}}$ with constant $\Delta\bar h_{\mathrm{sub}}$:

   ```{math}
   \begin{aligned}
   \ln\frac{P_{\mathrm{tp}}}{P_1}&=\frac{\Delta\bar h_{\mathrm{sub}}}{R}\left(\frac{1}{T_1}-\frac{1}{T_{\mathrm{tp}}}\right)\\
   &=(2802.3\ \mathrm{K})\left(4.7778\times10^{-3}-4.4731\times10^{-3}\right)\mathrm{K^{-1}}\\
   &=0.8540,\\
   P_{\mathrm{tp}}&=(1.01325\ \mathrm{bar})\,e^{0.8540}=2.38\ \mathrm{bar}.
   \end{aligned}
   ```

   The estimate is $2.9\ \%$ above the measured $2.314\ \mathrm{bar}$.

4. **Liquid at atmospheric pressure**

   No. On the phase diagram the liquid region begins at the triple point, so the liquid is a stable phase only at pressures above $P_{\mathrm{tp}}=2.314\ \mathrm{bar}$. Atmospheric pressure lies below $P_{\mathrm{tp}}$. A path of rising temperature at $1.01325\ \mathrm{bar}$ crosses only the sublimation curve, at $209.3\ \mathrm{K}$, and the solid passes directly into the vapor.

**Result.** $dP/dT=6.48\ \mathrm{kPa\,K^{-1}}$ at $1.01325\ \mathrm{bar}$, and $P_{\mathrm{tp}}\approx2.38\ \mathrm{bar}$ against the measured $2.314\ \mathrm{bar}$. Liquid sulfur hexafluoride is not stable at atmospheric pressure.
````

### Problem 4. Vapor pressure of propane from an equation

Draws on [Section 6.2](section-02.md).

Between $231\ \mathrm{K}$ and $321\ \mathrm{K}$, the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/) represents the vapor pressure of propane by

```{math}
\log_{10}\left(\frac{P_{\mathrm{sat}}}{\mathrm{bar}}\right)=A-\frac{B}{T+C},
```

with $A=3.98292$, $B=819.296\ \mathrm{K}$, and $C=-24.417\ \mathrm{K}$. Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}=0.083145\ \mathrm{L\,bar\,mol^{-1}\,K^{-1}}$.

**(a)** Find the normal boiling point, where $P_{\mathrm{sat}}=1.01325\ \mathrm{bar}$, from the equation. Compare it with the measured $231.04\ \mathrm{K}$.

**(b)** Differentiate the equation to find $d\ln P_{\mathrm{sat}}/dT$. Use the Clausius–Clapeyron equation to estimate $\Delta\bar h_{\mathrm{vap}}$ at $231.04\ \mathrm{K}$, and compare the estimate with the calorimetric $18.77\ \mathrm{kJ\,mol^{-1}}$.

**(c)** At the normal boiling point, $\bar v_{\mathrm g}=18.25\ \mathrm{L\,mol^{-1}}$ and $\bar v_{\mathrm l}=0.076\ \mathrm{L\,mol^{-1}}$. Find $\Delta Z$, and repeat the estimate of part (b) with the exact relation between the slope and $\Delta\bar h_{\mathrm{vap}}$.

**(d)** A cylinder that holds liquid propane and its vapor stands in the sun at $318.15\ \mathrm{K}$ ($45\ ^\circ\mathrm{C}$). Find the pressure in the cylinder, and compare it with $15.34\ \mathrm{bar}$ from the NIST Chemistry WebBook.

````{admonition} Solution
:class: dropdown

1. **Normal boiling point**

   Solve the equation for the temperature. With $\log_{10}(1.01325)=0.00572$,

   ```{math}
   \begin{aligned}
   T&=\frac{B}{A-\log_{10}(P_{\mathrm{sat}}/\mathrm{bar})}-C\\
   &=\frac{819.296\ \mathrm{K}}{3.98292-0.00572}+24.417\ \mathrm{K}\\
   &=206.00\ \mathrm{K}+24.42\ \mathrm{K}=230.41\ \mathrm{K}.
   \end{aligned}
   ```

   The equation places the boiling point $0.63\ \mathrm{K}$ below the measured $231.04\ \mathrm{K}$.

2. **Slope and the Clausius–Clapeyron estimate**

   Because $\ln P_{\mathrm{sat}}=\ln(10)\,\log_{10}P_{\mathrm{sat}}$,

   ```{math}
   \frac{d\ln P_{\mathrm{sat}}}{dT}=\frac{\ln(10)\,B}{(T+C)^2}.
   ```

   The Clausius–Clapeyron equation, $d\ln P_{\mathrm{sat}}/dT=\Delta\bar h_{\mathrm{vap}}/(RT^2)$, then gives

   ```{math}
   \begin{aligned}
   \Delta\bar h_{\mathrm{vap}}&=RT^2\,\frac{d\ln P_{\mathrm{sat}}}{dT}=R\,\ln(10)\,B\left(\frac{T}{T+C}\right)^2\\
   &=(8.3145)(2.3026)(819.296)\left(\frac{231.04}{206.623}\right)^2\ \mathrm{J\,mol^{-1}}\\
   &=(15\,685)(1.2503)\ \mathrm{J\,mol^{-1}}=19.61\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

   The estimate is $4.5\ \%$ above the calorimetric $18.77\ \mathrm{kJ\,mol^{-1}}$.

3. **Estimate with $\Delta Z$**

   ```{math}
   \begin{aligned}
   \Delta Z&=\frac{P_{\mathrm{sat}}\,\Delta\bar v_{\mathrm{vap}}}{RT}\\
   &=\frac{(1.01325\ \mathrm{bar})(18.25-0.076)\ \mathrm{L\,mol^{-1}}}{(0.083145\ \mathrm{L\,bar\,mol^{-1}\,K^{-1}})(231.04\ \mathrm{K})}\\
   &=\frac{18.415}{19.210}=0.9586.
   \end{aligned}
   ```

   The exact relation is $RT^2\,d\ln P_{\mathrm{sat}}/dT=\Delta\bar h_{\mathrm{vap}}/\Delta Z$, so

   ```{math}
   \Delta\bar h_{\mathrm{vap}}=\Delta Z\;RT^2\,\frac{d\ln P_{\mathrm{sat}}}{dT}=(0.9586)(19.61)=18.80\ \mathrm{kJ\,mol^{-1}}.
   ```

   This value is within $0.2\ \%$ of the calorimetric one.

4. **Pressure in the cylinder**

   While liquid and vapor are both present, the pressure is the vapor pressure at the temperature of the cylinder. At $T=318.15\ \mathrm{K}$, $T+C=293.733\ \mathrm{K}$ and

   ```{math}
   \begin{aligned}
   \log_{10}\left(\frac{P_{\mathrm{sat}}}{\mathrm{bar}}\right)&=3.98292-\frac{819.296}{293.733}=1.19367,\\
   P_{\mathrm{sat}}&=10^{1.19367}\ \mathrm{bar}=15.62\ \mathrm{bar}.
   \end{aligned}
   ```

   The result is $1.8\ \%$ above $15.34\ \mathrm{bar}$.

**Result.** The equation gives a normal boiling point of $230.41\ \mathrm{K}$ and $P_{\mathrm{sat}}=15.62\ \mathrm{bar}$ at $318.15\ \mathrm{K}$. At the boiling point, $\Delta\bar h_{\mathrm{vap}}=19.61\ \mathrm{kJ\,mol^{-1}}$ from the Clausius–Clapeyron equation and $18.80\ \mathrm{kJ\,mol^{-1}}$ with $\Delta Z=0.9586$, against the calorimetric $18.77\ \mathrm{kJ\,mol^{-1}}$.
````

## 4. Further Reading

**Textbooks**

- @mcquarrie1997, Chapter 23, Sections 23.1 to 23.5 (phase equilibria).
- @dill2011, Chapter 14 (equilibria between liquids, solids, and gases) and Chapter 25 (phase transitions).

**Original papers**

- @clapeyron1834 relates the slope of the vapor-pressure curve to the heat of vaporization, the first form of the Clapeyron equation of [Section 6.2](section-02.md). The English translation by R. Taylor is in *Scientific Memoirs* **1**, 347–376 (1837).
- @clausius1850 rederives Clapeyron's relation from the equivalence of heat and work ([Section 6.2](section-02.md)). An English translation appeared in *Philosophical Magazine* **2**, 1–21 and 102–119 (1851).
- @trouton1884 observes that the molar heat of vaporization divided by the boiling temperature is nearly the same for many liquids, which is Trouton's rule of [Section 6.1](section-01.md).
