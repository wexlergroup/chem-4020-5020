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

# 5.3. Ammonia Formation

## Overview

In this section we use the **Haber–Bosch ammonia formation reaction** as a case study for how
thermodynamics and kinetics both matter in chemical engineering:

```{math}
\mathrm{N_2(g) + 3H_2(g) \rightarrow 2NH_3(g)}
```

- At room temperature the reaction is **thermodynamically favorable** ($\Delta G_r^\circ < 0$),
  but it is **kinetically slow** because the $\mathrm{N\equiv N}$ bond is very strong.
- Raising the temperature helps overcome the **activation barrier** (faster kinetics), but it can also
  change the **reaction spontaneity** because $G = H - TS$.
- For gas-phase reactions, we can also “turn the pressure knob.” Pressure can strongly affect $\Delta G$
  when the reaction changes the **number of gas molecules**.

The goal is to understand (and quantify) these tradeoffs.

---

## Thermodynamics vs. kinetics

### Standard thermodynamic data at 298.15 K

From the lecture notes, at $T = 298.15\ \mathrm{K}$ ($25^\circ\mathrm{C}$):

```{math}
\Delta H_r^\circ(298.15\,\mathrm{K}) \approx -92\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_r^\circ(298.15\,\mathrm{K}) \approx -198\ \mathrm{J\,(K\,mol)^{-1}}
```

```{math}
\Delta G_r^\circ(298.15\,\mathrm{K}) \approx -33\ \mathrm{kJ\,mol^{-1}}
```

So, **thermodynamically**, ammonia formation is spontaneous at 298 K because $\Delta G_r^\circ < 0$.

:::{admonition} Data skills: where do these $\Delta H_r^\circ$, $\Delta S_r^\circ$, and $\Delta G_r^\circ$ numbers come from?
:class: dropdown

In practice, standard reaction values are usually computed from **tabulated species data**:

1. Pull $\Delta_f H^\circ(298.15\ \mathrm{K})$ and $S^\circ(298.15\ \mathrm{K})$ for each species from a data source (NIST WebBook, JANAF, ATcT).
2. Compute $\Delta_r H^\circ$ and $\Delta_r S^\circ$ using stoichiometry (products minus reactants).
3. Compute $\Delta_r G^\circ = \Delta_r H^\circ - T\Delta_r S^\circ$.
4. Convert $\Delta_r G^\circ$ to an equilibrium constant via $K=\exp(-\Delta_r G^\circ/RT)$.

A step-by-step workflow (including a worked ammonia example) is given in **Section 7.1.7** and the main data sources are collected in **Appendix Y**.
:::

:::{admonition} Key idea
:class: tip
Thermodynamics ($\Delta G$) tells us *which direction is favored at equilibrium*.
Kinetics (activation barrier / $E_a$) tells us *how fast the system gets there*.
:::

### Why it is still “too slow” at room temperature

Even if $\Delta G_r^\circ < 0$, the reaction can be extremely slow if the activation barrier is large.
For ammonia formation, the strong $\mathrm{N\equiv N}$ bond leads to a very high barrier for the uncatalyzed pathway.

A catalyst (e.g., an Fe-based catalyst) provides a different reaction pathway with a lower barrier.
Raising temperature further increases rates (Arrhenius behavior), so at $500^\circ\mathrm{C}$ the barrier can be
“surmountable” (fast kinetics).

---

## Temperature dependence of $G$

We want to know: **if we increase temperature from 25°C to 500°C, do we change spontaneity?**

A practical strategy is to compute $G^\circ(T_f) - G^\circ(T_i)$ in two pieces:

1. **Enthalpy correction**: $H(T_f) - H(T_i)$
2. **Entropy correction**: $S(T_f) - S(T_i)$

since $G(T)=H(T)-TS(T)$.

### Enthalpy correction from $C_P$

At constant pressure,

```{math}
\left(\frac{\partial H}{\partial T}\right)_P = C_P
\qquad\Rightarrow\qquad
dH = C_P\,dT
```

so

```{math}
H(T_f) - H(T_i) = \int_{T_i}^{T_f} C_P(T)\,dT
```

### Entropy correction from $C_P$

Using the equilibrium identity (constant composition)

```{math}
dH = T\,dS + V\,dP
```

At constant pressure ($dP=0$),

```{math}
dS = \frac{1}{T}\,dH = \frac{C_P}{T}\,dT
```

so

```{math}
S(T_f) - S(T_i) = \int_{T_i}^{T_f} \frac{C_P(T)}{T}\,dT
```

---

## Ideal-gas heat capacities from degrees of freedom

To get a simple estimate, the notes assume $\mathrm{N_2}$, $\mathrm{H_2}$, and $\mathrm{NH_3}$ behave as ideal gases
and use a constant-$C_P$ model based on the number of (active) degrees of freedom $f$:

```{math}
C_P = \left(\frac{f}{2} + 1\right)R
```

- For linear molecules ($\mathrm{N_2}$, $\mathrm{H_2}$):
  $f = 3$ translational $+\ 2$ rotational $=5$, so

  ```{math}
  C_P = \frac{7}{2}R
  ```

- For non-linear $\mathrm{NH_3}$:
  $f = 3$ translational $+\ 3$ rotational $=6$, so

  ```{math}
  C_P = 4R
  ```

With constant $C_P$, the integrals become:

```{math}
H(T_f)-H(T_i) = C_P(T_f-T_i)
```

```{math}
S(T_f)-S(T_i) = C_P\ln\!\left(\frac{T_f}{T_i}\right)
```

---

## Worked example: temperature effect from 298.15 K to 773.15 K

Let

- $T_i = 298.15\ \mathrm{K}$ (25°C)
- $T_f = 773.15\ \mathrm{K}$ (500°C)

### Species enthalpy/entropy shifts (ideal-gas approximation)

**For $\mathrm{N_2}$ and $\mathrm{H_2}$** ($C_P = \tfrac{7}{2}R$):

```{math}
H(T_f) \approx H(T_i) + \frac{7}{2}R(T_f-T_i) \approx H(T_i) + 13.8\ \mathrm{kJ\,mol^{-1}}
```

```{math}
S(T_f) \approx S(T_i) + \frac{7}{2}R\ln\!\left(\frac{T_f}{T_i}\right)
\approx S(T_i) + 27.7\ \mathrm{J\,(K\,mol)^{-1}}
```

**For $\mathrm{NH_3}$** ($C_P = 4R$):

```{math}
H(T_f) \approx H(T_i) + 4R(T_f-T_i) \approx H(T_i) + 15.8\ \mathrm{kJ\,mol^{-1}}
```

```{math}
S(T_f) \approx S(T_i) + 4R\ln\!\left(\frac{T_f}{T_i}\right)
\approx S(T_i) + 31.7\ \mathrm{J\,(K\,mol)^{-1}}
```

### Reaction $\Delta G_r^\circ$ at 773.15 K

Start from

```{math}
\Delta G_r^\circ(T_f)=
2G_{\mathrm{NH_3}}^\circ(T_f)-G_{\mathrm{N_2}}^\circ(T_f)-3G_{\mathrm{H_2}}^\circ(T_f)
```

and use $G=H-TS$:

```{math}
\Delta G_r^\circ(T_f)=\Delta H_r^\circ(T_f)-T_f\,\Delta S_r^\circ(T_f)
```

Compute the *corrections* from $T_i\to T_f$:

- Enthalpy correction:

  ```{math}
  \Delta H_r^\circ(T_f)
  \approx \Delta H_r^\circ(T_i) +
  \Big[2(15.8)-\big(13.8\big)-3\big(13.8\big)\Big]\ \mathrm{kJ\,mol^{-1}}
  = \Delta H_r^\circ(T_i) - 23.6\ \mathrm{kJ\,mol^{-1}}
  ```

- Entropy correction:

  ```{math}
  \Delta S_r^\circ(T_f)
  \approx \Delta S_r^\circ(T_i) +
  \Big[2(31.7)-\big(27.7\big)-3\big(27.7\big)\Big]\ \mathrm{J\,(K\,mol)^{-1}}
  = \Delta S_r^\circ(T_i) - 47.4\ \mathrm{J\,(K\,mol)^{-1}}
  ```

Equivalently (as written in the notes),

```{math}
\Delta G_r^\circ(T_f)
\approx \Delta H_r^\circ(T_i)-T_f\,\Delta S_r^\circ(T_i) - 23.6 + (0.0474)T_f
```

(with $T_f$ in K and energy in kJ/mol).

Using $\Delta H_r^\circ(298.15\,\mathrm{K})\approx-92\ \mathrm{kJ\,mol^{-1}}$ and
$\Delta S_r^\circ(298.15\,\mathrm{K})\approx-0.198\ \mathrm{kJ\,(K\,mol)^{-1}}$:

```{math}
\Delta G_r^\circ(773.15\,\mathrm{K}) \approx 74.1\ \mathrm{kJ\,mol^{-1}}
```

So at $500^\circ\mathrm{C}$ (and near 1 bar), ammonia formation becomes **non-spontaneous** ($\Delta G_r^\circ > 0$),
even though it would be much faster kinetically.

---

## Pressure dependence of $G$ for an ideal gas

Temperature isn’t the only knob; for gases, pressure matters.

From $dG = -S\,dT + V\,dP$, at **constant temperature**:

```{math}
dG = V\,dP
```

For an ideal gas, $V = \dfrac{nRT}{P}$. Integrate from $P_i\to P_f$:

```{math}
G(P_f)-G(P_i)=\int_{P_i}^{P_f} \frac{nRT}{P}\,dP
= nRT\ln\!\left(\frac{P_f}{P_i}\right)
```

### Apply to ammonia formation

At fixed $T=T_f$,

```{math}
\Delta G_r(P_f,T_f) \approx \Delta G_r^\circ(T_f) + \Delta\nu\,RT_f\ln\!\left(\frac{P_f}{P_i}\right)
```

where $\Delta\nu$ is the change in gas moles:

```{math}
\Delta\nu = (\text{moles products})-(\text{moles reactants})
=2-(1+3)=-2
```

At the pressure where the reaction becomes just spontaneous, $\Delta G_r(P_f,T_f)=0$:

```{math}
0 = \Delta G_r^\circ(T_f) + \Delta\nu\,RT_f\ln\!\left(\frac{P_f}{P_i}\right)
```

Solve for $P_f$:

```{math}
P_f = P_i\exp\!\left(-\frac{\Delta G_r^\circ(T_f)}{\Delta\nu\,RT_f}\right)
```

Using $\Delta G_r^\circ(773.15\,\mathrm{K})\approx 74.1\ \mathrm{kJ\,mol^{-1}}$, $\Delta\nu=-2$,
and taking $P_i=1\ \mathrm{bar}$:

```{math}
P_f \approx 3.19\times 10^2\ \mathrm{bar} \;\approx\; 319\ \mathrm{bar}
```

**Interpretation:** because $\Delta\nu<0$, increasing pressure favors the side with *fewer* gas molecules
(Le Châtelier’s principle), helping make $\Delta G$ negative again.

---

## Summary of the “engineering knobs”

- **Catalyst:** lowers the activation barrier $\Rightarrow$ faster kinetics (does *not* change $\Delta G$).
- **Temperature:** higher $T$ $\Rightarrow$ faster kinetics, but can push $\Delta G$ in an unfavorable direction
  when $\Delta S_r^\circ < 0$.
- **Pressure:** for $\Delta\nu<0$ gas reactions, higher pressure can restore thermodynamic favorability.

This is the core thermodynamics/kinetics tradeoff that motivates the operating conditions for the Haber–Bosch process.

---

## (Preview) Phase equilibria

**Phase equilibria** applies the thermodynamic equilibrium condition to decide what phase (solid/liquid/gas)
a system exists in.

A typical single-component $P$-$T$ phase diagram contains:

- **Coexistence curves** (lines): solid–liquid, liquid–gas, solid–gas
- **Single-phase regions** (areas)
- **Triple point (tp):** where all three phases coexist
- **Critical point (cp):** end of the liquid–gas coexistence curve

The lecture also previewed that phase diagrams can describe *solid–solid* transitions (e.g., different solid
forms of Sn) and even motivate historical anecdotes (“phase diagrams & the downfall of Napoleon”).
