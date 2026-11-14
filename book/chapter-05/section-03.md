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

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Ammonia synthesis is one of the most important industrial chemical processes: the Haber–Bosch reaction produces the nitrogen feedstock for most synthetic fertilizers, supporting roughly half the world's food production. It is also a textbook example of the tension between thermodynamics and kinetics. At room temperature the reaction is thermodynamically favorable ($\Delta_r G^\circ < 0$) but kinetically inert. Raising the temperature makes the kinetics manageable but can push $\Delta_r G^\circ$ positive. Increasing pressure can restore favorability for reactions that decrease the number of gas molecules. Understanding and quantifying these tradeoffs is the goal of this section.

We use the free-energy machinery from [Section 5.1](section-01.md) and the absolute-entropy framework from [Section 5.2](section-02.md) to track how $\Delta_r G^\circ$ changes with temperature and pressure for the reaction

```{math}
\mathrm{N_2(g) + 3\,H_2(g) \rightarrow 2\,NH_3(g)}.
```

The analysis relies on two tools: heat-capacity-based temperature corrections to $\Delta H$ and $\Delta S$, and the ideal-gas pressure dependence of $G$.

---

Learning objectives:

- Articulate the difference between thermodynamic favorability ($\Delta G$) and kinetic accessibility (activation barrier, catalysis).
- Use $\Delta G = \Delta H - T\Delta S$ to reason qualitatively about temperature effects when $\Delta_r S^\circ < 0$.
- Estimate temperature corrections to $\Delta_r H^\circ$ and $\Delta_r S^\circ$ from constant-$C_P$ models based on active degrees of freedom.
- Derive the ideal-gas pressure dependence of the molar Gibbs energy, $G(T,P) - G^\circ(T) = RT\ln(P/P^\circ)$, and apply it to reactions with $\Delta\nu \neq 0$.
- Explain how catalyst, temperature, and pressure serve as complementary "engineering knobs" in the Haber–Bosch process.

## Core Ideas and Derivations

### Thermodynamics vs. Kinetics

#### Standard Thermodynamic Data at 298.15 K

The standard reaction quantities for ammonia synthesis at $T = 298.15\ \mathrm{K}$ are:

```{math}
\Delta_r H^\circ(298.15\,\mathrm{K}) \approx -92\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta_r S^\circ(298.15\,\mathrm{K}) \approx -198\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

```{math}
\Delta_r G^\circ(298.15\,\mathrm{K}) \approx -33\ \mathrm{kJ\,mol^{-1}}
```

These are consistent via $\Delta_r G^\circ = \Delta_r H^\circ - T\,\Delta_r S^\circ$. Since $\Delta_r G^\circ < 0$, the reaction is **thermodynamically spontaneous** at 298 K and 1 bar.

:::{admonition} Data Skills: Where Do These Numbers Come From?
:class: dropdown

Standard reaction quantities are computed from **tabulated species data**:

1. Look up $\Delta_f H^\circ(298.15\ \mathrm{K})$ and $S^\circ(298.15\ \mathrm{K})$ for each species from a reference source (NIST WebBook, NIST-JANAF, ATcT).
2. Compute $\Delta_r H^\circ$ and $\Delta_r S^\circ$ by stoichiometry (products minus reactants). The absolute entropies $S^\circ$ used here are third-law values of the type discussed in [Section 5.2](section-02.md).
3. Combine via $\Delta_r G^\circ = \Delta_r H^\circ - T\,\Delta_r S^\circ$.

<!-- A step-by-step workflow (including a worked ammonia example) is given in Section 7.1.7 and the main data sources are collected in Appendix Y. -->
:::

:::{admonition} Key Idea
:class: tip
**Thermodynamics** ($\Delta G$) tells us *which direction is favored at equilibrium*.
**Kinetics** (activation barrier, $E_a$) tells us *how fast the system gets there*.
These are independent questions: a reaction can be thermodynamically favorable yet kinetically inert, or vice versa.
:::

#### Why It Is Too Slow at Room Temperature

Even with $\Delta_r G^\circ < 0$, the reaction can be extremely slow if the activation barrier is large. For ammonia formation, the $\mathrm{N \equiv N}$ triple bond (bond dissociation energy $\approx 945\ \mathrm{kJ\,mol^{-1}}$) makes the uncatalyzed barrier prohibitively high at 298 K.

A **catalyst** provides an alternative reaction pathway with a lower activation barrier. Fritz Haber demonstrated the synthesis in 1909 with an osmium catalyst; the iron-based catalyst of the industrial process was found by Alwin Mittasch at BASF, where Carl Bosch scaled the process up (Nobel Prizes: Haber 1918, Bosch 1931). The catalyst accelerates both the forward and reverse reactions equally, so it **does not change $\Delta_r G^\circ$ or the equilibrium constant** — it only changes how fast equilibrium is reached. Raising the temperature further increases the rate (Arrhenius behavior), so industrial reactors operate at $400\text{–}500\,\text{°C}$ where the barrier is surmountable.

But raising $T$ introduces a thermodynamic penalty, which we now quantify.

---

### Qualitative Effect of Temperature

Since $\Delta_r S^\circ < 0$ for this reaction, the $-T\,\Delta_r S^\circ$ term in

```{math}
\Delta_r G^\circ = \Delta_r H^\circ - T\,\Delta_r S^\circ
```

is **positive** and grows with $T$. At sufficiently high temperature, this positive term overwhelms the negative $\Delta_r H^\circ$, and $\Delta_r G^\circ$ changes sign — the reaction becomes non-spontaneous at 1 bar.

The qualitative rule is: for reactions with $\Delta_r S^\circ < 0$, higher temperatures work *against* thermodynamic favorability.

---

### Temperature Corrections from $C_P$

To quantify the temperature effect, we need $\Delta_r H^\circ(T)$ and $\Delta_r S^\circ(T)$ at the operating temperature. The strategy is to correct each species' enthalpy and entropy from the reference temperature $T_i = 298.15\ \mathrm{K}$ to the target temperature $T_f$ using heat-capacity data.

#### Enthalpy Correction

At constant pressure ([Section 3.3](../chapter-03/section-03.md)):

```{math}
H(T_f) - H(T_i) = \int_{T_i}^{T_f} C_P(T)\,dT
```

#### Entropy Correction

From the enthalpy differential $dH = T\,dS + V\,dP$ at constant $P$ ($dP = 0$):

```{math}
S(T_f) - S(T_i) = \int_{T_i}^{T_f} \frac{C_P(T)}{T}\,dT
```

This is the same integral that appears in the absolute-entropy formula ([Section 5.2](section-02.md), Eq. {eq}`eq:absolute-entropy`), applied over a finite temperature range rather than from $0$.

---

### Ideal-Gas Heat Capacities from Degrees of Freedom

For a simple estimate we treat $\mathrm{N_2}$, $\mathrm{H_2}$, and $\mathrm{NH_3}$ as ideal gases and approximate $C_P$ as temperature-independent, using the number of active (translational + rotational) degrees of freedom $f$:

```{math}
:label: eq:Cp-dof
C_P = \left(\frac{f}{2} + 1\right)R
```

For **linear molecules** ($\mathrm{N_2}$, $\mathrm{H_2}$): $f = 3\ (\mathrm{trans}) + 2\ (\mathrm{rot}) = 5$, giving

```{math}
C_P = \tfrac{7}{2}R = 29.1\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

For **nonlinear** $\mathrm{NH_3}$: $f = 3\ (\mathrm{trans}) + 3\ (\mathrm{rot}) = 6$, giving

```{math}
C_P = 4R = 33.3\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

:::{admonition} Approximations in This Model
:class: note

Equation {eq}`eq:Cp-dof` counts only translational and rotational degrees of freedom — it is the "rigid rotor" approximation. Vibrational modes are neglected. At 773 K, some vibrational modes (particularly the lower-frequency bending modes of $\mathrm{NH_3}$) are partially excited, so the true $C_P$ values are somewhat larger than these estimates. This is the same "freeze-out" phenomenon discussed in [Section 2.6](../chapter-02/section-06.md): vibrational contributions to $C_P$ become significant only when $T$ is comparable to or larger than the vibrational temperature $\Theta_{\mathrm{vib}}$.

Despite this limitation, the constant-$C_P$ model captures the correct qualitative behavior. Its errors in $\Delta_r H^\circ$ and in $T\,\Delta_r S^\circ$ largely cancel, so that $\Delta_r G^\circ(T)$ stays within about $3\ \mathrm{kJ\,mol^{-1}}$ of the value computed from tabulated heat capacities up to $500\,\text{°C}$ (see the figure in the Summary below).
:::

With constant $C_P$, the integrals simplify to:

```{math}
H(T_f) - H(T_i) = C_P(T_f - T_i)
```

```{math}
S(T_f) - S(T_i) = C_P\ln\!\left(\frac{T_f}{T_i}\right)
```

---

### Temperature Effect: 298.15 K to 773.15 K (500°C)

Set $T_i = 298.15\ \mathrm{K}$ and $T_f = 773.15\ \mathrm{K}$.

#### Per-Species Corrections

**$\mathrm{N_2}$ and $\mathrm{H_2}$** ($C_P = \tfrac{7}{2}R$):

```{math}
\Delta H_{\mathrm{species}} = \tfrac{7}{2}R\,(773.15 - 298.15) \approx 13.82\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_{\mathrm{species}} = \tfrac{7}{2}R\,\ln\!\left(\frac{773.15}{298.15}\right) \approx 27.73\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

**$\mathrm{NH_3}$** ($C_P = 4R$):

```{math}
\Delta H_{\mathrm{species}} = 4R\,(773.15 - 298.15) \approx 15.80\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_{\mathrm{species}} = 4R\,\ln\!\left(\frac{773.15}{298.15}\right) \approx 31.69\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

#### Reaction-Level Corrections

Apply stoichiometric weights ($2 \times \mathrm{NH_3} - 1 \times \mathrm{N_2} - 3 \times \mathrm{H_2}$):

```{math}
\begin{aligned}
\Delta_r H^\circ(T_f) &\approx \Delta_r H^\circ(T_i) + \big[2(15.80) - (13.82) - 3(13.82)\big] \\
&= \Delta_r H^\circ(T_i) - 23.7\ \mathrm{kJ\,mol^{-1}}
\end{aligned}
```

```{math}
\begin{aligned}
\Delta_r S^\circ(T_f) &\approx \Delta_r S^\circ(T_i) + \big[2(31.69) - (27.73) - 3(27.73)\big] \\
&= \Delta_r S^\circ(T_i) - 47.5\ \mathrm{J\,K^{-1}\,mol^{-1}}
\end{aligned}
```

The reaction-level $\Delta_r C_P$ is negative (products have fewer degrees of freedom than reactants), so both $\Delta_r H^\circ$ and $\Delta_r S^\circ$ become more negative with increasing temperature.

#### $\Delta_r G^\circ$ at 773.15 K

Combining the corrected values:

```{math}
\Delta_r H^\circ(773) = -92.0 + (-23.7) = -115.7\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta_r S^\circ(773) = -198 + (-47.5) = -245.5\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

```{math}
:label: eq:dG-773
\Delta_r G^\circ(773) = -115.7 - (773.15)(-0.2455) \approx +74.1\ \mathrm{kJ\,mol^{-1}}
```

At $500\,\text{°C}$ and 1 bar, ammonia formation is **non-spontaneous** ($\Delta_r G^\circ > 0$), even though it would be much faster kinetically than at room temperature. This is the central dilemma of the Haber–Bosch process.

---

### Pressure Dependence of $G$ for an Ideal Gas

Temperature is not the only knob available. For gas-phase reactions, **pressure** can strongly affect $\Delta G$ when the reaction changes the number of gas molecules.

From the Gibbs differential $dG = -S\,dT + V\,dP$ ([Section 5.1](section-01.md)) at constant temperature:

```{math}
dG\big|_T = V\,dP
```

For one mole of an ideal gas, $V = RT/P$. Integrating from $P^\circ$ to $P$:

```{math}
:label: eq:G-pressure
G(T, P) - G^\circ(T) = RT\ln\!\left(\frac{P}{P^\circ}\right)
```

where $G^\circ(T)$ is the standard-state molar Gibbs energy at pressure $P^\circ = 1\ \mathrm{bar}$.

:::{admonition} From One Gas to a Reacting Mixture
:class: note
Equation {eq}`eq:G-pressure` describes a pure gas. For a reacting *mixture*, the relevant quantity is the reaction quotient $Q_p$ (a ratio of partial pressures weighted by stoichiometric coefficients), and the general expression is $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$. The simplified form used below, with each species at the same pressure $P$, captures the correct pressure *scale* for restoring spontaneity but does not describe the equilibrium composition of a mixture. [Section 7.1](../chapter-07/section-01.md) develops the general treatment.
:::

#### Applying to Ammonia Formation

If each species is compressed from $P^\circ$ to $P$ at fixed $T$, the pressure correction to the reaction Gibbs energy is

```{math}
\Delta_r G(T, P) \approx \Delta_r G^\circ(T) + \Delta\nu\,RT\ln\!\left(\frac{P}{P^\circ}\right),
```

where $\Delta\nu$ is the change in the number of moles of gas:

```{math}
\Delta\nu = 2 - (1 + 3) = -2.
```

Setting $\Delta_r G = 0$ to find the pressure at which the reaction is just spontaneous:

```{math}
:label: eq:P-crossover
P = P^\circ \exp\!\left(-\frac{\Delta_r G^\circ(T)}{\Delta\nu\,RT}\right)
```

Inserting $\Delta_r G^\circ(773) \approx +74.1\ \mathrm{kJ\,mol^{-1}}$, $\Delta\nu = -2$, and $T = 773.15\ \mathrm{K}$:

```{math}
\ln\!\left(\frac{P}{P^\circ}\right)
= -\frac{74.1 \times 10^3}{(-2)(8.314)(773.15)} = 5.76
\quad\Rightarrow\quad
P \approx 320\ \mathrm{bar}
```

Because $\Delta\nu < 0$, increasing pressure favors the side with fewer gas molecules (Le Châtelier's principle). At $500\,\text{°C}$, pressures on the order of several hundred bar can restore thermodynamic favorability.

---

### Summary: The Engineering Knobs

The figure below collects the two results of this section: the standard reaction Gibbs energy at 1 bar as a function of temperature, and the pressure at which the reaction Gibbs energy vanishes.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import R

# Constant-C_P model of this section (rigid molecules, no vibration)
T_i = 298.15                      # K
dH_i, dS_i = -92.0e3, -198.0      # J/mol and J/(mol K) at 298.15 K
dCp = -6 * R                      # 2(4R) - (7R/2) - 3(7R/2)

def dG_model(T):
    dH = dH_i + dCp * (T - T_i)
    dS = dS_i + dCp * np.log(T / T_i)
    return dH - T * dS

# Reference: Shomate equations from the NIST Chemistry WebBook (Chase, 1998).
# Coefficients A-H; each tuple is (T_min, T_max, coefficients).
shomate = {
    "NH3": [(298.0, 1400.0, (19.99563, 49.77119, -15.37599, 1.921168, 0.189174, -53.30667, 203.8591, -45.89806))],
    "N2": [(100.0, 500.0, (28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0.0)),
           (500.0, 2000.0, (19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0))],
    "H2": [(298.0, 1000.0, (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0))],
}
dfH = {"NH3": -45.90e3, "N2": 0.0, "H2": 0.0}   # J/mol at 298.15 K
nu = {"NH3": 2, "N2": -1, "H2": -3}

def H_and_S(species, T):
    """Molar enthalpy (J/mol, relative to the elements at 298.15 K) and entropy (J/(mol K))."""
    for T_min, T_max, (A, B, C, D, E, F, G, H) in shomate[species]:
        if T_min <= T <= T_max:
            t = T / 1000
            dH = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H      # kJ/mol
            S = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
            return dfH[species] + 1000 * dH, S

def dG_nist(T):
    dH = sum(nu[s] * H_and_S(s, T)[0] for s in nu)
    dS = sum(nu[s] * H_and_S(s, T)[1] for s in nu)
    return dH - T * dS

T = np.linspace(298.15, 873.15, 231)
theta = T - 273.15                               # temperature in degrees Celsius
g_model = dG_model(T) / 1000                     # kJ/mol
g_nist = np.array([dG_nist(Tk) for Tk in T]) / 1000
T_cross = T[np.argmin(np.abs(g_model))]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.6), constrained_layout=True)
color = "tab:blue"

# (a) Standard reaction Gibbs energy at 1 bar
ax1.axvspan(400, 500, color="0.88", lw=0)
ax1.axhline(0, color="0.5", lw=0.8)
ax1.plot(theta, g_model, "-", color=color, lw=2, label=r"constant-$C_P$ model")
ax1.plot(theta, g_nist, "--", color="0.25", lw=1.5, label="NIST heat capacities")
ax1.plot(T_cross - 273.15, 0, "o", color=color, ms=6)
ax1.annotate(f"{T_cross:.0f} K", (T_cross - 273.15, 0), xytext=(8, -14), textcoords="offset points", fontsize=9)
ax1.text(450, -27, "operating\nrange", ha="center", va="center", fontsize=9)
ax1.set_xlabel("Temperature (°C)")
ax1.set_ylabel(r"$\Delta_r G^\circ$ (kJ mol$^{-1}$)")
ax1.set_xlim(25, 600)
ax1.set_title("(a) At 1 bar", loc="left", fontsize=10)
ax1.legend(frameon=False, fontsize=9, loc="upper left")

# (b) Pressure at which the reaction Gibbs energy vanishes
mask = T > T_cross + 1
P_model = np.exp(g_model[mask] * 1000 / (2 * R * T[mask]))
P_nist = np.exp(g_nist[mask] * 1000 / (2 * R * T[mask]))
ax2.fill_between([400, 500], 150, 350, color="0.88", lw=0)
ax2.plot(theta[mask], P_model, "-", color=color, lw=2)
ax2.plot(theta[mask], P_nist, "--", color="0.25", lw=1.5)
for th in (400, 500):
    P_th = np.exp(dG_model(th + 273.15) / (2 * R * (th + 273.15)))
    ax2.plot(th, P_th, "o", color=color, ms=6, zorder=3)
ax2.text(450, 420, "operating\nwindow", ha="center", va="bottom", fontsize=9)
ax2.text(215, 200, r"$\Delta_r G < 0$", fontsize=10)
ax2.text(480, 4, r"$\Delta_r G > 0$", fontsize=10)
ax2.set_yscale("log")
ax2.set_xlabel("Temperature (°C)")
ax2.set_ylabel(r"Pressure for $\Delta_r G = 0$ (bar)")
ax2.set_xlim(180, 600)
ax2.set_ylim(1, 2000)
ax2.set_title("(b) Pressure needed", loc="left", fontsize=10)

plt.show()
plt.close(fig)
```

Ammonia synthesis as a function of temperature. (a) $\Delta_r G^\circ$ at 1 bar from the constant-$C_P$ model of this section (solid line) and from the Shomate heat-capacity equations of the NIST Chemistry WebBook (dashed line). Both change sign at 456 K (183 °C). (b) The pressure at which $\Delta_r G = 0$ when each species is at the pressure $P$, from Eq. {eq}`eq:P-crossover`; the points mark the estimates of this section at 400 °C and 500 °C. The shaded areas mark the operating range of the Haber–Bosch process.

The Haber–Bosch process operates at $400\text{–}500\,\text{°C}$ and $150\text{–}350\ \mathrm{bar}$ over an iron-based catalyst. Each operating parameter addresses a different aspect of the thermodynamics/kinetics tradeoff:

- **Catalyst:** lowers the activation barrier, increasing the reaction rate. Does *not* change $\Delta_r G^\circ$ or the equilibrium constant — it only determines how fast equilibrium is reached.
- **Temperature:** higher $T$ gives faster kinetics (Arrhenius), but for this reaction ($\Delta_r S^\circ < 0$) it pushes $\Delta_r G^\circ$ in the unfavorable (positive) direction.
- **Pressure:** for reactions with $\Delta\nu < 0$, higher pressure favors the product side, compensating for the unfavorable $\Delta_r G^\circ(T)$ at elevated temperature.

Industrial conditions represent a compromise: high enough temperature for acceptable rates, high enough pressure for acceptable equilibrium yield, with a catalyst to make the whole process feasible.

:::{admonition} Looking Ahead: Phase Equilibria
:class: seealso
The Gibbs free energy also governs **phase equilibria** — the conditions under which a substance exists as solid, liquid, or gas. The minimization principle $dG \le 0$ at constant $T$ and $P$ ([Section 5.1](section-01.md)) determines which phase is stable under given conditions, generating the familiar $P$–$T$ phase diagrams with coexistence curves, triple points, and critical points. [Section 6.1](../chapter-06/section-01.md) develops this application.
:::

---

## Worked Examples

### Example 1: Crossover Temperature at 1 Bar

**Problem.** At what temperature does $\Delta_r G^\circ$ for ammonia formation change sign (from negative to positive) at 1 bar? Estimate this (a) neglecting $C_P$ corrections, and (b) including the constant-$C_P$ corrections from this section.

**Solution.**

**(a) Without $C_P$ corrections.** If we treat $\Delta_r H^\circ$ and $\Delta_r S^\circ$ as temperature-independent, the crossover occurs when $\Delta_r G^\circ = \Delta_r H^\circ - T\,\Delta_r S^\circ = 0$, giving

```{math}
T_{\mathrm{cross}} = \frac{\Delta_r H^\circ}{\Delta_r S^\circ}
= \frac{-92{,}000\ \mathrm{J\,mol^{-1}}}{-198\ \mathrm{J\,K^{-1}\,mol^{-1}}}
= 465\ \mathrm{K}
\approx 192\,\text{°C}.
```

**(b) With $C_P$ corrections.** Including the temperature-dependent corrections from this section, the reaction-level $\Delta_r C_P = 2C_{P,\mathrm{NH_3}} - C_{P,\mathrm{N_2}} - 3C_{P,\mathrm{H_2}} = 2(4R) - (7R/2) - 3(7R/2) = -6R \approx -49.9\ \mathrm{J\,K^{-1}\,mol^{-1}}$. We need to solve

```{math}
\Delta_r H^\circ(T_i) + \Delta_r C_P(T - T_i) - T\!\left[\Delta_r S^\circ(T_i) + \Delta_r C_P\ln\!\left(\frac{T}{T_i}\right)\right] = 0,
```

which cannot be solved in closed form but can be solved numerically. The result is

```{math}
T_{\mathrm{cross}} \approx 456\ \mathrm{K} \approx 183\,\text{°C}.
```

**Result.** The reaction becomes non-spontaneous (at 1 bar) around $456\text{–}465\ \mathrm{K}$, depending on whether $C_P$ corrections are included. Either way, the crossover is well below the $400\text{–}500\,\text{°C}$ operating range of Haber–Bosch, confirming that **pressure is essential** to restore thermodynamic favorability at operating temperatures.

---

### Example 2: Pressure Needed at 400°C

**Problem.** Repeat the pressure-for-spontaneity estimate at $T = 673.15\ \mathrm{K}$ ($400\,\text{°C}$) instead of $500\,\text{°C}$.

**Solution.**

First, compute $\Delta_r G^\circ(673\ \mathrm{K})$ using the same constant-$C_P$ framework. The temperature shift is $\Delta T = 673.15 - 298.15 = 375.0\ \mathrm{K}$:

```{math}
\Delta_r H^\circ(673) = -92.0 + (-6R)(375.0) = -92.0 - 18.7 = -110.7\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta_r S^\circ(673) = -198 + (-6R)\ln\!\left(\frac{673.15}{298.15}\right) = -198 - 40.6 = -238.6\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

```{math}
\Delta_r G^\circ(673) = -110.7 - (673.15)(-0.2386) = -110.7 + 160.6 = +49.9\ \mathrm{kJ\,mol^{-1}}
```

Now use Eq. {eq}`eq:P-crossover`:

```{math}
\ln\!\left(\frac{P}{1\ \mathrm{bar}}\right)
= -\frac{49.9 \times 10^3}{(-2)(8.314)(673.15)} = 4.46
\quad\Rightarrow\quad
P \approx 86\ \mathrm{bar}.
```

**Result.** At $400\,\text{°C}$, only about 86 bar is needed for spontaneity, compared with $\sim$320 bar at $500\,\text{°C}$. This illustrates the tradeoff: a lower operating temperature demands less pressure for thermodynamic favorability but gives a slower reaction rate. The industrial compromise ($400\text{–}500\,\text{°C}$, $150\text{–}350\ \mathrm{bar}$) sits squarely in the range our estimates predict.

## Concept Checks

1. Why does $\Delta_r S^\circ < 0$ make higher temperatures less favorable for ammonia formation at equilibrium?
2. How does the sign of $\Delta\nu$ determine whether increasing pressure helps or hurts product formation?
3. Which approximations are embedded in using $G(P) - G^\circ = RT\ln(P/P^\circ)$ for each species? How would the analysis change for a non-ideal gas?
4. A catalyst accelerates both the forward and reverse reactions. Why does this not shift the equilibrium?
5. An engineer proposes running at 150°C (where $\Delta_r G^\circ < 0$ at 1 bar) to avoid the need for high pressure. What practical problem does this create?

## Key Takeaways

- Thermodynamics ($\Delta G$) and kinetics (activation barriers, catalysts) are independent: $\Delta G < 0$ does not guarantee a useful rate, and a catalyst does not change $\Delta G$.
- For reactions with $\Delta_r S^\circ < 0$, increasing $T$ drives $\Delta_r G^\circ$ upward, potentially reversing spontaneity.
- Heat-capacity corrections ($\Delta H$ from $\int C_P\,dT$, $\Delta S$ from $\int C_P/T\,dT$) quantify how reaction thermodynamics shift with temperature.
- The ideal-gas pressure dependence of $G$ introduces $\Delta\nu\,RT\ln(P/P^\circ)$ terms; for $\Delta\nu < 0$, high pressure favors products.
- The Haber–Bosch operating conditions ($400\text{–}500\,\text{°C}$, $150\text{–}350\ \mathrm{bar}$, Fe catalyst) represent an engineering compromise between kinetic rate and thermodynamic yield.
