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

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Ammonia synthesis is one of the most important industrial chemical processes: the Haber–Bosch reaction produces the nitrogen feedstock for most synthetic fertilizers, supporting roughly half the world's food production. It is also a textbook example of the tension between thermodynamics and kinetics. At room temperature the reaction is thermodynamically favorable ($\Delta G_r^\circ < 0$) but kinetically inert. Raising the temperature makes the kinetics manageable but can push $\Delta G_r^\circ$ positive. Increasing pressure can restore favorability for reactions that decrease the number of gas molecules. Understanding and quantifying these tradeoffs is the goal of this section.

We use the free-energy machinery from Section 5.1 and the absolute-entropy framework from Section 5.2 to track how $\Delta G_r^\circ$ changes with temperature and pressure for the reaction

```{math}
\mathrm{N_2(g) + 3\,H_2(g) \rightarrow 2\,NH_3(g)}.
```

The analysis relies on two tools: heat-capacity-based temperature corrections to $\Delta H$ and $\Delta S$, and the ideal-gas pressure dependence of $G$.

---

Learning objectives:

- Articulate the difference between thermodynamic favorability ($\Delta G$) and kinetic accessibility (activation barrier, catalysis).
- Use $\Delta G = \Delta H - T\Delta S$ to reason qualitatively about temperature effects when $\Delta S_r^\circ < 0$.
- Estimate temperature corrections to $\Delta H_r^\circ$ and $\Delta S_r^\circ$ from constant-$C_P$ models based on active degrees of freedom.
- Derive the ideal-gas pressure dependence $G(P_f) - G(P_i) = nRT\ln(P_f/P_i)$ and apply it to reactions with $\Delta\nu \neq 0$.
- Explain how catalyst, temperature, and pressure serve as complementary "engineering knobs" in the Haber–Bosch process.

## Core Ideas and Derivations

### Thermodynamics vs. kinetics

#### Standard thermodynamic data at 298.15 K

The standard reaction quantities for ammonia synthesis at $T = 298.15\ \mathrm{K}$ are:

```{math}
\Delta H_r^\circ(298.15\,\mathrm{K}) \approx -92\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_r^\circ(298.15\,\mathrm{K}) \approx -198\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

```{math}
\Delta G_r^\circ(298.15\,\mathrm{K}) \approx -33\ \mathrm{kJ\,mol^{-1}}
```

These are consistent via $\Delta G_r^\circ = \Delta H_r^\circ - T\,\Delta S_r^\circ$. Since $\Delta G_r^\circ < 0$, the reaction is **thermodynamically spontaneous** at 298 K and 1 bar.

:::{admonition} Data skills: where do these numbers come from?
:class: dropdown

Standard reaction quantities are computed from **tabulated species data**:

1. Look up $\Delta_f H^\circ(298.15\ \mathrm{K})$ and $S^\circ(298.15\ \mathrm{K})$ for each species from a reference source (NIST WebBook, NIST-JANAF, ATcT).
2. Compute $\Delta_r H^\circ$ and $\Delta_r S^\circ$ by stoichiometry (products minus reactants). The absolute entropies $S^\circ$ used here are third-law values of the type discussed in Section 5.2.
3. Combine via $\Delta_r G^\circ = \Delta_r H^\circ - T\,\Delta_r S^\circ$.

<!-- A step-by-step workflow (including a worked ammonia example) is given in Section 7.1.7 and the main data sources are collected in Appendix Y. -->
:::

:::{admonition} Key idea
:class: tip
**Thermodynamics** ($\Delta G$) tells us *which direction is favored at equilibrium*.
**Kinetics** (activation barrier, $E_a$) tells us *how fast the system gets there*.
These are independent questions: a reaction can be thermodynamically favorable yet kinetically inert, or vice versa.
:::

#### Why it is too slow at room temperature

Even with $\Delta G_r^\circ < 0$, the reaction can be extremely slow if the activation barrier is large. For ammonia formation, the $\mathrm{N \equiv N}$ triple bond (bond dissociation energy $\approx 945\ \mathrm{kJ\,mol^{-1}}$) makes the uncatalyzed barrier prohibitively high at 298 K.

A **catalyst** (e.g., an Fe-based catalyst developed by Haber and Bosch, ca. 1909; Nobel Prize 1918) provides an alternative reaction pathway with a lower activation barrier. The catalyst accelerates both the forward and reverse reactions equally, so it **does not change $\Delta G_r^\circ$ or the equilibrium constant** — it only changes how fast equilibrium is reached. Raising the temperature further increases the rate (Arrhenius behavior), so industrial reactors operate at $400\text{–}500\,\text{°C}$ where the barrier is surmountable.

But raising $T$ introduces a thermodynamic penalty, which we now quantify.

---

### Qualitative effect of temperature

Since $\Delta S_r^\circ < 0$ for this reaction, the $-T\,\Delta S_r^\circ$ term in

$$
\Delta G_r^\circ = \Delta H_r^\circ - T\,\Delta S_r^\circ
$$

is **positive** and grows with $T$. At sufficiently high temperature, this positive term overwhelms the negative $\Delta H_r^\circ$, and $\Delta G_r^\circ$ changes sign — the reaction becomes non-spontaneous at 1 bar.

The qualitative rule is: for reactions with $\Delta S_r^\circ < 0$, higher temperatures work *against* thermodynamic favorability.

---

### Temperature corrections from $C_P$

To quantify the temperature effect, we need $\Delta H_r^\circ(T)$ and $\Delta S_r^\circ(T)$ at the operating temperature. The strategy is to correct each species' enthalpy and entropy from the reference temperature $T_i = 298.15\ \mathrm{K}$ to the target temperature $T_f$ using heat-capacity data.

#### Enthalpy correction

At constant pressure (Section 3.3):

```{math}
H(T_f) - H(T_i) = \int_{T_i}^{T_f} C_P(T)\,dT
```

#### Entropy correction

From the enthalpy differential $dH = T\,dS + V\,dP$ at constant $P$ ($dP = 0$):

```{math}
S(T_f) - S(T_i) = \int_{T_i}^{T_f} \frac{C_P(T)}{T}\,dT
```

This is the same integral that appears in the absolute-entropy formula (Section 5.2, Eq. {eq}`eq:absolute-entropy`), applied over a finite temperature range rather than from $0$.

---

### Ideal-gas heat capacities from degrees of freedom

For a simple estimate we treat N$_2$, H$_2$, and NH$_3$ as ideal gases and approximate $C_P$ as temperature-independent, using the number of active (translational + rotational) degrees of freedom $f$:

```{math}
:label: eq:Cp-dof
C_P = \left(\frac{f}{2} + 1\right)R
```

For **linear molecules** (N$_2$, H$_2$): $f = 3\ (\text{trans}) + 2\ (\text{rot}) = 5$, giving

```{math}
C_P = \tfrac{7}{2}R = 29.1\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

For **nonlinear** NH$_3$: $f = 3\ (\text{trans}) + 3\ (\text{rot}) = 6$, giving

```{math}
C_P = 4R = 33.3\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

:::{admonition} Approximations in this model
:class: note

Equation {eq}`eq:Cp-dof` counts only translational and rotational degrees of freedom — it is the "rigid rotor" approximation. Vibrational modes are neglected. At 773 K, some vibrational modes (particularly the lower-frequency bending modes of NH$_3$) are partially excited, so the true $C_P$ values are somewhat larger than these estimates. This is the same "freeze-out" phenomenon discussed in Section 2.6: vibrational contributions to $C_P$ become significant only when $T$ is comparable to or larger than the vibrational temperature $\Theta_{\mathrm{vib}}$.

Despite this limitation, the constant-$C_P$ model captures the correct qualitative behavior and gives order-of-magnitude estimates for $\Delta G_r^\circ(T)$.
:::

With constant $C_P$, the integrals simplify to:

```{math}
H(T_f) - H(T_i) = C_P(T_f - T_i)
```

```{math}
S(T_f) - S(T_i) = C_P\ln\!\left(\frac{T_f}{T_i}\right)
```

---

### Temperature effect: 298.15 K to 773.15 K (500°C)

Set $T_i = 298.15\ \mathrm{K}$ and $T_f = 773.15\ \mathrm{K}$.

#### Per-species corrections

**N$_2$ and H$_2$** ($C_P = \tfrac{7}{2}R$):

```{math}
\Delta H_{\text{species}} = \tfrac{7}{2}R\,(773.15 - 298.15) \approx 13.8\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_{\text{species}} = \tfrac{7}{2}R\,\ln\!\left(\frac{773.15}{298.15}\right) \approx 27.7\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

**NH$_3$** ($C_P = 4R$):

```{math}
\Delta H_{\text{species}} = 4R\,(773.15 - 298.15) \approx 15.8\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_{\text{species}} = 4R\,\ln\!\left(\frac{773.15}{298.15}\right) \approx 31.7\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

#### Reaction-level corrections

Apply stoichiometric weights ($2 \times \text{NH}_3 - 1 \times \text{N}_2 - 3 \times \text{H}_2$):

```{math}
\Delta H_r^\circ(T_f) \approx \Delta H_r^\circ(T_i) + \big[2(15.8) - (13.8) - 3(13.8)\big]
= \Delta H_r^\circ(T_i) - 23.6\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_r^\circ(T_f) \approx \Delta S_r^\circ(T_i) + \big[2(31.7) - (27.7) - 3(27.7)\big]
= \Delta S_r^\circ(T_i) - 47.4\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

The reaction-level $\Delta C_P$ is negative (products have fewer degrees of freedom than reactants), so both $\Delta H_r^\circ$ and $\Delta S_r^\circ$ become more negative with increasing temperature.

#### $\Delta G_r^\circ$ at 773.15 K

Combining the corrected values:

```{math}
\Delta H_r^\circ(773) = -92.0 + (-23.6) = -115.6\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_r^\circ(773) = -198 + (-47.4) = -245.4\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

```{math}
:label: eq:dG-773
\Delta G_r^\circ(773) = -115.6 - (773.15)(-0.2454) \approx +74.1\ \mathrm{kJ\,mol^{-1}}
```

At $500\,\text{°C}$ and 1 bar, ammonia formation is **non-spontaneous** ($\Delta G_r^\circ > 0$), even though it would be much faster kinetically than at room temperature. This is the central dilemma of the Haber–Bosch process.

---

### Pressure dependence of $G$ for an ideal gas

Temperature is not the only knob available. For gas-phase reactions, **pressure** can strongly affect $\Delta G$ when the reaction changes the number of gas molecules.

From the Gibbs differential $dG = -S\,dT + V\,dP$ (Section 5.1) at constant temperature:

```{math}
dG\big|_T = V\,dP
```

For one mole of an ideal gas, $V = RT/P$. Integrating from $P^\circ$ to $P$:

```{math}
:label: eq:G-pressure
G(T, P) - G^\circ(T) = RT\ln\!\left(\frac{P}{P^\circ}\right)
```

where $G^\circ(T)$ is the standard-state molar Gibbs energy at pressure $P^\circ = 1\ \mathrm{bar}$.

:::{admonition} Assumptions behind Eq. {eq}`eq:G-pressure`
:class: note
This expression assumes (i) ideal-gas behavior ($PV = nRT$), (ii) constant temperature, and (iii) a pure substance (fixed composition). For a reacting *mixture*, the relevant quantity is the reaction quotient $Q_p$ (a ratio of partial pressures weighted by stoichiometric coefficients), and the general expression is $\Delta G_r = \Delta G_r^\circ + RT\ln Q_p$. The simplified form used below — with all species at the same total pressure — captures the correct pressure *scale* for restoring spontaneity but does not describe the equilibrium composition of a mixture. We will develop the general treatment when we introduce the equilibrium constant.
:::

#### Applying to ammonia formation

If each species is compressed from $P^\circ$ to $P$ at fixed $T$, the pressure correction to the reaction Gibbs energy is

```{math}
\Delta G_r(T, P) \approx \Delta G_r^\circ(T) + \Delta\nu\,RT\ln\!\left(\frac{P}{P^\circ}\right),
```

where $\Delta\nu$ is the change in the number of moles of gas:

```{math}
\Delta\nu = 2 - (1 + 3) = -2.
```

Setting $\Delta G_r = 0$ to find the pressure at which the reaction is just spontaneous:

```{math}
:label: eq:P-crossover
P = P^\circ \exp\!\left(-\frac{\Delta G_r^\circ(T)}{\Delta\nu\,RT}\right)
```

Inserting $\Delta G_r^\circ(773) \approx +74.1\ \mathrm{kJ\,mol^{-1}}$, $\Delta\nu = -2$, and $T = 773.15\ \mathrm{K}$:

```{math}
\ln\!\left(\frac{P}{P^\circ}\right)
= -\frac{74.1 \times 10^3}{(-2)(8.314)(773.15)} = 5.77
\quad\Rightarrow\quad
P \approx 320\ \mathrm{bar}
```

Because $\Delta\nu < 0$, increasing pressure favors the side with fewer gas molecules (Le Châtelier's principle). At $500\,\text{°C}$, pressures on the order of several hundred bar can restore thermodynamic favorability.

---

### Summary: the engineering knobs

The Haber–Bosch process operates at $400\text{–}500\,\text{°C}$ and $150\text{–}350\ \mathrm{bar}$ over an iron-based catalyst. Each operating parameter addresses a different aspect of the thermodynamics/kinetics tradeoff:

- **Catalyst:** lowers the activation barrier, increasing the reaction rate. Does *not* change $\Delta G_r^\circ$ or the equilibrium constant — it only determines how fast equilibrium is reached.
- **Temperature:** higher $T$ gives faster kinetics (Arrhenius), but for this reaction ($\Delta S_r^\circ < 0$) it pushes $\Delta G_r^\circ$ in the unfavorable (positive) direction.
- **Pressure:** for reactions with $\Delta\nu < 0$, higher pressure favors the product side, compensating for the unfavorable $\Delta G_r^\circ(T)$ at elevated temperature.

Industrial conditions represent a compromise: high enough temperature for acceptable rates, high enough pressure for acceptable equilibrium yield, with a catalyst to make the whole process feasible.

:::{admonition} Looking ahead: phase equilibria
:class: seealso
The Gibbs free energy also governs **phase equilibria** — the conditions under which a substance exists as solid, liquid, or gas. The minimization principle $dG \le 0$ at constant $T$ and $P$ (Section 5.1) determines which phase is stable under given conditions, generating the familiar $P$–$T$ phase diagrams with coexistence curves, triple points, and critical points. We will develop this application in the next chapter.
:::

---

## Worked Examples

### Example 1: Crossover temperature at 1 bar

**Problem.** At what temperature does $\Delta G_r^\circ$ for ammonia formation change sign (from negative to positive) at 1 bar? Estimate this (a) neglecting $C_P$ corrections, and (b) including the constant-$C_P$ corrections from this section.

**Solution.**

**(a) Without $C_P$ corrections.** If we treat $\Delta H_r^\circ$ and $\Delta S_r^\circ$ as temperature-independent, the crossover occurs when $\Delta G_r^\circ = \Delta H_r^\circ - T\,\Delta S_r^\circ = 0$, giving

```{math}
T_{\mathrm{cross}} = \frac{\Delta H_r^\circ}{\Delta S_r^\circ}
= \frac{-92{,}000\ \mathrm{J\,mol^{-1}}}{-198\ \mathrm{J\,K^{-1}\,mol^{-1}}}
= 465\ \mathrm{K}
\approx 192\,\text{°C}.
```

**(b) With $C_P$ corrections.** Including the temperature-dependent corrections from this section, the reaction-level $\Delta C_P = 2C_{P,\mathrm{NH_3}} - C_{P,\mathrm{N_2}} - 3C_{P,\mathrm{H_2}} = 2(4R) - (7R/2) - 3(7R/2) = -6R \approx -49.9\ \mathrm{J\,K^{-1}\,mol^{-1}}$. We need to solve

$$
\Delta H_r^\circ(T_i) + \Delta C_P(T - T_i) - T\!\left[\Delta S_r^\circ(T_i) + \Delta C_P\ln\!\left(\frac{T}{T_i}\right)\right] = 0,
$$

which cannot be solved in closed form but can be solved numerically. The result is

```{math}
T_{\mathrm{cross}} \approx 456\ \mathrm{K} \approx 183\,\text{°C}.
```

**Result.** The reaction becomes non-spontaneous (at 1 bar) around $456\text{–}465\ \mathrm{K}$, depending on whether $C_P$ corrections are included. Either way, the crossover is well below the $400\text{–}500\,\text{°C}$ operating range of Haber–Bosch, confirming that **pressure is essential** to restore thermodynamic favorability at operating temperatures.

---

### Example 2: Pressure needed at 400°C

**Problem.** Repeat the pressure-for-spontaneity estimate at $T = 673.15\ \mathrm{K}$ ($400\,\text{°C}$) instead of $500\,\text{°C}$.

**Solution.**

First, compute $\Delta G_r^\circ(673\ \mathrm{K})$ using the same constant-$C_P$ framework. The temperature shift is $\Delta T = 673.15 - 298.15 = 375.0\ \mathrm{K}$:

```{math}
\Delta H_r^\circ(673) = -92.0 + (-6R)(375.0) = -92.0 - 18.7 = -110.7\ \mathrm{kJ\,mol^{-1}}
```

```{math}
\Delta S_r^\circ(673) = -198 + (-6R)\ln\!\left(\frac{673.15}{298.15}\right) = -198 - 40.7 = -238.7\ \mathrm{J\,K^{-1}\,mol^{-1}}
```

```{math}
\Delta G_r^\circ(673) = -110.7 - (673.15)(-0.2387) = -110.7 + 160.7 = +50.0\ \mathrm{kJ\,mol^{-1}}
```

Now use Eq. {eq}`eq:P-crossover`:

```{math}
\ln\!\left(\frac{P}{1\ \mathrm{bar}}\right)
= -\frac{50.0 \times 10^3}{(-2)(8.314)(673.15)} = 4.47
\quad\Rightarrow\quad
P \approx 87\ \mathrm{bar}.
```

**Result.** At $400\,\text{°C}$, only about 87 bar is needed for spontaneity, compared with $\sim$320 bar at $500\,\text{°C}$. This illustrates the tradeoff: a lower operating temperature demands less pressure for thermodynamic favorability but gives a slower reaction rate. The industrial compromise ($400\text{–}500\,\text{°C}$, $150\text{–}350\ \mathrm{bar}$) sits squarely in the range our estimates predict.

## Concept Checks

1. Why does $\Delta S_r^\circ < 0$ make higher temperatures less favorable for ammonia formation at equilibrium?
2. How does the sign of $\Delta\nu$ determine whether increasing pressure helps or hurts product formation?
3. Which approximations are embedded in using $G(P) - G^\circ = RT\ln(P/P^\circ)$ for each species? How would the analysis change for a non-ideal gas?
4. A catalyst accelerates both the forward and reverse reactions. Why does this not shift the equilibrium?
5. An engineer proposes running at 300°C (where $\Delta G_r^\circ < 0$ at 1 bar) to avoid the need for high pressure. What practical problem does this create?

## Key Takeaways

- Thermodynamics ($\Delta G$) and kinetics (activation barriers, catalysts) are independent: $\Delta G < 0$ does not guarantee a useful rate, and a catalyst does not change $\Delta G$.
- For reactions with $\Delta S_r^\circ < 0$, increasing $T$ drives $\Delta G_r^\circ$ upward, potentially reversing spontaneity.
- Heat-capacity corrections ($\Delta H$ from $\int C_P\,dT$, $\Delta S$ from $\int C_P/T\,dT$) quantify how reaction thermodynamics shift with temperature.
- The ideal-gas pressure dependence of $G$ introduces $\Delta\nu\,RT\ln(P/P^\circ)$ terms; for $\Delta\nu < 0$, high pressure favors products.
- The Haber–Bosch operating conditions ($400\text{–}500\,\text{°C}$, $150\text{–}350\ \mathrm{bar}$, Fe catalyst) represent an engineering compromise between kinetic rate and thermodynamic yield.
