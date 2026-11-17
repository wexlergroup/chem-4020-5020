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


# 6.1. Phase Diagrams

See also: [Course-wide Conventions & Notation](../notation.md)

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=6&amp;here=6.1"
  title="Course map, opened at Module 6"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-6.png
:class: hidden
:alt: Course map with Module 6 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6, which has no microscopic step. Below the map, one card for each of Sections 6.1 and 6.2 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## Overview and Learning Objectives

In [Section 5.1](../chapter-05/section-01.md) we introduced the Gibbs free energy $G$ and showed that at constant $T$ and $P$, a closed system moves spontaneously in the direction of decreasing $G$ and reaches equilibrium when $G$ is minimized. [Section 5.3](../chapter-05/section-03.md) closed by pointing out that this minimization principle also governs **phase equilibria** — the conditions under which a substance exists as solid, liquid, or gas, and the $(T,P)$ curves along which two phases coexist. That is the application we develop in this module.

A one-component substance at fixed $(T,P)$ can in principle exist in any of several phases, but only one is thermodynamically stable except on special curves in the $(T,P)$ plane where two (or at isolated points, three) phases coexist. These curves generate the familiar **phase diagram**. Our task in this section is to:

1. introduce the **chemical potential** $\mu$ as the open-system generalization of the Gibbs differential,
2. derive the **phase-equilibrium condition** $\mu_\alpha = \mu_\beta$ from the requirement that $G$ be minimized at fixed $T$ and $P$, and
3. interpret that condition graphically on plots of $\bar g$ versus $T$ at fixed $P$ (and $\bar g$ versus $P$ at fixed $T$), and use it to identify the latent heat and entropy jump at a first-order phase transition.

The *quantitative* consequence of the equilibrium condition — the slope $dP/dT$ of a coexistence curve and its relation to the latent heat — is the Clapeyron equation, which we derive in [Section 6.2](section-02.md).

Learning objectives:

- Identify the key features of a one-component $P$–$T$ phase diagram: single-phase regions, coexistence curves, triple point, and critical point.
- Define the chemical potential $\mu_i$ as the open-system extension of $dG = -S\,dT + V\,dP$, and recognize that for a pure substance $\mu = \bar g$, the molar Gibbs free energy.
- Derive the phase-equilibrium condition $\mu_\alpha = \mu_\beta$ from Gibbs minimization at fixed $T,P$ and interpret it as "no driving force for matter transfer."
- Read off the stable phase from intersecting $\bar g$–$T$ or $\bar g$–$P$ curves and relate the slopes to $-\bar s$ and $+\bar v$.
- Relate the latent heat and entropy jump at a first-order transition via $\Delta \bar h_{\mathrm{trs}} = T_{\mathrm{trs}}\,\Delta \bar s_{\mathrm{trs}}$, and use it together with Trouton's rule to estimate $\Delta \bar s_{\mathrm{vap}}$ from tabulated enthalpies of vaporization.

## Core Ideas and Derivations

(section:phase-diagram-tour)=
### A Quick Tour of a One-Component $P$–$T$ Phase Diagram

A schematic $P$–$T$ phase diagram for a single component contains regions where one phase is stable (solid, liquid, gas) separated by **coexistence curves** on which two phases are simultaneously in equilibrium.

<div style="width: 100%; border: 1px solid #cbd5e1; border-radius: 12px; overflow: hidden; box-shadow: 0 4px 12px rgba(0,0,0,0.08);">
  <iframe
    src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/phase-diagram-schematic.html"
    title="Phase Diagram Schematic"
    style="width: 100%; height: 850px; border: 0;"
    loading="lazy"
  ></iframe>
</div>

Two special points deserve named attention:

- **Triple point (tp):** the unique $(T_{\mathrm{tp}},P_{\mathrm{tp}})$ at which solid, liquid, and gas coexist simultaneously.
- **Critical point (cp):** the point at which the liquid–gas coexistence curve *terminates*. Beyond the critical point, liquid and gas become indistinguishable and the system is called a **supercritical fluid**.

The critical point exists because liquid and gas differ only *quantitatively* — they are both disordered, fluid phases whose molar volumes can be made arbitrarily close to one another by heating the liquid and compressing the gas. Solid and liquid, by contrast, differ *qualitatively* (long-range translational order is either present or absent), so the solid–liquid line generally continues indefinitely as pressure increases rather than ending at a critical point.

Moving across a coexistence curve at fixed $P$ (or fixed $T$) changes which phase is thermodynamically stable — the defining feature of a **phase transition**. Our goal for the rest of this section is to characterize these transitions in terms of $G$, and in particular to identify the condition that picks out the coexistence curves in the $(T,P)$ plane.

(section:chemical-potential)=
### Open Systems and the Chemical Potential

For a closed system at fixed composition, we derived in [Section 5.1](../chapter-05/section-01.md) the Gibbs differential

```{math}
:label: eq:dG-closed
dG = -S\,dT + V\,dP.
```

A closed system cannot, however, describe two coexisting phases on its own, because at a phase transition matter is transferred from one phase to the other. To handle this, we regard each phase as a subsystem whose amount of matter can change. Each phase is then an **open system** — one that can exchange energy *and* matter with its surroundings (in this case, the other phase).

For an open, multicomponent system, the Gibbs differential generalizes to

```{math}
:label: eq:dG-open
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i,
```

where the sum runs over the chemical components $i$ and the **chemical potential** $\mu_i$ is defined by the partial derivative

```{math}
:label: eq:mu-definition
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,\,n_{j\neq i}}.
```

In words: $\mu_i$ is how much the Gibbs free energy changes when one mole of component $i$ is added while holding $T$, $P$, and the amounts of all other components fixed. A *component* is one of the independent chemical species needed to specify the composition of every phase in the system, and the number of components is the minimum number of such species; for a pure substance that number is one and there is a single chemical potential $\mu$.

:::{admonition} Picking up a thread from Module 5
:class: note
[Section 5.1](../chapter-05/section-01.md) stated the Euler relation $U = TS - PV + \mu N$ and its consequence $G = \mu N$ for a one-component system, before the chemical potential was defined. In terms of the amount $n$ in moles and the chemical potential per mole defined above, the relation reads

```{math}
:label: eq:G-equals-mu-n
G = \mu\,n,
```

so that the molar Gibbs free energy and the chemical potential are the same object:

```{math}
:label: eq:mu-is-gbar
\mu = \frac{G}{n} \equiv \bar g.
```

We use an overbar for *molar* quantities throughout this module ($\bar g, \bar s, \bar v, \bar h$). For a pure substance, "chemical potential," "molar Gibbs free energy," and "$G$ per mole" all refer to the same quantity, and we will use whichever name is most convenient.
:::

```{admonition} Doubling the Amount of a Phase
:class: important
Suppose that the amount of a phase is doubled at fixed $T$ and $P$. By what factor does its Gibbs free energy $G$ change, and by what factor does its chemical potential $\mu$ change? Use Eqs. {eq}`eq:G-equals-mu-n` and {eq}`eq:mu-is-gbar`. Then decide which of the two quantities can settle the question "which of two phases is stable at this $T$ and $P$?"
```

(section:phase-equilibrium)=
### Phase Equilibrium as "No Driving Force for Matter Transfer"

Consider a closed overall system containing two phases of the *same* pure substance, labeled $\alpha$ and $\beta$, separated by a boundary across which matter can flow. Let the two phases hold $n_\alpha$ and $n_\beta$ moles respectively, with $n_\alpha + n_\beta = n_{\mathrm{tot}}$ held fixed by overall matter conservation. Let the amount of phase $\alpha$ change by an infinitesimal $dn_\alpha$ at fixed $T$ and $P$, through transfer of matter to or from phase $\beta$; then

```{math}
dn_\beta = -dn_\alpha.
```

Applying the open-system Gibbs differential {eq}`eq:dG-open` to each phase and summing, and recognizing that at fixed $T$ and $P$ only the composition terms contribute,

```{math}
dG = \mu_\alpha\,dn_\alpha + \mu_\beta\,dn_\beta
    = (\mu_\alpha - \mu_\beta)\,dn_\alpha.
```

At fixed $T$ and $P$, the equilibrium condition is that $G$ be a minimum with respect to every internal variable, including the amount transferred. Requiring $dG = 0$ for an *arbitrary* small $dn_\alpha$ of either sign gives

```{math}
:label: eq:phase-equilibrium
\boxed{\mu_\alpha(T,P) = \mu_\beta(T,P)}.
```

This is the central result of Section 6.1: **two phases of the same pure substance are in equilibrium at a given $(T,P)$ if and only if they have the same chemical potential.**

The physical reading is that $\mu$ is the "price per mole" of adding substance to a phase. Matter spontaneously flows from higher $\mu$ to lower $\mu$ — from "expensive" to "cheap" — and stops flowing exactly when the two prices are equal. Away from equilibrium, the sign of $(\mu_\alpha - \mu_\beta)$ tells us which direction the transition runs:

- If $\mu_\alpha > \mu_\beta$, then $dG < 0$ when matter moves from $\alpha$ to $\beta$ ($dn_\alpha < 0$), so the $\alpha \to \beta$ transition is spontaneous.
- If $\mu_\alpha < \mu_\beta$, the reverse transition is spontaneous.
- If $\mu_\alpha = \mu_\beta$, the two phases coexist at equilibrium.

Because $\mu$ depends on $T$ and $P$, the equation $\mu_\alpha(T,P) = \mu_\beta(T,P)$ is a single scalar constraint on two variables and therefore picks out a *curve* in the $(T,P)$ plane — precisely the coexistence curve we saw in the phase diagram of [A Quick Tour of a One-Component $P$–$T$ Phase Diagram](#section:phase-diagram-tour).

(section:competing-phases)=
### Competing Phases on $\bar g$-vs-$T$ and $\bar g$-vs-$P$ Plots

The condition $\mu_\alpha = \mu_\beta$ becomes much more intuitive when we draw it graphically. For a pure substance, $\mu = \bar g$, so a phase transition is the intersection of $\bar g$-versus-$T$ (or $\bar g$-versus-$P$) curves for the competing phases. The *slopes* of those curves encode entropy and volume via the Gibbs differential {eq}`eq:dG-closed`.

#### (A) $\bar g$ versus $T$ at Fixed Pressure

At constant pressure, Eq. {eq}`eq:dG-closed` gives

```{math}
\left(\frac{\partial \bar g}{\partial T}\right)_P = -\bar s.
```

On a $\bar g$-vs-$T$ plot at fixed $P$:

- All phases have *negative* slope (entropy is positive).
- The phase with the *larger* molar entropy has the *more negative* (steeper) slope.

Because typically $\bar s_{\mathrm{g}} > \bar s_{\mathrm{l}} > \bar s_{\mathrm{s}}$, the gas line is steepest, then liquid, then solid. Starting at low $T$ the solid curve is the lowest (most stable); as $T$ increases, the liquid curve eventually crosses the solid curve, and then the gas curve crosses the liquid. Each crossing marks a first-order phase transition at which the lowest-$\bar g$ phase — the stable one — changes. The figure below shows the three curves for water at 1 bar.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

# Water at 1 bar, with a constant heat capacity for each phase
T_fus, T_vap = 273.15, 372.76             # K
dh_fus, dh_vap = 6.01e3, 40.67e3          # J/mol
c_s, c_l, c_g = 37.8, 75.3, 36.0          # J/(mol K)
s_l_298 = 69.95                           # J/(mol K), liquid at 298.15 K

# Molar enthalpy (zero: liquid at T_fus) and entropy at the two transitions
s_l_fus = s_l_298 + c_l * np.log(T_fus / 298.15)
h_l_vap = c_l * (T_vap - T_fus)
s_l_vap = s_l_fus + c_l * np.log(T_vap / T_fus)

def g_phase(T, T_ref, h_ref, s_ref, c):
    # Molar Gibbs free energy in kJ/mol for a constant heat capacity c
    h = h_ref + c * (T - T_ref)
    s = s_ref + c * np.log(T / T_ref)
    return (h - T * s) / 1000

phases = {
    "ice": (T_fus, -dh_fus, s_l_fus - dh_fus / T_fus, c_s, "#0072B2"),
    "liquid": (T_fus, 0.0, s_l_fus, c_l, "#009E73"),
    "steam": (T_vap, h_l_vap + dh_vap, s_l_vap + dh_vap / T_vap, c_g, "#D55E00"),
}
stable = {"ice": (230, T_fus), "liquid": (T_fus, T_vap), "steam": (T_vap, 420)}
drawn = {"ice": (230, 330), "liquid": (230, 420), "steam": (318, 420)}

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
for name, (T_ref, h_ref, s_ref, c, color) in phases.items():
    T = np.linspace(*drawn[name], 200)
    ax.plot(T, g_phase(T, T_ref, h_ref, s_ref, c), "--", color=color, lw=1.2)
    T = np.linspace(*stable[name], 200)
    ax.plot(T, g_phase(T, T_ref, h_ref, s_ref, c), "-", color=color, lw=2.5)

for T_x, label in ((T_fus, "273.15 K"), (T_vap, "372.76 K")):
    g_x = g_phase(T_x, T_fus, 0.0, s_l_fus, c_l)
    ax.plot(T_x, g_x, "o", color="0.2", ms=5, zorder=3)
    ax.annotate(label, (T_x, g_x), xytext=(-8, -14), textcoords="offset points",
                ha="right", fontsize=9)

ax.text(240, -17.4, "ice", color="#0072B2", fontsize=10)
ax.text(318, -22.2, "liquid", color="#009E73", fontsize=10)
ax.text(402, -30.2, "steam", color="#D55E00", fontsize=10)
ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\bar g$ (kJ mol$^{-1}$)")
ax.set_xlim(230, 420)
ax.set_ylim(-35, -13)

plt.show()
plt.close(fig)
```

Molar Gibbs free energy of ice, liquid water, and steam at 1 bar. Each curve is solid where its phase is the stable one and dashed where it is not. The curves cross at the melting point ($273.15\ \mathrm{K}$) and at the boiling point ($372.76\ \mathrm{K}$ at 1 bar). The curves were computed with a constant heat capacity for each phase ($37.8$, $75.3$, and $36\ \mathrm{J\,mol^{-1}\,K^{-1}}$), the standard entropy of the liquid at $298.15\ \mathrm{K}$ ($69.95\ \mathrm{J\,mol^{-1}\,K^{-1}}$), and the enthalpies of fusion ($6.01\ \mathrm{kJ\,mol^{-1}}$) and vaporization ($40.67\ \mathrm{kJ\,mol^{-1}}$). The zero of enthalpy is that of the liquid at $273.15\ \mathrm{K}$.

#### (B) $\bar g$ versus $P$ at Fixed Temperature

At constant temperature, Eq. {eq}`eq:dG-closed` gives

```{math}
\left(\frac{\partial \bar g}{\partial P}\right)_T = \bar v.
```

On a $\bar g$-vs-$P$ plot at fixed $T$:

- All phases have *positive* slope (molar volumes are positive).
- The phase with the *larger* molar volume has the *steeper* slope.

Because typically $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}} \gtrsim \bar v_{\mathrm{s}}$, the gas line rises most steeply with pressure. Compressing the system therefore favors the phases with smaller molar volume (liquid or solid over gas) — the graphical statement of Le Châtelier's principle for phase transitions. Water is the familiar counter-example on the solid–liquid side: $\bar v_{\mathrm{s}} > \bar v_{\mathrm{l}}$ for water (ice is less dense than liquid water). Concept Check 5 asks what this does to the two curves.

In both views, the **stable phase at each point** is the one with the *lowest* $\bar g$, and a phase transition occurs at each intersection. The graphical picture is particularly useful for understanding *why* the coexistence curves have the shapes they do, for example why increasing pressure generally raises the melting and boiling points of normal substances.

(section:latent-heat)=
### Latent Heat and Entropy Jump at a First-Order Transition

At most phase transitions (fusion, vaporization, sublimation, solid–solid polymorphic transitions) the first derivatives of $G$ — namely $S$ and $V$ — are **discontinuous** across the coexistence curve. Transitions of this kind are called **first-order**. The macroscopic consequences are two quantities that every student of thermodynamics encounters: a **latent heat** and an **entropy jump**.

For a transition $\alpha \to \beta$ occurring at $(T_{\mathrm{trs}},P_{\mathrm{trs}})$, define the molar changes

```{math}
\Delta \bar h_{\mathrm{trs}} \equiv \bar h_\beta(T_{\mathrm{trs}}) - \bar h_\alpha(T_{\mathrm{trs}}),
\qquad
\Delta \bar s_{\mathrm{trs}} \equiv \bar s_\beta(T_{\mathrm{trs}}) - \bar s_\alpha(T_{\mathrm{trs}}).
```

At the coexistence curve, $\mu_\alpha = \mu_\beta$ implies $\bar g_\alpha = \bar g_\beta$ and therefore $\Delta \bar g_{\mathrm{trs}} = 0$. Combined with $\Delta \bar g = \Delta \bar h - T\,\Delta \bar s$, this gives the key identity

```{math}
:label: eq:latent-entropy
\boxed{\Delta \bar h_{\mathrm{trs}} = T_{\mathrm{trs}}\,\Delta \bar s_{\mathrm{trs}}}.
```

Equation {eq}`eq:latent-entropy` connects a quantity one can *measure calorimetrically* — the latent heat $\Delta \bar h_{\mathrm{trs}}$, released or absorbed when matter crosses the boundary at constant pressure — to a quantity that usually appears only in tables, the entropy jump $\Delta \bar s_{\mathrm{trs}}$. It is the reason $\Delta \bar s_{\mathrm{vap}}$ is rarely tabulated separately from $\Delta \bar h_{\mathrm{vap}}$: given one, the other follows by division.

**Common special cases.** At the melting temperature $T_{\mathrm{fus}}$ (at the chosen pressure),

```{math}
\Delta \bar h_{\mathrm{fus}} = \bar h_{\mathrm{l}}(T_{\mathrm{fus}}) - \bar h_{\mathrm{s}}(T_{\mathrm{fus}}), \qquad
\Delta \bar s_{\mathrm{fus}} = \bar s_{\mathrm{l}}(T_{\mathrm{fus}}) - \bar s_{\mathrm{s}}(T_{\mathrm{fus}}),
```

and at the boiling temperature $T_{\mathrm{vap}}$,

```{math}
\Delta \bar h_{\mathrm{vap}} = \bar h_{\mathrm{g}}(T_{\mathrm{vap}}) - \bar h_{\mathrm{l}}(T_{\mathrm{vap}}), \qquad
\Delta \bar s_{\mathrm{vap}} = \bar s_{\mathrm{g}}(T_{\mathrm{vap}}) - \bar s_{\mathrm{l}}(T_{\mathrm{vap}}).
```

Both are positive at ordinary first-order transitions: the higher-temperature phase is the one with greater molar enthalpy *and* greater molar entropy.

(section:trouton-rule)=
### Trouton's Rule and Water as an Exception

Frederick Trouton (1884) noticed that for many ordinary liquids at their *normal* boiling point,

```{math}
:label: eq:trouton
\Delta \bar s_{\mathrm{vap}} \approx 85.9\ \mathrm{J\,mol^{-1}\,K^{-1}} \approx 10.3\,R.
```

The rough constancy of $\Delta \bar s_{\mathrm{vap}}$ across chemically unrelated liquids suggests that the entropy gained when one mole of molecules leaves a liquid for a nearly ideal gas is governed mostly by the phase change itself (the volume increase, the loss of short-range order) rather than by the chemical identity of the liquid.

Not every liquid follows the rule. Hydrogen-bonded liquids lie well above Trouton's value: at the normal boiling point, $\Delta \bar s_{\mathrm{vap}}$ is about $97\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for ammonia, $104$ for methanol, $109$ for water, and $110$ for ethanol. Water is the familiar example: $\Delta \bar s_{\mathrm{vap}}(\mathrm{H_2O}) \approx 109\ \mathrm{J\,mol^{-1}\,K^{-1}}$ is roughly 25 % above Trouton's value. Liquids that boil at very low temperatures lie well below it: $\Delta \bar s_{\mathrm{vap}}$ is $64\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for neon, $44$ for hydrogen, and $19.5$ (only $2.3\,R$) for helium.

```{admonition} Reading a Deviation from Trouton's Rule
:class: important
At 1 atm the vapors of all these liquids are nearly ideal gases. If $\Delta \bar s_{\mathrm{vap}}$ of a liquid lies well above Trouton's value, what does that say about the liquid? Concept Check 4 asks this for water.
```

We will revisit the molecular interpretation of $\Delta \bar h_{\mathrm{vap}}$ for water in [Section 6.2](section-02.md) after introducing the Clausius–Clapeyron equation as a tool for extracting enthalpies of vaporization from vapor-pressure data.

:::{admonition} Looking Ahead: Slopes of Coexistence Curves
:class: seealso
Equation {eq}`eq:phase-equilibrium` defines each coexistence curve implicitly as the set of $(T,P)$ pairs for which $\mu_\alpha(T,P) = \mu_\beta(T,P)$. To obtain the *slope* $dP/dT$ of that curve explicitly, we differentiate the equilibrium condition along the curve and apply the Gibbs differential to each phase. That calculation produces the **Clapeyron equation** and — in the liquid–vapor limit with an ideal gas — its specialization the **Clausius–Clapeyron equation**. Both are developed in [Section 6.2](section-02.md).
:::

## Worked Example

### Heating Water from 298.15 K to 800 K at 1 bar

**Goal.** Estimate $\Delta \bar s$ and $\Delta \bar h$ when one mole of water is heated reversibly at constant pressure $P = 1$ bar from **298.15 K** (liquid) to **800 K** (superheated steam), crossing the liquid–vapor coexistence curve. This example combines the single-phase heating integrals developed in [Section 3.3](../chapter-03/section-03.md) and [Section 5.3](../chapter-05/section-03.md) with the first-order transition relation {eq}`eq:latent-entropy`, and illustrates how strongly the overall entropy change is dominated by the phase transition itself.

Within a single phase at constant pressure,

```{math}
\Delta \bar s = \int_{T_1}^{T_2} \frac{\bar c_P(T)}{T}\,dT,
\qquad
\Delta \bar h = \int_{T_1}^{T_2} \bar c_P(T)\,dT,
```

and at the coexistence temperature $T_{\mathrm{trs}}$ the entropy has a finite jump

```{math}
\Delta \bar s_{\mathrm{vap}}(T_{\mathrm{trs}}) = \frac{\Delta \bar h_{\mathrm{vap}}(T_{\mathrm{trs}})}{T_{\mathrm{trs}}}.
```

We use a coarse set of $\bar c_P^{\,\circ}(T)$ values from the **NIST–JANAF thermochemical tables** for $\mathrm{H_2O(l)}$ and $\mathrm{H_2O(g)}$ at $P^\circ = 0.1\ \mathrm{MPa} = 1\ \mathrm{bar}$. In that dataset, the liquid–vapor transition occurs at

```{math}
T_{\mathrm{trs}} = 372.780\ \mathrm{K}\quad\text{(the normal boiling point, at 1 atm, is 373.12 K).}
```

#### Data Used (Excerpt)

| Phase | $T$ (K) | $\bar c_P^{\,\circ}$ (J mol$^{-1}$ K$^{-1}$) |
| --- | ---: | ---: |
| $\mathrm{H_2O(l)}$ | 298.15 | 75.351 |
| | 300 | 75.349 |
| | 320 | 75.344 |
| | 340 | 75.388 |
| | 360 | 75.679 |
| | 372.780 | 75.962 |
| $\mathrm{H_2O(g)}$ | 300 | 33.596 |
| | 400 | 34.262 |
| | 500 | 35.226 |
| | 600 | 36.325 |
| | 700 | 37.495 |
| | 800 | 38.721 |

The gas-phase heat capacity at $T_{\mathrm{trs}}$ is obtained by linear interpolation of the 300–400 K JANAF entries. The latent heat at $T_{\mathrm{trs}}$ is obtained from the tabulated standard enthalpies of formation:

```{math}
\Delta \bar h_{\mathrm{vap}}(T_{\mathrm{trs}})
\approx
\Delta_f H^\circ_{\mathrm{g}}(T_{\mathrm{trs}}) - \Delta_f H^\circ_{\mathrm{l}}(T_{\mathrm{trs}}).
```

```{code-cell} ipython3
import numpy as np

# ---- Cp(T) data (J/mol/K) ----
T_tr = 372.780  # K (liquid <-> vapor at 1 bar in the JANAF tables)

T_liq = np.array([298.15, 300, 320, 340, 360, T_tr])
Cp_liq = np.array([75.351, 75.349, 75.344, 75.388, 75.679, 75.962])

T_gas_tab = np.array([300, 400, 500, 600, 700, 800])
Cp_gas_tab = np.array([33.596, 34.262, 35.226, 36.325, 37.495, 38.721])

# Interpolate Cp_gas at T_tr (between 300 and 400 K)
Cp_gas_tr = np.interp(T_tr, T_gas_tab, Cp_gas_tab)

T_gas = np.concatenate([[T_tr], T_gas_tab[1:]])
Cp_gas = np.concatenate([[Cp_gas_tr], Cp_gas_tab[1:]])

# ---- Heating contributions (trapezoid rule) ----
dS_liq = np.trapezoid(Cp_liq/T_liq, T_liq)          # J/mol/K
dH_liq = np.trapezoid(Cp_liq, T_liq)/1000           # kJ/mol

dS_gas = np.trapezoid(Cp_gas/T_gas, T_gas)          # J/mol/K
dH_gas = np.trapezoid(Cp_gas, T_gas)/1000           # kJ/mol

# ---- Vaporization jump from the JANAF enthalpies of formation ----
# (linear interpolation to T_tr; values in kJ/mol)
Hf_liq_tr = np.interp(T_tr, [360, 400], [-283.874, -282.591])   # liquid
Hf_gas_tr = np.interp(T_tr, [300, 400], [-241.844, -242.846])   # gas

dH_vap = Hf_gas_tr - Hf_liq_tr                 # kJ/mol
dS_vap = dH_vap*1000/T_tr                      # J/mol/K

dH_total = dH_liq + dH_vap + dH_gas
dS_total = dS_liq + dS_vap + dS_gas

rows = [
    (f"Heating the liquid, 298.15 K to {T_tr:.2f} K", dH_liq, dS_liq),
    (f"Vaporization at {T_tr:.2f} K", dH_vap, dS_vap),
    (f"Heating the gas, {T_tr:.2f} K to 800 K", dH_gas, dS_gas),
    ("Total", dH_total, dS_total),
]
print(f"{'Step':<42}{'ΔH (kJ/mol)':>12}{'ΔS (J/mol/K)':>15}")
for name, dH, dS in rows:
    print(f"{name:<42}{dH:>12.2f}{dS:>15.2f}")
```

**Results.** Heating the liquid contributes $\Delta \bar h = 5.63\ \mathrm{kJ\,mol^{-1}}$ and $\Delta \bar s = 16.87\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Vaporization contributes $40.89\ \mathrm{kJ\,mol^{-1}}$ and $109.69\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Heating the gas contributes $15.48\ \mathrm{kJ\,mol^{-1}}$ and $27.57\ \mathrm{J\,mol^{-1}\,K^{-1}}$. The totals are $\Delta \bar h \approx 62.01\ \mathrm{kJ\,mol^{-1}}$ and $\Delta \bar s \approx 154.1\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

The enthalpy of vaporization found here, $40.89\ \mathrm{kJ\,mol^{-1}}$, is that of the tables, in which the gas is in its ideal-gas standard state at 1 bar. For real steam the measured value at 1 bar is $40.67\ \mathrm{kJ\,mol^{-1}}$, which gives $\Delta \bar s_{\mathrm{vap}} = 109.1\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**What to notice.** The total entropy change along this path is dominated by the phase transition: $\Delta \bar s_{\mathrm{vap}} \approx 1.1\times 10^2\ \mathrm{J\,mol^{-1}\,K^{-1}}$, which alone is about 2.5 times the combined contribution of heating the liquid by about 75 K and the vapor by about 430 K.

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt

# Cumulative entropy relative to 298.15 K along the path (trapezoid rule)
def cumulative(T, Cp):
    steps = 0.5 * (Cp[1:] / T[1:] + Cp[:-1] / T[:-1]) * np.diff(T)
    return np.concatenate([[0.0], np.cumsum(steps)])

S_liq_cum = cumulative(T_liq, Cp_liq)
S_gas_cum = S_liq_cum[-1] + dS_vap + cumulative(T_gas, Cp_gas)

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.plot(T_liq, S_liq_cum, "-", color="#009E73", lw=2.5)
ax.plot(T_gas, S_gas_cum, "-", color="#D55E00", lw=2.5)
ax.plot([T_tr, T_tr], [S_liq_cum[-1], S_gas_cum[0]], "-", color="0.3", lw=2.5)

ax.text(325, 17, "liquid", color="#009E73", fontsize=10)
ax.text(560, 131, "steam", color="#D55E00", fontsize=10)
ax.text(T_tr + 10, 70, "vaporization\nat 372.78 K", color="0.3", fontsize=10,
        va="center")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\Delta \bar s$ from 298.15 K (J mol$^{-1}$ K$^{-1}$)")
ax.set_xlim(280, 820)
ax.set_ylim(-5, 165)

plt.show()
plt.close(fig)
```

Entropy of one mole of water at 1 bar, relative to the liquid at $298.15\ \mathrm{K}$, along the heating path of this example. The vertical segment at $T_{\mathrm{trs}} = 372.78\ \mathrm{K}$ is the entropy of vaporization.

The sharp vertical step at $T_{\mathrm{trs}}$ is the graphical signature of a first-order transition: a finite entropy jump at a single temperature, sitting atop a smooth continuous rise from heat-capacity integration through each single-phase region.

## Concept Checks

1. Why is the chemical potential, rather than the Gibbs free energy itself, the natural variable for analyzing phase equilibria?
2. On a $\bar g$-vs-$T$ plot at fixed pressure, what does the slope of each curve represent, and why do different phases have different slopes? Use your answer to predict the order of the solid, liquid, and gas slopes.
3. Starting from $\mu_\alpha = \mu_\beta$ at coexistence, argue in words why this condition defines a *curve* in the $(T,P)$ plane rather than an isolated point or a two-dimensional region. What happens when a *third* phase is added to the argument?
4. Water has $\Delta \bar s_{\mathrm{vap}} \approx 109\ \mathrm{J\,mol^{-1}\,K^{-1}}$, about 25 % higher than Trouton's value. What does this tell you about the liquid relative to an "average" liquid? (Nothing yet about the gas.)
5. Ice floats on water, so $\bar v_{\mathrm{s}} > \bar v_{\mathrm{l}}$ for water. On a $\bar g$-vs-$P$ plot at fixed $T$ just below the normal freezing point, sketch the ice and water curves and identify which has the steeper slope. What does your sketch predict for the effect of increasing $P$ on the stable phase at that temperature?

## Key Takeaways

- For an open system, the Gibbs differential acquires a composition term: $dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i$. The chemical potential $\mu_i \equiv (\partial G/\partial n_i)_{T,P,n_{j\neq i}}$ is the intensive quantity conjugate to amount of substance, and for a pure substance $\mu = \bar g$.
- Two phases of a pure substance are in equilibrium if and only if they have equal chemical potentials: $\mu_\alpha(T,P) = \mu_\beta(T,P)$. Because this is a single constraint on two variables, it defines the coexistence curve in the $(T,P)$ plane.
- The stable phase at a given $(T,P)$ is the one with the lowest $\bar g$. On $\bar g$–$T$ plots (fixed $P$) and $\bar g$–$P$ plots (fixed $T$), the slopes are $-\bar s$ and $+\bar v$ respectively; intersections mark first-order phase transitions.
- At a first-order transition, the latent heat and entropy jump are tied together by $\Delta \bar h_{\mathrm{trs}} = T_{\mathrm{trs}}\,\Delta \bar s_{\mathrm{trs}}$.
- Trouton's rule, $\Delta \bar s_{\mathrm{vap}} \approx 85.9\ \mathrm{J\,mol^{-1}\,K^{-1}}$ at the normal boiling point, is a useful benchmark; hydrogen-bonded liquids such as water lie well above it, and liquids that boil at very low temperatures lie well below it.
- The *slope* $dP/dT$ of a coexistence curve — the quantitative statement of how a phase boundary bends through the $(T,P)$ plane — is the subject of [Section 6.2](section-02.md).
