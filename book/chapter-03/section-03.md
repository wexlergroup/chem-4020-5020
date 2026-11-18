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
numbering:
  enumerator: 3.3.%s
---


# 3.3. Enthalpy

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 3.3](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-03/section-03.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

Most chemistry happens at constant pressure—reactions in open beakers, biological processes at atmospheric pressure, industrial reactors vented to the atmosphere. Under these conditions, the heat we measure is not $\Delta U$ but something slightly different: some of the energy goes into $PV$ work as the system expands or contracts against the atmosphere. Enthalpy $H = U + PV$ provides the natural accounting for constant-pressure heat flow, absorbing that $PV$ bookkeeping into a single state function.

This section derives enthalpy from the First Law at constant pressure, relates $q_P$ to $\Delta H$, introduces the constant-pressure heat capacity $C_P$, and develops the standard enthalpy framework (formation and reaction enthalpies, Hess's Law) used throughout thermochemistry.

---

Learning objectives:

- Derive the definition $H=U+PV$ and show that $q_P=\Delta H$ for constant-pressure processes with only $PV$ work.
- Distinguish conditions under which $q_P = \Delta H$ holds from those where it does not (e.g., non-$PV$ work at constant pressure).
- Define heat capacity at constant pressure $C_P$ and use it to compute enthalpy changes over temperature intervals.
- Explain standard states and the meaning of standard enthalpy of formation and reaction.
- Apply Hess's Law to compute $\Delta_r H^\circ$ from formation enthalpies.

## Core Ideas and Derivations

### Defining Enthalpy

#### The Problem: $\delta q_V = dU$ Is Nice—Can We Do the Same at Constant $P$?

In [Section 3.2](section-02.md), we saw that at constant volume ($dV = 0$), the First Law simplifies beautifully:

```{math}
\delta q_V = dU.
```

The heat absorbed at constant volume equals the change in a state function ($U$), which means $q_V$ is path-independent for any process between two fixed states at the same volume. This is experimentally powerful: measuring $q$ in a bomb calorimeter (constant $V$) directly gives $\Delta U$.

But most chemistry happens at constant *pressure*, not constant volume. At constant $P$, with only $PV$ work and the pressure of the system equal to the constant external pressure ($P = P_{\mathrm{ext}}$, so that $\delta w = -P\,dV$), the First Law gives

```{math}
\delta q_P = dU + P\,dV,
```

which involves *two* terms. Can we define a single state function that plays the same role at constant $P$ that $U$ plays at constant $V$? That is, can we find a function $H$ such that $\delta q_P = dH$?

```{admonition} A recurring strategy in thermodynamics
:class: note

The move we are about to make—defining a new state function to simplify expressions under a specific constraint—is a *general strategy* you will see again. When the natural variables for a problem don't match your existing state function, define a new one. We'll use this same approach to define the Helmholtz free energy $A = U - TS$ (convenient at constant $T$ and $V$) and the Gibbs free energy $G = H - TS$ (convenient at constant $T$ and $P$) in [Section 5.1](../chapter-05/section-01.md).
```

#### Derivation

Consider the First Law expression when $T$ and $P$ are the independent variables. Starting from $\delta q = dU + P\,dV$ (with $P = P_{\mathrm{ext}}$, as above) and writing total differentials of $U(T,P)$ and $V(T,P)$:

```{math}
:label: delta-q-TP
\delta q 
= \left[\left(\frac{\partial U}{\partial T}\right)_P + P \left(\frac{\partial V}{\partial T}\right)_P\right] dT 
\;+\; \left[\left(\frac{\partial U}{\partial P}\right)_T 
+ P \left(\frac{\partial V}{\partial P}\right)_T\right] dP.
```

At constant pressure ($dP = 0$), this becomes:

```{math}
:label: differential_first_law_constant_pressure
\delta q_P 
= \left[\left(\frac{\partial U}{\partial T}\right)_P + P \left(\frac{\partial V}{\partial T}\right)_P\right] dT.
```

The coefficient of $dT$ is what we call the **heat capacity at constant pressure**:

```{math}
:label: Cp-definition
C_P 
\;=\; \left(\frac{\partial U}{\partial T}\right)_P 
\;+\; P \left(\frac{\partial V}{\partial T}\right)_P.
```

Compare this with the heat capacity at constant volume, $C_V = \left(\frac{\partial U}{\partial T}\right)_V$, from [Section 3.2](section-02.md).

We seek a state function $H$ such that

```{math}
:label: differential_enthalpy
\delta q_P 
= \left(\frac{\partial H}{\partial T}\right)_P \, dT
= dH \quad\text{(at constant }P\text{)}.
```

Comparing Equations {eq}`differential_first_law_constant_pressure` and {eq}`differential_enthalpy`, we need $(\partial H/\partial T)_P = C_P$. This motivates the definition of **enthalpy**:

```{math}
:label: enthalpy
\boxed{H \;=\; U + PV.}
```

```{admonition} A direct check
:class: important

Write the total differential of $H = U + PV$. What remains when $dP = 0$? Compare the result with $\delta q_P$.
```

````{admonition} Verification: $H = U + PV$ yields $\delta q_P = dH$
:class: dropdown

Taking the partial derivative of $H = U + PV$ with respect to $T$ at constant $P$:

```{math}
\begin{align*}
\left(\frac{\partial H}{\partial T}\right)_P
&= \frac{\partial}{\partial T} \bigl(U + PV\bigr)\bigg|_P \\[6pt]
&= \left(\frac{\partial U}{\partial T}\right)_P 
  + P \left(\frac{\partial V}{\partial T}\right)_P 
  + V \underbrace{\left(\frac{\partial P}{\partial T}\right)_P}_{=\,0} \\[6pt]
&= \left(\frac{\partial U}{\partial T}\right)_P
  \;+\; P \left(\frac{\partial V}{\partial T}\right)_P
\;=\; C_P.
\end{align*}
```

This matches Equation {eq}`Cp-definition`, confirming that $\delta q_P = C_P\,dT = (\partial H/\partial T)_P\,dT = dH$ at constant pressure. ✓
````

````{admonition} **Quick check**
:class: important

In [Section 3.2](section-02.md), we showed that $\Delta U = 0$ for an isothermal ideal-gas process ($U$ depends only on $T$). Is $\Delta H$ also zero for an isothermal ideal-gas process? Why or why not?

```{admonition} Answer
:class: dropdown

**Yes, $\Delta H = 0$.** For an ideal gas, $PV = nRT$, so

$$\Delta H = \Delta U + \Delta(PV) = \Delta U + nR\,\Delta T = 0 + 0 = 0.$$

Both $\Delta U$ and $\Delta(PV)$ vanish when $\Delta T = 0$. This makes sense: $H$, like $U$, depends only on $T$ for an ideal gas.
```
````

```{admonition} Condition for $q_P = \Delta H$
:class: warning

The relation $q_P = \Delta H$ holds only when $PV$ work is the *sole* form of work. If the system also does another kind of work at constant pressure ({numref}`generalized-work-table` in [Section 3.1](section-01.md) lists several), then $q_P \neq \Delta H$.
```

---

### Measuring Enthalpy and Enthalpy Changes

At constant pressure (with $PV$-only work), the heat absorbed or released by a process equals the change in enthalpy:

```{math}
:label: qP-DeltaH
q_P = \Delta H.
```

When the process is a change of temperature from $T_i$ to $T_f$ with no reaction or phase change, $\Delta H = H(T_f) - H(T_i)$. If $C_P$ is approximately constant over the temperature range $\Delta T = T_f - T_i$, then

```{math}
:label: DeltaH-Cp
\Delta H 
= \int_{T_i}^{T_f} C_P(T)\,dT 
\;\approx\; C_P \,\Delta T.
```

#### Calorimetry

A typical way to measure $\Delta H$ experimentally is via **calorimetry**—measuring the temperature change of a known mass of material (often water) that exchanges heat with the process of interest. [](#fig:calorimeter) shows the parts of a constant-pressure calorimeter.

```{code-cell} ipython3
:label: fig:calorimeter
:caption: Schematic of a constant-pressure calorimeter. The reaction occurs in the sample container; the released or absorbed heat ($q_P$) flows into the surrounding water, producing a measurable temperature change. Insulating walls minimize heat loss to the room. Since the process occurs at atmospheric pressure, $q_P = \Delta H$.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs


def draw_calorimeter(c):
    fig, ax = plt.subplots(figsize=(5.3, 3.7), constrained_layout=True)

    # Insulating walls, inner container, and water
    ax.add_patch(Rectangle((0.3, 0.2), 3.4, 3.2, fc=c.fill, ec=c.muted, lw=1.6, ls=(0, (4, 3))))
    ax.add_patch(Rectangle((0.6, 0.4), 2.8, 2.6, fc=c.ground, ec="none"))
    ax.add_patch(Rectangle((0.6, 0.4), 2.8, 2.1, fc=c.tint(c.blue, 0.75), ec="none"))
    ax.add_patch(Rectangle((0.6, 0.4), 2.8, 2.6, fc="none", ec=c.ink, lw=1.6))
    ax.text(1.35, 0.98, "water\n(known mass,\nknown $C_P$)", ha="center", va="center")

    # Sample container and the heat that flows from it into the water
    ax.add_patch(Rectangle((2.2, 0.5), 0.95, 0.7, fc=c.tint(c.vermilion, 0.75), ec=c.vermilion, lw=1.6))
    ax.text(2.675, 0.85, "sample", ha="center", va="center")
    fs.arrow(ax, (2.675, 1.28), (2.675, 2.05), c.vermilion, lw=1.8)
    ax.text(2.82, 1.66, "$q_P$", color=c.ink, va="center")

    # Thermometer and stirrer, which pass through the top wall
    ax.plot([0.95, 0.95], [1.9, 3.65], color=c.ink, lw=2.6, solid_capstyle="round")
    ax.plot([0.95], [1.9], "o", color=c.ink, ms=8)
    ax.text(0.95, 3.75, "thermometer", ha="center", va="bottom")
    ax.plot([1.9, 1.9], [1.65, 3.65], color=c.ink, lw=1.6)
    ax.plot([1.68, 1.9, 2.12], [1.8, 1.65, 1.8], color=c.ink, lw=1.6)
    ax.text(1.9, 3.75, "stirrer", ha="center", va="bottom")

    # Label of the walls
    ax.plot([3.55, 3.95], [2.75, 2.75], color=c.muted, lw=1.0)
    ax.text(4.03, 2.75, "insulating walls", color=c.muted, va="center")

    ax.set(xlim=(0.2, 5.3), ylim=(0.1, 4.0), aspect="equal")
    ax.axis("off")
    return fig


fs.show(draw_calorimeter, alt="Schematic of a constant-pressure calorimeter. A container of water with a known mass "
                              "and heat capacity stands inside insulating walls. A sample container sits in the "
                              "water, and an arrow labeled q P points from the sample into the water. A thermometer "
                              "and a stirrer reach into the water from above.")
```

````{admonition} **In-class activity**
:class: important

Dissolving ammonium nitrate ($\mathrm{NH_4NO_3}$) in water is endothermic—this is the chemistry behind instant cold packs. Suppose you dissolve $10\,\mathrm{g}$ of $\mathrm{NH_4NO_3}$ ($M = 80.04\,\mathrm{g/mol}$) in $100\,\mathrm{mL}$ of water in a Styrofoam cup and measure a temperature drop of $5.3\,\text{°C}$.

Estimate $\Delta H_{\mathrm{soln}}$ per mole of $\mathrm{NH_4NO_3}$. Assume $C_{P,\mathrm{soln}} \approx C_{P,\mathrm{water}} = 4.18\,\mathrm{J\,g^{-1}\,K^{-1}}$ and that the total solution mass is $\approx 110\,\mathrm{g}$.

```{admonition} Solution
:class: dropdown

1. **Heat absorbed by the dissolution process:**

   The solution *cools*, so heat flows *from* the water *into* the dissolving salt. The heat absorbed by the dissolution is the heat lost by the water:

   $$q_{\mathrm{soln}} = -q_{\mathrm{water}} = -m_{\mathrm{soln}}\,C_P\,\Delta T = -(110)(4.18)(-5.3) = +2{,}437\,\mathrm{J}.$$

   The positive sign confirms an endothermic process (the system absorbs heat from the water).

2. **Moles of $\mathrm{NH_4NO_3}$:**

   $$n = \frac{10\,\mathrm{g}}{80.04\,\mathrm{g/mol}} = 0.125\,\mathrm{mol}.$$

3. **Molar enthalpy of solution:**

   $$\Delta H_{\mathrm{soln}} = \frac{q_P}{n} = \frac{2{,}437}{0.125}\,\mathrm{J/mol} \approx +19.5\,\mathrm{kJ/mol}.$$

The literature value is about $+25.7\,\mathrm{kJ/mol}$; the discrepancy comes from heat that leaks in from the warmer surroundings (Styrofoam is not a perfect insulator) and the approximation that $C_{P,\mathrm{soln}} = C_{P,\mathrm{water}}$.
```
````

---

### Defining Common Enthalpy Changes

#### Standard Conditions

Standard conditions are defined as $P^\circ = 1\text{ bar}$. The superscript $^\circ$ denotes a standard-state quantity (see [Notation](../notation.md)). Reference databases—e.g., the [NIST-JANAF Thermochemical Tables](https://janaf.nist.gov/) and [Active Thermochemical Tables](https://atct.anl.gov/)—typically tabulate properties at $T = 25\,\text{°C}$ (298.15 K), though standard-state quantities can be defined at any temperature.

#### Standard Enthalpy of Formation

The **standard enthalpy of formation**, $\Delta_f H^\circ$, is the enthalpy change when **1 mole** of a compound in its standard state is formed from its constituent elements **in their standard states**.

- The standard state of an element is its most stable form at $P^\circ = 1\text{ bar}$ (and a specified $T$).
- By convention, $\Delta_f H^\circ$ for an element **in its standard state** is **zero**.

```{list-table} Standard States of the Elements
:header-rows: 1
:name: standard-state-elements

* - Standard State
  - Elements
* - Monatomic ideal gas
  - He, Ne, Ar, Kr, Xe, Rn
* - Homonuclear diatomic ideal gas
  - H, N, O, F, Cl
* - Liquid
  - Br, Hg
* - Solid (most stable crystal structure)
  - All other elements
```

The key principle is that the standard state is the **most thermodynamically stable form** at $P^\circ$. For most elements this is a solid (e.g., C as graphite, not diamond; Fe as bcc iron). The noble gases are monatomic, the common nonmetals form diatomics, and Br and Hg are liquids at 298 K. Phosphorus is the usual exception: thermochemical tables take white phosphorus as the reference state, although black phosphorus is more stable.

````{admonition} Detailed solid standard states by crystal structure
:class: dropdown

For reference, the specific crystal structures defining the standard state of each solid element (phosphorus, whose reference state is white phosphorus, is not listed):

```{list-table}
:header-rows: 1
:name: solid-standard-states

* - Crystal Structure
  - Elements
* - Body-centered cubic
  - Alkali metals (Li, Na, K, Rb, Cs), Ba, group 5 transition metals (V, Nb, Ta), group 6 transition metals (Cr, Mo, W), Mn, Fe, & Eu
* - Hexagonal
  - Be, Mg, group 3 transition metals (Sc, Y, Lu), group 4 transition metals (Ti, Zr, Hf), Tc, Re, Ru, Os, Co, group 12 transition metals (Zn, Cd), Tl, C (graphite), Se, Te, most lanthanides (La, Ce, Pr, Nd, Pm, Gd, Tb, Dy, Ho, Er, Tm)
* - Face-centered cubic
  - Ca, Sr, Rh, Ir, group 10 transition metals (Ni, Pd, Pt), group 11 transition metals (Cu, Ag, Au), Al, Si (diamond cubic), Ge (diamond cubic), Pb, Yb
* - Rhombohedral
  - B, As, Sb, Bi, Sm
* - Orthorhombic
  - Ga, S, I, U, and others
* - Body-centered tetragonal
  - In, Sn ($\beta$, white)
* - Simple cubic
  - Po
```
````

#### Standard Enthalpy of Reaction

The **standard enthalpy of reaction**, $\Delta_r H^\circ$, is the enthalpy change when a reaction is carried out under standard conditions. Mathematically:

```{math}
:label: standard_enthalpy_reaction
\Delta_r H^\circ
\;=\;\sum_{\text{products}} \nu_p \,H_p^\circ \;-\; \sum_{\text{reactants}} \nu_r \,H_r^\circ,
```

where $\nu_p$ and $\nu_r$ are the (positive) stoichiometric coefficients of products and reactants, respectively, and $H_p^\circ$ and $H_r^\circ$ are their standard molar enthalpies. Only differences of enthalpy can be measured, so the individual values of $H^\circ$ are not known, and Equation {eq}`standard_enthalpy_reaction` cannot be evaluated term by term.

`````{admonition} Hess's Law & Standard Enthalpies of Formation
:class: note

**Hess's Law**: Because enthalpy is a state function, enthalpy changes are **additive** and **path independent**. The enthalpy change for a reaction depends only on the initial and final states, not on the route taken between them.

As a consequence, $\Delta_r H^\circ$ can be computed from tabulated formation enthalpies:

```{math}
:label: hess-law-unsigned
\Delta_r H^\circ
\;=\; \sum_{p} \nu_p \,\Delta_f H_p^\circ
\;-\;\sum_{r} \nu_r \,\Delta_f H_r^\circ.
```

```{admonition} Alternative notation with signed stoichiometric coefficients
:class: dropdown

Some texts use *signed* stoichiometric coefficients $\nu_i$, defined as positive for products and negative for reactants. In that convention, Hess's Law takes the compact form

$$\Delta_r H^\circ = \sum_{i} \nu_i \,\Delta_f H_i^\circ.$$

Both forms are equivalent. For example, in the reaction $\mathrm{CH_4}(g) + 2\,\mathrm{O_2}(g) \to \mathrm{CO_2}(g) + 2\,\mathrm{H_2O}(g)$, the signed coefficients are $\nu_{\mathrm{CH_4}} = -1$, $\nu_{\mathrm{O_2}} = -2$, $\nu_{\mathrm{CO_2}} = +1$, $\nu_{\mathrm{H_2O}} = +2$.
```

#### Example: Steam–Methane Reforming

Consider the steam–methane reforming reaction:

```{math}
\mathrm{CH_4}(g) + \mathrm{H_2O}(g) \;\;\longrightarrow\;\; \mathrm{CO}(g) + 3\,\mathrm{H_2}(g).
```

Applying Equation {eq}`hess-law-unsigned`:

```{math}
\Delta_r H^\circ
= \Bigl[\Delta_f H^\circ(\mathrm{CO}) + 3\,\underbrace{\Delta_f H^\circ(\mathrm{H_2})}_{=\,0}\Bigr]
- \Bigl[\Delta_f H^\circ(\mathrm{CH_4}) + \Delta_f H^\circ(\mathrm{H_2O})\Bigr].
```

Note that $\Delta_f H^\circ(\mathrm{H_2}) = 0$ because $\mathrm{H_2}(g)$ is hydrogen's standard state. With $\Delta_f H^\circ(\mathrm{CO}) = -110.5$, $\Delta_f H^\circ(\mathrm{CH_4}) = -74.8$, and $\Delta_f H^\circ(\mathrm{H_2O}) = -241.8\ \mathrm{kJ\,mol^{-1}}$ for the gases at 298.15 K,

```{math}
\Delta_r H^\circ
= \bigl[(-110.5) + 3(0)\bigr] - \bigl[(-74.8) + (-241.8)\bigr]
= +206.1\ \mathrm{kJ\,mol^{-1}}.
```

Alternatively, we can see *why* this works by decomposing the reaction into formation steps:

````{dropdown} Breakdown into formation reactions
```{math}
\begin{aligned}
&\mathrm{CH_4}(g) \longrightarrow \mathrm{C}(s, \mathrm{graphite}) + 2 \mathrm{H_2}(g) &\qquad &\Delta_r H^{\circ} = -\Delta_f H^{\circ}(\mathrm{CH_4}) \\
&\mathrm{H_2O}(g) \longrightarrow \mathrm{H_2}(g) + \tfrac{1}{2} \mathrm{O_2}(g) &\qquad &\Delta_r H^{\circ} = -\Delta_f H^{\circ}(\mathrm{H_2O}) \\
&\mathrm{C}(s, \mathrm{graphite}) + \tfrac{1}{2} \mathrm{O_2}(g) \longrightarrow \mathrm{CO}(g) &\qquad &\Delta_r H^{\circ} = \Delta_f H^{\circ}(\mathrm{CO})
\end{aligned}
```

Summing gives the net reaction, and the enthalpy changes add:

```{math}
\Delta_r H^\circ
= -\Delta_f H^\circ(\mathrm{CH_4})
- \Delta_f H^\circ(\mathrm{H_2O})
+ \Delta_f H^\circ(\mathrm{CO}).
```
````

Either approach gives the same result, thanks to Hess's Law.
`````

[](#fig:enthalpy-levels) visualizes the path independence that underlies Hess's Law. The direct route (green arrow) and the stepwise route through the elements (two purple arrows) connect the same initial and final states, so they produce the same $\Delta_r H^\circ$.

% Try this: Start from liquid water: change `'H2O': -241.8` to `'H2O': -285.8`, the standard enthalpy of formation of $\mathrm{H_2O}(l)$. By how much does $\Delta_r H^\circ$ change, and which process accounts for the difference? (The label still reads $\mathrm{H_2O(g)}$.)

```{code-cell} ipython3
:label: fig:enthalpy-levels
:caption: Enthalpy-level diagram for the steam–methane reforming reaction, drawn to scale in $\mathrm{kJ\,mol^{-1}}$ with the elements in their standard states as the zero. The direct path (green arrow) and the stepwise path through the elements (two purple arrows) both connect reactants to products. Because $H$ is a state function, both paths give the same $\Delta_r H^\circ$. This is Hess's Law.
:tags: [hide-input]

# Standard enthalpies of formation at 298.15 K (kJ/mol)
dHf = {'CH4': -74.8, 'H2O': -241.8, 'CO': -110.5, 'H2': 0.0}

H_elements = 0.0
H_reactants = dHf['CH4'] + dHf['H2O']
H_products = dHf['CO'] + 3 * dHf['H2']
dH_rxn = H_products - H_reactants


def signed(value):
    """Number with its sign and one decimal, with a true minus sign."""
    return f"{value:+.1f}".replace("-", "−")


def draw_enthalpy_levels(c):
    fig, ax = plt.subplots(figsize=(fs.WIDTH, 4.1), constrained_layout=True)

    # Enthalpy levels
    ax.plot([0.0, 2.0], [H_reactants, H_reactants], color=c.ink, solid_capstyle="butt")
    ax.plot([3.0, 5.0], [H_elements, H_elements], color=c.muted, ls=(0, (3, 1.5)), solid_capstyle="butt")
    ax.plot([6.0, 8.0], [H_products, H_products], color=c.ink, solid_capstyle="butt")
    ax.text(1.0, H_reactants - 12, "Reactants\n" + r"$\mathrm{CH_4(g) + H_2O(g)}$", ha="center", va="top")
    ax.text(4.0, H_elements + 12,
            "Elements in standard states\n" + r"$\mathrm{C(s,\,graphite) + 3\,H_2(g) + \frac{1}{2}\,O_2(g)}$",
            ha="center", va="bottom", color=c.muted)
    ax.text(7.0, H_products - 12, "Products\n" + r"$\mathrm{CO(g) + 3\,H_2(g)}$", ha="center", va="top")

    # Arrows: each starts and ends on an enthalpy level
    fs.arrow(ax, (2.0, H_reactants), (3.0, H_elements), c.purple, lw=1.8)     # reactants -> elements
    fs.arrow(ax, (5.0, H_elements), (6.0, H_products), c.purple, lw=1.8)      # elements -> products
    fs.arrow(ax, (2.0, H_reactants), (6.0, H_products), c.green, lw=1.8)      # direct path

    ax.text(2.33, H_reactants / 2 + 25, r"$-\sum \nu_r\,\Delta_f H_r^\circ$" + f"\n$= {signed(-H_reactants)}$",
            ha="right", va="center", color=c.purple)
    ax.text(5.67, H_products / 2 + 5, r"$\sum \nu_p\,\Delta_f H_p^\circ$" + f"\n$= {signed(H_products)}$",
            ha="left", va="center", color=c.purple)
    ax.text(4.3, (H_reactants + H_products) / 2 - 22, r"$\Delta_r H^\circ$" + f"$ = {signed(dH_rxn)}$",
            ha="left", va="top", color=c.green)

    ax.set(ylabel=r"Enthalpy relative to the elements (kJ mol$^{-1}$)", xlim=(-0.3, 8.3), ylim=(-390, 70), xticks=[])
    ax.spines["bottom"].set_visible(False)
    return fig


fs.show(draw_enthalpy_levels,
        alt="Enthalpy-level diagram with three levels. The reactants, methane and water vapor, lie at minus 316.6 "
            "kJ per mole. The elements in their standard states lie at zero. The products, carbon monoxide and "
            "hydrogen, lie at minus 110.5 kJ per mole. One arrow runs from the reactants up to the elements, plus "
            "316.6, and a second from the elements down to the products, minus 110.5. A third arrow runs directly "
            "from the reactants up to the products, plus 206.1.")
```

---

### Synthesis: The First Law Toolkit So Far

Enthalpy completes our First Law toolkit for the most common laboratory conditions:

| Constraint | Relevant state function | Key relation |
| :---: | :---: | :---: |
| Constant $V$ (bomb calorimeter) | $U$ | $q_V = \Delta U$ |
| Constant $P$ (open beaker, coffee-cup calorimeter) | $H = U + PV$ | $q_P = \Delta H$ |

Hess's Law lets us combine tabulated formation data to predict $\Delta_r H^\circ$ for *any* reaction without measuring it directly. In [Section 5.3](../chapter-05/section-03.md), we will see how the temperature dependence of $\Delta H$ (via $C_P$) allows us to extrapolate thermochemical data to other temperatures—an essential tool for real-world applications where reactions rarely occur at exactly 298 K.

---

````{admonition} **Try it: Hess's Law practice**
:class: important

Compute $\Delta_r H^\circ$ for the combustion of methane:

$$\mathrm{CH_4}(g) + 2\,\mathrm{O_2}(g) \longrightarrow \mathrm{CO_2}(g) + 2\,\mathrm{H_2O}(g).$$

Use the following standard enthalpies of formation:

| Species | $\Delta_f H^\circ$ (kJ/mol) |
|:-------:|:---------------------------:|
| $\mathrm{CH_4}(g)$ | −74.8 |
| $\mathrm{CO_2}(g)$ | −393.5 |
| $\mathrm{H_2O}(g)$ | −241.8 |
| $\mathrm{O_2}(g)$ | $0$ (element in standard state) |

Set up the Hess's Law calculation before checking the solution.

```{admonition} Solution
:class: dropdown

Applying Equation {eq}`hess-law-unsigned`:

$$\Delta_r H^\circ
= \bigl[\Delta_f H^\circ(\mathrm{CO_2}) + 2\,\Delta_f H^\circ(\mathrm{H_2O})\bigr]
- \bigl[\Delta_f H^\circ(\mathrm{CH_4}) + 2\,\underbrace{\Delta_f H^\circ(\mathrm{O_2})}_{=\,0}\bigr].$$

Substituting numerical values:

$$\Delta_r H^\circ
= \bigl[(-393.5) + 2(-241.8)\bigr] - \bigl[(-74.8) + 0\bigr]$$

$$= (-393.5 - 483.6) - (-74.8)$$

$$= -877.1 + 74.8 = \boxed{-802.3\,\mathrm{kJ/mol}}.$$

The large negative value confirms that methane combustion is strongly exothermic, consistent with its use as a fuel. Note that $\Delta_f H^\circ(\mathrm{O_2}) = 0$ because $\mathrm{O_2}(g)$ is oxygen's standard state.
```
````

## Worked Example

### Heating at Constant Pressure

A sample of liquid water ($n = 1\,\mathrm{mol}$) has constant-pressure heat capacity $C_P=75.0\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (approximately constant over the range). Find $\Delta H$ when it is heated from $T_i=298\ \mathrm{K}$ to $T_f=350\ \mathrm{K}$ at constant pressure.

**Setup.** At constant pressure with $PV$-only work, $q_P = \Delta H$ (Eq. {eq}`qP-DeltaH`), and the enthalpy change is given by Eq. {eq}`DeltaH-Cp`:

```{math}
\Delta H = \int_{T_i}^{T_f} C_P\,dT \approx C_P\,(T_f - T_i).
```

**Calculation.**

```{math}
\Delta T = 350 - 298 = 52\ \mathrm{K},
```

```{math}
\Delta H \approx (75.0\,\mathrm{J\,mol^{-1}\,K^{-1}})(52\,\mathrm{K}) = 3.90\times10^{3}\ \mathrm{J\,mol^{-1}} = 3.90\ \mathrm{kJ\,mol^{-1}}.
```

**Result.** $\Delta H \approx 3.9\ \mathrm{kJ\,mol^{-1}}$ for this heating step at constant pressure. Because this is a constant-pressure process, this is also the heat absorbed: $q_P = 3.9\,\mathrm{kJ/mol}$.

```{admonition} Sanity check
:class: tip

The specific heat of water is $4.18\,\mathrm{J\,g^{-1}\,K^{-1}}$. For $1\,\mathrm{mol}$ of water ($18\,\mathrm{g}$), $C_P = 18 \times 4.18 \approx 75\,\mathrm{J\,mol^{-1}\,K^{-1}}$—consistent with the given value. Heating $18\,\mathrm{g}$ of water by $52\,\mathrm{K}$ should require about $18 \times 4.18 \times 52 \approx 3{,}900\,\mathrm{J}$. ✓
```

## Concept Checks

1. Why does adding the $PV$ term make $H$ especially convenient at constant pressure?
2. When is it *not* valid to identify $q_P$ with $\Delta H$? Give a specific physical example.
3. Why are formation enthalpies of elements in their standard states defined as zero?
4. How does Hess's Law justify computing reaction enthalpies without specifying a mechanism?
5. Is $C_P$ larger or smaller than $C_V$ for an ideal gas? Why? (*Hint:* think about what happens to the volume at constant $P$ when you add heat.)

## Key Takeaways

- Enthalpy is defined by $H=U+PV$ and satisfies $q_P=\Delta H$ for constant-pressure processes with $PV$-only work.
- Defining $H$ is an example of a general thermodynamic strategy: when the natural variables for a problem don't match your existing state function, define a new one. We will use this strategy again for the Helmholtz and Gibbs free energies.
- Heat capacities relate temperature changes to enthalpy changes via $\Delta H=\int C_P\,dT$.
- Standard enthalpy changes (formation, reaction) provide a consistent bookkeeping framework. The key convention is $\Delta_f H^\circ = 0$ for elements in their standard states.
- Hess's Law lets you build $\Delta_r H^\circ$ from tabulated formation enthalpies.
