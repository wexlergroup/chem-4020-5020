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


# 3.3. Enthalpy

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Enthalpy $H=U+PV$ is the natural energy-like potential for constant-pressure processes, which are common in chemistry and calorimetry. This section defines enthalpy from the First Law at constant pressure, relates $q_P$ to $\Delta H$, and introduces standard enthalpy changes used in thermochemistry.

Enthalpy $H$ is a thermodynamic potential particularly convenient for processes at constant pressure, which is common in chemistry.

---

Learning objectives:

- Derive the definition $H=U+PV$ and show that $q_P=\Delta H$ for constant-pressure processes with only $PV$ work.
- Define heat capacity at constant pressure $C_P$ and use it to compute enthalpy changes over temperature intervals.
- Explain standard states and the meaning of standard enthalpy of formation and reaction.
- Apply Hess’s Law to compute $\Delta H_{\mathrm{rxn}}^\circ$ from formation enthalpies.

## Core Ideas and Derivations

### Defining Enthalpy

Consider the first-law expression when $T$ and $P$ are the independent variables:

```{math}
\delta q 
= \left[\left(\frac{\partial U}{\partial T}\right)_P + P \left(\frac{\partial V}{\partial T}\right)_P\right] dT 
\;+\; \left[\left(\frac{\partial U}{\partial P}\right)_T 
+ P \left(\frac{\partial V}{\partial P}\right)_T\right] dP.
```

At constant pressure, this becomes:

```{math}
:label: differential_first_law_constant_pressure
\delta q_P 
= \left[\left(\frac{\partial U}{\partial T}\right)_P + P \left(\frac{\partial V}{\partial T}\right)_P\right] dT.
```

Define the heat capacity at constant pressure, $C_P$:

```{math}
C_P 
\;=\; \left(\frac{\partial U}{\partial T}\right)_P 
\;+\; P \left(\frac{\partial V}{\partial T}\right)_P.
```

We seek a state function $H$ such that

```{math}
:label: differential_enthalpy
\delta q_P 
= \left(\frac{\partial H}{\partial T}\right)_P \, dT.
```

This motivates the definition of **enthalpy**:

```{math}
:label: enthalpy
H \;=\; U + PV.
```

````{admonition} Demonstrating That $H$ Yields $\delta q_P$
:class: dropdown

By substituting $H = U + PV$, Equation {eq}`enthalpy`, into $\delta q_P = \left(\partial H / \partial T\right)_P dT$, Equation {eq}`differential_enthalpy`:

```{math}
\begin{align*}
\delta q_P
&= \left[\frac{\partial}{\partial T} \bigl(U + PV\bigr)\right]_P dT \\[6pt]
&= \left[ \left(\frac{\partial U}{\partial T}\right)_P 
  + P \left(\frac{\partial V}{\partial T}\right)_P 
  + V \underbrace{\left(\frac{\partial P}{\partial T}\right)_P}_{=0}\right] dT \\[6pt]
&= \left[\left(\frac{\partial U}{\partial T}\right)_P
  \;+\; P \left(\frac{\partial V}{\partial T}\right)_P\right] dT.
\end{align*}
```

This matches Equation {eq}`differential_first_law_constant_pressure`, thus confirming that 
$\delta q_P = \mathrm{d}H$ at constant pressure.
````

### Measuring Enthalpy and Enthalpy Changes

At constant pressure, the heat absorbed or released by a process is equal to the change in enthalpy of the system:

```{math}
q_P 
= \Delta H 
= H(T_f) - H(T_i).
```

If $C_P$ is approximately constant over the temperature range $\Delta T = T_f - T_i$, then

```{math}
\Delta H 
= \int_{T_i}^{T_f} C_P(T)\,dT 
\;\longrightarrow\; C_P \,\Delta T.
```

A typical way to measure $\Delta H$ experimentally is via **calorimetry**. Below is a schematic of a simple calorimeter:

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV
from labellines import labelLines
from myst_nb import glue

fig = plt.figure(figsize=(4, 4))

# Container
plt.plot([0, 2], [0, 0], 'k-', lw=2, label='Container base')
plt.plot([0, 2], [2, 2], 'k-', lw=2)
plt.plot([0, 0], [0, 2], 'k-', lw=2)
plt.plot([2, 2], [0, 2], 'k-', lw=2)

# Water
plt.fill_between([0, 2], 0, 2, color='blue', alpha=0.3, label='Water')

# Thermometer
plt.plot([0.25, 0.25], [1.5, 2.5], 'r-', lw=2, label='Thermometer')

# Stirrer
plt.plot([0.5, 0.5], [0.5, 2.5], 'k-', lw=2, label='Stirrer')
plt.plot([0.25, 0.5, 0.75], [0.4, 0.5, 0.4], 'k-', lw=2)
plt.plot([0.25, 0.5, 0.75], [0.6, 0.5, 0.6], 'k-', lw=2, alpha=0.5)

# Sample container
plt.plot([1.5, 1.9], [0.1, 0.1], 'C7-', lw=2, label='Sample Container')
plt.plot([1.5, 1.9], [0.5, 0.5], 'C7-', lw=2)
plt.plot([1.5, 1.5], [0.1, 0.5], 'C7-', lw=2)
plt.plot([1.9, 1.9], [0.1, 0.5], 'C7-', lw=2)

plt.fill_between([1.5, 1.9], 0.1, 0.5, color='white')

# Sample inside container
np.random.seed(42)  # For reproducibility
x_sample = np.random.uniform(1.6, 1.8, 20)
y_sample = np.random.uniform(0.2, 0.3, 20)
plt.scatter(x_sample, y_sample, color='orange', alpha=0.5, label='Sample')

# Ignition source
plt.plot([1.7, 1.7], [0.3, 2.5], 'm-', lw=2, label='Ignition Source')

plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0., frameon=False)
plt.axis('off')
plt.tight_layout()
plt.show()
plt.close(fig)
```

A simplified calorimeter. The process occurs in the sample container, which transfers heat to or from the surrounding water at constant pressure. The thermometer and stirrer ensure accurate, uniform temperature readings.

````{admonition} Build Your Own Calorimeter
:class: tip

**Question**: How might you build a rudimentary calorimeter with a Styrofoam cup?

Think about a simple chemical reaction in your kitchen.

```{admonition} **Hint**
:class: dropdown

Consider a reaction that involves dissolving a solid in water. You could use a thermometer to measure the temperature change.
```
````

### Defining Common Enthalpy Changes

#### Standard Conditions

Standard conditions are defined as $P^\circ = 1\text{ bar}$. Reference databases---e.g., the [NIST-JANAF Thermochemical Tables](https://doi.org/10.18434/T42S31) and [Active Thermochemical Tables](https://atct.anl.gov/)---often document properties at $T = 25^\circ C$ (298.15 K).

#### Standard Enthalpy of Formation

The **standard enthalpy of formation**, $\Delta H_f^\circ$, is the enthalpy change when **1 mole** of a compound is formed from its constituent elements **in their standard states**.

- The standard state of an element is its most stable form at $P^\circ = 1\text{ bar}$ (and a specified $T$).
- By convention, $\Delta H_f^\circ$ for an element **in its standard state** is **zero**.

```{list-table} Standard State and Elements
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
* - Solid
  - All other elements
```

```{list-table} Solid Standard States
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
  - Ga, P (black), S, I, U, and others
* - Body-centered tetragonal
  - In, Sn ($\beta$, white)
* - Simple cubic
  - Po
```

#### Standard Enthalpy of Reaction

The **standard enthalpy of reaction**, $\Delta H_{\mathrm{rxn}}^\circ$, is the enthalpy change when a reaction is carried out under standard conditions. Mathematically:

```{math}
:label: standard_enthalpy_reaction
\Delta H_{\mathrm{rxn}}^\circ
\;=\;\sum_{\text{products}} \nu_p H_p^\circ \;-\; \sum_{\text{reactants}} \nu_r H_r^\circ,
```

where $\nu_p$ and $\nu_r$ are stoichiometric coefficients of products and reactants, respectively.

`````{admonition} Hess’s Law & Standard Enthalpies of Formation
:class: note

**Hess’s Law**: Enthalpy changes are **additive** and **path independent**. Consequently,  
```{math}
\Delta H_{\mathrm{rxn}}^\circ 
\;=\; \sum_{p} \nu_p \,\Delta H_{f,p}^\circ 
\;-\;\sum_{r} \nu_r \,\Delta H_{f,r}^\circ.
```

For example, consider the steam–methane reforming reaction:

```{math}
\text{CH}_4(g) + \text{H}_2\text{O}(g) \;\;\longrightarrow\;\; \text{CO}(g) + 3\,\text{H}_2(g).
```

Directly:

```{math}
\Delta H_{\mathrm{rxn}}^\circ 
= \nu_{\mathrm{CO}}\,H_{\mathrm{CO}}^\circ
+ 3\,\nu_{\mathrm{H}_2}\,H_{\mathrm{H}_2}^\circ
- \nu_{\mathrm{CH}_4}\,H_{\mathrm{CH}_4}^\circ 
- \nu_{\mathrm{H}_2\mathrm{O}}\,H_{\mathrm{H}_2\mathrm{O}}^\circ.
```

Alternatively, break each species into formation (or reverse formation) reactions from the elemental forms, then sum their enthalpy changes:

````{dropdown} Breakdown of the Steam–Methane Reforming Reaction
```{math}
\begin{align*}
\text{CH}_4(g) &\longrightarrow \text{C}(s, \text{graphite}) + 2 \text{H}_2(g) &\quad \Delta H^{\circ}_{\text{rxn}} &= -\Delta H^{\circ}_{f, \text{CH}_4} \\
\text{H}_2\text{O}(g) &\longrightarrow \text{H}_2(g) + \tfrac{1}{2} \text{O}_2(g) &\quad \Delta H^{\circ}_{\text{rxn}} &= -\Delta H^{\circ}_{f, \text{H}_2\text{O}} \\
\text{C}(s, \text{graphite}) + \tfrac{1}{2} \text{O}_2(g) &\longrightarrow \text{CO}(g) &\quad \Delta H^{\circ}_{\text{rxn}} &= \Delta H^{\circ}_{f, \text{CO}}
\end{align*}
```
````

```{math}
\Delta H_{\mathrm{rxn}}^\circ 
= \bigl[-\Delta H_{f}^\circ(\text{CH}_4)\bigr] 
+ \bigl[-\Delta H_{f}^\circ(\text{H}_2\mathrm{O})\bigr] 
+ \Delta H_{f}^\circ(\text{CO})
+ 3\,\Delta H_{f}^\circ(\text{H}_2).
```

Either approach gives the same result, thanks to Hess’s Law.
`````

## Worked Example

### Heating at constant pressure

A sample has constant-pressure heat capacity $C_P=75.0\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (approximately constant over the range). Find $\Delta H$ when it is heated from $T_i=298\ \mathrm{K}$ to $T_f=350\ \mathrm{K}$ at constant pressure.

At constant pressure,

```{math}
\Delta H=\int_{T_i}^{T_f} C_P\,dT \approx C_P(T_f-T_i).
```

Compute $\Delta T$:

```{math}
\Delta T = 350-298 = 52\ \mathrm{K}.
```

Then

```{math}
\Delta H \approx (75.0)(52)\ \mathrm{J\,mol^{-1}} = 3.90\times10^{3}\ \mathrm{J\,mol^{-1}}=3.90\ \mathrm{kJ\,mol^{-1}}.
```

**Result.** $\Delta H \approx 3.9\ \mathrm{kJ\,mol^{-1}}$ for this heating step at constant pressure.

## Concept Checks

1. Why does adding the $PV$ term make $H$ especially convenient at constant pressure?
2. When is it *not* valid to identify $q_P$ with $\Delta H$?
3. Why are formation enthalpies of elements in their standard states defined as zero?
4. How does Hess’s Law justify computing reaction enthalpies without specifying a mechanism?

## Key Takeaways

- Enthalpy is defined by $H=U+PV$ and satisfies $q_P=\Delta H$ for constant-pressure $PV$-only work.
- Heat capacities relate temperature changes to enthalpy changes via $\Delta H=\int C_P\,dT$.
- Standard enthalpy changes (formation, reaction) provide a consistent bookkeeping framework.
- Hess’s Law lets you build $\Delta H_{\mathrm{rxn}}^\circ$ from tabulated formation enthalpies.
