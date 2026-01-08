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


<!-- TODO: Add content from lecture-26-notes-jc.pdf -->

# 4.1. Entropy

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Entropy completes the thermodynamic description of spontaneity by accounting for the dispersal of energy and matter. This section motivates entropy through examples where enthalpy alone fails, defines entropy as a state function via reversible heat, and derives the fundamental relation $dU=T\,dS-P\,dV$ for simple compressible systems.

In Module 3, we established the foundations of thermodynamics and explored the concept of enthalpy. We will learn that **exothermicity** (a negative enthalpy change) often favors spontaneity but does not fully determine it. The missing piece is **entropy**, a state function that provides crucial insight into the direction of spontaneous change. In this sub-module, we build on those ideas by formally introducing entropy, examining its definition, and understanding its role in the **arrow of time** for spontaneous processes.

Learning objectives:

- Explain why exothermicity ($\Delta H<0$) does not by itself guarantee spontaneity.
- Define entropy $S$ as a state function using $dS=\delta q_{\mathrm{rev}}/T$.
- Use the First Law for a reversible $PV$-only process to derive $dU=T\,dS-P\,dV$.
- Compute entropy changes by integrating $\delta q_{\mathrm{rev}}/T$ along a convenient reversible path.

## Core Ideas and Derivations

### Exothermicity and Spontaneity

Exothermic reactions are often—but not always—spontaneous. An exothermic process releases heat to its surroundings, corresponding to a negative enthalpy change, $\Delta H < 0$. While this heat release can drive a process forward, enthalpy alone does not guarantee spontaneity; we must also consider entropy.

````{margin}
```{note}
$\Delta H < 0$ is **exothermic**  
$\Delta H > 0$ is **endothermic**
```
````

#### An Example: Formation of Liquid Water

Consider the formation of two moles of liquid water from hydrogen and oxygen gases at 298.15 K and 1 bar:

```{math}
\ce{2H2(g) + O2(g) -> 2H2O(l)}
```

The standard enthalpy change for this reaction is $-571.66$ kJ for two moles of H$_2$O(l), indicating an exothermic process. The reaction releases heat:

```{admonition} Exothermicity in Chemical Terms
:class: tip
In chemical terms, **exothermicity** is often associated with the formation of stronger chemical bonds in the products than those in the reactants. In this example, forming $\ce{H-O}$ bonds in water releases more energy than what was required to break $\ce{H-H}$ and $\ce{O=O}$ bonds.
```

It is fortunate that water formation is exothermic, as this makes the reaction more likely to proceed spontaneously under standard conditions.

#### A Counter Example: Mixing of Two Ideal Gases

If we mix two ideal gases, such as Ar and Kr, their intermolecular interactions are negligible (ideal-gas behavior), so $\Delta H = 0$. Nonetheless, mixing occurs spontaneously once the partition is removed.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV
from labellines import labelLines
from myst_nb import glue

rng = np.random.default_rng(8251991)

fig, axs = plt.subplot_mosaic([[0, 1]], figsize=(8, 4), constrained_layout=True, sharex=True, sharey=True)

# Before
axs[0].set_title('Before Opening the Stopper')
x_Ar = rng.uniform(-4, 4, 6)
y_Ar = rng.uniform(2, 5, 6)
axs[0].scatter(x_Ar, y_Ar, color='r', label='Ar', s=100)

x_Kr = rng.uniform(-4, 4, 6)
y_Kr = rng.uniform(-5, -2, 6)
axs[0].scatter(x_Kr, y_Kr, color='b', label='Kr', s=100)

# Ar region
axs[0].plot([1, 1], [0, 1], color='k')
axs[0].plot([1, 5], [1, 1], color='k')
axs[0].plot([5, 5], [1, 6], color='k')
axs[0].plot([5, -5], [6, 6], color='k')
axs[0].plot([-5, -5], [6, 1], color='k')
axs[0].plot([-5, -1], [1, 1], color='k')
axs[0].plot([-1, -1], [1, 0], color='k')

# Kr region
axs[0].plot([1, 1], [0, -1], color='k')
axs[0].plot([1, 5], [-1, -1], color='k')
axs[0].plot([5, 5], [-1, -6], color='k')
axs[0].plot([5, -5], [-6, -6], color='k')
axs[0].plot([-5, -5], [-6, -1], color='k')
axs[0].plot([-5, -1], [-1, -1], color='k')
axs[0].plot([-1, -1], [-1, 0], color='k')

# Closed stopper
axs[0].plot([-1.5, 1.5], [0, 0], color='m', lw=2)
axs[0].text(2, 0, "Stopper (closed)", fontsize=10, ha='left', va='center', color='m')
axs[0].set_aspect('equal')
axs[0].axis('off')

# After
axs[1].set_title('After Opening the Stopper')
x_Ar = rng.uniform(-4, 4, 6)
y_Ar_1 = rng.uniform(2, 5, 3)
y_Ar_2 = rng.uniform(-5, -2, 3)
y_Ar = np.concatenate((y_Ar_1, y_Ar_2))
axs[1].scatter(x_Ar, y_Ar, color='r', label='Ar', s=100)

x_Kr = rng.uniform(-4, 4, 6)
y_Kr_1 = rng.uniform(2, 5, 3)
y_Kr_2 = rng.uniform(-5, -2, 3)
y_Kr = np.concatenate((y_Kr_1, y_Kr_2))
axs[1].scatter(x_Kr, y_Kr, color='b', label='Kr', s=100)

axs[1].plot([1, 1], [0, 1], color='k')
axs[1].plot([1, 5], [1, 1], color='k')
axs[1].plot([5, 5], [1, 6], color='k')
axs[1].plot([5, -5], [6, 6], color='k')
axs[1].plot([-5, -5], [6, 1], color='k')
axs[1].plot([-5, -1], [1, 1], color='k')
axs[1].plot([-1, -1], [1, 0], color='k')
axs[1].plot([1, 1], [0, -1], color='k')
axs[1].plot([1, 5], [-1, -1], color='k')
axs[1].plot([5, 5], [-1, -6], color='k')
axs[1].plot([5, -5], [-6, -6], color='k')
axs[1].plot([-5, -5], [-6, -1], color='k')
axs[1].plot([-5, -1], [-1, -1], color='k')
axs[1].plot([-1, -1], [-1, 0], color='k')

# Open stopper
axs[1].plot([-5, -2], [0, 0], color='m', lw=2)
axs[1].text(2, 0, "Stopper (open)", fontsize=10, ha='left', va='center', color='m')
axs[1].axis('off')
axs[1].set_aspect('equal')
axs[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

fig.suptitle(r'$\Delta H = 0$')

plt.show()
plt.close(fig)
```

Mixing of two ideal gases before and after opening the stopper that separates them.

Even though $\Delta H = 0$, the process is still **spontaneous**. If you were shown the “before” and “after” snapshots without additional labels, you would know the natural direction of mixing. This spontaneous behavior underscores that **enthalpy** alone cannot capture whether a process will occur without additional driving forces. That driving force is related to **entropy**.

### Definition of Entropy

```{glossary}
Entropy
: A state function denoted by $S$, quantifying the degree of *dispersal* or *spread* of energy and matter, thereby predicting the direction of spontaneous change.
```

#### First Law for a Reversible Process

Recall the differential form of the **First Law of Thermodynamics** for a reversible process in which only $PV$ work is done on an **ideal monatomic gas** ($PV = Nk_{\mathrm{B}}T$ and $C_V = \tfrac{3}{2}Nk_{\mathrm{B}}$):

```{math}
\begin{align*}
\delta q_{\text{rev}} &= dU - \delta w_{\text{rev}} \\
&= dU + P\,dV \\
&= C_V\,dT + \left[ \left( \frac{\partial U}{\partial V} \right)_T + P \right] dV \\
&= \frac{3}{2} N k_{\mathrm{B}}\,dT + \frac{N k_{\mathrm{B}} T}{V}\,dV.
\end{align*}
```

Because $\delta q_{\text{rev}}$ depends on the path, it is an **inexact** differential.

##### Dividing by $T$ to Obtain an Exact Differential

Observe that if we divide $\delta q_{\mathrm{rev}}$ by $T$, the resulting expression becomes exact:

```{math}
\frac{\delta q_{\text{rev}}}{T} \;=\; \frac{3}{2}\,\frac{N k_{\mathrm{B}}}{T} \, dT \;+\; \frac{N k_{\mathrm{B}}}{V}\,dV.
```

Hence, we define **entropy** $S$ through the exact differential:

```{math}
dS \;=\; \frac{\delta q_{\text{rev}}}{T}.
```

The change in entropy between two states $A$ and $B$ is then:

```{math}
\Delta S \;=\; \int_{A}^{B} \frac{\delta q_{\text{rev}}}{T}.
```

#### Fundamental Thermodynamic Relation

From the definition of entropy, we obtain:

```{math}
\delta q_{\text{rev}} \;=\; T\,dS.
```

For a simple closed system with only $PV$ work,

```{math}
\delta q_{\text{rev}} = dU + P\,dV \;=\; T\,dS \quad\Longrightarrow\quad dU = T\,dS - P\,dV.
```

Here, $U$ is naturally expressed as a function of the **extensive** variables $S$ and $V$. In systems with additional types of work (e.g., electrical, surface, magnetic), the more general relation becomes:

```{math}
dU \;=\; T\,dS \;+\; \sum_{i} \vec{F}_{i}\cdot d\vec{x}_{i}.
```

```{admonition} Natural Variables
:class: tip
A state function’s **natural variables** are those that appear as independent variables in its total differential. For internal energy $U$, the natural variables are $S$ and $V$ for a simple system with only $PV$ work.
```

## Worked Example

### Entropy change for isothermal ideal-gas expansion

For a reversible isothermal expansion of an ideal gas from $V_1$ to $V_2$ at temperature $T$,

```{math}
\Delta S = \int_{1}^{2}\frac{\delta q_{\mathrm{rev}}}{T}.
```

For a reversible isotherm of an ideal gas, $\delta q_{\mathrm{rev}} = nRT\,dV/V$, so

```{math}
\Delta S = nR\int_{V_1}^{V_2}\frac{dV}{V}=nR\ln\left(\frac{V_2}{V_1}\right).
```

Take $n=1\ \mathrm{mol}$ and $V_2=2V_1$:

```{math}
\Delta S = R\ln 2 = (8.314)\ln 2 = 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** Doubling the volume reversibly at constant $T$ increases the entropy by $R\ln 2$ per mole.

## Concept Checks

1. In the ideal-gas mixing example, why can $\Delta H=0$ but the process still be spontaneous?
2. Why must $\Delta S$ be computed along a reversible path even for an irreversible process between the same states?
3. What are the natural variables of $U$ for a simple compressible system, and how do they appear in $dU$?
4. How does the sign of $\Delta S$ relate to “directionality” (arrow of time) in spontaneous processes?

## Key Takeaways

- Spontaneity depends on more than enthalpy; entropy accounts for dispersal and multiplicity of accessible states.
- Entropy is defined by $dS=\delta q_{\mathrm{rev}}/T$ and is a **state function**.
- For reversible $PV$-only processes, $dU=T\,dS-P\,dV$ is the fundamental thermodynamic relation.
- Entropy changes can be computed by integrating $\delta q_{\mathrm{rev}}/T$ along any convenient reversible path.
