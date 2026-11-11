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


# 2.1. Introduction to Statistical Mechanics

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Statistical mechanics bridges microscopic states and macroscopic thermodynamics by treating measurable quantities as ensemble averages. This section introduces expected values, defines the common ensembles (microcanonical, canonical, and grand canonical), and states the fundamental postulate for an isolated system.

```{mermaid}
%%{init: {"theme": "neutral", "look": "handDrawn", "layout": "elk"}}%%
flowchart LR
  %% Beginnings and ending
  CM([Classical mechanics])
  QM([Quantum mechanics])
  Th([Thermodynamics])

  %% Processes
  KT[[Kinetic theory]]
  SM[[Statistical mechanics]]

  %% Inputs/Outputs
  StateCM[/"<i>r<sup>N</sup></i>, <i>p<sup>N</sup></i>"/]
  StateQM[/"&Psi;(<i>r<sup>N</sup></i>)"/]
  
  %% Decision
  LimitCM{"Classical limit?"}

  subgraph Microscopic World
    QM
    CM
    StateCM
    StateQM
    LimitCM
  end
  
  subgraph Bridges
    KT
    SM
  end
  
  subgraph Macroscopic World
    Th
  end
  
  CM --> StateCM
  StateCM --> KT --> Th
  StateCM --> SM --> Th
  
  QM --> LimitCM
  LimitCM -- Yes --> CM
  LimitCM -- No --> StateQM --> SM
```

Learning objectives:

- Distinguish an arithmetic average (sample mean) from an ensemble (expected) average $\langle X\rangle$.
- Define microstates, macrostates, and an ensemble, and identify the constraints defining the microcanonical, canonical, and grand canonical ensembles.
- State the fundamental postulate of statistical mechanics for the microcanonical ensemble.
- Compute microcanonical probabilities $p_i=1/M$ and simple expected values from a discrete distribution.

## Core Ideas and Derivations

### Macroscopic Properties as Expected Values of Microscopic Properties

A core principle of statistical mechanics is that **macroscopic thermodynamic properties** are **statistical averages** (expected values) of microscopic properties.

#### Arithmetic Average vs. Expected Value

In basic statistics, the **arithmetic average** (sample mean) $\bar{X}$ of a dataset $X=\{X_1, X_2, \ldots, X_M\}$ is:

```{math}
\bar{X} = \frac{1}{M} \sum_{i=1}^M X_i.
```

````{admonition} Example of an Arithmetic Average
:class: tip
The arithmetic average of $X = \{ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4 \}$ is

```{math}
\bar{X} = \frac{1}{10} \left(1 + 2 + 2 + 3 + 3 + 3 + 4 + 4 + 4 + 4\right) = 3.
```
````

In **statistical mechanics**, we typically work with an **expected value**, which weights each outcome by its probability. For a discrete random variable $X$ that takes values $\{X_i\}$ with probabilities $\{p_i\}$, the expected value is

```{math}
:label: expected-value
\langle X \rangle = \sum_{i=1}^M X_i \, p_i,
```

where $p_i$ is the probability of the $i$-th outcome (microstate), and the sum runs over the microstates included in the ensemble.

````{admonition} Example of an Expected Value
:class: tip
If $X$ takes the values $\{1, 2, 3, 4\}$ with probabilities $p = \{0.1, 0.2, 0.3, 0.4\}$, then

```{math}
\langle X \rangle = 1 \times 0.1 + 2 \times 0.2 + 3 \times 0.3 + 4 \times 0.4 = 3.
```
````

````{admonition} Expected Value of the Number of Tails in 100 Fair Coin Flips
:class: tip
Let $X=1$ for tails and $X=0$ for heads. For one flip,

```{math}
\langle X \rangle_1 = 0 \times 0.5 + 1 \times 0.5 = 0.5.
```

By linearity of expectation, the expected number of tails in 100 flips is

```{math}
\langle N_{\text{tails}}\rangle = 100\,\langle X\rangle_1 = 100 \times 0.5 = 50.
```
````

```{list-table} Statistical Variables and Their Definitions
:header-rows: 1
:name: statistical-variables

* - Symbol
  - Meaning
* - $M$
  - Number of possible microstates (outcomes)
* - $i$
  - Index of a microstate
* - $X_i$
  - Value of a microscopic property in microstate $i$
* - $\langle X \rangle$
  - Expected (ensemble) average of $X$
* - $p_i$
  - Probability of finding the system in microstate $i$
```

In **thermodynamics**, typical choices of $X$ include the **internal energy**, **enthalpy**, or other measurable properties. We will compute such quantities by specifying the probabilities $\{p_i\}$ appropriate to the ensemble of interest.

### Ensembles of Microstates

```{glossary}
Ensemble
: The set of microstates consistent with specified macroscopic constraints.

Microcanonical ensemble
: Microstates are sampled at fixed $\left(N, V, E\right)$.

Canonical ensemble
: Microstates are sampled at fixed $\left(N, V, T\right)$.

Grand canonical ensemble
: Microstates are sampled at fixed $\left(\mu, V, T\right)$.
```

### Probability of a Microstate in the Microcanonical Ensemble

For an **isolated system**—one that exchanges neither energy nor matter with its surroundings—the appropriate statistical description is the **microcanonical ensemble**.

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from myst_nb import glue

# Helper function to plot a system
def plot_system(ax, title, annotations, boundary_color='b'):
    box = mpatches.FancyBboxPatch((0, 0), 1, 1, boxstyle='roundtooth', ec=boundary_color, fc='w')
    ax.add_patch(box)
    ax.set_title(title, fontsize=14)
    ax.text(0.5, 0.5, 'System', ha='center', va='center', fontsize=12)
    ax.text(0.5, -0.65, 'Surroundings', ha='center', va='center', fontsize=12)
    ax.text(0.5, 1.3, 'Boundary', ha='center', va='bottom', fontsize=12, color=boundary_color)
    for annotation in annotations:
        if "arrowprops" in annotation:  # Arrow annotations
            ax.annotate('', **annotation)
        else:  # Text annotations
            ax.text(**annotation)
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 2)
    ax.set_aspect('equal')
    ax.axis('off')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plot_system(ax, "", [])

plt.show()
plt.close(fig)
```

An isolated system exchanges neither energy nor matter with its surroundings.

#### Fundamental Postulate of Statistical Mechanics

```{admonition} Fundamental Postulate
:class: tip
**For an isolated system (microcanonical ensemble), each accessible microstate is equally probable.**
```

Therefore, the probability of finding the system in the $i$-th microstate is

```{math}
:label: microcanonical-probability
p_i = \frac{1}{M},
```

where $M$ is the total number of accessible microstates compatible with $\left(N, V, E\right)$.

````{admonition} Isolated Spin-Up Electron in an f Orbital
:class: tip
Consider an isolated spin-up electron in an $f$ orbital. The possible quantum states have magnetic quantum numbers $m_\ell = \{-3, -2, -1, 0, 1, 2, 3\}$, giving $M=7$. By the fundamental postulate, the probability of each microstate is

```{math}
p_i = \frac{1}{7}.
```
````

## Worked Example

### Expected value vs. arithmetic mean (die rolls)

A fair die is rolled $n=300$ times. Let $X_j=1$ if the $j$-th outcome is a six and $X_j=0$ otherwise.

1. **One roll**

   ```{math}
   \langle X_j\rangle = 1\cdot \frac{1}{6} + 0\cdot \frac{5}{6} = \frac{1}{6}.
   ```

2. **Linearity of expectation**

   For independent, identical trials, the expected total number of sixes is the sum of the individual expectations:

   ```{math}
   \langle N_{\text{sixes}}\rangle = \sum_{j=1}^{300}\langle X_j\rangle = 300\left(\frac{1}{6}\right)=50.
   ```

**Result.** The expected number of sixes in 300 rolls is $50$.

## Concept Checks

1. In what sense is “temperature” absent from the microcanonical ensemble description?
2. Why is the expected value $\langle X\rangle$ more informative than the arithmetic mean of a single dataset when predicting thermodynamic properties?
3. What physical meaning does “accessible microstate” carry in the fundamental postulate?
4. How would you modify the probability assignment if only a subset of the nominal $M$ microstates were accessible?

## Key Takeaways

- Macroscopic properties are computed as **expected values** over microstates.
- An **ensemble** is a probability model for microstates consistent with specified macroscopic constraints.
- In the microcanonical ensemble, each accessible microstate has equal probability: $p_i=1/M$.
- Choosing the right ensemble is choosing the right constraints: $(N,V,E)$, $(N,V,T)$, or $(\mu,V,T)$.
