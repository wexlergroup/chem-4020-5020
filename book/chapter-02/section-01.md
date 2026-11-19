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
  enumerator: 2.1.%s
---


# 2.1. Introduction to Statistical Mechanics

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 2.1](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-02/section-01.ipynb), in which you can change and rerun the Python cells of this section in your browser.

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=2&amp;here=2.1"
  title="Course map, opened at Module 2"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-2.png
:class: hidden
:alt: Course map with Module 2 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6. Module 2 is the microscopic strand itself. Below the map, one card for each of Sections 2.1 to 2.8 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

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

In basic statistics, the **arithmetic average** (sample mean) $\bar{X}$ of a dataset of $n$ values, $X=\{X_1, X_2, \ldots, X_n\}$, is:

```{math}
\bar{X} = \frac{1}{n} \sum_{j=1}^n X_j.
```

````{admonition} Example of an Arithmetic Average
:class: tip
The arithmetic average of the $n = 10$ values $X = \{ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4 \}$ is

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
Microstate
: A complete specification of the microscopic state of a system: the positions and momenta of all particles in classical mechanics, or the quantum state of the system in quantum mechanics.

Macrostate
: The thermodynamic state of a system, specified by a small number of macroscopic variables such as $\left(N, V, E\right)$. Many microstates are usually consistent with one macrostate.

Ensemble
: The set of microstates consistent with specified macroscopic constraints. An ensemble can be pictured as a very large collection of copies of the system, each in one of these microstates; $p_i$ is then the fraction of copies in microstate $i$.

Microcanonical ensemble
: Microstates are sampled at fixed $\left(N, V, E\right)$.

Canonical ensemble
: Microstates are sampled at fixed $\left(N, V, T\right)$.

Grand canonical ensemble
: Microstates are sampled at fixed $\left(\mu, V, T\right)$.
```

[](#fig:ensembles) shows the system and its surroundings in each of the three ensembles.

```{code-cell} ipython3
:tags: [hide-input]
:label: fig:ensembles
:caption: The three ensembles. (a) Microcanonical ensemble: an isolated system, with $N$, $V$, and $E$ fixed. (b) Canonical ensemble: a closed system that exchanges energy with a heat bath at temperature $T$, with $N$, $V$, and $T$ fixed. (c) Grand canonical ensemble: an open system that exchanges energy and matter with a reservoir at temperature $T$ and chemical potential $\mu$, with $\mu$, $V$, and $T$ fixed.

import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs


def draw_ensemble(ax, c, fixed, bath=None, energy=False, matter=False):
    """Sketch of a system in its surroundings for one ensemble, with the exchanges that its boundary allows."""
    if bath:
        ax.add_patch(FancyBboxPatch((0.05, 0.17), 0.90, 0.72, boxstyle="round,pad=0.03,rounding_size=0.06",
                                    fc=c.fill, ec="none"))
    ax.text(0.5, 0.83, bath or "Surroundings", ha="center", va="center", color=c.note)
    ax.add_patch(FancyBboxPatch((0.36, 0.34), 0.28, 0.32, boxstyle="round,pad=0.03,rounding_size=0.06",
                                fc=c.ground, ec=c.ink, lw=1.6))
    ax.text(0.5, 0.59, "System", ha="center", va="center")
    if energy:
        ax.annotate("", xy=(0.09, 0.42), xytext=(0.45, 0.42),
                    arrowprops=dict(arrowstyle="<|-|>", color=c.vermilion, lw=1.6, mutation_scale=12))
        ax.text(0.20, 0.46, "Energy", ha="center", va="bottom", color=c.vermilion)
    if matter:
        ax.annotate("", xy=(0.55, 0.42), xytext=(0.91, 0.42),
                    arrowprops=dict(arrowstyle="<|-|>", color=c.blue, lw=1.6, mutation_scale=12))
        ax.text(0.80, 0.46, "Matter", ha="center", va="bottom", color=c.blue)
    ax.text(0.5, 0.05, fixed, ha="center", va="center")
    ax.set(xlim=(0, 1), ylim=(0, 0.95), aspect="equal")
    ax.axis("off")


def draw(c):
    fig, axes = plt.subplots(1, 3, figsize=(fs.WIDTH, 2.55), constrained_layout=True)
    draw_ensemble(axes[0], c, "$N$, $V$, $E$ fixed")
    draw_ensemble(axes[1], c, "$N$, $V$, $T$ fixed", bath="Heat bath at $T$", energy=True)
    draw_ensemble(axes[2], c, r"$\mu$, $V$, $T$ fixed", bath=r"Reservoir at $T$ and $\mu$", energy=True, matter=True)
    for ax, letter, title in zip(axes, "abc", ["Microcanonical", "Canonical", "Grand canonical"]):
        fs.panel(ax, letter, title)
    return fig


fs.show(draw, alt="Three sketches of a system, drawn as a box. In panel a, the microcanonical ensemble, nothing crosses "
                  "the boundary, and N, V, and E are fixed. In panel b, the canonical ensemble, the system sits in a "
                  "heat bath at temperature T, a two-headed arrow labeled energy crosses the boundary, and N, V, and T "
                  "are fixed. In panel c, the grand canonical ensemble, the system sits in a reservoir at temperature "
                  "T and chemical potential mu, two-headed arrows labeled energy and matter cross the boundary, "
                  "and mu, V, and T are fixed.")
```

### Probability of a Microstate in the Microcanonical Ensemble

For an **isolated system**—one that exchanges neither energy nor matter with its surroundings—the appropriate statistical description is the **microcanonical ensemble** (panel (a) of [](#fig:ensembles)).

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

   The expected total number of sixes is the sum of the individual expected values:

   ```{math}
   \langle N_{\text{sixes}}\rangle = \sum_{j=1}^{300}\langle X_j\rangle = 300\left(\frac{1}{6}\right)=50.
   ```

3. **Arithmetic means from simulated rolls**

   In five simulated sets of 300 rolls, the number of sixes was 56, 42, 49, 43, and 41. The arithmetic mean of $X_j$ in each set, $\bar{X} = N_{\text{sixes}}/300$, was therefore 0.187, 0.140, 0.163, 0.143, and 0.137, compared with $\langle X_j\rangle = 1/6 \approx 0.167$.

**Result.** The expected number of sixes in 300 rolls is $50$; the five simulated sets gave between 41 and 56.

## Concept Checks

1. In what sense is “temperature” absent from the microcanonical ensemble description?
2. Why is the expected value $\langle X\rangle$ more informative than the arithmetic mean of a single dataset when predicting thermodynamic properties?
3. What physical meaning does “accessible microstate” carry in the fundamental postulate?
4. How would you modify the probability assignment if only a subset of the nominal $M$ microstates were accessible?

## Key Takeaways

- Macroscopic properties are computed as **expected values** over microstates.
- An **ensemble** is the set of microstates consistent with specified macroscopic constraints, together with their probabilities $p_i$.
- In the microcanonical ensemble, each accessible microstate has equal probability: $p_i=1/M$.
- Choosing the right ensemble is choosing the right constraints: $(N,V,E)$, $(N,V,T)$, or $(\mu,V,T)$.
