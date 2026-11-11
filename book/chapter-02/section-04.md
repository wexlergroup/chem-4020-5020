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


# 2.4. Molecular Partition Functions

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Molecular partition functions extend the canonical framework from a single particle to many particles, and then to molecules with internal degrees of freedom. This section shows how distinguishability affects counting (and therefore $Q$), introduces the $1/N!$ correction for identical particles, and motivates factorization into translational, rotational, vibrational, and electronic contributions.

In Sections 2.2 and 2.3, we derived the canonical partition function and ensemble averages for single-particle systems. Here we extend these ideas to closed systems of **many identical, independent particles with internal degrees of freedom**. We also introduce **indistinguishability**, which is central to quantum statistical mechanics.

```{note}
**Assumption (Single Occupancy):**
Throughout the microstate-counting examples in this section, we assume that each one-particle microstate can be occupied by **at most one** particle. Under this assumption, the total number of ways to place $N$ distinguishable particles in $M$ one-particle microstates is $\tfrac{M!}{(M-N)!}$ rather than $M^N$. If multiple occupancy were allowed, the counting would differ accordingly (e.g., $M^N$ for unlimited occupancy).

**Flag:** The factorized derivation $Q=q^N$ below uses independent sums over one-particle states (allowing repeated indices). If you intend to enforce single occupancy in the partition-function derivations as well, those sums would need an $i\neq j\neq\cdots$ restriction and $Q$ would not factorize exactly.
```

Learning objectives:

* Count many-particle microstates for distinguishable vs. indistinguishable particles under the stated single-occupancy assumption.
* Derive the factorized form $Q=q^N$ for independent (distinguishable) particles and $Q=q^N/N!$ for identical indistinguishable particles.
* Explain why indistinguishability is essential for obtaining extensive thermodynamics (avoiding the Gibbs paradox).
* Write the molecular partition function as $q=q_{\mathrm{trans}}q_{\mathrm{rot}}q_{\mathrm{vib}}q_{\mathrm{elec}}$ within the Born–Oppenheimer approximation.

## Core Ideas and Derivations

### Partition Function for Distinguishable Particles

````{admonition} Two Distinguishable Particles in a Four-State System
:class: dropdown

The $N$-particle microstates of a system of two distinguishable particles in a four-state system are listed in the table below:

```{list-table}
:header-rows: 1

* - $N$-Particle Microstate
  - One-Particle Microstate 1
  - One-Particle Microstate 2
  - One-Particle Microstate 3
  - One-Particle Microstate 4
* - I
  - a
  - b
  - 
  - 
* - II
  - a
  - 
  - b
  - 
* - III
  - a
  - 
  - 
  - b
* - IV
  - 
  - a
  - b
  - 
* - V
  - 
  - a
  - 
  - b
* - VI
  - 
  - 
  - a
  - b
* - VII
  - b
  - a
  - 
  - 
* - VIII
  - b
  - 
  - a
  - 
* - IX
  - b
  - 
  - 
  - a
* - X
  - 
  - b
  - a
  - 
* - XI
  - 
  - b
  - 
  - a
* - XII
  - 
  - 
  - b
  - a
```

There are $12$ $N$-particle microstates in this example.

Here, we treat each of the four one-particle states as an exclusive "slot," so two distinguishable particles must occupy different states, yielding $4\times 3=12$ arrangements rather than $4\times 4=16$.

---

More generally, for $N$ distinguishable particles in a system with $M$ one-particle microstates (single occupancy), the number of $N$-particle microstates is the number of permutations of $N$ particles in $M$ states:

```{math}
\frac{M!}{(M - N)!}
```

This assumes no more than one particle per one-particle state.
````

The energy $E$ of an $N$-particle microstate is the sum of the one-particle energies for the occupied states:

```{math}
E_l = \varepsilon_i^a + \varepsilon_j^b + \varepsilon_k^c + \cdots
```

Here $i, j, k, \ldots$ index the one-particle microstates occupied by particles $a, b, c, \ldots$, and $l$ indexes the corresponding $N$-particle microstate of the system.

The canonical partition function $Q$ for $N$ distinguishable, independent particles is given by the sum over all possible microstates:

```{math}
\begin{aligned}
Q &= \sum_l e^{-\beta E_l} \\
  &= \sum_{i = 1}^M \sum_{j = 1}^M \sum_{k = 1}^M \cdots
     \exp\!\left[-\beta\left(\varepsilon_i^a + \varepsilon_j^b + \varepsilon_k^c + \cdots\right)\right] \\
  &= \left(\sum_{i = 1}^M e^{-\beta \varepsilon_i^a}\right)
     \left(\sum_{j = 1}^M e^{-\beta \varepsilon_j^b}\right)
     \left(\sum_{k = 1}^M e^{-\beta \varepsilon_k^c}\right)\cdots \\
  &= \boxed{q_a q_b q_c \cdots}
\end{aligned}
```

where $q_a$, $q_b$, and $q_c$ are the canonical partition functions for the particles $a$, $b$, and $c$, respectively.

If the particles share the same one-particle spectrum so that $q_a=q_b=q_c=\cdots=q$, then

```{math}
Q = q^N.
```

### Partition Function for Indistinguishable Particles

````{admonition} Two Indistinguishable Particles in a Four-State System
:class: dropdown

The $N$-particle microstates of a system of two indistinguishable particles in a four-state system are listed in the table below:

```{list-table}
:header-rows: 1

* - $N$-Particle Microstate
  - One-Particle Microstate 1
  - One-Particle Microstate 2
  - One-Particle Microstate 3
  - One-Particle Microstate 4
* - I
  - x
  - x
  - 
  - 
* - II
  - x
  - 
  - x
  - 
* - III
  - x
  - 
  - 
  - x
* - IV
  - 
  - x
  - x
  - 
* - V
  - 
  - x
  - 
  - x
* - VI
  - 
  - 
  - x
  - x
```

There are $6$ $N$-particle microstates in this example.

---

In general, for $N$ indistinguishable particles in a system with $M$ one-particle microstates (single occupancy), the number of $N$-particle microstates is the number of combinations of $N$ particles in $M$ states:

```{math}
\frac{1}{N!} \frac{M!}{(M - N)!} = \binom{M}{N}.
```

This assumes no more than one particle per one-particle state.
````

If the particles are indistinguishable and identical, we can write the canonical partition function as

```{math}
Q = \frac{q^N}{N!}.
```

### Partition Function for Molecules

Within the Born–Oppenheimer approximation, the total energy of a molecule is given by the sum of the translational, rotational, vibrational, and electronic energies:

```{math}
\varepsilon_{\lambda}
= \varepsilon_i^{\mathrm{trans}}
+ \varepsilon_j^{\mathrm{rot}}
+ \varepsilon_k^{\mathrm{vib}}
+ \varepsilon_m^{\mathrm{elec}}
```

where $i$, $j$, $k$, and $m$ are the indices of the (one-degree-of-freedom) microstates of the translational, rotational, vibrational, and electronic energies, respectively, and $\lambda$ is the index of the (many-degree-of-freedom) microstate of the molecule.

The canonical partition function $q$ for a molecule is given by the sum over all possible microstates:

```{math}
\begin{aligned}
q
&= \sum_{\lambda} e^{-\beta \varepsilon_{\lambda}} \\
&= \sum_{i} \sum_{j} \sum_{k} \sum_{m}
\exp\!\left[-\beta\left(\varepsilon_i^{\mathrm{trans}} + \varepsilon_j^{\mathrm{rot}} + \varepsilon_k^{\mathrm{vib}} + \varepsilon_m^{\mathrm{elec}}\right)\right] \\
&= \sum_{i} e^{-\beta \varepsilon_i^{\mathrm{trans}}}
   \sum_{j} e^{-\beta \varepsilon_j^{\mathrm{rot}}}
   \sum_{k} e^{-\beta \varepsilon_k^{\mathrm{vib}}}
   \sum_{m} e^{-\beta \varepsilon_m^{\mathrm{elec}}} \\
&= q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}.
\end{aligned}
```

where $q_{\mathrm{trans}}$, $q_{\mathrm{rot}}$, $q_{\mathrm{vib}}$, and $q_{\mathrm{elec}}$ are the canonical partition functions for the translational, rotational, vibrational, and electronic energies, respectively.

## Worked Example

### Counting microstates and the $1/N!$ factor (single occupancy)

Two particles occupy $M=4$ one-particle microstates, with **at most one particle per microstate**.

1. **Distinguishable particles**

   The number of distinct arrangements (permutations) is

   ```{math}
   W_{\mathrm{dist}}=\frac{M!}{(M-N)!}=\frac{4!}{2!}=12.
   ```

2. **Indistinguishable particles**

   Exchanging labels does not create a new microstate, so

   ```{math}
   W_{\mathrm{indist}}=\frac{1}{N!}\frac{M!}{(M-N)!}=\frac{1}{2}\cdot 12=6.
   ```

3. **Connection to partition functions**

   If each allowed arrangement has the same energy $E_0$, then the canonical partition function is proportional to the count:

   ```{math}
   Q_{\mathrm{dist}} = W_{\mathrm{dist}}\,e^{-\beta E_0},\qquad
   Q_{\mathrm{indist}} = W_{\mathrm{indist}}\,e^{-\beta E_0}=\frac{Q_{\mathrm{dist}}}{2!}.
   ```

**Result.** The $1/N!$ factor corrects overcounting for identical particles and changes the thermodynamics derived from $\ln Q$.

## Concept Checks

1. What physical experiment would fail to distinguish two microstates that differ only by swapping particle labels?
2. Why does $Q$ factorize for independent particles but not for interacting particles?
3. What assumption about occupancy is built into the combinatorial formulas in this section?
4. How would the counting change if multiple occupancy were allowed?

## Key Takeaways

* Counting microstates depends on whether particles are **distinguishable**.
* For identical particles, $Q=q^N/N!$ corrects overcounting and yields consistent extensive properties.
* Molecular energies often separate into translation/rotation/vibration/electronic parts, motivating $q=q_{\mathrm{trans}}q_{\mathrm{rot}}q_{\mathrm{vib}}q_{\mathrm{elec}}$.
* Partition functions are bookkeeping devices: changes in counting change $\ln Q$ and therefore the thermodynamics.
