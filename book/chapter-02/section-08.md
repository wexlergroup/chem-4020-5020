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
  enumerator: 2.8.%s
---

# 2.8. Molecular Statistical Mechanics

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Molecular statistical mechanics often builds the total molecular partition function as a product of translational, rotational, vibrational, and electronic factors ([Section 2.4](section-04.md)). The translational factor is derived in [Section 2.5](section-05.md), the vibrational factor of a single oscillator in [Section 2.6](section-06.md), and the rotational factor of a linear molecule in [Section 2.7](section-07.md). This section collects high-temperature rigid-rotor formulas for different molecular shapes, reviews vibrational and electronic partition functions, and summarizes standard textbook approximations used to connect molecular constants to thermodynamic properties. We use *rotational temperatures* $\Theta_{\mathrm{rot}}$ and *vibrational temperatures* $\Theta_{\mathrm{vib}}$, which package molecular constants into energy scales expressed in Kelvin (i.e., divided by Boltzmann’s constant $k_{\mathrm{B}}$).

Learning objectives:

* Identify how molecular shape (linear, spherical top, symmetric top, asymmetric top) affects the rotational partition function.
* Use the rotational symmetry number $\sigma$ and rotational temperatures $\Theta_{\mathrm{rot}}$ to write high-$T$ expressions for $q_{\mathrm{rot}}$.
* Write the harmonic-oscillator vibrational partition function for diatomic and polyatomic molecules in terms of $\Theta_{\mathrm{vib}}$.
* Interpret the electronic partition function and justify the ground-state approximation $q_{\mathrm{elec}}\approx g_1$ when excited states are high in energy.

## Core Ideas and Derivations

### Rotational Symmetry

A nonlinear rigid molecule can rotate about three perpendicular axes that pass through its center of mass. For one choice of these axes, called the principal axes $a$, $b$, and $c$, the classical rotational energy is

```{math}
E = \frac{L_a^2}{2I_A} + \frac{L_b^2}{2I_B} + \frac{L_c^2}{2I_C},
```

where $L_a$, $L_b$, and $L_c$ are the components of the angular momentum along the principal axes, and $I_A$, $I_B$, and $I_C$ are the moments of inertia about these axes, which are called the principal moments of inertia. The corresponding expression for a linear molecule in [Section 2.7](section-07.md) has two terms and one moment of inertia. As in Section 2.7, each moment of inertia defines a rotational temperature:

```{math}
\Theta_{\mathrm{rot},i} = \frac{\hbar^2}{2 k_{\mathrm{B}} I_i}, \qquad i = A, B, C.
```

Nonlinear molecules are classified by their principal moments of inertia:

* In a **spherical top**, all three are equal, so the molecule has a single rotational temperature $\Theta_{\mathrm{rot}}$. Examples are $\mathrm{CH_4}$ and $\mathrm{SF_6}$.
* In a **symmetric top**, two are equal. Examples are $\mathrm{NH_3}$ and $\mathrm{C_6H_6}$.
* In an **asymmetric top**, all three are different. An example is $\mathrm{H_2O}$.

[](#fig:tops) shows one molecule of each kind with its principal axes.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/molecular-tops.html
:label: fig:tops
:width: 100%
:title: Rotatable models of methane, ammonia, and water with their principal axes
:placeholder: ../_static/molecular-tops/molecular-tops.png

Principal axes $a$, $b$, and $c$ of (a) a spherical top, CH$_4$, (b) a symmetric top, NH$_3$, and (c) an asymmetric top, H$_2$O. The axes pass through the center of mass, and the relation below each molecule compares its principal moments of inertia. In CH$_4$, any three perpendicular axes through the carbon atom are principal axes. In NH$_3$, the $c$ axis is the threefold axis. In H$_2$O, the $b$ axis is the twofold axis and the $c$ axis is perpendicular to the plane of the molecule.
```

The *rotational partition function* $q_{\mathrm{rot}}$ depends on:

1. The rotational constants (encoded as rotational temperatures $\Theta_{\mathrm{rot},A}$, $\Theta_{\mathrm{rot},B}$, $\Theta_{\mathrm{rot},C}$, etc.).
2. The rotational symmetry number $\sigma$.
3. The temperature $T$.

The symmetry number $\sigma$ is the number of orientations of the molecule that can be reached by rotating it as a rigid body and that cannot be distinguished from one another, because they differ only by an exchange of identical nuclei. The count includes the original orientation. For a linear molecule, $\sigma$ is 1 or 2 (see the tip on the symmetry number in [Section 2.7](section-07.md)). For $\mathrm{NH_3}$, rotations by 0°, 120°, and 240° about the axis that passes through the nitrogen atom and is perpendicular to the plane of the three hydrogen atoms give indistinguishable orientations, so $\sigma=3$.

Below is a summary table of approximate expressions for rigid rotors in the high-temperature limit $\bigl(T \gg \Theta_{\mathrm{rot}}\bigr)$. The expressions for linear molecules are derived in [Section 2.7](section-07.md). Here, $\Theta_{\mathrm{rot}}$ and $T$ are both expressed in Kelvin, so ratios like $T/\Theta_{\mathrm{rot}}$ are dimensionless. In the symmetric-top expressions below, we take $\Theta_{\mathrm{rot},A}=\Theta_{\mathrm{rot},B}$.

```{list-table} Rotational Partition Functions of Rigid Rotors
:header-rows: 1
:name: rot_partition_functions

* - Molecular Class
  - Symmetry Number ($\sigma$)
  - Rotational Partition Function ($q_{\mathrm{rot}}$)
  - Examples ($\sigma$)
* - **Heteronuclear Diatomic (Linear)**
  - 1
  - $\frac{T}{\Theta_{\mathrm{rot}}}$
  - CO, NO, HF, HCl, etc.
* - **Unsymmetric Polyatomic (Linear)**
  - 1
  - $\frac{T}{\Theta_{\mathrm{rot}}}$
  - HCN, N₂O, FCN, etc.
* - **Homonuclear Diatomic (Linear)**
  - 2
  - $\frac{T}{2\,\Theta_{\mathrm{rot}}}$
  - H₂, O₂, N₂, etc.
* - **Symmetric Polyatomic (Linear)**
  - 2
  - $\frac{T}{2\,\Theta_{\mathrm{rot}}}$
  - CO₂, CS₂, XeF₂
* - **Spherical Top**
  - 12 ($T_d$) or 24 ($O_h$)
  - $\frac{\sqrt{\pi}}{\sigma}\,\biggl(\frac{T}{\Theta_{\mathrm{rot}}}\biggr)^{3/2}$
  - CH₄, SiH₄, SF₆, etc.
* - **Symmetric Top**
  - Varies with symmetry
  - $\frac{\sqrt{\pi}}{\sigma}\,\sqrt{\frac{T^3}{\Theta_{\mathrm{rot},A}^2\,\Theta_{\mathrm{rot},C}}}$
  - C₆H₆ ($\sigma=12$), NH₃ ($\sigma=3$), XeF₄ ($\sigma=8$), etc.
* - **Asymmetric Top**
  - Varies with symmetry
  - $\frac{\sqrt{\pi}}{\sigma}\,\sqrt{\frac{T^3}{\Theta_{\mathrm{rot},A}\,\Theta_{\mathrm{rot},B}\,\Theta_{\mathrm{rot},C}}}$
  - H₂O ($\sigma=2$), NO₂ ($\sigma=2$), SO₂ ($\sigma=2$), most large molecules (often $\sigma=1$)
```

```{admonition} Rotational energy and heat capacity of a nonlinear molecule
:class: important
Use the high-temperature $q_{\mathrm{rot}}$ of a nonlinear molecule and the derivative in [Section 2.3](section-03.md) to obtain $U_{\mathrm{rot}}$ and $C_V^{(\mathrm{rot})}$ per molecule. Compare the results with those for the linear rotor in [Section 2.7](section-07.md) and with the prediction of the equipartition theorem ([Section 2.6](section-06.md)) for the classical rotational energy above.
```

### Nonlinear Rigid Rotors

{numref}`rotational-temperatures-nonlinear` lists the rotational constants, rotational temperatures, and symmetry numbers of one molecule of each kind. The rotational temperatures are calculated from the rotational constants $\tilde{A}$, $\tilde{B}$, and $\tilde{C}$ as $\Theta_{\mathrm{rot},A} = hc\tilde{A}/k_{\mathrm{B}}$ (and likewise for $B$ and $C$), with $hc/k_{\mathrm{B}} = 1.43878\ \mathrm{cm\,K}$, as in [Section 2.7](section-07.md). The rotational constants of $\mathrm{CH_4}$ and $\mathrm{H_2O}$ are from the [NIST Computational Chemistry Comparison and Benchmark Database](https://cccbdb.nist.gov/), and those of $\mathrm{NH_3}$ are from S. Urban et al., *J. Mol. Spectrosc.* **101**, 1 (1983).

```{list-table} Rotational constants, rotational temperatures, and symmetry numbers of three nonlinear molecules
:header-rows: 1
:name: rotational-temperatures-nonlinear

* - Molecule
  - Kind of top
  - $\tilde{A}$, $\tilde{B}$, $\tilde{C}$ ($\mathrm{cm^{-1}}$)
  - $\Theta_{\mathrm{rot},A}$, $\Theta_{\mathrm{rot},B}$, $\Theta_{\mathrm{rot},C}$ (K)
  - $\sigma$
* - $\mathrm{CH_4}$
  - Spherical
  - 5.241, 5.241, 5.241
  - 7.541, 7.541, 7.541
  - 12
* - $\mathrm{NH_3}$
  - Symmetric
  - 9.947, 9.947, 6.227
  - 14.31, 14.31, 8.959
  - 3
* - $\mathrm{H_2O}$
  - Asymmetric
  - 27.88, 14.51, 9.285
  - 40.11, 20.88, 13.36
  - 2
```

### Vibrations of Polyatomic Molecules

For a diatomic molecule approximated as a *harmonic oscillator* with characteristic vibrational temperature $\Theta_{\mathrm{vib}}$, the vibrational partition function is ([Section 2.6](section-06.md))

```{math}
q_{\mathrm{vib}} \;=\; \frac{e^{-\Theta_{\mathrm{vib}}/(2T)}}{\,1 - e^{-\Theta_{\mathrm{vib}}/T}\,}.
```

The vibrations of a polyatomic molecule are described by normal modes. A **normal mode** is a vibration in which all of the atoms oscillate at the same frequency. In the harmonic approximation, each normal mode behaves as a one-dimensional harmonic oscillator with its own vibrational temperature, and the normal modes are independent of one another.

The number of normal modes follows from the number of atoms. In this subsection, $N$ denotes the number of atoms in one molecule. The positions of the atoms are specified by $3N$ coordinates. Three combinations of these coordinates describe the translation of the center of mass, and three more describe the rotation of a nonlinear molecule (two for a linear molecule; see [Section 2.7](section-07.md)). The remaining combinations describe vibrations, so the number of normal modes is $\alpha = 3N - 6$ for a nonlinear molecule and $\alpha = 3N - 5$ for a linear molecule. For example, $\mathrm{H_2O}$ has three normal modes, and $\mathrm{CO_2}$ has four.

The vibrational partition function is a product of the contributions of the normal modes:

```{math}
q_{\mathrm{vib}} \;=\; \prod_{j=1}^\alpha \; \frac{\,e^{-\Theta_{\mathrm{vib},j}/(2T)}\,}{\,1 - e^{-\Theta_{\mathrm{vib},j}/T}\,}.
```

Here, $\Theta_{\mathrm{vib},j}$ is the vibrational temperature of the $j$-th normal mode. Compare this product with the partition function of a system of independent, distinguishable particles in [Section 2.4](section-04.md) and with the Einstein model of a crystal in [Section 2.6](section-06.md).

### Electronic Partition Function

The *electronic partition function* is

```{math}
q_{\mathrm{elec}} \;=\; \sum_{i=1}^{\infty} g_i \, e^{-\beta\,\varepsilon_i}
\;\;\underset{\varepsilon_1 = 0}{\longrightarrow}\;\;
g_1 \;+\; g_2\,e^{-\beta\,(\varepsilon_2 - \varepsilon_1)} \;+\;\cdots,
```

where $g_i$ is the degeneracy of the $i$-th electronic level, $\varepsilon_i$ is the energy of that level, and $\beta = 1/(k_{\mathrm{B}} T)$. The energies are measured from the ground level, so $\varepsilon_1 = 0$. If every term after the first is negligible, then $q_{\mathrm{elec}}\approx g_1$.

The size of the second term depends on the temperature and on the energy of the first excited level. Two examples at $T = 300\ \mathrm{K}$, where $k_{\mathrm{B}}T = 0.0259\ \mathrm{eV}$, show the range:

* For a level with $\varepsilon_2 = 1\ \mathrm{eV}$, $e^{-\beta\varepsilon_2} = e^{-38.7} = 1.6\times10^{-17}$.
* For $\mathrm{NO}$, $g_1 = g_2 = 2$ and $\varepsilon_2/k_{\mathrm{B}} = 172.4\ \mathrm{K}$, so $e^{-\beta\varepsilon_2} = e^{-172.4/300} = 0.563$ and $q_{\mathrm{elec}} \approx 2 + 2(0.563) = 3.13$.

{numref}`electronic-levels` lists $g_1$ for several atoms and molecules, together with the energy of the first excited electronic level, given as a wavenumber, $\varepsilon_2/hc$, and as a temperature, $\varepsilon_2/k_{\mathrm{B}}$. The values are from the [NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database) (atoms) and the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/) (molecules).

```{list-table} Degeneracy of the ground electronic level and energy of the first excited electronic level of selected atoms and molecules
:header-rows: 1
:name: electronic-levels

* - Species
  - $g_1$
  - $\varepsilon_2/hc$ ($\mathrm{cm^{-1}}$)
  - $\varepsilon_2/k_{\mathrm{B}}$ (K)
* - $\mathrm{H}$
  - 2
  - 82 258.9
  - $1.184\times10^{5}$
* - $\mathrm{O}$
  - 5
  - 158.27
  - 227.7
* - $\mathrm{N_2}$
  - 1
  - 50 203.6
  - $7.223\times10^{4}$
* - $\mathrm{O_2}$
  - 3
  - 7918.1
  - $1.139\times10^{4}$
* - $\mathrm{NO}$
  - 2
  - 119.82
  - 172.4
```

```{admonition} Validity of $q_{\mathrm{elec}}\approx g_1$
:class: important
Use the two examples above to state a condition, in terms of $\varepsilon_2$ and $k_{\mathrm{B}}T$, under which $q_{\mathrm{elec}}\approx g_1$. For which species in {numref}`electronic-levels` does the condition hold at 300 K?
```

## Summary

The table below collects the main partition function formulas for atoms and molecules under typical textbook approximations (high-$T$ rotation, harmonic vibration, and negligible population of excited electronic states). In the table, $N$ is the number of atoms in the molecule, and "none" means that an atom has no rotational or vibrational factor.

```{list-table} Partition functions for various systems
:header-rows: 1
:name: tab:partfunc

* - **System**
  - $q_{\mathrm{trans}}$
  - $q_{\mathrm{rot}}$
  - $q_{\mathrm{vib}}$
  - $q_{\mathrm{elec}}$
* - Atom
  - $\frac{V}{\Lambda^3}$
  - none
  - none
  - $g_1$
* - Diatomic molecule
  - $\frac{V}{\Lambda^3}$
  - $\frac{T}{\sigma\,\Theta_{\mathrm{rot}}}$
  - $\frac{e^{-\Theta_{\mathrm{vib}}/(2T)}}{\,1 - e^{-\Theta_{\mathrm{vib}}/T}\,}$
  - $g_1$
* - Linear polyatomic molecule
  - $\frac{V}{\Lambda^3}$
  - $\frac{T}{\sigma\,\Theta_{\mathrm{rot}}}$
  - $\prod_{j=1}^{3N-5}\,\frac{e^{-\Theta_{\mathrm{vib},j}/(2T)}}{1-e^{-\Theta_{\mathrm{vib},j}/T}}$
  - $g_1$
* - Nonlinear molecule, spherical top
  - $\frac{V}{\Lambda^3}$
  - $\frac{\sqrt{\pi}}{\sigma}\,\bigl(\tfrac{T}{\Theta_{\mathrm{rot}}}\bigr)^{3/2}$
  - $\prod_{j=1}^{3N-6}\,\frac{e^{-\Theta_{\mathrm{vib},j}/(2T)}}{1-e^{-\Theta_{\mathrm{vib},j}/T}}$
  - $g_1$
* - Nonlinear molecule, symmetric top
  - $\frac{V}{\Lambda^3}$
  - $\frac{\sqrt{\pi}}{\sigma}\,\sqrt{\frac{T^3}{\Theta_{\mathrm{rot},A}^2\,\Theta_{\mathrm{rot},C}}}$
  - $\prod_{j=1}^{3N-6}\,\frac{e^{-\Theta_{\mathrm{vib},j}/(2T)}}{1-e^{-\Theta_{\mathrm{vib},j}/T}}$
  - $g_1$
* - Nonlinear molecule, asymmetric top
  - $\frac{V}{\Lambda^3}$
  - $\frac{\sqrt{\pi}}{\sigma}\,\sqrt{\frac{T^3}{\Theta_{\mathrm{rot},A}\,\Theta_{\mathrm{rot},B}\,\Theta_{\mathrm{rot},C}}}$
  - $\prod_{j=1}^{3N-6}\,\frac{e^{-\Theta_{\mathrm{vib},j}/(2T)}}{1-e^{-\Theta_{\mathrm{vib},j}/T}}$
  - $g_1$
```

The translational partition function $q_{\mathrm{trans}} = \frac{V}{\Lambda^3}$ (with $\Lambda$ the thermal de Broglie wavelength; [Section 2.5](section-05.md)) applies to all gas-phase systems. For many molecules, the total partition function is well approximated as the product

```{math}
q \;\approx\; q_{\mathrm{trans}}\;q_{\mathrm{rot}}\;q_{\mathrm{vib}}\;q_{\mathrm{elec}},
```

with each factor evaluated using the formulas above. The product holds when the energy of a molecule is a sum of translational, rotational, vibrational, and electronic parts ([Section 2.4](section-04.md)).

### Characteristic Temperatures and the Heat Capacity

Because $q$ is a product of factors, $\ln q$ is a sum, and the heat capacity of a gas is a sum of translational, rotational, vibrational, and electronic contributions. The rotational and vibrational temperatures mark the temperatures near which the rotational and vibrational contributions appear. [](#fig:heat-capacity-steps) shows the molar heat capacity at constant volume, $C_{V,m}$, of a gas of $\mathrm{CO}$ molecules in the rigid-rotor and harmonic-oscillator models, with $\Theta_{\mathrm{rot}} = 2.779\ \mathrm{K}$ ([Section 2.7](section-07.md)) and $\Theta_{\mathrm{vib}} = 3122\ \mathrm{K}$ ([Section 2.6](section-06.md)).

Translation contributes $\tfrac{3}{2}R$ at all temperatures in the figure ([Section 2.5](section-05.md)). Rotation adds $R$ when $T$ exceeds $\Theta_{\mathrm{rot}}$, and vibration adds up to $R$ more when $T$ approaches and exceeds $\Theta_{\mathrm{vib}}$. At $300\ \mathrm{K}$, which lies far above $\Theta_{\mathrm{rot}}$ and far below $\Theta_{\mathrm{vib}}$, $C_{V,m} = 2.50\,R$. The figure leaves out the electronic contribution, which is negligible for $\mathrm{CO}$ below several thousand kelvin.

```{code-cell} ipython3
:label: fig:heat-capacity-steps
:caption: Molar heat capacity at constant volume of a gas of $\mathrm{CO}$ molecules as a function of temperature on a logarithmic axis, calculated as the sum of a translational contribution ($\tfrac{3}{2}R$), the rotational contribution of a rigid rotor with $\Theta_{\mathrm{rot}} = 2.779\ \mathrm{K}$ (sum over the levels $J$), and the vibrational contribution of a harmonic oscillator with $\Theta_{\mathrm{vib}} = 3122\ \mathrm{K}$. The dashed lines mark the two characteristic temperatures, and the open point marks $300\ \mathrm{K}$. The models are applied at all temperatures, although $\mathrm{CO}$ condenses at $82\ \mathrm{K}$ at a pressure of 1 atm, and at the highest temperatures anharmonic vibration, excited electronic levels, and dissociation change the heat capacity of the real gas.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

theta_rot, theta_vib = 2.779, 3122.0             # K, for CO (Sections 2.7 and 2.6)
T = np.logspace(np.log10(0.3), np.log10(4e4), 700)

# Rotation: C/R = (<e^2> - <e>^2) / T^2 with e_J = J (J + 1) theta_rot in K and the weights (2J + 1) exp(-e_J / T)
J = np.arange(0, 1500)
e_J = J * (J + 1) * theta_rot
weight = (2 * J + 1) * np.exp(-e_J / T[:, None])
mean = (weight * e_J).sum(axis=1) / weight.sum(axis=1)
mean_square = (weight * e_J**2).sum(axis=1) / weight.sum(axis=1)
C_rot = (mean_square - mean**2) / T**2

# Vibration: C/R = x^2 e^(-x) / (1 - e^(-x))^2 with x = theta_vib / T
x = theta_vib / T
C_vib = x**2 * np.exp(-x) / (1 - np.exp(-x)) ** 2

C_trans = np.full_like(T, 1.5)
C_total = C_trans + C_rot + C_vib


def draw_motion(ax, c, kind, color, position):
    """Small sketch of a diatomic molecule with arrows for one kind of motion, at a position in axes fractions."""
    icon = ax.inset_axes([position[0], position[1], 0.095, 0.145])
    angle = np.radians(35 if kind == "rotation" else 0)
    dx, dy = 0.42 * np.cos(angle), 0.42 * np.sin(angle)
    icon.plot([-dx, dx], [-dy, dy], color=c.ink, lw=1.6, solid_capstyle="round", zorder=1)
    for sign in (-1, 1):
        icon.add_patch(Circle((sign * dx, sign * dy), 0.24, fc=color, ec="none", zorder=2))
    arrows = dict(arrowstyle="-|>", mutation_scale=8, lw=1.3, color=c.ink, shrinkA=0, shrinkB=0, zorder=3)
    if kind == "translation":
        icon.add_patch(FancyArrowPatch((0.8, 0), (1.45, 0), **arrows))
    elif kind == "rotation":
        icon.add_patch(FancyArrowPatch((0.95, -0.35), (0.35, 0.9), connectionstyle="arc3,rad=0.45", **arrows))
        icon.add_patch(FancyArrowPatch((-0.95, 0.35), (-0.35, -0.9), connectionstyle="arc3,rad=0.45", **arrows))
    else:
        icon.add_patch(FancyArrowPatch((0.78, 0), (1.4, 0), **arrows))
        icon.add_patch(FancyArrowPatch((-0.78, 0), (-1.4, 0), **arrows))
    icon.set(xlim=(-1.5, 1.5), ylim=(-1.05, 1.05), aspect="equal")
    icon.axis("off")


def draw(c):
    fig, ax = plt.subplots(figsize=(fs.WIDTH, 3.7), constrained_layout=True)
    for level in (1.5, 2.5, 3.5):
        ax.axhline(level, color=c.rule, lw=0.9, zorder=0)
    ax.fill_between(T, 0, C_trans, color=c.tint(c.blue, 0.75), lw=0)
    ax.fill_between(T, C_trans, C_trans + C_rot, color=c.tint(c.vermilion, 0.75), lw=0)
    ax.fill_between(T, C_trans + C_rot, C_total, color=c.tint(c.green, 0.75), lw=0)
    ax.plot(T, C_total, color=c.ink)
    for theta in (theta_rot, theta_vib):
        ax.plot([theta, theta], [0, 3.72], color=c.muted, lw=1.2, ls=(0, (4, 3)))
    ax.text(theta_rot, 3.9, rf"$\Theta_{{\mathrm{{rot}}}}$ = {theta_rot:.3f} K", ha="center", va="center")
    ax.text(theta_vib, 3.9, rf"$\Theta_{{\mathrm{{vib}}}}$ = {theta_vib:.0f} K", ha="center", va="center")
    ax.plot(300, np.interp(300, T, C_total), "o", mfc=c.ground, mec=c.ink, mew=1.6, ms=6.5, zorder=4)
    ax.text(300, 2.72, "300 K", ha="center", va="center")

    # the three kinds of motion, each inside its contribution
    draw_motion(ax, c, "translation", c.blue, (0.36, 0.11))
    ax.text(0.47, 0.1825, "Translation", transform=ax.transAxes, ha="left", va="center")
    draw_motion(ax, c, "rotation", c.vermilion, (0.36, 0.415))
    ax.text(0.47, 0.4875, "Rotation", transform=ax.transAxes, ha="left", va="center")
    draw_motion(ax, c, "vibration", c.green, (0.795, 0.66))
    ax.text(0.9, 0.7325, "Vibration", transform=ax.transAxes, ha="left", va="center")

    ax.set(xscale="log", xlabel="Temperature $T$ (K)", ylabel="Molar heat capacity $C_{V,m}/R$", xlim=(0.3, 4e4),
           ylim=(0, 4.1))
    ax.set_xticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1000", "10 000"])
    ax.minorticks_off()
    ax.set_yticks([0, 1.5, 2.5, 3.5], ["0", "3/2", "5/2", "7/2"])
    return fig


fs.show(draw, alt="Heat capacity of carbon monoxide gas divided by R against temperature from 0.3 to 40 000 kelvin on a "
                  "logarithmic axis. Three shaded bands are stacked. The translational band has a constant height of "
                  "3/2. The rotational band rises near the rotational temperature of 2.779 kelvin, passes through a "
                  "shallow maximum, and then has a height of 1, which brings the total to 5/2 from about 10 kelvin "
                  "to several hundred kelvin, including 300 kelvin. The vibrational band rises around the "
                  "vibrational temperature of 3122 kelvin, and the total approaches 7/2 at the highest temperatures. "
                  "Small sketches of a diatomic molecule that moves, rotates, and vibrates label the three bands.")
```

## Worked Example

### High-$T$ rotational partition function for a spherical top (CH$_4$)

For a spherical top,

```{math}
q_{\mathrm{rot}} \approx \frac{\sqrt{\pi}}{\sigma}\left(\frac{T}{\Theta_{\mathrm{rot}}}\right)^{3/2}.
```

Using the rotational temperature of $\mathrm{CH_4}$ in {numref}`rotational-temperatures-nonlinear`, $\Theta_{\mathrm{rot}}=7.541\ \mathrm{K}$, and its symmetry number, $\sigma=12$:

1. **Compute $T/\Theta_{\mathrm{rot}}$ at $T=300\ \mathrm{K}$**

   ```{math}
   \frac{T}{\Theta_{\mathrm{rot}}}=\frac{300}{7.541}=39.78.
   ```

2. **Evaluate $q_{\mathrm{rot}}$**

   ```{math}
   q_{\mathrm{rot}}
   \approx \frac{\sqrt{\pi}}{12}(39.78)^{3/2}
   =\frac{1.7725}{12}\,(39.78)\sqrt{39.78}
   \approx 37.1.
   ```

**Result.** At room temperature, the rotational partition function of CH$_4$ is large, indicating many thermally accessible rotational states.

## Concept Checks

1. Why do nonlinear molecules typically have $q_{\mathrm{rot}}\propto T^{3/2}$ in the high-$T$ limit rather than $q_{\mathrm{rot}}\propto T$?
2. What is the physical meaning of the symmetry number $\sigma$, and why does it reduce $q_{\mathrm{rot}}$?
3. When is it reasonable to approximate $q_{\mathrm{elec}}\approx g_1$?
4. For a polyatomic molecule, why is $q_{\mathrm{vib}}$ a product over normal modes?

## Key Takeaways

* Total molecular partition functions are built as products of translational, rotational, vibrational, and electronic factors.
* Rigid-rotor high-$T$ formulas depend on molecular symmetry and rotational temperatures.
* Harmonic-oscillator vibration introduces $\Theta_{\mathrm{vib}}$ and freezes out at $T\ll\Theta_{\mathrm{vib}}$.
* The electronic partition function often reduces to the degeneracy of the ground level, $q_{\mathrm{elec}}\approx g_1$.
