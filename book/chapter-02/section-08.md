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

Below is a visualization comparing *spherical top*, *symmetric top*, and *asymmetric top* molecules. The ellipsoid semi-axes are scaled by the three principal rotational temperatures $\Theta_{\mathrm{rot},i}$, highlighting how anisotropy in the rotational constants maps onto the 3D “rotational profile.”

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from myst_nb import glue

u = np.linspace(0, 2*np.pi, 100)
v = np.linspace(0, np.pi, 100)
u, v = np.meshgrid(u, v)

def surface_points(a, b, c):
    """Return x, y, z for an ellipsoid with semi-axes a, b, c."""
    x = a * np.cos(u) * np.sin(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(v)
    return x, y, z

# Rotational temperatures in K (see the table above), listed as the
# semi-axes along x, y, and z
theta_ch4 = (7.54, 7.54, 7.54)     # Spherical top
theta_nh3 = (14.31, 14.31, 8.96)   # Symmetric top: the c axis (threefold axis) is z
theta_h2o = (40.1, 13.4, 20.9)     # Asymmetric top: a is x, c is y (normal to the
                                   # molecular plane), and b is z (twofold axis)

fig = plt.figure(figsize=(9, 3.4))

# 1) CH4: Spherical top
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
x1, y1, z1 = surface_points(*theta_ch4)
surf1 = ax1.plot_surface(x1, y1, z1, color='skyblue', alpha=0.4, edgecolor='none')
ax1.set_title("Spherical Top", fontsize=11)
ax1.set_box_aspect((1,1,1), zoom=0.85)
ax1.set_xlim(-7.54, 7.54)
ax1.set_ylim(-7.54, 7.54)
ax1.set_zlim(-7.54, 7.54)
ax1.view_init(elev=20, azim=35)
ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$')
ax1.set_zlabel('$z$')
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_zticks([])

# Add a CH4 molecule to the plot (schematic coordinates)
positions = np.array([
    [ 0.0000,  0.0000,  0.0000],
    [ 0.6276,  0.6276,  0.6276],
    [ 0.6276, -0.6276, -0.6276],
    [-0.6276,  0.6276, -0.6276],
    [-0.6276, -0.6276,  0.6276]
]) * 7.54/2

for i, pos in enumerate(positions):
    if i == 0:
        ax1.plot(*pos, 'o', markersize=6, markeredgecolor='black', markerfacecolor='black')
    else:
        ax1.plot(*pos, 'o', markersize=3.5, markeredgecolor='black', markerfacecolor='white')
    if i > 0:
        ax1.plot([positions[0][0], pos[0]],
                 [positions[0][1], pos[1]],
                 [positions[0][2], pos[2]], color='black', linewidth=1)

ax1.text(7.54/4, -7.54/4, 7.54/16, 'CH$_4$', color='black', fontsize=10)

# 2) NH3: Symmetric top
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
x2, y2, z2 = surface_points(*theta_nh3)
surf2 = ax2.plot_surface(x2, y2, z2, color='salmon', alpha=0.4, edgecolor='none')
ax2.set_title("Symmetric Top", fontsize=11)
ax2.set_box_aspect((1,1,1), zoom=0.85)
ax2.set_xlim(-14.31, 14.31)
ax2.set_ylim(-14.31, 14.31)
ax2.set_zlim(-14.31, 14.31)
ax2.view_init(elev=20, azim=35)
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$')
ax2.set_zlabel('$z$')
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_zticks([])

positions = np.array([
    [ 0.0000,  0.0000,  0.0000],
    [ 0.0000, -0.9377, -0.3816],
    [ 0.8121,  0.4689, -0.3816],
    [-0.8121,  0.4689, -0.3816]
]) * 14.31/2

for i, pos in enumerate(positions):
    if i == 0:
        ax2.plot(*pos, 'o', markersize=6, markeredgecolor='black', markerfacecolor='skyblue')
    else:
        ax2.plot(*pos, 'o', markersize=3.5, markeredgecolor='black', markerfacecolor='white')
    if i > 0:
        ax2.plot([positions[0][0], pos[0]],
                 [positions[0][1], pos[1]],
                 [positions[0][2], pos[2]], color='black', linewidth=1)

ax2.text(14.31/4, -14.31/4, 14.31/16, 'NH$_3$', color='black', fontsize=10)

# 3) H2O: Asymmetric top
ax3 = fig.add_subplot(1, 3, 3, projection='3d')
x3, y3, z3 = surface_points(*theta_h2o)
surf3 = ax3.plot_surface(x3, y3, z3, color='plum', alpha=0.4, edgecolor='none')
ax3.set_title("Asymmetric Top", fontsize=11)
ax3.set_box_aspect((1,1,1), zoom=0.85)
ax3.set_xlim(-40.1, 40.1)
ax3.set_ylim(-40.1, 40.1)
ax3.set_zlim(-40.1, 40.1)
ax3.view_init(elev=20, azim=35)
ax3.set_xlabel('$x$')
ax3.set_ylabel('$y$')
ax3.set_zlabel('$z$')
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_zticks([])

positions = np.array([
    [ 0.0000,  0.0000,  0.1173],
    [ 0.7572,  0.0000, -0.4692],
    [-0.7572,  0.0000, -0.4692]
]) * 40.1/2

for i, pos in enumerate(positions):
    if i == 0:
        ax3.plot(*pos, 'o', markersize=6, markeredgecolor='black', markerfacecolor='red')
    else:
        ax3.plot(*pos, 'o', markersize=3.5, markeredgecolor='black', markerfacecolor='white')
    if i > 0:
        ax3.plot([positions[0][0], pos[0]],
                 [positions[0][1], pos[1]],
                 [positions[0][2], pos[2]], color='black', linewidth=1)

ax3.text(40.1/4, -40.1/4, 40.1/16, 'H$_2$O', color='black', fontsize=10)

fig.subplots_adjust(left=0.0, right=0.97, bottom=0.0, top=0.95, wspace=0.0)
plt.show()
plt.close(fig)
```

Approximate visualization of rotational shapes: (left) a spherical top (CH₄), (middle) a symmetric top (NH₃), and (right) an asymmetric top (H₂O). The semi-axes of each ellipsoid are proportional to the rotational temperatures in {numref}`rotational-temperatures-nonlinear` and lie along the corresponding principal axes. Each panel has its own scale.

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
