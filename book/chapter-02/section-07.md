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
  enumerator: 2.7.%s
---


# 2.7. Linear Rigid Rotor

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The linear rigid rotor models molecular rotation (especially for diatomics) and provides the rotational partition function used in molecular thermodynamics. We derive the quantized rotational levels and their degeneracies, develop the high-$T$ approximation $q_{\mathrm{rot}}\approx T/(\sigma\Theta_{\mathrm{rot}})$, and obtain the corresponding rotational contributions to the internal energy and heat capacity.

Learning objectives:

* State the rigid-rotor energy levels $E_J=\hbar^2J(J+1)/(2I)$ and degeneracy $g_J=2J+1$.
* Write the rotational partition function as a sum over $J$, including degeneracy.
* Derive the high-$T$ approximation and define the rotational temperature $\Theta_{\mathrm{rot}}=\hbar^2/(2k_{\mathrm{B}}I)$.
* Compute rotational contributions to $U$ and $C_V$ in the classical limit, and interpret the symmetry number $\sigma$.

## Core Ideas and Derivations

### Review of the Linear Rigid Rotor

[](#fig:rotor-levels) shows the first four energy levels of the linear rigid rotor that is reviewed in this subsection.

```{code-cell} ipython3
:label: fig:rotor-levels
:caption: Energy levels of a linear rigid rotor in units of $\hbar^2/2I$. Each short line is one microstate. The level $E_J$ has degeneracy $g_J = 2J + 1$, because the magnetic quantum number $m$ takes the $2J+1$ integer values from $-J$ to $J$. For example, the $J=2$ level contains five degenerate microstates, with $m=-2,-1,0,1,2$ (the numbers below the lines). The spacing between adjacent levels increases with $J$.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs


def draw(c):
    fig, ax = plt.subplots(figsize=(5.6, 3.7), constrained_layout=True)
    for J in range(4):
        energy = J * (J + 1)                 # in units of hbar^2 / (2 I)
        for m in range(-J, J + 1):
            ax.plot([m - 0.36, m + 0.36], [energy, energy], color=c.ink, solid_capstyle="butt")
            ax.text(m, energy - 0.28, str(m).replace("-", "−"), color=c.muted, ha="center", va="top",
                    fontsize=fs.BASE - 1.5)
        ax.text(-J - 0.62, energy - 0.28, "$m =$", color=c.muted, ha="right", va="top", fontsize=fs.BASE - 1.5)
        ax.text(4.0, energy, rf"$J = {J}$,  $g_J = {2 * J + 1}$", va="center", clip_on=False)
    ax.set(ylabel=r"Energy $E_J$ in units of $\hbar^2/2I$", xlim=(-4.6, 3.8), ylim=(-1.4, 13), xticks=[],
           yticks=[0, 2, 6, 12])
    ax.spines["bottom"].set_visible(False)
    return fig


fs.show(draw, alt="Energy level diagram of a linear rigid rotor. The levels J equal to 0, 1, 2, and 3 lie at 0, 2, "
                  "6, and 12 in units of h-bar squared over 2 I. Each level is drawn as a row of short lines, one "
                  "for each value of m from minus J to J, so the rows have 1, 3, 5, and 7 lines.")
```

A linear rigid rotor is a linear molecule whose bond lengths are fixed. The molecule can rotate about two perpendicular axes, $a$ and $b$, that pass through its center of mass and are perpendicular to the molecular axis. Rotation about the molecular axis is not counted, because the nuclei lie on that axis. In classical mechanics, the rotational energy is

```{math}
E = \frac{L_a^2}{2I} + \frac{L_b^2}{2I},
```

where $L_a$ and $L_b$ are the components of the angular momentum along the two axes and $I$ is the moment of inertia, which is the same for both axes.

In quantum mechanics, the rotational energy is quantized. For a linear rigid rotor with moment of inertia $I$, the energy levels are:

```{math}
E_J \;=\; \frac{\hbar^2}{2I}\,J\bigl(J+1\bigr)
\quad\text{for}\quad 
J = 0, 1, 2, \dots
```

Here, $\hbar$ is the reduced Planck constant, $J$ is the rotational quantum number, and $g_J=2J+1$ is the degeneracy of the level $E_J$. The quantum number $m$ specifies the component of the angular momentum along a chosen axis in space, $m\hbar$, and takes the values $m=-J,\ldots,J$. The energy $E_J$ does not depend on $m$.

````{admonition} Physical Significance of the Moment of Inertia
:class: note
The moment of inertia $I$ measures how mass is distributed about the rotation axis. For a diatomic molecule with atoms A and B,
```{math}
  I \;=\; \mu\,r^2, 
  \quad
  \mu = \frac{m_{\mathrm{A}}\,m_{\mathrm{B}}}{m_{\mathrm{A}} + m_{\mathrm{B}}},
```
where $r$ is the bond length and $\mu$ is the reduced mass.
````

### Partition Function for a Linear Rigid Rotor

In the canonical ensemble, the rotational partition function is

```{math}
q_{\mathrm{rot}} 
\;=\; 
\sum_{i=0}^\infty e^{-\beta E_i}
\;=\;
\sum_{J=0}^\infty g_J \, e^{-\beta\,E_J}
\;=\;
\sum_{J=0}^\infty (2J+1)\,\exp\!\Bigl[-\beta \,\frac{\hbar^2}{2I}\,J(J+1)\Bigr].
```

```{admonition} Summation vs. Levels
:class: dropdown
- **Sum over $i$**: sums over individual microstates (each degenerate sub-level).
- **Sum over $J$**: reorganizes the sum by energy level $E_J$, factoring in the degeneracy $g_J = 2J+1$.
```

#### High-Temperature Approximation

When $k_{\mathrm{B}}T \gg \frac{\hbar^2}{2I}$, we can approximate the discrete sum by treating $J$ as continuous and converting the sum to an integral. Define

```{math}
x \;=\; J(J+1),
\quad
dx \;=\; (2J+1)\,dJ.
```

Then

```{math}
q_{\mathrm{rot}}
\;\approx\;
\int_{0}^{\infty} (2J+1)\,\exp\Bigl[-\beta \,\frac{\hbar^2}{2I}\,J(J+1)\Bigr]\;dJ
\;=\;
\int_{x=0}^{\infty} 
\exp\Bigl[-\beta \,\tfrac{\hbar^2}{2I}\,x\Bigr]\;dx.
```

Evaluating the integral gives

```{math}
\int_{0}^{\infty} 
\exp\Bigl[-\beta \,\tfrac{\hbar^2}{2I}\,x\Bigr]\;dx
\;=\;
\frac{1}{\beta \,\frac{\hbar^2}{2I}}
\;=\;
\frac{2I\,k_{\mathrm{B}} T}{\hbar^2}.
```

We define the **rotational temperature** $\Theta_{\mathrm{rot}}$ by

```{math}
\Theta_{\mathrm{rot}} 
\;=\;
\frac{\hbar^2}{2k_{\mathrm{B}}\,I}.
```

In terms of $\Theta_{\mathrm{rot}}$, the energy levels are $E_J = k_{\mathrm{B}}\Theta_{\mathrm{rot}}\,J(J+1)$, so $\Theta_{\mathrm{rot}}$ sets the temperature scale of the level spacing, and the condition for the high-temperature approximation is $T \gg \Theta_{\mathrm{rot}}$.

Thus, for a heteronuclear diatomic rotor (symmetry number $\sigma=1$),

```{math}
q_{\mathrm{rot}}
\;\approx\;
\frac{2I\,k_{\mathrm{B}} T}{\hbar^2}
\;=\;
\frac{T}{\Theta_{\mathrm{rot}}}.
```

````{admonition} Symmetry Number $\sigma$
:class: tip
For a **homonuclear** diatomic molecule (and, more generally, a symmetric linear molecule such as $\mathrm{CO_2}$), a rotation by 180° about an axis perpendicular to the molecular axis exchanges identical nuclei and gives an orientation that cannot be distinguished from the original. The high-temperature result above counts these two orientations as different states, so it must be divided by the **symmetry number** $\sigma=2$ (compare the factor of $1/N!$ in [Section 2.4](section-04.md)):
```{math}
q_{\mathrm{rot}} \;\approx\; \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}.
```
For a heteronuclear diatomic molecule, $\sigma=1$. The symmetry number is required whenever the value of $q_{\mathrm{rot}}$ matters, as in the entropy ([Section 4.3](../chapter-04/section-03.md)) and in equilibrium constants ([Section 7.2](../chapter-07/section-02.md)). It has no effect on $U_{\mathrm{rot}}$ or $C_V^{(\mathrm{rot})}$, which depend only on the temperature derivative of $\ln q_{\mathrm{rot}}$.
````

{numref}`rotational-temperatures` lists $\Theta_{\mathrm{rot}}$ and $\sigma$ for the diatomic molecules whose vibrational temperatures are given in [Section 2.6](section-06.md). The values of $\Theta_{\mathrm{rot}}$ are calculated from the rotational constant $\tilde{B} = h/(8\pi^2 c I)$ as $\Theta_{\mathrm{rot}} = hc\tilde{B}/k_{\mathrm{B}}$, where $c$ is the speed of light and $hc/k_{\mathrm{B}} = 1.43878\ \mathrm{cm\,K}$. The rotational constants and bond lengths are equilibrium values from the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/); the rotational constant of $\mathrm{O_2}$ is derived from the constants listed there for its lowest vibrational level.

```{list-table} Rotational constants, bond lengths, rotational temperatures, and symmetry numbers of selected diatomic molecules
:header-rows: 1
:name: rotational-temperatures

* - Molecule
  - $\tilde{B}$ ($\mathrm{cm^{-1}}$)
  - $r$ (Å)
  - $\Theta_{\mathrm{rot}}$ (K)
  - $\sigma$
* - $\mathrm{H_2}$
  - 60.853
  - 0.7414
  - 87.55
  - 2
* - $\mathrm{HCl}$
  - 10.593
  - 1.275
  - 15.24
  - 1
* - $\mathrm{N_2}$
  - 1.9982
  - 1.098
  - 2.875
  - 2
* - $\mathrm{CO}$
  - 1.9313
  - 1.128
  - 2.779
  - 1
* - $\mathrm{O_2}$
  - 1.4456
  - 1.208
  - 2.080
  - 2
* - $\mathrm{Cl_2}$
  - 0.24399
  - 1.988
  - 0.3510
  - 2
* - $\mathrm{Br_2}$
  - 0.082107
  - 2.281
  - 0.1181
  - 2
* - $\mathrm{I_2}$
  - 0.037372
  - 2.666
  - 0.05377
  - 2
```

### Ensemble Averages

#### Natural Logarithm of the Partition Function

From the high-$T$ approximation (with $\sigma=1$; a different value of $\sigma$ adds only a constant to $\ln q_{\mathrm{rot}}$ and does not change the results below),

```{math}
\ln q_{\mathrm{rot}}
\;=\;
\ln \Bigl(\tfrac{T}{\Theta_{\mathrm{rot}}}\Bigr)
\;=\;
\ln T \;-\;\ln \Theta_{\mathrm{rot}}.
```

#### Internal Energy

The rotational internal energy (per rotor) follows from the derivative in [Section 2.3](section-03.md):

```{math}
U_{\mathrm{rot}} 
\;=\;
- \left(\frac{\partial \ln q_{\mathrm{rot}}}{\partial \beta}\right)_{N,V}
\;=\;
k_{\mathrm{B}}\,T^2 \left(\frac{\partial \ln q_{\mathrm{rot}}}{\partial T}\right)_{N,V}.
```

Since $\ln q_{\mathrm{rot}}=\ln T-\ln \Theta_{\mathrm{rot}}$, we have

```{math}
\frac{\partial \ln q_{\mathrm{rot}}}{\partial T} 
\;=\; 
\frac{\partial}{\partial T}\bigl(\ln T\bigr) 
\;=\; 
\frac{1}{T}.
```

Hence,

```{math}
U_{\mathrm{rot}}
\;=\;
k_{\mathrm{B}} \, T.
```

#### Heat Capacity at Constant Volume

The rotational contribution to the heat capacity (per rotor) is

```{math}
C_V^{(\mathrm{rot})}
\;=\;
\left(\frac{\partial U_{\mathrm{rot}}}{\partial T}\right)_{N,V}
\;=\;
k_{\mathrm{B}}.
```

A single linear rotor therefore contributes $k_{\mathrm{B}}$ to the heat capacity in the classical (high-$T$) limit.

```{admonition} Comparing with equipartition
:class: important
Count the quadratic terms in the classical rotational energy of the linear rigid rotor, and use the equipartition theorem ([Section 2.6](section-06.md)) to predict $U_{\mathrm{rot}}$ and $C_V^{(\mathrm{rot})}$. Compare the predictions with the results above.
```

```{note}
In a more rigorous quantum treatment (and at lower temperatures), the partition function and the resulting averages must use the full sum over $J$. At sufficiently low $T$, only the $J=0$ and $J=1$ levels are significantly populated, which reduces the effective heat capacity below $k_{\mathrm{B}}$.
```

### Computational Studio: Linear Rigid Rotor

In [](#fig:rigid-rotor-studio), you can choose the masses, the bond length, and the temperature of a diatomic molecule, or one of the molecules of {numref}`rotational-temperatures`, and follow $\Theta_{\mathrm{rot}}$, $q_{\mathrm{rot}}$, and the probabilities of the rotational levels.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/rigid-rotor-studio.html
:label: fig:rigid-rotor-studio
:width: 100%
:title: Linear rigid rotor
:placeholder: ../_static/rigid-rotor-studio/rigid-rotor-studio.png

A diatomic molecule as a linear rigid rotor. Panel (a) shows the molecule to scale, turning about its center of mass (the cross). Panel (b) lists the reduced mass, the moment of inertia, the rotational constant, the rotational temperature, the partition function from the sum over $J$ and from the high-temperature approximation, and the molar entropy of rotation, which [Section 4.3](../chapter-04/section-03.md) obtains from the partition function. Panel (c) shows the probabilities $p_J$ of the levels, and panel (d) shows their energies beside the temperature (dashed line). For two identical atoms, panel (c) does not include the effect of the nuclear spins on the levels with even and odd $J$.
```

You can also open the studio in a new tab: [Linear Rigid Rotor Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/rigid-rotor-studio.html).

## Worked Example

### Rotational temperature and $q_{\mathrm{rot}}$ for CO

Approximate CO as a rigid rotor with bond length $r = 1.128\ \text{\AA} = 1.128\times10^{-10}\ \mathrm{m}$. Use $m_{\mathrm{C}}=12u$, $m_{\mathrm{O}}=16u$, and $u=1.66054\times10^{-27}\ \mathrm{kg}$.

1. **Reduced mass**

   ```{math}
   \mu=\frac{m_{\mathrm{C}}m_{\mathrm{O}}}{m_{\mathrm{C}}+m_{\mathrm{O}}}
   =\frac{(12u)(16u)}{28u}=6.857u
   =1.14\times10^{-26}\ \mathrm{kg}.
   ```

2. **Moment of inertia**

   ```{math}
   I=\mu r^2=(1.14\times10^{-26})(1.128\times10^{-10})^2
   =1.45\times10^{-46}\ \mathrm{kg\,m^2}.
   ```

3. **Rotational temperature**

   ```{math}
   \Theta_{\mathrm{rot}}=\frac{\hbar^2}{2k_{\mathrm{B}}I}
   =\frac{(1.055\times10^{-34})^2}{2(1.381\times10^{-23})(1.45\times10^{-46})}
   \approx 2.78\ \mathrm{K}.
   ```

4. **High-$T$ partition function (heteronuclear, $\sigma=1$) at $T=300\ \mathrm{K}$**

   ```{math}
   q_{\mathrm{rot}}\approx \frac{T}{\Theta_{\mathrm{rot}}}=\frac{300}{2.78}=107.9.
   ```

5. **Direct sum for comparison**

   Adding the terms of the sum over $J$ numerically (the terms beyond $J=60$ are negligible) gives

   ```{math}
   q_{\mathrm{rot}}=\sum_{J=0}^{\infty}(2J+1)\,e^{-\Theta_{\mathrm{rot}}J(J+1)/T}=108.2.
   ```

**Result.** For CO at room temperature, $T\gg \Theta_{\mathrm{rot}}$, and the high-$T$ approximation agrees with the direct sum to within 0.3%.

## Concept Checks

1. Why does each $J$ level have degeneracy $2J+1$? What symmetry is responsible?
2. What changes in the partition function when the molecule is homonuclear rather than heteronuclear?
3. Why does the classical (high-$T$) rotor have $U_{\mathrm{rot}}=k_{\mathrm{B}}T$ per molecule?
4. What physical parameter(s) of the molecule increase $\Theta_{\mathrm{rot}}$?

## Key Takeaways

* Rigid-rotor levels scale as $J(J+1)$ with degeneracy $2J+1$.
* In the high-$T$ limit, $q_{\mathrm{rot}}\approx T/(\sigma\Theta_{\mathrm{rot}})$.
* $\Theta_{\mathrm{rot}}$ is set by the moment of inertia $I=\mu r^2$.
* Rotational contributions approach equipartition values at sufficiently high temperature.
