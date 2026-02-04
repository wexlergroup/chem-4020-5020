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


# 2.5. Particle in a Box

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

A particle in a box is the simplest quantized model of translation and provides a clean route to the translational partition function. We review the 1D and 3D energy levels, derive $q_{\mathrm{trans}}=V/\Lambda^3$ in the high-temperature/large-volume limit, and use $Q=q^N/N!$ to recover ideal-gas thermodynamic relations.

Learning objectives:

- Write the quantized energy levels for a particle in a 1D box and generalize to a 3D rectangular box and cube.
- Derive the single-particle partition function and its continuum approximation in the limit of small level spacing.
- Define the thermal de Broglie wavelength $\Lambda$ and express $q_{\mathrm{trans}}=V/\Lambda^3$.
- Use $Q=(V/\Lambda^3)^N/N!$ to obtain $U=\tfrac{3}{2}Nk_{\mathrm B}T$, $C_V=\tfrac{3}{2}Nk_{\mathrm B}$, and $P=Nk_{\mathrm B}T/V$.

## Core Ideas and Derivations

### Review of the Particle in a Box

#### Particle in a One-Dimensional Box

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.constants import k, eV
from labellines import labelLines
from myst_nb import glue

fig, axs = plt.subplot_mosaic([[0]], figsize=(4, 4))

# Plot a one-dimensional box
axs[0].plot([-0.5, -0.5], [18, 0], color='black', zorder=2.5, lw=4)
axs[0].plot([-0.5, 0.5], [0, 0], color='black', zorder=2.5, lw=4)
axs[0].plot([0.5, 0.5], [0, 18], color='black', zorder=2.5, lw=4)

# Plot the energy levels (blue lines)
for n in range(1, 5):
    energy = axs[0].plot([-0.5, 0.5], [n**2, n**2], color='blue', label=r'$E_{%d}$' % n)
    labelLines(energy, xvals=[-1/6], zorder=2.5)

# Plot the wavefunctions (red curves)
x = np.linspace(-0.5, 0.5, 100)
for n in range(1, 5):
    wavefunction = axs[0].plot(x, np.sin(n * np.pi * (x + 0.5)) + n**2, color='red', label=r'$\psi_{%d}$' % n)
    labelLines(wavefunction, xvals=[1/6], zorder=2.5)

axs[0].set_xlabel('Position')
axs[0].set_ylabel('Energy')
axs[0].set_xticks([-0.5, 0, 0.5])
axs[0].set_xticklabels([r'$-L/2$', r'$0$', r'$L/2$'])
axs[0].set_yticks([])
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[0].spines['left'].set_visible(False)

plt.tight_layout()
plt.show()
plt.close(fig)
```

Energy levels $E_n$ (blue) and wavefunctions $\psi_n$ (red) for a particle in a one-dimensional box with boundaries at $x=\pm L/2$.

For a particle in a one-dimensional box of length $L$, the quantized energy levels are

```{math}
E_n = \frac{h^2}{8 m L^2} \, n^2 \quad \text{for} \; n = 1, 2, 3, \ldots,
```

where $h$ is [Planck’s constant](https://physics.nist.gov/cgi-bin/cuu/Value?h), $m$ is the particle mass, and $n$ is a positive integer.

#### Particle in a Three-Dimensional Box

In three dimensions, the energy levels are

```{math}
E_{n_x, n_y, n_z} = \frac{h^2}{8 m}
\left( \frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} + \frac{n_z^2}{L_z^2} \right),
```

where $n_x$, $n_y$, and $n_z$ are positive integers, and $L_x$, $L_y$, and $L_z$ are the respective side lengths of a rectangular box.

##### Particle in a Cube

For a cube of side $L$, $L_x = L_y = L_z = L$. The energy simplifies to

```{math}
E_{n_x, n_y, n_z} 
= \frac{h^2}{8mL^2} \left( n_x^2 + n_y^2 + n_z^2 \right)
= \frac{h^2}{8 m V^{2/3}}
\left( n_x^2 + n_y^2 + n_z^2 \right),
```

where $V = L^3$ is the cube volume.

### Partition Function for a Particle in a Cube

The single-particle partition function for a three-dimensional cube is

```{math}
\begin{aligned}
q 
&= \sum_{n_x, n_y, n_z} 
   \exp \left[-\beta E_{n_x, n_y, n_z}\right] \\
&= \sum_{n_x=1}^{\infty} \sum_{n_y=1}^{\infty} \sum_{n_z=1}^{\infty} 
   \exp \left(-\beta \frac{h^2}{8 m L^2} (n_x^2 + n_y^2 + n_z^2)\right).
\end{aligned}
```

Define $\alpha = \beta h^2 / (8 m V^{2/3})$. Then

```{math}
:label: eq:partition-function-cube
q 
= \left( \sum_{n=1}^{\infty} e^{-\alpha n^2} \right)^3.
```

#### Approximation of the Sum

When the level spacing is small (large $T$ or large $L$), the sum $\sum_{n=1}^{\infty} e^{-\alpha n^2}$ is well approximated by an integral:

```{math}
\sum_{n = 1}^{\infty} e^{-\alpha n^2} 
\approx \int_0^{\infty} e^{-\alpha x^2} \, dx 
= \left(\frac{\pi}{4 \alpha}\right)^{\!\!1/2}.
```

Substituting $\alpha = \beta h^2/(8 m V^{2/3})$ gives

```{math}
\left(\frac{\pi}{4 \alpha}\right)^{\!\!1/2}
= \left(\frac{\pi}{4} \cdot \frac{8 m V^{2/3}}{\beta h^2}\right)^{\!\!1/2}
= \left(\frac{2 \pi m}{\beta h^2}\right)^{\!1/2} V^{1/3}.
```

We often define the thermal de Broglie wavelength $\Lambda$ as

```{math}
\Lambda 
= \left(\frac{h^2}{2 \pi m k_{\mathrm B} T}\right)^{\!1/2}.
```

````{admonition} What is the thermal de Broglie wavelength?
:class: tip

The thermal de Broglie wavelength, $\Lambda$, is a characteristic length scale for the quantum mechanical delocalization of a particle at temperature $T$. When the typical distance between particles is much larger than $\Lambda$, the system behaves classically and particles can be treated as effectively distinguishable. Conversely, when the interparticle spacing is on the order of $\Lambda$ (or smaller), quantum effects become significant and particle indistinguishability must be accounted for.
````

Hence,

```{math}
\sum_{n = 1}^{\infty} e^{-\alpha n^2} 
\approx \frac{V^{1/3}}{\Lambda}.
```

So the single-particle partition function becomes

```{math}
q 
= \frac{V}{\Lambda^3}.
```

### Partition Function for $N$ Particles in a Cube

For $N$ identical, non-interacting, indistinguishable particles, the total partition function $Q$ is

```{math}
:label: eq:partition-function-many-particles-cube
Q 
= \frac{q^N}{N!} 
= \frac{1}{N!}\left(\frac{V}{\Lambda^3}\right)^{\!N}.
```

The factor of $N!$ corrects for particle indistinguishability: permutations of identical particles do not produce new states.

### Ensemble Averages for $N$ Particles in a Cube

#### Natural Logarithm of the Partition Function

Taking the logarithm of $Q$ gives

```{math}
\ln Q 
= N \ln \left(\frac{V}{\Lambda^3}\right) - \ln N!.
```

Often, the Stirling approximation ($\ln N! \approx N \ln N - N$) is used for large $N$.

#### Internal Energy

The internal energy $U$ follows from

```{math}
:label: eq:internal-energy-cube
U 
= -\left(\frac{\partial \ln Q}{\partial \beta}\right)_{N,V}
= \frac{3}{2} \frac{N}{\beta}
= \frac{3}{2} N k_{\mathrm B} T.
```

```{important}
This result is consistent with equipartition: each of the three translational degrees of freedom contributes $\frac{1}{2}k_{\mathrm B}T$ per particle. Note that $U$ is independent of $V$ because the particles do not interact.
```

#### Heat Capacity at Constant Volume

From $U = \frac{3}{2}N k_{\mathrm B}T$, it follows that

```{math}
:label: eq:heat-capacity-cube
C_V 
= \left(\frac{\partial U}{\partial T}\right)_{N,V}
= \frac{3}{2} N k_{\mathrm B}.
```

```{important}
This temperature-independent $C_V$ is a hallmark of an ideal (monatomic) gas.
```

#### Pressure

To find the pressure, we use

```{math}
:label: eq:pressure-cube
P 
= k_{\mathrm B} T \left(\frac{\partial \ln Q}{\partial V}\right)_{N,T}
= \frac{N k_{\mathrm B} T}{V},
```

recovering the ideal gas law.

### Translational Partition Function

Because the particle-in-a-box spectrum describes translational motion, the single-particle partition function

```{math}
:label: eq:translational-partition-function
q_{\mathrm{trans}} 
= \frac{V}{\Lambda^3}
```

is called the translational partition function. It underpins the classical description of dilute gases when $\Lambda \ll$ (typical particle spacing), so quantum effects are negligible at ordinary densities and temperatures.

### Computational Studio: Ideal Gas

Use the interactive studio below (or open it in a new tab) to vary particle count ($N$), volume ($V$), and temperature ($T$). The studio visualizes how these parameters affect the gas while showing $P$, $U$, and $S$ derived from the partition function.

<!-- <div style="width: 100%; border: 1px solid #cbd5e1; border-radius: 12px; overflow: hidden; box-shadow: 0 4px 12px rgba(0,0,0,0.08);">
  <iframe
    src="https://chem-4020-5020-s763.vercel.app/"
    title="Ideal Gas Computational Studio"
    style="width: 100%; height: 900px; border: 0;"
    loading="lazy"
  ></iframe>
</div> -->

You can open the studio in a new tab: [Ideal Gas Computational Studio](https://chem-4020-5020-s763.vercel.app/).

## Worked Example

### Estimating $\Lambda$ and $q_{\mathrm{trans}}$ (helium, 300 K)

For a particle of mass $m$, the thermal de Broglie wavelength is

```{math}
\Lambda=\frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}.
```

Take helium: $m=4.0026\,u=6.65\times10^{-27}\ \mathrm{kg}$, $T=300\ \mathrm{K}$.

1. **Compute $\Lambda$**

   ```{math}
   \Lambda=\frac{6.626\times10^{-34}}{\sqrt{2\pi(6.65\times10^{-27})(1.381\times10^{-23})(300)}}
   \approx 5.0\times10^{-11}\ \mathrm{m}.
   ```

2. **Compute the single-particle translational partition function in $V=1.00\ \mathrm{L}=10^{-3}\ \mathrm{m^3}$**

   ```{math}
   q_{\mathrm{trans}}=\frac{V}{\Lambda^3}
   =\frac{10^{-3}}{(5.0\times10^{-11})^3}
   \approx 7.8\times10^{27}.
   ```

**Result.** $q_{\mathrm{trans}}\gg 1$ in macroscopic volumes, consistent with the continuum approximation used to derive $V/\Lambda^3$.

## Concept Checks

1. Why does the sum over quantum numbers become well-approximated by an integral at high $T$ or large $V$?
2. What does the condition "interparticle spacing $\gg \Lambda$" mean physically?
3. Why does translation give $U\propto T$ but not $U\propto V$ for an ideal gas?
4. Where does the $1/N!$ factor enter when connecting microscopic states to macroscopic entropy?

## Key Takeaways

- Quantized translation in a box leads, in the continuum limit, to $q_{\mathrm{trans}}=V/\Lambda^3$.
- $\Lambda$ sets the scale of quantum wavepacket "size"; classical behavior emerges when $\Lambda$ is small compared to typical spacings.
- For an ideal monatomic gas, $Q=(V/\Lambda^3)^N/N!$ reproduces $PV=Nk_{\mathrm B}T$ and $U=\tfrac{3}{2}Nk_{\mathrm B}T$.
- Partition functions connect microscopic spectra to thermodynamic equations of state via derivatives of $\ln Q$.
