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
  enumerator: 1.2.%s
---

# 1.2. Kinetic Theory

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 1.2](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-01/section-02.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

Kinetic theory provides a microscopic route to macroscopic gas behavior by modeling a gas as a large number of rapidly moving particles that undergo elastic collisions. In this section, we derive the pressure of an ideal gas from particle–wall momentum transfer, connect temperature to average kinetic energy, and interpret molecular speed distributions.

---

Learning objectives:

- List the assumptions of kinetic theory used to model an ideal gas.
- Derive the relation $P=\tfrac{1}{3}(N/V)m\langle v^{2}\rangle$ from particle–wall momentum transfer.
- Use $\langle E_{\mathrm{kin}}\rangle = \tfrac{3}{2}k_{\mathrm{B}}T$ to relate temperature to molecular motion.
- Compute characteristic speeds ($v_{\mathrm{mp}},\ \langle v\rangle,\ v_{\mathrm{rms}}$) and interpret the Maxwell–Boltzmann distribution.

## Core Ideas and Derivations

### Foundational Assumptions of Kinetic Theory

1. **Large Number of Particles**:  
   A gas contains a very large number of identical particles moving randomly in all directions.

2. **Point Particles**:  
   Each particle’s size is negligible compared to the average distance between particles.

3. **Elastic Collisions**:  
   Collisions between particles and between particles and the container walls conserve both momentum and kinetic energy.

4. **No Long-Range Interparticle Forces**:  
   Particles exert no forces on one another except during collisions (i.e., there are no long-range attractive or repulsive forces).

5. **Classical Mechanics Applies**:  
   Particle motion follows Newton’s second law:

   ```{math}
   \vec{F} = \frac{d\vec{p}}{dt},
   ```

   where $\vec{F}$ is the net force on a particle, $\vec{p}$ is its linear momentum, and $t$ is time.

Assumptions 2 and 4 are approximations, because real molecules have a finite size and exert forces on one another. [Section 1.3](section-03.md) and [Section 1.4](section-04.md) examine the conditions under which these approximations hold and the models that replace them when they fail.

---

### Deriving Pressure from Particle-Wall Collisions

```{code-cell} ipython3
:label: fig:particle-in-container
:caption: Two-dimensional schematic of a single gas particle of mass $m$ in a cuboid container. The particle moves toward the right wall with velocity component $v_x$, crosses the container with $-v_x$ after the collision, and moves to the right again after the collision with the left wall. The length $L_y$ is perpendicular to the plane of the drawing.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs


def draw(c):
    fig, ax = plt.subplots(figsize=(6.4, 2.1), constrained_layout=True)
    Lx, Lz = 10.0, 2.4
    ax.add_patch(Rectangle((0, 0), Lx, Lz, fc=c.fill, ec=c.ink, lw=1.6))

    # The particle and its three legs: toward the right wall, back across the box, and forward again.
    y_out, y_back, y_again = 0.78 * Lz, 0.50 * Lz, 0.22 * Lz
    ax.add_patch(Circle((0.5 * Lx, y_out), 0.24, fc=c.blue, ec="none", zorder=5))
    ax.text(0.5 * Lx, y_out, "$m$", color=fs.on(c.blue), ha="center", va="center", zorder=6)
    fs.arrow(ax, (0.5 * Lx + 0.35, y_out), (Lx - 0.08, y_out), c.vermilion)
    fs.arrow(ax, (Lx - 0.08, y_back), (0.08, y_back), c.vermilion)
    fs.arrow(ax, (0.08, y_again), (0.5 * Lx, y_again), c.vermilion)
    for x, y, label in [(0.76 * Lx, y_out, "$v_x$"), (0.76 * Lx, y_back, "$-v_x$"), (0.24 * Lx, y_again, "$v_x$")]:
        ax.text(x, y + 0.12, label, color=c.vermilion, ha="center", va="bottom")

    # Dimensions of the container
    kw = dict(arrowstyle="<->", color=c.muted, lw=1.0, shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(0, -0.35), xytext=(Lx, -0.35), arrowprops=kw)
    ax.text(Lx / 2, -0.55, "$L_x$", ha="center", va="top", color=c.muted)
    ax.annotate("", xy=(-0.3, 0), xytext=(-0.3, Lz), arrowprops=kw)
    ax.text(-0.5, Lz / 2, "$L_z$", ha="right", va="center", color=c.muted)

    ax.set(xlim=(-1.2, Lx + 0.3), ylim=(-1.1, Lz + 0.2), aspect="equal")
    ax.axis("off")
    return fig


fs.show(draw, alt="A rectangular container of length L x and height L z that holds one particle of mass m. Arrows show "
                  "the particle moving to the right wall with velocity v x, crossing the whole box to the left wall "
                  "with velocity minus v x, and moving to the right again with velocity v x.")
```

#### Microscopic Picture of Pressure

Pressure is the force exerted per unit area on the container walls. Microscopically, it arises from momentum transfer during particle–wall collisions.

##### Particle Momentum Change

Consider an elastic collision of a particle of mass $m$ with a wall perpendicular to the $x$-axis. The $x$-component of the velocity reverses ($v_x \to -v_x$). If we take $v_x>0$ to denote the *magnitude* of the $x$-component, then the magnitude of the particle’s momentum change is

```{math}
\Delta p_x = 2 m v_x.
```

##### Time Between Collisions

If the container has length $L_x$ in the $x$-direction, the time between successive collisions of the same particle with that wall is

```{math}
\Delta t = \frac{2 L_x}{v_x}.
```

##### Force on the Wall

A single particle’s average force on the wall (in the $x$-direction) is then

```{math}
F_{\text{p}, x} = \frac{\Delta p_x}{\Delta t} = \frac{m v_x^2}{L_x}.
```

#### Total Pressure

The wall perpendicular to the $x$-axis has area $L_y L_z$. Summing the single-particle forces over $N$ identical particles and dividing by this area gives the pressure on the wall,

```{math}
P = \frac{1}{L_y L_z} \sum_{i=1}^N \frac{m v_{x,i}^2}{L_x} = \frac{N m \langle v_x^2 \rangle}{V},
```

where $v_{x,i}$ is the $x$-component of the velocity of the $i$-th particle, $V = L_x L_y L_z$ is the volume of the container, and

```{math}
\langle v_x^2 \rangle = \frac{1}{N} \sum_{i=1}^N v_{x,i}^2
```

is the mean of the squared $x$-components.

The speed $v_i$ of the $i$-th particle satisfies $v_i^2 = v_{x,i}^2 + v_{y,i}^2 + v_{z,i}^2$, so that $\langle v^2 \rangle = \langle v_x^2 \rangle + \langle v_y^2 \rangle + \langle v_z^2 \rangle$. Because the particles move randomly in all directions (isotropic motion), no direction is preferred, and

```{math}
\langle v_x^2 \rangle = \langle v_y^2 \rangle = \langle v_z^2 \rangle = \frac{1}{3} \langle v^2 \rangle.
```

For $N$ identical particles with isotropic motion in a volume $V$, the total pressure $P$ is therefore

```{math}
P = \frac{1}{V} \sum_{i=1}^N \frac{1}{3} m v_i^2,
```

where $v_i$ is the speed of the $i$-th particle. Using the mean-square speed $\langle v^2 \rangle$, we obtain

```{math}
:label: pressure-kinetic-theory
P = \frac{N m \langle v^2 \rangle}{3 V}.
```

This equation shows how macroscopic pressure depends on the microscopic particle speeds.

---

### Kinetic Energy and Temperature

The average translational kinetic energy per particle is

```{math}
\langle E_{\mathrm{kin}} \rangle = \frac{1}{2} m \langle v^2 \rangle.
```

Equating Eq. {eq}`pressure-kinetic-theory` with the ideal-gas equation of state, $PV = N k_{\mathrm{B}} T$ (discussed in [Section 1.3](section-03.md)), gives

```{math}
:label: equipartition-theorem
\frac{1}{2} m \langle v^2 \rangle = \frac{3}{2} k_{\mathrm{B}} T,
```

where $k_{\mathrm{B}}$ is the Boltzmann constant and $T$ is the absolute temperature. This result—often presented as an application of equipartition—shows that temperature is directly proportional to the average translational kinetic energy of the particles.

The square root of the mean-square speed is the root-mean-square (rms) speed. From Eq. {eq}`equipartition-theorem`,

```{math}
:label: rms-speed
v_{\mathrm{rms}} = \sqrt{\langle v^2\rangle} = \sqrt{\frac{3k_{\mathrm{B}}T}{m}}.
```

`````{admonition} Complete Derivation of the Relationship Between Kinetic Energy and Temperature
:class: dropdown

**1. Kinetic Energy of a Single Particle**  
Consider a single particle with mass $m$ and speed $v$. Its translational kinetic energy is

```{math}
E_{\text{p},\,\mathrm{kin}} \;=\; \frac{1}{2}\,m\,v^{2}.
```

**2. Total Kinetic Energy of $N$ Particles**  
For $N$ particles with masses $m_1, m_2, \ldots, m_N$ and respective speeds $v_1, v_2, \ldots, v_N$, the total translational kinetic energy is

```{math}
E_{\mathrm{kin}}\left(m_1,\ldots,m_N\right) \;=\; \frac{1}{2}\,\sum_{i=1}^{N} m_i\,v_i^2.
```

If all particles are identical with mass $m$, this simplifies to

```{math}
E_{\mathrm{kin}} \;=\; \frac{1}{2}\,m\,\sum_{i=1}^{N} v_i^2.
```

**3. Defining the Average of the Speed Squared**  
Define the mean-square speed

```{math}
\langle v^2 \rangle \;=\; \frac{1}{N}\,\sum_{i=1}^{N} v_i^2,
```

so that

```{math}
\sum_{i=1}^{N} v_i^2 \;=\; N\,\langle v^2 \rangle.
```

Substituting back gives

````{important}
```{math}
E_{\mathrm{kin}} \;=\; \frac{1}{2}\,N\,m\,\langle v^2 \rangle.
```
````

**4. Relating Pressure to Kinetic Energy**  
From Eq. {eq}`pressure-kinetic-theory`, the pressure $P$ in a volume $V$ can be written as

```{math}
P \;=\; \frac{N\,m\,\langle v^2 \rangle}{3\,V}.
```

Recognizing that $N m \langle v^2 \rangle / 2 = E_{\mathrm{kin}}$, we obtain

```{math}
P \;=\; \frac{2\,E_{\mathrm{kin}}}{3\,V}.
```

**5. Equating to the Ideal-Gas Equation of State**  
For an ideal gas (discussed in [Section 1.3](section-03.md)),

```{math}
P \;=\; \frac{N\,k_{\mathrm{B}}\,T}{V}.
```

Equating the two expressions for $P$,

```{math}
\frac{2\,E_{\mathrm{kin}}}{3\,V} \;=\; \frac{N\,k_{\mathrm{B}}\,T}{V},
```

and solving for $E_{\mathrm{kin}}$ gives

```{math}
E_{\mathrm{kin}} \;=\; \frac{3}{2}\,N\,k_{\mathrm{B}}\,T.
```

**6. Average Kinetic Energy Per Particle (Equipartition Theorem)**  
Dividing by $N$ yields the average kinetic energy per particle:

````{important}
```{math}
\langle E_{\mathrm{kin}} \rangle 
\;=\; \frac{E_{\mathrm{kin}}}{N}
\;=\; \frac{1}{2}\,m\,\langle v^2 \rangle
\;=\; \frac{3}{2}\,k_{\mathrm{B}}\,T.
```
````

Each translational degree of freedom contributes $k_{\mathrm{B}}T/2$ to the average kinetic energy.
`````

---

### Maxwell–Boltzmann Speed Distribution

The rms speed $v_{\mathrm{rms}}$ is a useful *single-number* summary, but in thermal equilibrium a gas has a **distribution** of particle speeds.

For an ideal gas in three dimensions, the **Maxwell–Boltzmann speed distribution** gives the probability density $f(v)$ of molecular speeds, where $f(v)\,dv$ is the probability of finding a molecule with speed between $v$ and $v+dv$:

```{math}
:label: maxwell-boltzmann-speed-distribution-eq
f(v)
=
4\pi\left(\frac{m}{2\pi k_{\mathrm{B}}T}\right)^{3/2} v^2\,\exp\!\left(-\frac{m v^2}{2k_{\mathrm{B}}T}\right),
\qquad v\ge 0.
```

It is normalized so that $\int_0^\infty f(v)\,dv = 1$.

#### Most probable speed and mean speed

From $f(v)$ we can define several “typical” speeds:

```{math}
:label: most-probable-speed
v_{\mathrm{mp}} = \sqrt{\frac{2k_{\mathrm{B}}T}{m}} \qquad \text{(speed at the peak of } f(v)\text{)}
```

```{math}
:label: mean-speed
\langle v\rangle = \int_0^\infty v\,f(v)\,dv = \sqrt{\frac{8k_{\mathrm{B}}T}{\pi m}}.
```

The rms speed (from Eq. {eq}`rms-speed`) is

```{math}
v_{\mathrm{rms}} = \sqrt{\langle v^2\rangle} = \sqrt{\frac{3k_{\mathrm{B}}T}{m}}.
```

For any Maxwell–Boltzmann distribution, these satisfy

```{math}
v_{\mathrm{mp}} < \langle v\rangle < v_{\mathrm{rms}}.
```

#### Comparing gases and temperatures

The Maxwell–Boltzmann curves below illustrate two trends:

- **Lighter gases** (smaller $m$) have distributions that are shifted to higher speeds and are broader.
- **Higher temperatures** shift the distribution to higher speeds and make it broader.

% Try this: In the list `cases`, change the temperature of the dashed N$_2$ curve from `600` to `1200`, and its label to match. How do the position and the height of the maximum compare with those of the curve at 300 K?

```{code-cell} ipython3
:label: fig:maxwell-boltzmann
:caption: Maxwell–Boltzmann speed distributions of three gases at 300 K (solid) and of N$_2$ at 600 K (dashed).
:tags: [hide-input]

import numpy as np
from scipy.constants import k as k_B, N_A


def f_MB(v, M, T):
    """Maxwell–Boltzmann distribution of speeds (s/m) for molar mass M (kg/mol) at temperature T (K)."""
    m = M / N_A
    return 4 * np.pi * (m / (2 * np.pi * k_B * T)) ** 1.5 * v**2 * np.exp(-m * v**2 / (2 * k_B * T))


def draw(c):
    # Color names the gas. The dashed curve is the same gas at the higher temperature.
    cases = [(44.0095e-3, 300, c.vermilion, "-", r"CO$_2$, 300 K", (520, 2.25)),
             (28.0134e-3, 300, c.blue, "-", r"N$_2$, 300 K", (640, 1.82)),
             (28.0134e-3, 600, c.blue, "--", r"N$_2$, 600 K", (930, 1.20)),
             (4.002602e-3, 300, c.green, "-", "He, 300 K", (1680, 0.60))]
    v = np.linspace(0, 3000, 2000)
    fig, ax = plt.subplots(figsize=(6, 3.8), constrained_layout=True)
    for M, T, color, style, label, (x, y) in cases:
        ax.plot(v, 1e3 * f_MB(v, M, T), color=color, ls=style)
        ax.text(x, y, label, color=color, ha="left", va="bottom")
    ax.set(xlabel="Speed $v$ (m/s)", ylabel=r"Probability density $f(v)$ ($10^{-3}$ s/m)", xlim=(0, 3000), ylim=(0, 2.5))
    return fig


fs.show(draw, alt="Maxwell–Boltzmann distributions of molecular speeds from 0 to 3000 meters per second. Carbon dioxide "
                  "at 300 K has the narrowest curve, with its peak near 340 meters per second. Nitrogen at 300 K "
                  "peaks near 420, nitrogen at 600 K near 600, and helium at 300 K has the broadest curve, with its "
                  "peak near 1100 meters per second.")
```

### Gas in a Box

The simulation in [](#fig:gas-in-a-box) brings the results of this section together. Two hundred particles move in a cubic box and collide elastically with the walls and with one another. The pressure is calculated from the momentum that the particles deliver to the six walls, as in the derivation of Eq. {eq}`pressure-kinetic-theory`, and it is shown beside the ideal-gas value $N k_{\mathrm{B}} T / V$. The histogram of the particle speeds is drawn with the Maxwell–Boltzmann distribution of Eq. {eq}`maxwell-boltzmann-speed-distribution-eq`.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/gas-in-a-box.html
:label: fig:gas-in-a-box
:width: 100%
:title: Gas in a box: an interactive simulation of particles in a cubic box
:placeholder: ../_static/gas-in-a-box/gas-in-a-box.png

Particles of a gas in a cubic box, seen from one side. The pressure from the collisions with the walls is shown beside the ideal-gas value $N k_{\mathrm{B}} T / V$, as an average and for successive intervals of 10 ps. The bars are the speeds of the particles, and the curve is the Maxwell–Boltzmann distribution of Eq. [](#maxwell-boltzmann-speed-distribution-eq).
```

Three things to try:

1. Reduce the number of particles to 20, and compare the pressure in successive intervals with the pressure at 300 particles.
2. Give every particle the same speed, and watch how the histogram of the speeds changes with time.
3. Change the gas from N$_2$ to He at the same temperature, and compare the pressure and the speeds before and after the change.

## Worked Example

### Root-mean-square speed at room temperature

Estimate the rms speed of $\mathrm{N_2}$ molecules at $T=300\ \mathrm{K}$.

**Assumptions.** Ideal-gas kinetic theory and equipartition; $v_{\mathrm{rms}}=\sqrt{3k_{\mathrm{B}}T/m}$.  
Take $m(\mathrm{N_2}) = 28.0\,u = 28.0(1.66054\times10^{-27})\ \mathrm{kg} = 4.65\times10^{-26}\ \mathrm{kg}$.

1. **Insert numbers**

   ```{math}
   v_{\mathrm{rms}}=\sqrt{\frac{3k_{\mathrm{B}}T}{m}}
   =\sqrt{\frac{3(1.38065\times10^{-23}\ \mathrm{J/K})(300\ \mathrm{K})}{4.65\times10^{-26}\ \mathrm{kg}}}.
   ```

2. **Evaluate**

   ```{math}
   v_{\mathrm{rms}}=\sqrt{2.67\times10^{5}\ \mathrm{m^2/s^2}}
   \approx 5.17\times10^{2}\ \mathrm{m/s}.
   ```

**Result.** $v_{\mathrm{rms}}(\mathrm{N_2},300\ \mathrm{K})\approx 5.2\times10^{2}\ \mathrm{m/s}$.

## Concept Checks

1. Which kinetic-theory assumption is most directly violated at high pressures or low temperatures?
2. Why does pressure depend on $\langle v_x^2\rangle$ (or $\langle v^2\rangle$) rather than on $\langle v_x\rangle$?
3. How would doubling the absolute temperature change $v_{\mathrm{rms}}$?
4. What physical information is encoded in the *width* of the Maxwell–Boltzmann speed distribution?

## Key Takeaways

- Pressure arises from **momentum transfer** during particle–wall collisions.
- For an ideal gas, temperature measures **average translational kinetic energy**.
- Characteristic molecular speeds scale as $\sqrt{T/m}$.
- The Maxwell–Boltzmann distribution explains why gases contain a range of speeds even at fixed $T$.
