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
  enumerator: 3.1.%s
---


# 3.1. Conservation of Energy

See also: [Course-wide Conventions & Notation](../notation.md)

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=3&amp;here=3.1"
  title="Course map, opened at Module 3"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-3.png
:class: hidden
:alt: Course map with Module 3 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6, and for Module 3 it is the internal energy as the ensemble average of the energy. Below the map, one card for each of Sections 3.1 to 3.3 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## Overview and Learning Objectives

Every chemical reaction either absorbs or releases energy—but where does that energy go, and how do we keep track of it? The First Law of Thermodynamics provides the bookkeeping framework. Its origins lie in 19th-century efforts to understand the relationship between mechanical work and heat, and its consequences reach into every area of chemistry: reaction energetics, calorimetry, materials processing, and beyond.

This section motivates the First Law from historical measurements, introduces sign conventions, and catalogs common forms of work as generalized force–displacement pairs.

Learning objectives:

- State the First Law in finite and differential forms and recognize that $q$ and $w$ are path functions while $U$ is a state function.
- Apply sign conventions to determine whether heat/work is done *on* or *by* the system.
- Compute $PV$ work for expansion/compression under specified external pressure conditions.
- Identify the generalized force and generalized displacement for a given work mode and write the corresponding $\delta w$ expression.

## Core Ideas and Derivations

### Conservation of Mechanical Energy

Between 1732 and 1736, Bernoulli and Euler combined the discoveries of Newton (laws of motion) and Leibniz (the connection between weight × vertical displacement and weight × velocity squared) into an early form of the law of conservation of mechanical energy. A simple example is the interchange between potential and kinetic energy:

```{math}
\frac{1}{2} m v^2 = m g h.
```

When an object of mass $m$ is dropped from a height $h$, its gravitational potential energy $mgh$ is converted into kinetic energy $\tfrac{1}{2}mv^2$, assuming negligible air resistance.

````{admonition} Apollo 15
:class: tip

During the Apollo 15 mission, astronaut David Scott famously dropped a hammer and a feather on the Moon. With negligible air resistance, both hit the ground simultaneously. This illustrates that in free fall, all objects (no matter their mass) accelerate at the same rate, and thus gain the same velocity over the same drop distance:

```{math}
\frac{1}{2}\,\cancel{m} v^2 = \cancel{m} g h
\quad\Rightarrow\quad
v = \sqrt{2 g h}.
```

```{figure} ../_static/chapter-03/section-01/apollo-15-flag.jpg
---
name: apollo-15
width: 50%
align: center
alt: An astronaut in a white space suit salutes beside the flag of the United States on the gray surface of the Moon. The lunar module stands behind the flag, the lunar roving vehicle is parked on the right, and a mountain rises in the background under a black sky.
---
Astronaut James Irwin salutes the flag at the Apollo 15 landing site, beside the lunar module and the lunar roving vehicle. Photograph by David Scott, NASA (AS15-88-11866, public domain).
```
````

### Mechanical Equivalent of Heat

```{figure} ../_static/chapter-03/section-01/joule-apparatus.png
---
name: joule-apparatus
width: 50%
align: center
---
Joule's apparatus for measuring the mechanical equivalent of heat. The falling weights turn a paddle wheel immersed in water, converting a known amount of mechanical work into a measurable temperature rise. Modern calorimeters ([Section 3.3](section-03.md)) use the same principle—measuring temperature changes to quantify energy transfers—but with electrical heating in place of falling weights. Engraving from *Harper's New Monthly Magazine* (1869), via [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Joule%27s_Apparatus_(Harper%27s_Scan).png) (public domain).
```

In 1847, James Prescott Joule measured how mechanical work converts into heat. He famously found that $778\,\mathrm{ft}\cdot\mathrm{lbf}$ of mechanical work (done, for example, by falling weights that turn a paddle in water) raised the temperature of $1$ lb of water by $1\,\text{°F}$. Converting to SI units:

```{math}
\mathrm{MEH}
= \frac{778 \;\mathrm{ft}\cdot\mathrm{lbf}}{1\;\mathrm{lb}\,\text{°F}}
\approx \frac{1{,}055 \;\mathrm{J}}{(453.6 \;\mathrm{g}) \cdot (5/9\,\text{°C})}
\approx 4.19 \;\mathrm{J}\,\mathrm{g}^{-1}\,\text{°C}^{-1}.
```

Here, $4.19\,\mathrm{J}\,\mathrm{g}^{-1}\,\text{°C}^{-1}$ is the **specific heat** of water—i.e., the heat capacity per gram. The specific heat of water varies slightly with temperature; at $25\,\text{°C}$ it is $4.18\,\mathrm{J}\,\mathrm{g}^{-1}\,\text{°C}^{-1}$, the value used in [Section 3.3](section-03.md).

### First Law of Thermodynamics

Joule's result showed that a definite quantity of mechanical work always produces the same quantity of heat. This interconvertibility means we need a single energy-conservation statement that accounts for both. To include heat ($q$) and work ($w$) in one statement of energy conservation for a closed system (one that exchanges energy but not matter with its surroundings), we write

```{math}
:label: first-law-finite
\Delta U = q + w,
```

where

- $\Delta U$ is the **change in internal energy**,
- $q$ is the heat **absorbed by** the system,
- $w$ is the work **done on** the system.

Some textbooks instead define $w$ as the work done **by** the system, so the sign convention of a source should be checked before one of its formulas is used.

#### Differential Form

In differential form,

```{math}
:label: first-law-differential
dU = \delta q + \delta w.
```

Rearranging,

```{math}
:label: first-law-rearranged
\boxed{\delta q = dU - \delta w.}
```

This rearranged form is useful because we often want to find the heat exchanged in a process. Since $U$ is a state function (calculable from state variables alone) and we can often evaluate $w$ from the process path, we can determine $q$ by difference.

```{admonition} Why $dU$ but $\delta q$ and $\delta w$?
:class: note

We write $dU$ (with the standard "$d$") because $U$ is a **state function**—its differential is *exact*, meaning $\oint dU = 0$ around any closed cycle. The value of $\Delta U$ depends only on the initial and final states, not on the path taken between them.

In contrast, $\delta q$ and $\delta w$ use "$\delta$" to signal that these are **inexact differentials**: their values depend on the path, and integrating them around a cycle need not give zero. Two processes connecting the same initial and final states can involve different amounts of heat and work, even though they produce the same $\Delta U$.

This distinction becomes central when we define entropy in [Section 4.1](../chapter-04/section-01.md), where dividing the reversible heat $\delta q_{\mathrm{rev}}$ by $T$ produces the exact differential $dS$.

Some texts use the crossed-$d$ symbol $\text{đ}$ instead of $\delta$ for inexact differentials. We use $\delta$ throughout this course.
```

````{admonition} **Quick check**
:class: important

A gas is compressed ($\Delta V < 0$) and simultaneously heated ($q > 0$). Is $\Delta U$ positive, negative, or indeterminate?

```{admonition} Answer
:class: dropdown

**Positive.** Compression means the surroundings do work on the gas, so $w > 0$. The gas also absorbs heat, so $q > 0$. Therefore $\Delta U = q + w > 0$. Both contributions increase the internal energy.
```
````

### Types of Work

In [Section 1.1](../chapter-01/section-01.md), we defined work as "energy transferred when a force acts over a distance." Mathematically:

```{math}
\delta w = \vec{F}\cdot d\vec{r} 
= (F_x,\,F_y,\,F_z)\cdot(dx,\,dy,\,dz) 
= F_x\,dx + F_y\,dy + F_z\,dz.
```

Many physical processes fit a "generalized force" $\times$ "generalized displacement" pattern:

```{list-table} Common forms of work as generalized force–displacement pairs
:header-rows: 1
:name: generalized-work-table

* - Generalized "Force"
  - Generalized "Displacement"
  - $\delta w$
  - Example
* - Mechanical $F$
  - $x$
  - $F \,dx$
  - Lifting a weight
* - Tension $f$
  - Extension $l = x - x_0$
  - $f \,dl$
  - Stretching a spring ($f = kl$ for a Hookean spring)
* - Surface tension $\gamma$
  - $A$
  - $\gamma \, dA$
  - Blowing a soap bubble
* - Pressure $P_{\mathrm{ext}}$
  - $V$
  - $-P_{\mathrm{ext}}\,dV$
  - Compressing a gas
* - Chemical potential $\mu$ ([Section 6.1](../chapter-06/section-01.md))
  - Number of molecules $N$
  - $\mu\,dN$
  - Adding molecules to the system
* - Electrical $\mathcal{E}$
  - Charge $q_{\text{el}}$
  - $\mathcal{E}\,dq_{\text{el}}$
  - Moving an electric charge
```

```{admonition} Signs in the work table
:class: note

All entries in the table above give $\delta w > 0$ when the surroundings do work on the system, consistent with our sign convention ($\Delta U = q + w$). For $PV$ work, the negative sign in $-P_{\mathrm{ext}}\,dV$ ensures that *compression* ($dV < 0$) gives $w > 0$ (energy flows into the system). For surface tension, *increasing* the area ($dA > 0$) requires the surroundings to do work on the system ($w > 0$), so the sign is positive: $\delta w = +\gamma\,dA$.

Note: some texts define surface work from the system's perspective (e.g., "work done by the system to expand the surface"), which flips the sign. Always check the convention.
```

In this course, $PV$ work will be our workhorse, but the generalized structure shows up again when we discuss surface phenomena, electrochemistry, and chemical potentials. The chemical term applies when molecules enter or leave the system, which is then open rather than closed.

The expression for $PV$ work follows from the definition of work. Consider a gas in a cylinder that is closed by a piston of area $A$. The surroundings exert the external pressure $P_{\mathrm{ext}}$ on the piston, which is a force $F = P_{\mathrm{ext}}A$ directed into the gas. When the piston moves inward by a distance $dx$, the surroundings do the work $\delta w = F\,dx = P_{\mathrm{ext}}A\,dx$ on the gas. The volume of the gas changes by $dV = -A\,dx$, so

```{math}
:label: pv-work
\delta w = -P_{\mathrm{ext}}\,dV.
```

The external pressure equals the pressure $P$ of the gas only when the gas remains in mechanical equilibrium with its surroundings throughout the process. This is the case for the reversible processes defined in [Section 3.2](section-02.md), for which $\delta w = -P\,dV$.

---

### Calculating Work

The examples below apply the generalized force–displacement framework from {numref}`generalized-work-table` to three different work modes: $PV$ work, surface-tension work, and elastic (spring) work.

```{admonition} Try it first
:class: important

Before looking at the Example 1 solution, try setting up the integral for work during a constant-pressure expansion. What are your integration limits? What is the integrand?
```

#### Example 1. Expanding Against Constant Pressure

Calculate the work when an ideal gas expands from $20\,\mathrm{L}$ to $85\,\mathrm{L}$ against a constant external pressure of $2.5\,\mathrm{bar}$.  

````{admonition} Solution
:class: dropdown

1. **Convert Units to SI**  
   ```{math}
   P_{\mathrm{ext}} = 2.5\,\mathrm{bar} = 2.5\times10^5\,\mathrm{Pa},\quad
   V_i = 20\,\mathrm{L} = 20\times10^{-3}\,\mathrm{m^3},\quad
   V_f = 85\,\mathrm{L} = 85\times10^{-3}\,\mathrm{m^3}.
   ```

2. **Set Up the Work Integral**  
   For a constant external pressure,
   ```{math}
   \delta w = -P_{\mathrm{ext}}\,dV 
   \quad \Longrightarrow\quad
   w = -\int_{V_i}^{V_f} P_{\mathrm{ext}}\,dV.
   ```

3. **Perform the Integration**  
   ```{math}
   w 
   = -P_{\mathrm{ext}}\int_{V_i}^{V_f} dV
   = -P_{\mathrm{ext}}\bigl(V_f - V_i\bigr).
   ```

4. **Numerical Result**  
   ```{math}
   w 
   = -\bigl(2.5\times10^5\,\mathrm{Pa}\bigr)
     \Bigl[\bigl(85\times10^{-3}\bigr) - \bigl(20\times10^{-3}\bigr)\Bigr]\mathrm{m^3}
   = -16{,}250\,\mathrm{J}
   = \boxed{-16\,\mathrm{kJ}}.
   ```

5. **Interpretation**  
   - The **negative sign** indicates the system (gas) does $16\,\mathrm{kJ}$ of work **on the surroundings**.
   - According to the convention $\Delta U = q + w$, if no heat is exchanged ($q=0$), the gas would lose $16\,\mathrm{kJ}$ in internal energy because it expands.
````

[](#fig:work-constant-pressure) shows this process on a $P$–$V$ diagram. The magnitude of the work $|w|$ equals the shaded area under the constant-pressure line.

```{code-cell} ipython3
:label: fig:work-constant-pressure
:caption: Expansion of a gas from $V_i = 20\,\mathrm{L}$ to $V_f = 85\,\mathrm{L}$ against a constant external pressure $P_{\mathrm{ext}} = 2.5\,\mathrm{bar}$ (Example 1). The magnitude of the work, $|w| = P_{\mathrm{ext}}\,\Delta V$, is the shaded area under the line.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

V_i = 20       # L
V_f = 85       # L
P_ext = 2.5    # bar


def draw(c):
    fig, ax = plt.subplots(figsize=(5.4, 3.3), constrained_layout=True)

    # Shaded region, whose area is |w|
    ax.fill_between([V_i, V_f], 0, P_ext, color=c.tint(c.blue, 0.75), lw=0)
    ax.text((V_i + V_f) / 2, P_ext / 2, r"$|w| = P_{\mathrm{ext}}\,\Delta V$", ha="center", va="center")

    # Initial and final volumes
    for V, name, ha, dx in [(V_i, "$V_i$", "left", 1.5), (V_f, "$V_f$", "right", -1.5)]:
        ax.plot([V, V], [0, P_ext], color=c.muted, lw=1.2, ls=(0, (4, 3)))
        ax.text(V + dx, 0.12, name, ha=ha, va="bottom")

    # Constant-pressure path from the initial state to the final state
    ax.plot([V_i, V_f], [P_ext, P_ext], color=c.blue)
    ax.plot([V_i, V_f], [P_ext, P_ext], "o", color=c.blue, ms=6.5, zorder=5)
    ax.text(V_i - 2.5, P_ext, r"$P_{\mathrm{ext}}$", ha="right", va="center")
    fs.arrow(ax, (V_i + 14, P_ext + 0.3), (V_f - 14, P_ext + 0.3), c.ink)
    ax.text((V_i + V_f) / 2, P_ext + 0.42, "expansion", ha="center", va="bottom")

    ax.set(xlabel="Volume $V$ (L)", ylabel="Pressure $P$ (bar)", xlim=(0, 100), ylim=(0, 3.5))
    return fig


fs.show(draw, alt="Pressure against volume. A horizontal line at an external pressure of 2.5 bar runs from the initial "
                  "volume of 20 L to the final volume of 85 L, with an arrow that marks the direction of the "
                  "expansion. The rectangle under the line, between the two volumes, is shaded, and its area is the "
                  "magnitude of the work.")
```

{numref}`Figure %s <joule-apparatus>` showed how Joule connected mechanical work to heat. Here, the $P$–$V$ diagram provides the geometric interpretation of work that we will use throughout this chapter: **work equals the area under the pressure curve on a $P$–$V$ diagram**. Different paths between the same initial and final volumes can enclose different areas, which is precisely why work is path-dependent.

---

#### Example 2. Expanding a Soap Bubble

Calculate the work necessary to expand a soap bubble (two surfaces) from a radius of $1\,\mathrm{cm}$ to $3.25\,\mathrm{cm}$. Treat the film as pure water, with surface tension $\gamma = 0.072\,\mathrm{J/m^2}$ ($= 72\,\mathrm{mN/m}$, the value for a water–air interface at room temperature). Soap lowers the surface tension of water, so the work for a real soap film is smaller.

```{code-cell} ipython3
:label: fig:soap-bubble
:caption: Cross-section of a soap bubble showing the inner and outer surfaces. Because a bubble is a thin film enclosing air, its total surface area is $2 \times 4\pi r^2 = 8\pi r^2$, and the differential area change upon inflation is $dA = 16\pi r\,dr$. The thickness of the film is exaggerated.
:tags: [hide-input]

from matplotlib.patches import Circle


def draw_bubble(c):
    fig, ax = plt.subplots(figsize=(5.8, 3.3), constrained_layout=True)
    r_outer, r_inner = 2.0, 1.8              # the film between the two surfaces, not to scale

    # Soap film between the outer surface and the inner surface
    ax.add_patch(Circle((0, 0), r_outer, fc=c.fill, ec=c.blue, lw=fs.LINE))
    ax.add_patch(Circle((0, 0), r_inner, fc=c.ground, ec=c.vermilion, lw=fs.LINE))

    # Radius
    angle = np.radians(35)
    fs.arrow(ax, (0, 0), (r_inner * np.cos(angle), r_inner * np.sin(angle)), c.ink)
    ax.text(0.95 * np.cos(angle) - 0.1, 0.95 * np.sin(angle) + 0.1, "$r$", ha="center", va="bottom")

    ax.text(0, -0.62, "air inside", ha="center", va="center", color=c.muted)
    ax.text(0, -2.45, "air outside", ha="center", va="center", color=c.muted)

    # Labels of the two surfaces and of the film, with leaders
    for y, radius, label, color in [(1.15, r_outer, r"outer surface, area $4\pi r^2$", c.blue),
                                    (0.35, r_inner, r"inner surface, area $4\pi r^2$", c.vermilion)]:
        x = np.sqrt(radius**2 - y**2)
        ax.plot([x, 2.55], [y, y], color=color, lw=1.0)
        ax.text(2.65, y, label, color=color, va="center")
    ax.plot([(r_outer + r_inner) / 2 * np.cos(np.radians(-14)), 2.55], [-0.46, -0.46], color=c.muted, lw=1.0)
    ax.text(2.65, -0.46, "soap film", color=c.muted, va="center")
    ax.text(2.65, -1.3, r"total area $A = 2 \times 4\pi r^2 = 8\pi r^2$", va="center")

    ax.set(xlim=(-2.2, 6.0), ylim=(-2.7, 2.2), aspect="equal")
    ax.axis("off")
    return fig


fs.show(draw_bubble, alt="Cross-section of a soap bubble, drawn as two concentric circles with the soap film between "
                         "them and air inside and outside. An arrow from the center marks the radius r. The outer "
                         "surface and the inner surface each have the area 4 pi r squared, so the total area of the "
                         "bubble is 8 pi r squared.")
```

````{admonition} Solution
:class: dropdown

1. **Convert Units**  
   ```{math}
   \gamma = 0.072\,\mathrm{J/m^2}, \quad
   r_i = 1\times10^{-2}\,\mathrm{m}, \quad
   r_f = 3.25\times10^{-2}\,\mathrm{m}.
   ```

2. **Identify the Surface Area Change**  
   - For a **single** spherical surface, $A = 4\pi r^2$, so $dA = 8\pi r\,dr$.  
   - A bubble has **two** surfaces (inner and outer), so its total area is $2 \times 4\pi r^2 = 8\pi r^2$. Hence,
     ```{math}
     \text{(bubble area)} = 8\pi r^2 
     \quad\Longrightarrow\quad
     dA = 16\pi r\,dr.
     ```

3. **Write the Work Expression**  
   From {numref}`generalized-work-table`, $\delta w = \gamma\,dA$. Here, the surroundings do work on the bubble film to increase its area:
   ```{math}
   w = \int_{A_i}^{A_f} \gamma \,dA 
     = \gamma \bigl(A_f - A_i\bigr).
   ```

4. **Perform the Integral**  
   Substituting $A = 8\pi r^2$ for the bubble:
   ```{math}
   \Delta A = A_f - A_i
   = 8\pi \bigl(r_f^2 - r_i^2\bigr).
   ```
   Thus,
   ```{math}
   w = \gamma\,\Delta A 
      = \gamma \,8\pi\bigl(r_f^2 - r_i^2\bigr).
   ```

5. **Calculate Numerically**  
   ```{math}
   w 
   = 0.072\,\mathrm{J/m^2}
     \,\times\,8\pi
     \Bigl[\bigl(3.25\times10^{-2}\bigr)^2 
           - \bigl(1\times10^{-2}\bigr)^2\Bigr]\mathrm{m^2}.
   ```

   Evaluating step by step:
   ```{math}
   r_f^2 - r_i^2 = (1.056 \times 10^{-3} - 1.000 \times 10^{-4})\,\mathrm{m^2}
   = 9.56\times10^{-4}\,\mathrm{m^2},
   ```
   ```{math}
   \Delta A = 8\pi \times 9.56\times10^{-4}\,\mathrm{m^2}
   = 2.40\times10^{-2}\,\mathrm{m^2},
   ```
   ```{math}
   w = 0.072 \times 2.40\times10^{-2}\,\mathrm{J}
   = \boxed{1.7\times10^{-3}\,\mathrm{J} \approx 1.7\,\mathrm{mJ}}.
   ```

6. **Interpretation**  
   - The positive sign means the surroundings do work **on** the system (the bubble film) to increase its surface area. This is consistent with everyday experience: you must blow air (do work) to inflate a soap bubble.
   - The magnitude is very small ($\sim$millijoules) because the surface tension of water is modest ($72\,\mathrm{mN/m}$) and the area change, while geometrically large ($\sim 240\,\mathrm{cm^2}$), is small in SI units.
````

---

#### Example 3. Stretching a Hookean Fiber

Calculate the work required to stretch a fiber obeying Hooke's law, with $k=100\,\mathrm{N/cm}$, by $0.15\,\mathrm{cm}$.

````{admonition} Solution
:class: dropdown

1. **Convert Units**  
   ```{math}
   k = 100\,\frac{\mathrm{N}}{\mathrm{cm}}
     = 100\times10^2\,\frac{\mathrm{N}}{\mathrm{m}}
     = 10^4\,\frac{\mathrm{N}}{\mathrm{m}}, 
     \quad
   l_f = 0.15\,\mathrm{cm} = 0.15\times10^{-2}\,\mathrm{m},
   ```
   where $l_f$ is the final extension.

2. **Write the Work Expression**  
   Hooke's law gives the tension at the extension $l = x - x_0$ as $f = kl$. From {numref}`generalized-work-table`, the infinitesimal work is:
   ```{math}
   \delta w = f\,dl = kl\,dl.
   ```

3. **Integrate**  
   If we stretch from $l = 0$ to $l = l_f$:
   ```{math}
   w
   = \int_{0}^{l_f} kl\,dl
   = \frac{1}{2} \,k\,l_f^2.
   ```

4. **Substitute Numbers**  
   ```{math}
   w 
   = \tfrac{1}{2}\,(10^4\,\mathrm{N/m})\,
     \bigl(1.5\times10^{-3}\,\mathrm{m}\bigr)^2 
   = \tfrac{1}{2}\times 10^4
     \times 2.25\times10^{-6}\,\mathrm{J}
   = 1.1\times10^{-2}\,\mathrm{J}.
   ```

5. **Interpretation**  
   - This is the energy required (work done **on** the system) to stretch the fiber by $0.15\,\mathrm{cm}$.
   - The positive sign indicates the system absorbs energy (an external force pulls on the fiber).
   - Note that unlike constant-pressure expansion, the Hookean force varies with displacement, producing the factor of $\tfrac{1}{2}$.
````

## Worked Example

### Free Expansion of an Ideal Gas

A rigid, insulated container is divided in half by a membrane. One side contains $n$ moles of an ideal gas at temperature $T$; the other side is evacuated. The membrane is punctured. Determine $w$, $q$, and $\Delta U$.

1. **Work**

   The gas expands into vacuum, so $P_{\mathrm{ext}} = 0$:

   ```{math}
   w = -\int_{V_i}^{V_f} P_{\mathrm{ext}}\,dV = 0.
   ```

2. **Heat**

   The container is insulated (adiabatic walls), so $q = 0$.

3. **Internal energy change**

   From the First Law:

   ```{math}
   \Delta U = q + w = 0 + 0 = 0.
   ```

**Result.** In free expansion, $w = 0$, $q = 0$, and $\Delta U = 0$. For an ideal gas (where $U$ depends only on $T$; see [Section 2.5](../chapter-02/section-05.md)), this means the temperature does not change either. Free expansion is the extreme case of an irreversible process ([Section 3.2](section-02.md)): the gas does no work and exchanges no heat, yet it undergoes a dramatic change of state (its volume doubles).

```{admonition} Comparison
:class: tip

Consider two ways to take an ideal gas from $20\,\mathrm{L}$ to $85\,\mathrm{L}$ with the same temperature at the start and at the end. In a free expansion into an evacuated volume of $65\,\mathrm{L}$, $w = 0$. In the expansion of Example 1 against $P_{\mathrm{ext}} = 2.5\,\mathrm{bar}$, carried out in thermal contact with surroundings that return the gas to its initial temperature, $w = -16\,\mathrm{kJ}$. The two processes connect the same initial and final states, so $\Delta U$ is the same for both, but the work differs because $P_{\mathrm{ext}}$ differs. This is a concrete illustration of why work depends on the path.
```

## Concept Checks

1. Why are $q$ and $w$ called *path functions* while $U$ is a *state function*?
2. In which sign convention does expansion work appear as positive? How would the First Law be written then?
3. Why does $PV$ work depend on the **external** pressure for irreversible expansions?
4. What is the physical meaning of the chemical work term $\mu\,dN$?
5. Two processes connect the same initial and final states. Process A involves more work than Process B. Which process involves more heat? (Assume only $PV$ work.)

## Key Takeaways

- The First Law ($\Delta U=q+w$) enforces energy conservation by bookkeeping heat and work.
- Heat and work are **path-dependent** ($\delta q$, $\delta w$ are inexact differentials); internal energy is a **state function** ($dU$ is an exact differential).
- Many work modes share a generalized force–displacement structure, including $PV$ work $(-P_{\mathrm{ext}}\,dV)$.
- The magnitude of work equals the area under the pressure curve on a $P$–$V$ diagram. Different paths between the same states enclose different areas.
- Careful sign conventions prevent systematic mistakes in energy balances.
