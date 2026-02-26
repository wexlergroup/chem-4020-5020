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


# 3.1. Conservation of Energy

[Course-wide Conventions & Notation](../notation.md)

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

When an object (mass $m$) of height $h$ is dropped, its gravitational potential energy $mgh$ is converted into kinetic energy $\tfrac{1}{2}mv^2$, assuming negligible air resistance.

````{admonition} Apollo 15
:class: tip

During the Apollo 15 mission, astronaut David Scott famously dropped a hammer and a feather on the Moon. With negligible air resistance, both hit the ground simultaneously. This illustrates that in free fall, all objects (no matter their mass) accelerate at the same rate, and thus gain the same velocity over the same drop distance:

```{math}
\frac{1}{2}\,\cancel{m} v^2 = \cancel{m} g h
\quad\Rightarrow\quad
v = \sqrt{2 g h}.
```

```{image} ../_static/chapter-03/section-01/apollo-15-flag.jpg
:width: 50%
:align: center
```
````

### Mechanical Equivalent of Heat

```{figure} https://upload.wikimedia.org/wikipedia/commons/c/c3/Joule%27s_Apparatus_%28Harper%27s_Scan%29.png
---
name: joule-apparatus
width: 50%
align: center
---
Joule's apparatus for measuring the mechanical equivalent of heat. The falling weights turn a paddle wheel immersed in water, converting a known amount of mechanical work into a measurable temperature rise. Modern calorimeters (Section 3.3) use the same principle—measuring temperature changes to quantify energy transfers—but with electrical heating in place of falling weights.
```

In 1847, James Prescott Joule measured how mechanical work converts into heat. He famously found that dropping a total of $778$ lb$\cdot$ft of weight (e.g., by turning a paddle in water) raised the temperature of $1$ lb of water by $1\,\text{°F}$. Converting to SI units:

```{math}
\mathrm{MEH}
= \frac{778 \;\mathrm{lb}\cdot\mathrm{ft}}{1\;\mathrm{lb}\,\text{°F}}
\approx \frac{1{,}055 \;\mathrm{J}}{(453.6 \;\mathrm{g}) \cdot (5/9\,\text{°C})}
\approx 4.18 \;\mathrm{J}\,\mathrm{g}^{-1}\,\text{°C}^{-1}.
```

Here, $4.18\,\mathrm{J}\,\mathrm{g}^{-1}\,\text{°C}^{-1}$ is the **specific heat** of water—i.e., the heat capacity per gram.

### First Law of Thermodynamics

Joule's result showed that a definite quantity of mechanical work always produces the same quantity of heat. This interconvertibility means we need a single energy-conservation statement that accounts for both. To include heat ($q$) and work ($w$) in one statement of energy conservation, we write

```{math}
:label: first-law-finite
\Delta U = q + w,
```

where

- $\Delta U$ is the **change in internal energy**,
- $q$ is the heat **absorbed by** the system,
- $w$ is the work **done on** the system.

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

This distinction becomes central when we define entropy in Chapter 4, where dividing the inexact differential $\delta q$ by $T$ produces the exact differential $dS$.

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

In Section 1.1, we defined work as "energy transferred when a force acts over a distance." Mathematically:

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
* - Linear tension $k$
  - $l=x - x_0$
  - $k \,dl$
  - Stretching a spring
* - Surface tension $\gamma$
  - $A$
  - $\gamma \, dA$
  - Blowing a soap bubble
* - Pressure $P$
  - $V$
  - $-P\,dV$
  - Compressing a gas
* - Chemical $\mu$
  - $N$ or $n$
  - $\mu\,dN$
  - Forming a molecule
* - Electrical $\mathcal{E}$
  - Charge $q_{\text{el}}$
  - $\mathcal{E}\,dq_{\text{el}}$
  - Moving an electric charge
```

```{admonition} Signs in the work table
:class: note

All entries in the table above give $\delta w > 0$ when the surroundings do work on the system, consistent with our sign convention ($\Delta U = q + w$). For $PV$ work, the negative sign in $-P\,dV$ ensures that *compression* ($dV < 0$) gives $w > 0$ (energy flows into the system). For surface tension, *increasing* the area ($dA > 0$) requires the surroundings to do work on the system ($w > 0$), so the sign is positive: $\delta w = +\gamma\,dA$.

Note: some texts define surface work from the system's perspective (e.g., "work done by the system to expand the surface"), which flips the sign. Always check the convention.
```

In this course, $PV$ work ($-P\,dV$) will be our workhorse, but the generalized structure shows up again when we discuss surface phenomena, electrochemistry, and chemical potentials.

---

### Worked Examples

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

The $P$–$V$ diagram below illustrates this process geometrically. The magnitude of the work $|w|$ equals the shaded area under the constant-pressure line.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5, 3.5))

V_i = 20  # L
V_f = 85  # L
P_ext = 2.5  # bar

# Shaded region: area = |w|
ax.fill_between([V_i, V_f], 0, P_ext, color='C0', alpha=0.25, label=r'$|w| = P_{\mathrm{ext}}\,\Delta V$')

# Constant-pressure line
ax.plot([V_i, V_f], [P_ext, P_ext], 'C0-', lw=2)

# Vertical dashed lines at V_i and V_f
ax.plot([V_i, V_i], [0, P_ext], 'k--', lw=1, alpha=0.5)
ax.plot([V_f, V_f], [0, P_ext], 'k--', lw=1, alpha=0.5)

# Markers for initial and final states
ax.plot(V_i, P_ext, 'ko', ms=6, zorder=5)
ax.plot(V_f, P_ext, 'ko', ms=6, zorder=5)

# Arrow showing direction of expansion
ax.annotate('', xy=(V_f - 2, P_ext + 0.3), xytext=(V_i + 2, P_ext + 0.3),
            arrowprops=dict(arrowstyle='->', color='C3', lw=1.5))
ax.text((V_i + V_f) / 2, P_ext + 0.45, 'expansion', ha='center', fontsize=10, color='C3')

# Labels
ax.text(V_i, -0.25, r'$V_i$', ha='center', fontsize=11)
ax.text(V_f, -0.25, r'$V_f$', ha='center', fontsize=11)
ax.text(2, P_ext, r'$P_{\mathrm{ext}}$', ha='left', va='center', fontsize=11)

ax.set_xlabel('Volume (L)', fontsize=11)
ax.set_ylabel('Pressure (bar)', fontsize=11)
ax.set_xlim(0, 100)
ax.set_ylim(0, 4)
ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
ax.set_title('Constant-Pressure Expansion', fontsize=12)

plt.tight_layout()
plt.show()
plt.close(fig)
```

{numref}`Figure %s <joule-apparatus>` showed how Joule connected mechanical work to heat. Here, the $P$–$V$ diagram provides the geometric interpretation of work that we will use throughout this chapter: **work equals the area under the pressure curve on a $P$–$V$ diagram**. Different paths between the same initial and final volumes can enclose different areas, which is precisely why work is path-dependent.

---

#### Example 2. Expanding a Soap Bubble

Calculate the work necessary to expand a soap bubble (two surfaces) with surface tension $\gamma = 0.072\,\mathrm{J/m^2}$ ($= 72\,\mathrm{mN/m}$, the value for a water–air interface at room temperature) from a radius of $1\,\mathrm{cm}$ to $3.25\,\mathrm{cm}$.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch

fig, ax = plt.subplots(figsize=(5, 4))

# Bubble parameters
center = (0, 0)
r_outer = 2.0
r_inner = 1.8  # thin film
film_thickness = r_outer - r_inner

# Draw the outer surface
theta = np.linspace(0, 2 * np.pi, 300)
ax.plot(r_outer * np.cos(theta), r_outer * np.sin(theta), 'C0-', lw=2.5, label='Outer surface')

# Draw the inner surface
ax.plot(r_inner * np.cos(theta), r_inner * np.sin(theta), 'C1-', lw=2.5, label='Inner surface')

# Shade the film region
theta_fill = np.linspace(0, 2 * np.pi, 300)
ax.fill_between(
    r_outer * np.cos(theta_fill),
    r_outer * np.sin(theta_fill),
    r_inner * np.sin(theta_fill),
    alpha=0.15, color='C0'
)
# Fill the other half of the annulus
ax.fill(
    np.concatenate([r_outer * np.cos(theta_fill), r_inner * np.cos(theta_fill[::-1])]),
    np.concatenate([r_outer * np.sin(theta_fill), r_inner * np.sin(theta_fill[::-1])]),
    alpha=0.2, color='C4', label='Soap film'
)

# Radius arrow to outer surface
angle_deg = 35
angle_rad = np.radians(angle_deg)
ax.annotate(
    '', xy=(r_outer * np.cos(angle_rad), r_outer * np.sin(angle_rad)),
    xytext=center,
    arrowprops=dict(arrowstyle='->', color='k', lw=1.5)
)
# Label r
r_mid = (r_outer * 0.55)
ax.text(r_mid * np.cos(angle_rad) - 0.05, r_mid * np.sin(angle_rad) + 0.12,
        r'$r$', fontsize=14, ha='center', va='bottom')

# Zoom inset showing film thickness
# Draw a magnified bracket on the right side
inset_angle = 0  # radians (right side)
x_outer = r_outer * np.cos(inset_angle)
y_outer = r_outer * np.sin(inset_angle)
x_inner = r_inner * np.cos(inset_angle)
y_inner = r_inner * np.sin(inset_angle)

# Bracket lines for film thickness
bracket_len = 0.4
ax.plot([x_outer, x_outer + bracket_len], [y_outer + 0.05, y_outer + 0.05], 'k-', lw=1)
ax.plot([x_inner, x_inner + bracket_len], [y_inner - 0.05, y_inner - 0.05], 'k-', lw=1)
ax.annotate(
    '', xy=(x_outer + bracket_len - 0.05, y_outer + 0.05),
    xytext=(x_inner + bracket_len - 0.05, y_inner - 0.05),
    arrowprops=dict(arrowstyle='<->', color='k', lw=1.2)
)
ax.text(x_outer + bracket_len + 0.08, (y_outer + y_inner) / 2,
        'film', fontsize=10, ha='left', va='center', style='italic')

# Label "air inside" and "air outside"
ax.text(0, 0, 'air\ninside', fontsize=11, ha='center', va='center', color='0.4')
ax.text(0, -2.55, 'air outside', fontsize=11, ha='center', va='center', color='0.4')

# Surface labels
ax.text(-r_inner - 0.15, 0.8, 'inner\nsurface', fontsize=9, ha='right', va='center', color='C1')
ax.text(-r_outer - 0.15, -0.8, 'outer\nsurface', fontsize=9, ha='right', va='center', color='C0')

# Annotations for areas
ax.text(r_outer + 0.15, 1.4, r'$A_{\mathrm{outer}} = 4\pi r^2$', fontsize=11,
        ha='left', va='center', color='C0',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='C0', alpha=0.8))
ax.text(r_outer + 0.15, -1.4, r'$A_{\mathrm{total}} = 2 \times 4\pi r^2 = 8\pi r^2$', fontsize=11,
        ha='left', va='center', color='C4',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='C4', alpha=0.8))

ax.set_xlim(-3.2, 5.2)
ax.set_ylim(-3.2, 3.2)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Cross-section of a soap bubble', fontsize=12)

plt.tight_layout()
plt.show()
plt.close(fig)
```

Cross-section of a soap bubble showing the inner and outer surfaces. Because a bubble is a thin film enclosing air, its total surface area is $2 \times 4\pi r^2 = 8\pi r^2$, and the differential area change upon inflation is $dA = 16\pi r\,dr$.

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
   l = 0.15\,\mathrm{cm} = 0.15\times10^{-2}\,\mathrm{m}.
   ```

2. **Write the Work Expression**  
   Hooke's law for tension: $F = k\,\Delta l$. The infinitesimal work is:
   ```{math}
   \delta w = F\,dx = k(x - x_0)\,dx.
   ```

3. **Integrate**  
   If we stretch from $x_0$ to $x = x_0 + l$:
   ```{math}
   w 
   = \int_{x_0}^{x_0 + l} k\,(x - x_0)\,dx
   = \frac{1}{2} \,k\,l^2.
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

**Result.** In free expansion, $w = 0$, $q = 0$, and $\Delta U = 0$. For an ideal gas (where $U$ depends only on $T$), this means the temperature does not change either. Free expansion is the extreme case of an irreversible process: the gas does no work and exchanges no heat, yet it undergoes a dramatic change of state (its volume doubles).

```{admonition} Comparison
:class: tip

Contrast this with the constant-pressure expansion in Example 1, where the same volume change against $P_{\mathrm{ext}} = 2.5\,\mathrm{bar}$ required $16\,\mathrm{kJ}$ of work. Both processes have the same $\Delta V$, but the work differs because $P_{\mathrm{ext}}$ differs—a concrete illustration of why work is path-dependent.
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
- Many work modes share a generalized force–displacement structure, including $PV$ work $(-P\,dV)$.
- The magnitude of work equals the area under the pressure curve on a $P$–$V$ diagram. Different paths between the same states enclose different areas.
- Careful sign conventions prevent systematic mistakes in energy balances.
