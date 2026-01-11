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

# Course-Wide Conventions & Notation

(notation)=
This page is a **one-page reference** for the conventions used throughout CHEM 4020/5020.

---

## Sign conventions (heat $q$ and work $w$)

We will use the **chemistry sign convention** throughout:

- **Heat** $q$ (or $\delta q$):
  - $q>0$ when heat is **absorbed by the system**.
  - $q<0$ when heat is **released by the system**.

- **Work** $w$ (or $\delta w$):
  - $w>0$ when work is **done on the system**.
  - $w<0$ when work is **done by the system**.

With this convention, the First Law reads

```{math}
\Delta U = q + w
\qquad\text{and}\qquad
dU = \delta q + \delta w.
```

```{admonition} Common alternative convention
:class: note
Some physics texts define $w$ as "work done by the system," yielding $\Delta U = q - w$.
Both conventions are equivalent as long as you are consistent.
```

### $P\,dV$ work (most common case)

For pressure–volume work against an external pressure $P_{\mathrm{ext}}$,

```{math}
\delta w = -P_{\mathrm{ext}}\,dV.
```

Therefore:

- Expansion ($dV>0$) gives $\delta w<0$ (system does work on surroundings).
- Compression ($dV<0$) gives $\delta w>0$ (surroundings do work on system).

For a **reversible** $PV$ process, $P_{\mathrm{ext}} = P$, so we often write
$\delta w_{\mathrm{rev}} = -P\,dV$.

---

## Standard state and the "$\circ$" symbol

A superscript $\circ$ denotes a **standard-state** quantity (e.g., $\mu^\circ$, $\Delta G^\circ$, $\Delta H^\circ$, $S^\circ$, $K^\circ$).

In this course, for gases, we take the standard-state pressure to be

```{math}
P^\circ = 1\ \mathrm{bar}.
```

Notes:

- Standard-state quantities can still depend on temperature (e.g., $\Delta G^\circ(T)$).
- When we use activities, the standard state corresponds to unit activity ($a=1$).

---

## Per-molecule vs. per-mole quantities ($k_{\mathrm{B}}$ vs. $R$)

We will be explicit about how we count:

- $N$ = number of molecules/particles
- $n$ = number of moles

The ideal-gas equation is written in either form:

```{math}
PV = N k_{\mathrm{B}} T = nRT.
```

The constants are related by

```{math}
R = N_{\mathrm{A}}\,k_{\mathrm{B}}.
```

Rule of thumb:

- Use $k_{\mathrm{B}}$ for **per-particle** (stat mech / microscopic) quantities.
- Use $R$ for **per-mole** (thermo / macroscopic) quantities.

---

## Reserved symbols and notation choices

Some symbols vary across textbooks. To minimize collisions, we adopt the following defaults.

### (A) $Z$ is reserved for compressibility

- **Compressibility factor**:

  ```{math}
  Z \equiv \frac{PV}{nRT}.
  ```

  - $Z=1$ for an ideal gas.
  - $Z\neq 1$ indicates non-ideal behavior.

We **do not** use $Z$ to denote a partition function.

### (B) Partition functions use $q$ and $Q$

- $q$ = **one-particle / one-molecule** partition function
- $Q$ = **many-particle (canonical)** partition function

(Many texts use $Z$ for a partition function; we will not.)

### (C) Chemical equilibrium uses $Q$ and $K$

- $Q$ (often written $Q_p$, $Q_c$, etc.) = **reaction quotient**
- $K$ = **equilibrium constant**

When the same letter appears in different contexts (e.g., $Q$ as a partition function vs. $Q$ as a reaction quotient), we will rely on **context and/or subscripts** ($Q(T,V,N)$ vs. $Q_p$, $Q_c$) to keep the meaning clear.

---

## Quick lookup table

| Concept | Symbol(s) | Meaning ("positive means…") / definition |
| - | - | - |
| Heat into system | $q$ | system absorbs heat |
| Work on system | $w$ | surroundings do work on system |
| $PV$ work | $\delta w=-P_{\mathrm{ext}}\,dV$ | expansion gives negative work |
| Standard-state pressure | $P^\circ$ | reference pressure (1 bar) |
| Per-particle constant | $k_{\mathrm{B}}$ | J K$^{-1}$ per particle |
| Per-mole constant | $R$ | J mol$^{-1}$ K$^{-1}$ |
| Compressibility factor | $Z$ | deviation from ideal gas ($Z=1$ ideal) |
| Partition functions | $q,\ Q$ | state counting in stat mech |
| Reaction quotient / equilibrium constant | $Q,\ K$ | composition measure / equilibrium value |
