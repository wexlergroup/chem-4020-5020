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


# 5.1. Free Energy

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Free energies package the second law into practical extremum principles: under common laboratory constraints, certain potentials decrease spontaneously until equilibrium. This section reviews the differentials of $U,H,A,G$, identifies which potential is minimized for given constraints, and connects Helmholtz free energy to the canonical partition function via $A=-k_{\mathrm B}T\ln Q$.

These notes bring together (i) the *thermodynamic* idea of **equilibrium as an extremum principle** and (ii) the *statistical-mechanical* idea that macroscopic quantities can be computed from the **partition function**.

- Recognize which thermodynamic potential is “natural” for a set of constraints (e.g. constant $T,V$ vs. constant $T,P$).
- Use the differentials of $U$, $H$, $A$ (Helmholtz), and $G$ (Gibbs) to read off measurable quantities like $P$ and $S$.
- Connect macroscopic free energies to microscopic counting via
  $$
  A(T,V,N) = -k_B T\ln Q(T,V,N)
  $$
  and then compute properties such as $P$ from derivatives of $\ln Q$.

---

Learning objectives:

- Identify the natural variables and differentials of $U$, $H$, $A$, and $G$.
- State which thermodynamic potential is minimized at equilibrium under fixed $(S,V)$, $(S,P)$, $(T,V)$, and $(T,P)$.
- Derive $A=-k_{\mathrm B}T\ln Q$ from the canonical entropy expression and define the partition function notation used.
- Use $P = k_{\mathrm B}T(\partial\ln Q/\partial V)_T$ to connect microscopic counting to an equation of state.

## Core Ideas and Derivations

### Where we are at: “tools” on the stat mech side and the thermo side

#### Statistical mechanics (canonical ensemble)

For a system in contact with a heat bath at temperature $T$ (canonical ensemble):

- **Microstate probabilities**
  $$
  p_j = \frac{e^{-\beta E_j}}{Q},\qquad \beta\equiv \frac{1}{k_B T}
  $$
- **Ensemble averages**
  $$
  \langle A\rangle = \sum_j p_j A_j
  $$
- **Internal energy**
  $$
  U \equiv \langle E\rangle = \sum_j p_j E_j
  = -\left(\frac{\partial \ln Q}{\partial \beta}\right)_{N,V}
  $$
- **Constant-volume heat capacity**
  $$
  C_V \equiv \left(\frac{\partial U}{\partial T}\right)_{N,V}
  $$
- **Entropy (canonical)**
  $$
  S = \frac{U}{T} + k_B\ln Q
  $$

> In these notes the canonical partition function is written as $Q$ (many texts use $Z$).

#### Thermodynamics

- **First law (energy conservation)**
  $$
  dU = \delta q + \delta w
  $$
- **Entropy and reversibility**
  - reversible: $dS = \dfrac{\delta q_{\text{rev}}}{T}$
  - irreversible (non-equilibrium): the equality becomes an inequality (schematically, “$=$ becomes $<$”).

A key inequality you’ll use repeatedly (for a simple compressible system doing only $PV$ work) is

```{math}
:label: eq:fundamental-ineq
dU \le T\,dS - P\,dV
```

Equality holds for reversible changes; the inequality direction encodes spontaneity for irreversible processes.

---

### Thermodynamic equilibrium as an extremum principle

Equation {eq}`eq:fundamental-ineq` gives an “equilibrium test” once you specify what is held fixed.

#### Constant $S$ and $V$  ⇒  minimize $U$

If $S$ and $V$ are fixed, then {eq}`eq:fundamental-ineq` reduces to

$$
dU\le 0\quad (S,V\;\text{constant})
$$

So **at constant $S$ and $V$, spontaneous evolution drives $U$ downward** until equilibrium is reached at a minimum of $U$.

#### Identifying $T$ and $P$ from derivatives of $U(S,V)$

For a simple compressible system:

```{math}
T = \left(\frac{\partial U}{\partial S}\right)_V,\qquad
P = -\left(\frac{\partial U}{\partial V}\right)_S
```

These derivative definitions are consistent with the differential form
$dU = T\,dS - P\,dV$ (at equilibrium / reversible).

---

### Euler’s theorem and the Euler relation

Many thermodynamic functions are **homogeneous** in their extensive variables. If a function
$f(x_1,\dots,x_n)$ is homogeneous of degree $k$ (i.e., $f(sx_1,\dots,sx_n)=s^k f(x_1,\dots,x_n)$),
Euler’s theorem states

$$
k f = \sum_{i=1}^n x_i\left(\frac{\partial f}{\partial x_i}\right)_{\text{others}}
$$

For a simple system with fixed composition ($N$ fixed), $U(S,V)$ is homogeneous of degree 1, so

```{math}
:label: eq:euler-simple
U = S\left(\frac{\partial U}{\partial S}\right)_V
  + V\left(\frac{\partial U}{\partial V}\right)_S
  = TS - PV
```

Equation {eq}`eq:euler-simple` is the **Euler relation** for this simplified case (fixed $N$).

---

### Thermodynamic potentials (“free energies”)

Thermodynamic potentials are obtained from $U$ by Legendre transforms that swap “awkward” variables
for variables that are experimentally easier to control.

#### Summary table

| potential (name) | symbol | definition | differential (equilibrium) | natural variables |
| --- | ---: | --- | --- | --- |
| internal energy | $U$ | — | $dU = T\,dS - P\,dV$ | $U=U(S,V)$ |
| enthalpy | $H$ | $H\equiv U+PV$ | $dH = T\,dS + V\,dP$ | $H=H(S,P)$ |
| Helmholtz free energy | $A$ | $A\equiv U-TS$ | $dA = -S\,dT - P\,dV$ | $A=A(T,V)$ |
| Gibbs free energy | $G$ | $G\equiv H-TS = U+PV-TS$ | $dG = -S\,dT + V\,dP$ | $G=G(T,P)$ |

**Non-equilibrium reminder:** for spontaneous (irreversible) changes, the equalities above turn into
inequalities of the form “$d(\text{potential})\le \cdots$”, which leads directly to minimization
principles.

#### Which potential is minimized?

From the inequalities derived in lecture:

- **Constant $S,V$:** $dU\le 0$  ⇒  $U$ decreases to a minimum.
- **Constant $S,P$:** $dH\le 0$  ⇒  $H$ decreases to a minimum.
- **Constant $T,V$:** $dA\le 0$  ⇒  $A$ decreases to a minimum.
- **Constant $T,P$:** $dG\le 0$  ⇒  $G$ decreases to a minimum.

These are extremely useful because many laboratory conditions are naturally “constant $T$” and either
“constant $V$” (rigid container) or “constant $P$” (open to the atmosphere).

---

### Helmholtz free energy from the partition function

A key bridge between **macro** thermodynamics and **micro** statistical mechanics comes from combining

- $A\equiv U-TS$ (definition of Helmholtz free energy)
- $S = \dfrac{U}{T} + k_B\ln Q$ (canonical entropy expression)

Substitute $S$ into $A$:

```{math}
:label: eq:A-from-Q
A = U - T\left(\frac{U}{T} + k_B\ln Q\right)
  = -k_B T\ln Q
```

So once you know $Q(T,V,N)$, you can compute the thermodynamics.

---

### Pressure from the partition function

From the Helmholtz differential

$$
dA = -S\,dT - P\,dV
$$

we can read off

```{math}
:label: eq:P-from-A
P = -\left(\frac{\partial A}{\partial V}\right)_T
```

Using {eq}`eq:A-from-Q`:

```{math}
:label: eq:P-from-Q
P
= -\left(\frac{\partial (-k_B T\ln Q)}{\partial V}\right)_T
= k_B T\left(\frac{\partial \ln Q}{\partial V}\right)_T
```

**Takeaway:** knowing $Q$ gives access to $U$, $C_V$, $S$, $A$, and also mechanical equations of state like $P(V,T)$.

---

### Example: monatomic ideal gas

For a monatomic ideal gas (canonical ensemble), the partition function factorizes and the result can be written as

```{math}
:label: eq:ideal-gas-Q
Q(T,V,N) = \frac{1}{N!}\,\frac{V^N}{\Lambda^{3N}}
```

where $\Lambda$ is the **thermal de Broglie wavelength**

```{math}
:label: eq:thermal-Lambda
\Lambda = \frac{h}{\sqrt{2\pi m k_B T}}
```

#### Helmholtz free energy

From {eq}`eq:A-from-Q`,

$$
A = -k_B T\ln Q
    = -k_B T\left[-\ln N! + N\ln V - 3N\ln\Lambda\right]
$$

Using Stirling’s approximation ($\ln N!\approx N\ln N - N$), we get the common form

```{math}
:label: eq:ideal-gas-A
A \approx -N k_B T\left[\ln\left(\frac{V}{N\Lambda^3}\right)+1\right]
```

#### Pressure (recovering the ideal gas law)

From {eq}`eq:P-from-Q` and the fact that $\ln Q$ depends on $V$ as $N\ln V + \text{(no }$V$\text{)}$:

$$
P = k_B T\left(\frac{\partial \ln Q}{\partial V}\right)_T
  = k_B T\left(\frac{N}{V}\right)
$$

so

$$
PV = N k_B T
$$

---

### “Energy free to do work”: Gibbs free energy and non-$PV$ work

The Gibbs free energy is

$$
G = U + PV - TS
$$

Differentiate:

$$
dG = dU + P\,dV + V\,dP - T\,dS - S\,dT
$$

For a **reversible** change, write the first law as

$$
dU = T\,dS + \delta w_{\text{tot,rev}}
= T\,dS + \delta w_{PV,\text{rev}} + \delta w_{\text{non-}PV,\text{rev}}
$$

and for $PV$ work, $\delta w_{PV,\text{rev}}=-P\,dV$. Substitute into $dG$:

$$
dG
= \delta w_{\text{non-}PV,\text{rev}} + V\,dP - S\,dT
$$

Therefore, at constant $T$ and $P$,

```{math}
:label: eq:dG-nonPV
dG\Big|_{T,P} = \delta w_{\text{non-}PV,\text{rev}}
```

So under typical lab constraints (constant $T,P$), **changes in $G$ track reversible “useful” work** (electrical, surface, shaft work, etc.) *other than* $PV$ expansion/compression.

#### Example: water splitting (electrolysis)

Consider the (standard-state) water splitting reaction:

$$
\mathrm{H_2O(\ell) \rightarrow H_2(g) + \tfrac{1}{2}O_2(g)}
$$

The notes quote a standard Gibbs free energy change of roughly

$$
\Delta G^\circ(\text{room }T) \approx +237\ \text{kJ mol}^{-1}
$$

Because $\Delta G^\circ>0$, the reaction is **non-spontaneous** under standard conditions and must be driven by supplying non-$PV$ work (electrical work in an electrolyzer). In the reversible limit, the required electrical work per mole is on the order of $\Delta G^\circ$.

---

## Worked Example

### Recovering the ideal-gas law from $Q$

For a monatomic ideal gas,

```{math}
Q(T,V,N)=\frac{1}{N!}\left(\frac{V}{\Lambda^3}\right)^N,
\qquad
\Lambda=\frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}.
```

Use $P=k_{\mathrm B}T\left(\frac{\partial\ln Q}{\partial V}\right)_{T,N}$.

1. **Compute $\ln Q$**

   ```{math}
   \ln Q = -\ln N! + N\ln V - 3N\ln\Lambda.
   ```

2. **Differentiate with respect to $V$**
   Only the $N\ln V$ term depends on $V$, so

   ```{math}
   \left(\frac{\partial\ln Q}{\partial V}\right)_{T,N}=\frac{N}{V}.
   ```

3. **Insert into the pressure formula**

   ```{math}
   P = k_{\mathrm B}T\frac{N}{V}
   \quad\Rightarrow\quad
   PV = Nk_{\mathrm B}T,
   ```

   which is the ideal gas law.

**Result.** The microscopic translational partition function reproduces the macroscopic equation of state via derivatives of $\ln Q$.

## Concept Checks

1. Why do extremum principles involve *minimization* of a potential rather than minimization of entropy in most lab settings?
2. What changes in the inequalities when a process is irreversible rather than reversible?
3. Why does $A$ (not $G$) naturally appear for systems at fixed $T,V$?
4. How does the interpretation of $\Delta G$ as non-$PV$ work depend on reversibility?

## Key Takeaways

- Thermodynamic potentials are Legendre transforms of $U$ tailored to controllable variables.
- Equilibrium corresponds to minimization of $U,H,A,$ or $G$ depending on constraints.
- In the canonical ensemble, $A=-k_{\mathrm B}T\ln Q$ connects free energy to microscopic counting.
- Derivatives of $\ln Q$ yield measurable properties such as pressure.
