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

In Section 4.3 we derived the canonical entropy $S = U/T + k_{\mathrm B}\ln Q$ and showed that rearranging it gives the Helmholtz free energy $A = U - TS = -k_{\mathrm B}T\ln Q$. We noted that this makes $A$ a "master potential" — if you can compute $Q$, you can derive all other thermodynamic quantities by taking appropriate derivatives. This section develops that machinery.

We begin by combining the Clausius inequality (Section 4.3) with the First Law to obtain a fundamental inequality, $dU \le T\,dS - P\,dV$, that governs the direction of spontaneous change. From this single inequality we derive **extremum principles**: rules for which thermodynamic potential ($U$, $H$, $A$, or $G$) is minimized at equilibrium under a given set of constraints. These potentials are related to one another by **Legendre transforms** that swap "awkward" natural variables (like $S$) for experimentally controllable ones (like $T$).

With the potentials and their differentials in hand, we can read off measurable quantities as partial derivatives. In particular, we derive the pressure formula $P = k_{\mathrm B}T(\partial\ln Q/\partial V)_T$ that was previewed in Section 2.3, and we show that changes in the Gibbs free energy $G$ track reversible non-$PV$ work under typical laboratory conditions.

---

Learning objectives:

- Derive the fundamental inequality $dU \le T\,dS - P\,dV$ from the Clausius inequality.
- Identify the natural variables and equilibrium differentials of $U$, $H$, $A$, and $G$.
- State which thermodynamic potential is minimized at equilibrium under fixed $(S,V)$, $(S,P)$, $(T,V)$, and $(T,P)$.
- Use $A = -k_{\mathrm B}T\ln Q$ (Section 4.3) to derive the pressure formula $P = k_{\mathrm B}T(\partial\ln Q/\partial V)_T$.
- Interpret $\Delta G\big|_{T,P}$ as reversible non-$PV$ work.

## Core Ideas and Derivations

### The fundamental inequality

In Section 4.1 we derived the fundamental thermodynamic relation for a reversible, $PV$-only process:

```{math}
dU = T\,dS - P\,dV.
```

This is an *equality* that holds when the process is reversible. For irreversible processes we need the **Clausius inequality** (Section 4.3),

```{math}
dS \ge \frac{\delta q}{T},
```

which says that entropy can be *produced* by irreversibility beyond the entropy carried in by heat flow. Rearranging gives $\delta q \le T\,dS$.

For a closed system doing only $PV$ work against an external pressure, the First Law reads $dU = \delta q + \delta w$ with $\delta w = -P_{\mathrm{ext}}\,dV$. Since $P_{\mathrm{ext}} \le P$ during a spontaneous expansion and $P_{\mathrm{ext}} \ge P$ during a spontaneous compression, we have $\delta w \ge -P\,dV$ in general, with equality in the reversible limit. Combining the two inequalities:

```{math}
dU = \delta q + \delta w \le T\,dS - P\,dV.
```

We therefore arrive at the **fundamental inequality** for a simple compressible system:

```{math}
:label: eq:fundamental-ineq
dU \le T\,dS - P\,dV
```

Equality holds for reversible changes; the inequality encodes the direction of spontaneous evolution for irreversible processes.

```{admonition} Scope of the inequality
:class: note
Equation {eq}`eq:fundamental-ineq` assumes a closed system ($N$ fixed) with only $PV$ work. When non-$PV$ work is present (electrical, surface, etc.), additional terms appear on the right-hand side. We will use Eq. {eq}`eq:fundamental-ineq` to derive extremum principles below, and then treat non-$PV$ work separately when we discuss the Gibbs free energy.
```

---

### Thermodynamic equilibrium as an extremum principle

Equation {eq}`eq:fundamental-ineq` gives an "equilibrium test" once you specify what is held fixed. The idea is simple: if certain variables are constant, the inequality constrains what happens to the remaining ones.

#### Constant $S$ and $V$  ⇒  minimize $U$

If $S$ and $V$ are both fixed, then $dS = 0$ and $dV = 0$, so {eq}`eq:fundamental-ineq` reduces to

$$
dU \le 0 \quad (S, V\ \text{constant}).
$$

**At constant $S$ and $V$, spontaneous evolution drives $U$ downward** until equilibrium is reached at a minimum of $U$.

#### Identifying $T$ and $P$ from derivatives of $U(S,V)$

Since the fundamental relation at equilibrium is $dU = T\,dS - P\,dV$, we can read off

```{math}
T = \left(\frac{\partial U}{\partial S}\right)_V,\qquad
P = -\left(\frac{\partial U}{\partial V}\right)_S.
```

These partial-derivative identities are the natural-variable counterparts of the differential form. In Section 4.1 we noted that a state function's **natural variables** are those that appear as independent variables in its total differential; for $U$, these are $S$ and $V$. The expressions above show that the *conjugate* intensive variables ($T$ and $P$) are encoded as slopes of $U$ along its natural-variable axes.

---

### Euler's theorem and the Euler relation

Many thermodynamic functions are **homogeneous** in their extensive variables. If a function $f(x_1,\dots,x_n)$ is homogeneous of degree $k$ — meaning $f(sx_1,\dots,sx_n) = s^k f(x_1,\dots,x_n)$ for all $s > 0$ — then Euler's theorem states

$$
k\,f = \sum_{i=1}^{n} x_i \left(\frac{\partial f}{\partial x_i}\right)_{\text{others}}.
$$

For a simple system at fixed composition ($N$ fixed), the internal energy $U(S,V)$ is homogeneous of degree 1 in its extensive arguments. Applying Euler's theorem with $k=1$:

```{math}
:label: eq:euler-simple
U = S\left(\frac{\partial U}{\partial S}\right)_V
  + V\left(\frac{\partial U}{\partial V}\right)_S
  = TS - PV
```

Equation {eq}`eq:euler-simple` is the **Euler relation** for this simplified case (fixed $N$).

````{admonition} A curious consequence for $G$
:class: dropdown

At fixed $N$, the Euler relation $U = TS - PV$ combined with the definition $G = U + PV - TS$ gives $G = 0$ for a single-component system. This seemingly paradoxical result is an artifact of holding $N$ fixed and not including the $\mu\,dN$ term. Once we allow the number of particles to vary, the Euler relation becomes $U = TS - PV + \mu N$ and the Gibbs free energy becomes $G = \mu N$, which is the physically meaningful result. We will return to this when we introduce the chemical potential.
````

---

### Thermodynamic potentials

Holding $S$ or $V$ constant is experimentally awkward — most laboratory work is done at controlled $T$ and/or controlled $P$. The solution is to define new state functions whose natural variables match the experimentally convenient ones. This is the same strategy we used in Section 3.3 when we defined enthalpy $H = U + PV$ to simplify the First Law at constant pressure; now we extend it systematically.

Each new potential is obtained from $U$ by a **Legendre transform**: a mathematical operation that swaps an extensive variable for its conjugate intensive variable while preserving all thermodynamic information. For example, replacing the natural variable $S$ in $U(S,V)$ with its conjugate $T = (\partial U/\partial S)_V$ gives the Helmholtz free energy $A(T,V) = U - TS$.

#### Summary table

| potential | symbol | definition | differential (equilibrium) | natural variables |
| --- | ---: | --- | --- | --- |
| internal energy | $U$ | — | $dU = T\,dS - P\,dV$ | $S, V$ |
| enthalpy | $H$ | $U + PV$ | $dH = T\,dS + V\,dP$ | $S, P$ |
| Helmholtz free energy | $A$ | $U - TS$ | $dA = -S\,dT - P\,dV$ | $T, V$ |
| Gibbs free energy | $G$ | $H - TS = U + PV - TS$ | $dG = -S\,dT + V\,dP$ | $T, P$ |

Each differential can be verified by direct computation. For example, for the Helmholtz free energy:

$$
dA = d(U - TS) = dU - T\,dS - S\,dT = (T\,dS - P\,dV) - T\,dS - S\,dT = -S\,dT - P\,dV.
$$

The other differentials follow by the same approach.

````{margin}
```{note}
These differentials assume fixed $N$ and only $PV$ work. When we allow particle exchange, each differential acquires a $+\mu\,dN$ term.
```
````

#### Which potential is minimized?

Each thermodynamic potential satisfies an inequality analogous to Eq. {eq}`eq:fundamental-ineq`, derived by applying the same fundamental inequality to the appropriate Legendre-transformed quantity. The results are:

- **Constant $S, V$:** $dU \le 0$ ⇒ $U$ decreases to a minimum.
- **Constant $S, P$:** $dH \le 0$ ⇒ $H$ decreases to a minimum.
- **Constant $T, V$:** $dA \le 0$ ⇒ $A$ decreases to a minimum.
- **Constant $T, P$:** $dG \le 0$ ⇒ $G$ decreases to a minimum.

```{admonition} Deriving the $A$ and $G$ inequalities
:class: dropdown

To see why $dA \le 0$ at constant $T$ and $V$, start from {eq}`eq:fundamental-ineq`:

$$
dU \le T\,dS - P\,dV.
$$

At constant $T$, we can write $T\,dS = d(TS) - S\,dT = d(TS)$ (since $dT = 0$). Rearranging:

$$
dU - d(TS) \le -P\,dV
\quad\Rightarrow\quad
d(U - TS) \le -P\,dV
\quad\Rightarrow\quad
dA \le -P\,dV.
$$

At constant $V$ ($dV = 0$), this gives $dA \le 0$.

The Gibbs case is similar. At constant $T$ and $P$, start from {eq}`eq:fundamental-ineq`, use $T\,dS = d(TS)$ and $-P\,dV = -d(PV)$ (since both $T$ and $P$ are constant), and combine:

$$
dU - d(TS) + d(PV) \le 0
\quad\Rightarrow\quad
d(U - TS + PV) \le 0
\quad\Rightarrow\quad
dG \le 0.
$$
```

These extremum principles are immensely practical. Most chemistry and biology occurs at constant temperature (thermostatted lab or regulated body temperature) and either constant volume (rigid container, computational simulation) or constant pressure (open to the atmosphere). Under these conditions:

- **Constant $T, V$:** the Helmholtz free energy $A$ is minimized at equilibrium.
- **Constant $T, P$:** the Gibbs free energy $G$ is minimized at equilibrium.

---

### The Helmholtz free energy as a master potential

Section 4.3 established the central connection between the canonical partition function and thermodynamics:

```{math}
:label: eq:A-from-Q
A = -k_{\mathrm B}T\ln Q.
```

From this single equation, combined with the Helmholtz differential $dA = -S\,dT - P\,dV$, we can extract every equilibrium property of the system by differentiation.

#### Entropy

Reading off the coefficient of $dT$ in the Helmholtz differential:

```{math}
:label: eq:S-from-A
S = -\left(\frac{\partial A}{\partial T}\right)_V
```

This is consistent with the canonical entropy expression $S = U/T + k_{\mathrm B}\ln Q$ derived in Section 4.3, as can be verified by differentiating Eq. {eq}`eq:A-from-Q` directly.

#### Pressure

Reading off the coefficient of $dV$:

```{math}
:label: eq:P-from-A
P = -\left(\frac{\partial A}{\partial V}\right)_T
```

Substituting $A = -k_{\mathrm B}T\ln Q$:

```{math}
:label: eq:P-from-Q
P = k_{\mathrm B}T\left(\frac{\partial\ln Q}{\partial V}\right)_T
```

This is the result previewed (without derivation) in Section 2.3. The route is now clear: the Helmholtz differential tells us that pressure is $-(\partial A/\partial V)_T$, and the partition-function bridge $A = -k_{\mathrm B}T\ln Q$ converts that into a derivative of $\ln Q$.

#### Internal energy

From $A = U - TS$ and the expressions above:

```{math}
U = A + TS = -k_{\mathrm B}T\ln Q + T\left[-\left(\frac{\partial A}{\partial T}\right)_V\right]
```

After simplification (or using the Section 2.3 result directly):

```{math}
U = -\left(\frac{\partial\ln Q}{\partial\beta}\right)_{N,V}, \qquad \beta = \frac{1}{k_{\mathrm B}T}
```

#### Heat capacity

```{math}
C_V = \left(\frac{\partial U}{\partial T}\right)_{N,V}
```

**Takeaway.** Knowing $Q(T,V,N)$ gives access to $A$, $S$, $U$, $C_V$, and the equation of state $P(T,V)$ — all from derivatives of a single function.

---

### Example: monatomic ideal gas

To see the derivative machinery in action, we return to a system whose partition function we derived in Section 2.5. For a monatomic ideal gas:

```{math}
:label: eq:ideal-gas-Q
Q(T,V,N) = \frac{1}{N!}\,\frac{V^N}{\Lambda^{3N}},
\qquad
\Lambda = \frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}.
```

#### Helmholtz free energy

Using Eq. {eq}`eq:A-from-Q`:

$$
A = -k_{\mathrm B}T\ln Q
  = -k_{\mathrm B}T\!\left[-\ln N! + N\ln V - 3N\ln\Lambda\right].
$$

Applying Stirling's approximation ($\ln N! \approx N\ln N - N$):

```{math}
:label: eq:ideal-gas-A
A \approx -Nk_{\mathrm B}T\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + 1\right]
```

#### Pressure (recovering the ideal gas law)

From Eq. {eq}`eq:P-from-Q`, and noting that $\ln Q$ depends on $V$ only through the $N\ln V$ term:

$$
P = k_{\mathrm B}T\!\left(\frac{\partial\ln Q}{\partial V}\right)_T
  = k_{\mathrm B}T\!\left(\frac{N}{V}\right)
\quad\Rightarrow\quad
PV = Nk_{\mathrm B}T.
$$

The microscopic translational partition function reproduces the macroscopic equation of state. This was also shown in Section 2.5 using the pressure formula stated there; the difference is that we can now trace its origin through the Helmholtz free energy.

#### Entropy (the Sackur–Tetrode equation)

From Eq. {eq}`eq:S-from-A`, noting that $\Lambda \propto T^{-1/2}$ so $\partial(3N\ln\Lambda)/\partial T = -3N/(2T)$:

```{math}
:label: eq:sackur-tetrode
S = Nk_{\mathrm B}\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + \frac{5}{2}\right]
```

This is the **Sackur–Tetrode equation** for the entropy of an ideal monatomic gas — a result that connects statistical mechanics (through $\Lambda$ and $Q$) to a measurable thermodynamic quantity.

---

### Gibbs free energy and non-$PV$ work

The Gibbs free energy is defined as

$$
G = H - TS = U + PV - TS.
$$

Its equilibrium differential is

$$
dG = -S\,dT + V\,dP.
$$

But what if the system can do work *other than* $PV$ expansion/compression — for example, electrical work in a battery or an electrolyzer? To see how $G$ handles this, consider a reversible process that includes both $PV$ work and non-$PV$ work. The First Law gives

$$
dU = T\,dS + \delta w_{PV,\text{rev}} + \delta w_{\text{non-}PV,\text{rev}}
   = T\,dS - P\,dV + \delta w_{\text{non-}PV,\text{rev}}.
$$

Substituting into $dG = dU + P\,dV + V\,dP - T\,dS - S\,dT$:

$$
dG = \delta w_{\text{non-}PV,\text{rev}} + V\,dP - S\,dT.
$$

At constant $T$ and $P$:

```{math}
:label: eq:dG-nonPV
dG\Big|_{T,P} = \delta w_{\text{non-}PV,\text{rev}}
```

Under the most common laboratory conditions (constant $T$ and $P$), **changes in $G$ equal the reversible non-$PV$ work**. This is why $G$ is sometimes called the "free energy" — it measures the energy "free" to do useful work beyond unavoidable $PV$ work against the atmosphere.

```{admonition} Irreversible non-$PV$ work
:class: dropdown

Equation {eq}`eq:dG-nonPV` applies to the reversible limit. For an irreversible process at constant $T$ and $P$, the inequality $dG \le 0$ (from the extremum principle) combines with the presence of non-$PV$ work to give

$$
dG \le \delta w_{\text{non-}PV}.
$$

Since the chemistry convention assigns $w > 0$ for work done *on* the system and $w < 0$ for work done *by* the system, this says that the magnitude of useful work *output* from a spontaneous process is bounded above by $|\Delta G|$. Irreversibilities reduce the actual work output below this bound.
```

#### Example: water splitting (electrolysis)

The standard-state water-splitting reaction at room temperature is

$$
\mathrm{H_2O(\ell) \rightarrow H_2(g) + \tfrac{1}{2}O_2(g)},
\qquad
\Delta G^\circ \approx +237\ \mathrm{kJ\,mol^{-1}}.
$$

Because $\Delta G^\circ > 0$, the reaction is non-spontaneous under standard conditions. It must be driven by supplying non-$PV$ work — in this case, electrical work in an electrolyzer. In the reversible limit, the minimum electrical work required per mole equals $\Delta G^\circ$.

---

## Worked Examples

### Example 1: Recovering the ideal-gas law from $Q$

**Problem.** Use the pressure formula Eq. {eq}`eq:P-from-Q` and the monatomic ideal-gas partition function Eq. {eq}`eq:ideal-gas-Q` to derive the equation of state.

**Solution.**

1. **Compute $\ln Q$.**

   ```{math}
   \ln Q = -\ln N! + N\ln V - 3N\ln\Lambda.
   ```

2. **Differentiate with respect to $V$.**
   Only the $N\ln V$ term depends on $V$, so

   ```{math}
   \left(\frac{\partial\ln Q}{\partial V}\right)_{T,N} = \frac{N}{V}.
   ```

3. **Insert into the pressure formula.**

   ```{math}
   P = k_{\mathrm B}T\,\frac{N}{V}
   \quad\Rightarrow\quad
   PV = Nk_{\mathrm B}T.
   ```

**Result.** The microscopic translational partition function reproduces the ideal gas law via derivatives of $\ln Q$.

---

### Example 2: Entropy and internal energy of the ideal gas from $A$

**Problem.** Starting from the Helmholtz free energy of the monatomic ideal gas, Eq. {eq}`eq:ideal-gas-A`, derive expressions for $S$ and $U$.

**Solution.**

The Helmholtz free energy is

$$
A = -Nk_{\mathrm B}T\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + 1\right].
$$

Since $\Lambda = h/\sqrt{2\pi m k_{\mathrm B}T}$, we have $\Lambda \propto T^{-1/2}$ and $\ln\Lambda = -\tfrac{1}{2}\ln T + \text{const}$. Therefore

$$
A = -Nk_{\mathrm B}T\!\left[\ln V - \ln N + \tfrac{3}{2}\ln T + \text{const} + 1\right],
$$

where "const" collects terms independent of $T$ and $V$.

1. **Entropy.** Using $S = -(\partial A/\partial T)_V$:

   $$
   S = Nk_{\mathrm B}\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + \frac{5}{2}\right].
   $$

   This is the Sackur–Tetrode equation.

2. **Internal energy.** From $U = A + TS$:

   $$
   U = -Nk_{\mathrm B}T\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + 1\right]
     + T \cdot Nk_{\mathrm B}\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + \frac{5}{2}\right]
   = \frac{3}{2}Nk_{\mathrm B}T.
   $$

**Result.** Differentiating $A$ with respect to $T$ gives the entropy (Sackur–Tetrode), and adding back $TS$ recovers the equipartition result $U = \frac{3}{2}Nk_{\mathrm B}T$. All thermodynamic properties of the ideal gas follow from a single function, $A(T,V,N)$.

## Concept Checks

1. Why do extremum principles involve *minimization* of a potential rather than maximization of entropy in most lab settings?
2. In the derivation of Eq. {eq}`eq:fundamental-ineq`, which inequality comes from the Clausius inequality and which comes from the work bound? What physical process does each represent?
3. Why does $A$ (not $G$) naturally appear for systems at fixed $T$ and $V$?
4. Equation {eq}`eq:dG-nonPV` gives $dG\big|_{T,P} = \delta w_{\text{non-}PV,\text{rev}}$ for a reversible process. How does this change for an irreversible process, and what does the change mean for the maximum useful work a spontaneous reaction can deliver?
5. Why does differentiating the *same* function $A(T,V)$ give both a thermal quantity ($S$) and a mechanical quantity ($P$)?

## Key Takeaways

- The fundamental inequality $dU \le T\,dS - P\,dV$ combines the First Law with the Clausius inequality and governs the direction of spontaneous change.
- Thermodynamic potentials ($U$, $H$, $A$, $G$) are related by Legendre transforms that match natural variables to experimental constraints.
- Equilibrium corresponds to minimization of $U$, $H$, $A$, or $G$ depending on which variables are held fixed; $A$ at constant $T,V$ and $G$ at constant $T,P$ are the most common laboratory cases.
- The connection $A = -k_{\mathrm B}T\ln Q$ (Section 4.3) makes $A$ a master potential: $S$, $U$, $P$, and $C_V$ all follow from partial derivatives.
- Under typical lab conditions (constant $T,P$), $\Delta G$ measures the reversible non-$PV$ work — the energy "free" to do useful work.
