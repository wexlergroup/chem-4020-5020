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

# 4.3. Microscopic View of Entropy

This section shifts from the **macroscopic** entropy statements used with Carnot cycles to a more **microscopic / statistical** interpretation: entropy as a measure of *how many ways* a system can realize a macrostate, and how this connects to probability distributions and the partition function.

## Roadmap (from the lecture)

**Review** (from earlier sections)

1. For an **isolated** system:  

   ```{math}
   dS \ge 0.
   ```

2. For a **reversible** process (definition of entropy change):  

   ```{math}
   dS = \frac{\delta q_{\mathrm{rev}}}{T}.
   ```

3. A useful microscopic expression for heat (when energy levels are fixed):  

   ```{math}
   \delta q = \sum_{i=1}^{M} E_i\, d p_i.
   ```

4. For a **closed** system at fixed $N,V,T$ (canonical ensemble), the microstate probabilities are

   ```{math}
   p_i = \frac{e^{-\beta E_i}}{Q},\qquad \beta \equiv \frac{1}{k_B T}.
   ```

**New topics** (Class 29)

- Entropy change and the **approach to equilibrium**
- Entropy and the **number of microstates** (Boltzmann)  
- Entropy in terms of **probabilities** (Gibbs/Shannon form)  
- Entropy in terms of the **partition function** $Q$

---

## Entropy increases until the system reaches equilibrium

For an **isolated system**, the lecture emphasized this qualitative picture:

- **Out of equilibrium:** entropy increases spontaneously  

  ```{math}
  dS > 0.
  ```

- **At equilibrium:** entropy stops changing (it “flattens out”)  

  ```{math}
  dS = 0.
  ```

Equilibrium here means: **subject to the constraints** (fixed total energy, volume, particle number, etc.), the system has reached a state where it cannot increase its entropy any further.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

# Schematic "approach to equilibrium" curve
t = np.linspace(0, 10, 400)
S_i = 1.0
S_f = 2.0

# Smooth approach to plateau + slight wiggle to mimic "spontaneous process" steps
S = S_f - (S_f - S_i) * np.exp(-t/2.0) + 0.03*np.sin(5*t)*np.exp(-t/2.0)

fig, ax = plt.subplots(figsize=(4,3))
ax.plot(t, S)
ax.axhline(S_f, linestyle="--")
ax.set_xlabel("time (schematic)")
ax.set_ylabel("entropy (schematic)")
ax.set_title("Entropy increases until equilibrium (isolated system)")
ax.text(0.4, S_i+0.05, r"$S_i$")
ax.text(7.5, S_f+0.03, r"$S_f$")
ax.set_ylim(S_i-0.1, S_f+0.2)
ax.grid(True, alpha=0.3)
plt.show()
plt.close(fig)
```

This plot is schematic, but it matches the lecture sketch: entropy rises from $S_i$ to $S_f$ during a spontaneous process and then becomes constant at equilibrium.

---

## A more generic second law statement: the Clausius inequality

The definition

```{math}
dS = \frac{\delta q_{\mathrm{rev}}}{T}
```

is specifically written for a **reversible** heat transfer.

The notes then give a more general statement that covers both reversible and irreversible processes:

```{math}
dS \ge \frac{\delta q}{T}.
```

- **Equality** holds for a **reversible** process:

  ```{math}
  dS = \frac{\delta q_{\mathrm{rev}}}{T}.
  ```

- **Strict inequality** holds for an **irreversible** process:

  ```{math}
  dS > \frac{\delta q}{T}.
  ```

Interpretation: irreversibility produces additional entropy beyond the entropy “carried in” by heat flow.

> Practical takeaway: entropy is a **state function**, but $\delta q$ is not.  
> To compute $\Delta S$, you can always choose a *convenient reversible path* between the same initial and final states and use $\int \delta q_{\mathrm{rev}}/T$.

---

## Entropy and the number of microstates: Boltzmann’s formula

For an **isolated system**, the lecture connects entropy to the number of accessible microstates:

```{math}
S = k_B \ln \Omega.
```

- $k_B$ is the Boltzmann constant.
- $\Omega$ is the number of microstates compatible with the constraints of the macrostate.

The notes emphasize that, for an isolated system, $\Omega$ can be viewed as the **degeneracy** at the system’s value of the internal energy $U$:

- Fix $U$ (and $N,V$, etc.)
- Count: “How many different states did I find?”

### Why the logarithm?

A key property of thermodynamic entropy is that it is **extensive**:

```{math}
S_{\text{total}} = S_A + S_B
```

for two independent subsystems $A$ and $B$.

But the number of microstates is multiplicative for independent systems:

```{math}
\Omega_{\text{total}} = \Omega_A\,\Omega_B.
```

Taking the logarithm turns products into sums:

```{math}
S_{\text{total}}
= k_B \ln(\Omega_A\Omega_B)
= k_B \ln \Omega_A + k_B \ln \Omega_B
= S_A + S_B.
```

That is the “why ln?” point made in the notes.

---

## Entropy in terms of probabilities: the Gibbs entropy

The lecture then generalizes from “counting microstates” to working with a probability distribution $\{p_i\}$ over microstates $i=1,\dots,M$:

```{math}
S = -k_B \sum_{i=1}^{M} p_i\,\ln p_i.
```

This formula has two important features:

1. It reduces to Boltzmann’s $S=k_B\ln\Omega$ for the **microcanonical** case (all accessible microstates equally likely).  
   If $p_i = 1/\Omega$ for each accessible microstate, then

   ```{math}
   S
   = -k_B \sum_{i=1}^{\Omega} \frac{1}{\Omega}\ln\left(\frac{1}{\Omega}\right)
   = -k_B\ln\left(\frac{1}{\Omega}\right)
   = k_B\ln\Omega.
   ```

2. It makes sense even when not all microstates are equally likely (e.g., when the system is in contact with a heat bath).

---

## Canonical ensemble and the partition function

For a **closed system** (fixed $N,V$) in thermal contact with a reservoir at temperature $T$, the canonical ensemble assigns probabilities

```{math}
p_i(N,V,\beta) = \frac{e^{-\beta E_i(N,V)}}{Q(N,V,\beta)},
\qquad
\beta \equiv \frac{1}{k_B T}.
```

The normalization factor

```{math}
Q(N,V,\beta) \equiv \sum_{i=1}^{M} e^{-\beta E_i}
```

is the **canonical partition function**.

### Entropy in terms of $U$ and $Q$

Starting from the Gibbs entropy,

```{math}
S = -k_B \sum_i p_i\ln p_i,
```

insert $p_i = e^{-\beta E_i}/Q$. Since

```{math}
\ln p_i = -\beta E_i - \ln Q,
```

we get

```{math}
\begin{aligned}
S
&= -k_B \sum_i p_i\,(-\beta E_i - \ln Q) \\
&= k_B\beta \sum_i p_i E_i + k_B(\ln Q)\sum_i p_i.
\end{aligned}
```

Now use $\sum_i p_i = 1$ and define the internal energy

```{math}
U \equiv \sum_i p_i E_i.
```

Then

```{math}
S = k_B\beta U + k_B\ln Q.
```

Finally, because $\beta = 1/(k_B T)$,

```{math}
\boxed{
S = \frac{U}{T} + k_B\ln Q
}
```

which is the final result on the last page of the lecture notes.

---

## Microscopic interpretation of heat (review connection)

The lecture’s “review” line

```{math}
\delta q = \sum_i E_i\, d p_i
```

fits naturally with the statistical definition of internal energy

```{math}
U = \sum_i p_i E_i.
```

Differentiating $U$ gives

```{math}
dU = \sum_i p_i\, dE_i + \sum_i E_i\, dp_i.
```

This is a helpful way to separate energy changes:

- $\sum_i p_i\, dE_i$: changes due to shifting energy levels (often associated with **work**, e.g., changing volume changes $E_i(V)$)
- $\sum_i E_i\, dp_i$: changes due to redistributing probabilities among fixed levels (associated with **heat**)

This is the microscopic counterpart to the macroscopic first law bookkeeping.

---

## Key takeaways

- **Entropy increase and equilibrium:** for an isolated system, spontaneous evolution drives $S$ upward until equilibrium where $dS=0$.
- **General second law (Clausius inequality):** $dS \ge \delta q/T$, with equality only for reversible processes.
- **Boltzmann entropy:** $S=k_B\ln\Omega$ connects entropy to the multiplicity of microstates in an isolated system.
- **Gibbs entropy:** $S=-k_B\sum p_i\ln p_i$ generalizes entropy to non-uniform probability distributions.
- **Canonical ensemble:** $p_i = e^{-\beta E_i}/Q$ and $S=U/T + k_B\ln Q$ connect thermodynamic entropy to the partition function.
