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

# Review

## 1. Checklist of Key Concepts

### Section 4.1: Entropy

See [Section 4.1](section-01.md).

1. **Exothermicity Does Not Guarantee Spontaneity**
   - $\Delta H < 0$ (exothermic) often favors spontaneity but does not determine it. Counter-example: mixing of two ideal gases is spontaneous even though $\Delta H = 0$.
   - A complete account of directionality requires **entropy**.

2. **Definition of Entropy**
   - Entropy is a state function defined by

   ```{math}
   dS = \frac{\delta q_{\mathrm{rev}}}{T}.
   ```

   - The finite change between states $A$ and $B$ is $\Delta S = \int_A^B \delta q_{\mathrm{rev}}/T$, evaluated along any reversible path connecting those states.

3. **Exactness of $\delta q_{\mathrm{rev}}/T$**
   - $\delta q_{\mathrm{rev}}$ is inexact, but dividing by $T$ produces an exact differential. Verification for an ideal gas: the cross-derivative test on $C_V/T$ and $Nk_{\mathrm{B}}/V$ gives matching mixed partials (both zero).

4. **Entropy of Irreversible Processes**
   - Because $S$ is a state function, $\Delta S$ depends only on the endpoints. To compute $\Delta S$ for an irreversible process, find any convenient reversible path between the same initial and final states and integrate $\delta q_{\mathrm{rev}}/T$ along that path.

5. **Fundamental Thermodynamic Relation**
   - Combining $\delta q_{\mathrm{rev}} = T\,dS$ and $\delta w_{\mathrm{rev}} = -P\,dV$ in the First Law:

   ```{math}
   dU = T\,dS - P\,dV.
   ```

   - The natural variables of $U$ are $S$ and $V$. This relation replaces two inexact differentials ($\delta q$, $\delta w$) with three exact ones ($dU$, $dS$, $dV$).
   - The relation is derived for a reversible process and holds for any change between equilibrium states, because it contains only state functions.

6. **Two Expansions Between the Same States** (Ideal Gas)
   - Reversible isothermal expansion: $\Delta S = nR\ln(V_2/V_1)$.
   - Free expansion into vacuum: $q = 0$ and $w = 0$, with the same $\Delta S$, because $S$ is a state function.

---

### Section 4.2: Carnot Cycle

See [Section 4.2](section-02.md).

1. **Structure of the Carnot Cycle**
   - A fully reversible cycle operating between a hot reservoir ($T_{\mathrm{hot}}$) and a cold reservoir ($T_{\mathrm{cold}}$), consisting of four steps:
     - $A \to B$: reversible isothermal expansion at $T_{\mathrm{hot}}$ (system absorbs heat $q_{AB} > 0$).
     - $B \to C$: reversible adiabatic expansion (system cools, $q = 0$).
     - $C \to D$: reversible isothermal compression at $T_{\mathrm{cold}}$ (system releases heat $q_{CD} < 0$).
     - $D \to A$: reversible adiabatic compression (system warms, $q = 0$).

2. **Entropy Bookkeeping Around the Cycle**
   - Adiabatic steps are isentropic: $\Delta S_{BC} = \Delta S_{DA} = 0$.
   - Isothermal steps carry entropy: $\Delta S_{AB} = q_{AB}/T_{\mathrm{hot}}$, $\Delta S_{CD} = q_{CD}/T_{\mathrm{cold}}$.
   - Over one complete cycle, $\Delta S_{\mathrm{cycle}} = 0$ (state function returns to initial value), so $\Delta S_{AB} = -\Delta S_{CD}$.

3. **Carnot Efficiency**
   - The efficiency of an engine is the fraction of the absorbed heat that becomes net work output: $\eta = W_{\mathrm{out}}/q_{AB} = (q_{AB} + q_{CD})/q_{AB}$.
   - For the Carnot cycle:

   ```{math}
   \eta_{\mathrm{Carnot}} = 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
   ```

   - This is the **maximum** possible efficiency for any engine operating between two temperatures. It depends only on the reservoir temperatures.

4. **Second Law for an Isolated System**
   - Spontaneous change requires $dS \ge 0$, with equality for a reversible change.

5. **Direction of Spontaneous Heat Flow**
   - Two subsystems at different temperatures inside an insulated boundary: $dS = dU_A(1/T_A - 1/T_B)$.
   - For $dS \ge 0$, heat must flow from hot to cold ($dU_A < 0$ when $T_A > T_B$).

---

### Section 4.3: Microscopic View of Entropy

See [Section 4.3](section-03.md).

1. **Entropy Increases Until Equilibrium** (Isolated Systems)
   - The second law for an isolated system: $dS \ge 0$, with equality at equilibrium.
   - Entropy rises during spontaneous processes and reaches its maximum at equilibrium.

2. **The Clausius Inequality**

   ```{math}
   dS \ge \frac{\delta q}{T}.
   ```

   - Equality for reversible processes; strict inequality for irreversible processes.
   - $T$ is the temperature at which the heat is transferred to the system (for a reservoir, the temperature of the reservoir).
   - Irreversibility **produces** additional entropy beyond the entropy carried by heat flow.
   - For an isolated system ($\delta q = 0$): reduces to $dS \ge 0$.
   - Summed over one cycle of any engine between two reservoirs, it gives $\eta \le \eta_{\mathrm{Carnot}}$.

3. **Boltzmann's Formula** (Microcanonical Ensemble)

   ```{math}
   S = k_{\mathrm{B}} \ln \Omega,
   ```

   where $\Omega$ is the number of accessible microstates at fixed $U$, $V$, $N$. The logarithm ensures entropy is extensive (additive for independent subsystems).

4. **Gibbs Entropy** (General Probability Distribution)

   ```{math}
   S = -k_{\mathrm{B}} \sum_i p_i \ln p_i.
   ```

   - Reduces to Boltzmann's formula when all $\Omega$ microstates are equally probable ($p_i = 1/\Omega$).
   - Applies to the canonical ensemble ($p_i = e^{-\beta E_i}/Q$) and any other ensemble.

5. **Canonical Entropy Identity**
   - Evaluating the Gibbs entropy with canonical probabilities yields

   ```{math}
   S = \frac{U}{T} + k_{\mathrm{B}} \ln Q.
   ```

   - Directly ties the macroscopic state function $S$ to the partition function $Q$ from Module 2.

6. **Helmholtz Free Energy**
   - Rearranging the canonical entropy identity:

   ```{math}
   A \equiv U - TS = -k_{\mathrm{B}}T\ln Q.
   ```

   - Natural variables of $A$ are $T$ and $V$.
   - All canonical thermodynamic quantities ($S$, $U$, $P$, $C_V$) can be derived from $A$ by differentiation.
   - Defining $A = U - TS$ is the same Legendre-transform strategy used to define $H = U + PV$ in [Section 3.3](../chapter-03/section-03.md): replace a hard-to-control variable ($S$) with an easy-to-control one ($T$).

7. **Microscopic Consistency**
   - For a reversible process, the Gibbs entropy with canonical probabilities gives $dS = \tfrac{1}{T}\sum_i E_i\,dp_i = \delta q_{\mathrm{rev}}/T$, the definition of entropy in Section 4.1.
   - This is consistent with the microscopic First Law from [Section 3.2](../chapter-03/section-02.md): $\delta q_{\mathrm{rev}} = \sum_i E_i\,dp_i$ and $\delta w_{\mathrm{rev}} = \sum_i p_i\,dE_i$.
   - The restriction to reversible processes matters: in a free expansion the entropy increases although no heat flows.

---

## 2. Checklist of Most Important Equations

Below is a unified list of the major equations from Sections 4.1–4.3.

A. **Entropy Definition**

```{math}
dS = \frac{\delta q_{\mathrm{rev}}}{T},
\qquad
\Delta S = \int_A^B \frac{\delta q_{\mathrm{rev}}}{T}.
```

- **Applicability**: any system. The integral must be evaluated along a *reversible* path, but $\Delta S$ itself depends only on the endpoints.

---

B. **Fundamental Thermodynamic Relation**

```{math}
dU = T\,dS - P\,dV.
```

- **Applicability**: simple closed system with $PV$-only work, for any change between equilibrium states. Natural variables of $U$ are $S$ and $V$.

---

C. **Entropy Change for an Ideal Gas** (from $V$ and $T$)

```{math}
\Delta S = \int_{T_1}^{T_2} \frac{C_V}{T}\,dT + Nk_{\mathrm{B}} \ln\!\left(\frac{V_2}{V_1}\right).
```

- For a monatomic ideal gas with constant $C_V = \tfrac{3}{2}Nk_{\mathrm{B}}$: $\Delta S = \tfrac{3}{2}Nk_{\mathrm{B}}\ln(T_2/T_1) + Nk_{\mathrm{B}}\ln(V_2/V_1)$.

---

D. **Carnot Efficiency**

```{math}
\eta_{\mathrm{Carnot}} = 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
```

- **Applicability**: maximum efficiency for any heat engine operating between two reservoir temperatures. The efficiency is defined as $\eta = W_{\mathrm{out}}/q_{AB}$.

---

E. **Clausius Inequality**

```{math}
dS \ge \frac{\delta q}{T}.
```

- Equality for reversible processes; strict inequality for irreversible processes. $T$ is the temperature at which the heat is transferred.
- For isolated systems ($\delta q = 0$): $dS \ge 0$.

---

F. **Boltzmann Entropy**

```{math}
S = k_{\mathrm{B}} \ln \Omega.
```

- **Applicability**: isolated systems (microcanonical ensemble) where all $\Omega$ accessible microstates are equally probable.

---

G. **Gibbs Entropy**

```{math}
S = -k_{\mathrm{B}} \sum_i p_i \ln p_i.
```

- **Applicability**: any probability distribution over microstates. Reduces to Boltzmann's formula in the microcanonical limit.

---

H. **Canonical Entropy Identity**

```{math}
S = \frac{U}{T} + k_{\mathrm{B}} \ln Q.
```

- **Applicability**: canonical ensemble (fixed $N$, $V$, $T$). Connects macroscopic entropy to the partition function.

---

I. **Helmholtz Free Energy**

```{math}
A = U - TS = -k_{\mathrm{B}}T \ln Q.
```

- **Applicability**: canonical ensemble. Natural variables are $T$ and $V$. All other canonical thermodynamic quantities can be obtained by differentiation.
