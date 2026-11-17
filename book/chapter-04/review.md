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

## 3. Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. Three reversible paths, one $\Delta S$

Draws on [Section 4.1](section-01.md).

One mole of oxygen, treated as an ideal gas with $C_{V,m}=\tfrac{5}{2}R$, is taken from state A ($T_A=300\ \mathrm{K}$, $V_A=10.0\ \mathrm{L}$) to state B ($T_B=420\ \mathrm{K}$, $V_B=40.0\ \mathrm{L}$). Three reversible paths connect the states:

- **Path 1.** The gas is heated at constant volume $V_A$ to $T_B$ and then expanded isothermally at $T_B$ to $V_B$.
- **Path 2.** The gas is expanded isothermally at $T_A$ to $V_B$ and then heated at constant volume $V_B$ to $T_B$.
- **Path 3.** The gas is compressed adiabatically until its temperature is $T_B$ and then expanded isothermally at $T_B$ to $V_B$.

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Along a reversible adiabat of this gas, $T\,V^{2/5}$ is constant ([Section 3.2](../chapter-03/section-02.md)).

**(a)** Find $q_{\mathrm{rev}}$ and $\Delta S$ for each step of Path 1.

**(b)** Find $q_{\mathrm{rev}}$ and $\Delta S$ for each step of Path 2.

**(c)** Add the steps to obtain $q_{\mathrm{rev}}$ and $\Delta S$ for Paths 1 and 2. Compare $\Delta S$ with the entropy change of an ideal gas from [Section 4.1](section-01.md).

**(d)** Find the volume at the end of the adiabatic step of Path 3. Then find $q_{\mathrm{rev}}$ and $\Delta S$ for Path 3.

````{admonition} Solution
:class: dropdown

1. **Relations for the two kinds of step**

   The heat capacity is $C_V=\tfrac{5}{2}(1\ \mathrm{mol})(8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}})=20.786\ \mathrm{J\,K^{-1}}$.

   - At constant volume, $\delta q_{\mathrm{rev}}=C_V\,dT$, so

     ```{math}
     q_{\mathrm{rev}}=C_V\,\Delta T,
     \qquad
     \Delta S=\int\frac{C_V}{T}\,dT=C_V\ln\left(\frac{T_{\mathrm{final}}}{T_{\mathrm{initial}}}\right).
     ```

   - At constant temperature, $\delta q_{\mathrm{rev}}=nRT\,dV/V$ (Example 1 of [Section 4.1](section-01.md)), so

     ```{math}
     \Delta S=nR\ln\left(\frac{V_{\mathrm{final}}}{V_{\mathrm{initial}}}\right),
     \qquad
     q_{\mathrm{rev}}=T\,\Delta S.
     ```

2. **Path 1, heating at constant volume ($300\ \mathrm{K}\to420\ \mathrm{K}$)**

   ```{math}
   \begin{aligned}
   q_{\mathrm{rev}}&=(20.786\ \mathrm{J\,K^{-1}})(120\ \mathrm{K})=2494\ \mathrm{J},\\
   \Delta S&=(20.786\ \mathrm{J\,K^{-1}})\ln\left(\frac{420}{300}\right)=6.994\ \mathrm{J\,K^{-1}}.
   \end{aligned}
   ```

3. **Path 1, isothermal expansion at $420\ \mathrm{K}$**

   ```{math}
   \begin{aligned}
   \Delta S&=(8.3145\ \mathrm{J\,K^{-1}})\ln\left(\frac{40.0}{10.0}\right)=11.526\ \mathrm{J\,K^{-1}},\\
   q_{\mathrm{rev}}&=(420\ \mathrm{K})(11.526\ \mathrm{J\,K^{-1}})=4841\ \mathrm{J}.
   \end{aligned}
   ```

4. **Path 2, isothermal expansion at $300\ \mathrm{K}$**

   The volume ratio is the same as in step 3, so $\Delta S=11.526\ \mathrm{J\,K^{-1}}$, and

   ```{math}
   q_{\mathrm{rev}}=(300\ \mathrm{K})(11.526\ \mathrm{J\,K^{-1}})=3458\ \mathrm{J}.
   ```

5. **Path 2, heating at constant volume ($300\ \mathrm{K}\to420\ \mathrm{K}$)**

   The heat capacity of an ideal gas does not depend on the volume, so this step has the values of step 2: $q_{\mathrm{rev}}=2494\ \mathrm{J}$ and $\Delta S=6.994\ \mathrm{J\,K^{-1}}$.

6. **Totals for Paths 1 and 2**

   | | $q_{\mathrm{rev}}$ (J) | $\Delta S$ ($\mathrm{J\,K^{-1}}$) |
   |---|---:|---:|
   | Path 1 | $2494+4841=7335$ | $6.994+11.526=18.52$ |
   | Path 2 | $3458+2494=5952$ | $11.526+6.994=18.52$ |

   For constant $C_V$, the entropy change of an ideal gas from [Section 4.1](section-01.md) is

   ```{math}
   \Delta S=C_V\ln\left(\frac{T_B}{T_A}\right)+nR\ln\left(\frac{V_B}{V_A}\right)=6.994+11.526=18.52\ \mathrm{J\,K^{-1}},
   ```

   which contains only the temperatures and volumes of the two states.

7. **Path 3, adiabatic compression**

   With $T_A V_A^{2/5}=T_B V'^{\,2/5}$, the volume $V'$ at the end of the adiabatic step is

   ```{math}
   V'=V_A\left(\frac{T_A}{T_B}\right)^{5/2}=(10.0\ \mathrm{L})\left(\frac{300}{420}\right)^{5/2}=4.312\ \mathrm{L}.
   ```

   On a reversible adiabat $\delta q_{\mathrm{rev}}=0$, so $q_{\mathrm{rev}}=0$ and $\Delta S=0$ for this step.

8. **Path 3, isothermal expansion at $420\ \mathrm{K}$**

   ```{math}
   \begin{aligned}
   \Delta S&=(8.3145\ \mathrm{J\,K^{-1}})\ln\left(\frac{40.0}{4.312}\right)=18.52\ \mathrm{J\,K^{-1}},\\
   q_{\mathrm{rev}}&=(420\ \mathrm{K})(18.52\ \mathrm{J\,K^{-1}})=7778\ \mathrm{J}.
   \end{aligned}
   ```

**Result.** The three paths have different heats and the same entropy change.

| | $q_{\mathrm{rev}}$ (J) | $\Delta S$ ($\mathrm{J\,K^{-1}}$) |
|---|---:|---:|
| Path 1 | 7335 | 18.52 |
| Path 2 | 5952 | 18.52 |
| Path 3 | 7778 | 18.52 |
````

### Problem 2. The four steps of a Carnot cycle

Draws on [Section 4.2](section-02.md).

A Carnot engine uses one mole of a diatomic ideal gas with $C_{V,m}=\tfrac{5}{2}R$ as its working substance. The reservoirs are at $T_{\mathrm{hot}}=600\ \mathrm{K}$ and $T_{\mathrm{cold}}=400\ \mathrm{K}$. The cycle has the four reversible steps of [Section 4.2](section-02.md): isothermal expansion $A\to B$ at $T_{\mathrm{hot}}$, adiabatic expansion $B\to C$, isothermal compression $C\to D$ at $T_{\mathrm{cold}}$, and adiabatic compression $D\to A$. The first two volumes are $V_A=4.00\ \mathrm{L}$ and $V_B=10.0\ \mathrm{L}$.

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Along a reversible adiabat of this gas, $T\,V^{2/5}$ is constant ([Section 3.2](../chapter-03/section-02.md)).

**(a)** Find $V_C$ and $V_D$.

**(b)** Find $q$ and $w$ for each of the four steps.

**(c)** Find the net work output and the efficiency from the heats of part (b). Compare the efficiency with $1-T_{\mathrm{cold}}/T_{\mathrm{hot}}$.

**(d)** Find $\Delta S$ of the gas for each step and for the whole cycle.

````{admonition} Solution
:class: dropdown

1. **Volumes at the ends of the adiabats**

   For step $B\to C$, $T_{\mathrm{hot}}V_B^{2/5}=T_{\mathrm{cold}}V_C^{2/5}$, so

   ```{math}
   V_C=V_B\left(\frac{T_{\mathrm{hot}}}{T_{\mathrm{cold}}}\right)^{5/2}=(10.0\ \mathrm{L})(1.5)^{5/2}=(10.0\ \mathrm{L})(2.756)=27.56\ \mathrm{L}.
   ```

   Step $D\to A$ connects the same two temperatures, so

   ```{math}
   V_D=V_A\left(\frac{T_{\mathrm{hot}}}{T_{\mathrm{cold}}}\right)^{5/2}=(4.00\ \mathrm{L})(2.756)=11.02\ \mathrm{L}.
   ```

   The two isotherms therefore have reciprocal volume ratios: $V_D/V_C=V_A/V_B=0.400$.

2. **Isothermal steps**

   For an ideal gas at constant temperature $\Delta U=0$, so $q=-w=nRT\ln(V_{\mathrm{final}}/V_{\mathrm{initial}})$.

   ```{math}
   \begin{aligned}
   q_{AB}&=(8.3145\ \mathrm{J\,K^{-1}})(600\ \mathrm{K})\ln(2.50)=+4571\ \mathrm{J},\\
   q_{CD}&=(8.3145\ \mathrm{J\,K^{-1}})(400\ \mathrm{K})\ln(0.400)=-3047\ \mathrm{J}.
   \end{aligned}
   ```

   The work is $w_{AB}=-4571\ \mathrm{J}$ and $w_{CD}=+3047\ \mathrm{J}$.

3. **Adiabatic steps**

   With $q=0$, $w=\Delta U=C_V\,\Delta T$, and $C_V=\tfrac{5}{2}(8.3145\ \mathrm{J\,K^{-1}})=20.786\ \mathrm{J\,K^{-1}}$.

   ```{math}
   \begin{aligned}
   w_{BC}&=(20.786\ \mathrm{J\,K^{-1}})(400-600)\ \mathrm{K}=-4157\ \mathrm{J},\\
   w_{DA}&=(20.786\ \mathrm{J\,K^{-1}})(600-400)\ \mathrm{K}=+4157\ \mathrm{J}.
   \end{aligned}
   ```

4. **Net work output and efficiency**

   The work of the two adiabatic steps cancels. In the notation of [Section 4.2](section-02.md),

   ```{math}
   W_{\mathrm{out}}=-w_{\mathrm{cycle}}=q_{AB}+q_{CD}=4571-3047=1524\ \mathrm{J},
   ```

   ```{math}
   \eta=\frac{W_{\mathrm{out}}}{q_{AB}}=\frac{1524}{4571}=0.333,
   \qquad
   1-\frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}=1-\frac{400}{600}=0.333.
   ```

5. **Entropy changes of the gas**

   On the isotherms, $\Delta S=q_{\mathrm{rev}}/T=nR\ln(V_{\mathrm{final}}/V_{\mathrm{initial}})$:

   ```{math}
   \begin{aligned}
   \Delta S_{AB}&=(8.3145\ \mathrm{J\,K^{-1}})\ln(2.50)=+7.618\ \mathrm{J\,K^{-1}},\\
   \Delta S_{CD}&=(8.3145\ \mathrm{J\,K^{-1}})\ln(0.400)=-7.618\ \mathrm{J\,K^{-1}}.
   \end{aligned}
   ```

   On the reversible adiabats, $\Delta S_{BC}=\Delta S_{DA}=0$. The sum over the cycle is $\Delta S_{\mathrm{cycle}}=0$.

**Result.** $V_C=27.56\ \mathrm{L}$ and $V_D=11.02\ \mathrm{L}$. The engine delivers $W_{\mathrm{out}}=1524\ \mathrm{J}$ per cycle with $\eta=0.333$, equal to $1-T_{\mathrm{cold}}/T_{\mathrm{hot}}$.

| Step | $q$ (J) | $w$ (J) | $\Delta U$ (J) | $\Delta S$ ($\mathrm{J\,K^{-1}}$) |
|---|---:|---:|---:|---:|
| $A\to B$ | +4571 | −4571 | 0 | +7.618 |
| $B\to C$ | 0 | −4157 | −4157 | 0 |
| $C\to D$ | −3047 | +3047 | 0 | −7.618 |
| $D\to A$ | 0 | +4157 | +4157 | 0 |
| Cycle | +1524 | −1524 | 0 | 0 |
````

### Problem 3. Heat flow between two bodies

Draws on [Section 4.1](section-01.md) and [Section 4.2](section-02.md).

A copper block of mass $150.0\ \mathrm{g}$ at $373.15\ \mathrm{K}$ is dropped into $250.0\ \mathrm{g}$ of water at $293.15\ \mathrm{K}$ inside an insulated container. The specific heat capacities are $0.385\ \mathrm{J\,g^{-1}\,K^{-1}}$ for copper and $4.184\ \mathrm{J\,g^{-1}\,K^{-1}}$ for water; treat them as constant. Neglect the thermal expansion of the copper and the water, so that the two bodies exchange only heat with each other, and neglect the heat capacity of the container.

**(a)** Find the final temperature.

**(b)** Find $\Delta S$ of the copper and $\Delta S$ of the water.

**(c)** Find the total entropy change, and compare its sign with the Second Law for an isolated system in [Section 4.2](section-02.md).

````{admonition} Solution
:class: dropdown

1. **Heat capacities**

   ```{math}
   \begin{aligned}
   C_{\mathrm{Cu}}&=(150.0\ \mathrm{g})(0.385\ \mathrm{J\,g^{-1}\,K^{-1}})=57.75\ \mathrm{J\,K^{-1}},\\
   C_{\mathrm{w}}&=(250.0\ \mathrm{g})(4.184\ \mathrm{J\,g^{-1}\,K^{-1}})=1046\ \mathrm{J\,K^{-1}}.
   \end{aligned}
   ```

2. **Final temperature**

   The container is insulated and no work is done, so the internal energy of the copper and the water together is constant:

   ```{math}
   C_{\mathrm{Cu}}\,(T_f-T_{\mathrm{Cu}})+C_{\mathrm{w}}\,(T_f-T_{\mathrm{w}})=0.
   ```

   Solving for $T_f$,

   ```{math}
   \begin{aligned}
   T_f&=\frac{C_{\mathrm{Cu}}T_{\mathrm{Cu}}+C_{\mathrm{w}}T_{\mathrm{w}}}{C_{\mathrm{Cu}}+C_{\mathrm{w}}}\\
   &=\frac{(57.75)(373.15)+(1046)(293.15)}{57.75+1046}\ \mathrm{K}=297.34\ \mathrm{K}.
   \end{aligned}
   ```

   The copper gives up $(57.75\ \mathrm{J\,K^{-1}})(373.15-297.34)\ \mathrm{K}=4378\ \mathrm{J}$ of heat to the water.

3. **Entropy change of each body**

   The flow of heat between the two bodies is irreversible. Entropy is a state function, so $\Delta S$ of each body is found on a reversible path between its initial and final states, on which the body is heated or cooled slowly. With no work, $\delta q_{\mathrm{rev}}=dU=C\,dT$, and

   ```{math}
   \Delta S=\int_{T_i}^{T_f}\frac{C}{T}\,dT=C\ln\left(\frac{T_f}{T_i}\right).
   ```

   Carrying one more figure in the final temperature, $T_f=297.336\ \mathrm{K}$,

   ```{math}
   \begin{aligned}
   \Delta S_{\mathrm{Cu}}&=(57.75\ \mathrm{J\,K^{-1}})\ln\left(\frac{297.336}{373.15}\right)=-13.12\ \mathrm{J\,K^{-1}},\\
   \Delta S_{\mathrm{w}}&=(1046\ \mathrm{J\,K^{-1}})\ln\left(\frac{297.336}{293.15}\right)=+14.83\ \mathrm{J\,K^{-1}}.
   \end{aligned}
   ```

4. **Total entropy change**

   ```{math}
   \Delta S=\Delta S_{\mathrm{Cu}}+\Delta S_{\mathrm{w}}=-13.12+14.83=+1.71\ \mathrm{J\,K^{-1}}.
   ```

   The copper and the water together form an isolated system, and the total is positive, in agreement with $dS\ge 0$ for a spontaneous change in an isolated system.

**Result.** $T_f=297.34\ \mathrm{K}$. The copper loses $13.12\ \mathrm{J\,K^{-1}}$ of entropy, the water gains $14.83\ \mathrm{J\,K^{-1}}$, and the total entropy increases by $1.71\ \mathrm{J\,K^{-1}}$.
````

### Problem 4. Electronic entropy of an ion with one $d$ electron

Draws on [Section 4.3](section-03.md) and [Section 2.2](../chapter-02/section-02.md).

The $\mathrm{Ti^{3+}}$ ion has a single $3d$ electron. Without spin–orbit coupling, its ten lowest states (five $d$ orbitals, each with two spin orientations) have the same energy. Spin–orbit coupling splits them into a lower level of four states and an upper level of six states, which lies $382.1\ \mathrm{cm^{-1}}$ above the lower level ([NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database)). The next level lies at $80\,389\ \mathrm{cm^{-1}}$ and is left out of this problem.

Treat one free ion at $T=300\ \mathrm{K}$ as the system, so that $Q$ is a sum over its ten states. Use $hc/k_{\mathrm B}=1.43878\ \mathrm{cm\,K}$, and give each entropy in units of $k_{\mathrm B}$.

**(a)** Find the entropy if the ten states have the same energy.

**(b)** With the splitting, find $Q$ and the probability of each state of the lower level and of the upper level.

**(c)** Find the Gibbs entropy from the probabilities of part (b).

**(d)** Use $S=U/T+k_{\mathrm B}\ln Q$ to find the entropy a second way, with $U$ measured from the lower level.

````{admonition} Solution
:class: dropdown

1. **Ten states of equal energy**

   The ten states are equally likely, $p_i=1/10$, so the Gibbs entropy reduces to Boltzmann's formula with $\Omega=10$:

   ```{math}
   S=k_{\mathrm B}\ln 10=2.303\,k_{\mathrm B}.
   ```

2. **The splitting in units of $k_{\mathrm B}T$**

   Let $\varepsilon$ be the energy of the upper level above the lower level.

   ```{math}
   \begin{aligned}
   \frac{\varepsilon}{k_{\mathrm B}}&=(1.43878\ \mathrm{cm\,K})(382.1\ \mathrm{cm^{-1}})=549.76\ \mathrm{K},\\
   \frac{\varepsilon}{k_{\mathrm B}T}&=\frac{549.76\ \mathrm{K}}{300\ \mathrm{K}}=1.8325.
   \end{aligned}
   ```

3. **Partition function and probabilities**

   Four states have the Boltzmann factor 1 and six states have $e^{-\varepsilon/k_{\mathrm B}T}=e^{-1.8325}=0.1600$.

   ```{math}
   Q=4+6\,e^{-\varepsilon/k_{\mathrm B}T}=4+6(0.1600)=4.960.
   ```

   Each state of the lower level has probability $p_a$, and each state of the upper level has probability $p_b$:

   ```{math}
   p_a=\frac{1}{4.960}=0.2016,
   \qquad
   p_b=\frac{0.1600}{4.960}=0.03226.
   ```

   The lower level holds $4p_a=0.8064$ of the probability and the upper level holds $6p_b=0.1936$.

4. **Gibbs entropy**

   ```{math}
   \begin{aligned}
   S&=-k_{\mathrm B}\left[4\,p_a\ln p_a+6\,p_b\ln p_b\right]\\
   &=-k_{\mathrm B}\left[0.8064\ln(0.2016)+0.1936\ln(0.03226)\right]\\
   &=k_{\mathrm B}\,(1.291+0.665)=1.956\,k_{\mathrm B}.
   \end{aligned}
   ```

5. **Entropy from the partition function**

   Only the six upper states contribute to the energy, so $U=6\,p_b\,\varepsilon=0.1936\,\varepsilon$ and

   ```{math}
   \frac{U}{k_{\mathrm B}T}=(0.1936)(1.8325)=0.355.
   ```

   Then

   ```{math}
   \frac{S}{k_{\mathrm B}}=\frac{U}{k_{\mathrm B}T}+\ln Q=0.355+\ln(4.960)=0.355+1.601=1.956.
   ```

**Result.** $S=2.303\,k_{\mathrm B}$ for ten states of equal energy. With the spin–orbit splitting, $Q=4.960$, $p_a=0.2016$, $p_b=0.03226$, and $S=1.956\,k_{\mathrm B}$ at $300\ \mathrm{K}$ by both routes.
````
