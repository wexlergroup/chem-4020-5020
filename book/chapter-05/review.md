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

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=5&amp;here=review"
  title="Course map, opened at Module 5"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-5.png
:class: hidden
:alt: Course map with Module 5 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6, and for Module 5 it is the Helmholtz free energy in terms of the canonical partition function. Below the map, one card for each of Sections 5.1 to 5.3 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## 1. Checklist of Key Concepts

### Section 5.1: Free Energy

See [Section 5.1](section-01.md).

1. **Fundamental Inequality**
   - For any change between equilibrium states of a closed system with $PV$-only work, $dU = T\,dS - P\,dV$ ([Section 4.1](../chapter-04/section-01.md)).
   - For a closed system with $PV$-only work that is in thermal and mechanical equilibrium with its surroundings, but not necessarily in equilibrium internally, the First Law and the Clausius inequality $dS \ge \delta q/T$ give

   ```{math}
   dU \le T\,dS - P\,dV,
   ```

   with equality when every change in the system is reversible. The strict inequality holds when an irreversible process takes place inside the system, and it gives the direction of spontaneous change.

2. **Extremum Principles for Thermodynamic Potentials**
   - Under different experimental constraints, a particular potential is minimized at equilibrium:
     - Constant $S, V$: $dU \le 0$ ⇒ $U$ minimized.
     - Constant $S, P$: $dH \le 0$ ⇒ $H$ minimized.
     - Constant $T, V$: $dA \le 0$ ⇒ $A$ minimized.
     - Constant $T, P$: $dG \le 0$ ⇒ $G$ minimized.
   - The $T,V$ and $T,P$ cases cover the vast majority of laboratory conditions.

3. **Legendre Transforms and Natural Variables**
   - Each potential is obtained from $U$ by a Legendre transform that swaps an extensive variable for its conjugate intensive variable, matching the natural variables to the experimentally controllable ones.
   - Summary of natural variables and equilibrium differentials:

   | potential | definition | differential | natural variables |
   | --- | --- | --- | --- |
   | $U$ | — | $dU = T\,dS - P\,dV$ | $S, V$ |
   | $H$ | $U + PV$ | $dH = T\,dS + V\,dP$ | $S, P$ |
   | $A$ | $U - TS$ | $dA = -S\,dT - P\,dV$ | $T, V$ |
   | $G$ | $U + PV - TS$ | $dG = -S\,dT + V\,dP$ | $T, P$ |

   - Defining $A$ and $G$ extends the strategy used in [Section 3.3](../chapter-03/section-03.md) to define $H$.

4. **Conjugate Variables as Derivatives**
   - At equilibrium, intensive variables are read off as slopes of a potential along its natural-variable axes. For example:

   ```{math}
   T = \left(\frac{\partial U}{\partial S}\right)_V,\qquad
   P = -\left(\frac{\partial U}{\partial V}\right)_S,
   ```

   with analogous identities for $H$, $A$, and $G$.

5. **Euler Relation**
   - The internal energy $U(S,V,N)$ is homogeneous of degree 1 in its extensive arguments. Euler's theorem gives

   ```{math}
   U = TS - PV + \mu N,
   ```

   where $\mu = (\partial U/\partial N)_{S,V}$ is the chemical potential ([Section 6.1](../chapter-06/section-01.md)). With $G = U + PV - TS$, it follows that $G = \mu N$.

6. **Helmholtz Free Energy as a Master Potential**
   - The bridge to statistical mechanics from [Section 4.3](../chapter-04/section-03.md) is

   ```{math}
   A = -k_{\mathrm B}T\ln Q.
   ```

   - Combined with $dA = -S\,dT - P\,dV$, every equilibrium property follows by differentiation:

   ```{math}
   S = -\left(\frac{\partial A}{\partial T}\right)_V,\quad
   P = -\left(\frac{\partial A}{\partial V}\right)_T = k_{\mathrm B}T\!\left(\frac{\partial\ln Q}{\partial V}\right)_T,
   ```

   with $U = A + TS$ and $C_V = (\partial U/\partial T)_{N,V}$ following in turn.

7. **Ideal Gas from $A$**
   - For a monatomic ideal gas, $Q = V^N/(N!\,\Lambda^{3N})$ and Stirling's approximation give

   ```{math}
   A \approx -Nk_{\mathrm B}T\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + 1\right].
   ```

   - Differentiating yields the **Sackur–Tetrode entropy** $S = Nk_{\mathrm B}[\ln(V/N\Lambda^3) + 5/2]$, the equation of state $PV = Nk_{\mathrm B}T$, and the internal energy $U = \tfrac{3}{2}Nk_{\mathrm B}T$ — all from a single function.

8. **Gibbs Free Energy and Non-$PV$ Work**
   - When non-$PV$ work is present, the First Law gains a $\delta w_{\text{non-}PV}$ term. At constant $T$ and $P$:

   ```{math}
   dG\big|_{T,P} = \delta w_{\text{non-}PV,\mathrm{rev}}.
   ```

   - $\Delta G$ measures the reversible non-$PV$ work available from a process — the "useful" work beyond unavoidable $PV$ work against the atmosphere.

---

### Section 5.2: Third Law

See [Section 5.2](section-02.md).

1. **Planck's Statement of the Third Law**
   - As $T \to 0\,\mathrm{K}$, the entropy of any pure crystalline substance approaches a constant. For an **ideal crystal** (no defects, impurities, or disorder), that constant is zero:

   ```{math}
   S(0\,\mathrm{K}) = 0 \quad\text{(ideal crystal)}.
   ```

   - This fixes the otherwise-undetermined integration constant in $\Delta S = \int \delta q_{\mathrm{rev}}/T$.

2. **Statistical-Mechanical Justification**
   - From Boltzmann's formula $S = k_{\mathrm B}\ln\Omega$ ([Section 4.3](../chapter-04/section-03.md)), a unique ground state has $\Omega_0 = 1$, so $S(0) = k_{\mathrm B}\ln 1 = 0$. The third law is the statement that a system with a nondegenerate ground state has zero entropy at absolute zero.

3. **Residual Entropy**
   - If the ground state is degenerate or the system has frozen-in disorder ($\Omega_0 > 1$), entropy approaches a nonzero constant:

   ```{math}
   S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0.
   ```

   - Examples:
     - **Solid CO** (two nearly isoenergetic orientations per molecule): $\Omega_0 = 2^N$, giving $S_{\mathrm{res,m}} = R\ln 2 \approx 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (experiment: $\approx 4.6$, indicating partial ordering).
     - **Ice** (Pauling's model with ice rules): $\Omega_0 \approx (3/2)^N$, giving $S_{\mathrm{res,m}} = R\ln(3/2) \approx 3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (experiment: $\approx 3.41$).

4. **Heat Capacities Must Vanish as $T \to 0$**
   - For the integral $S(T) = S(0) + \int_0^T (C_P/T')\,dT'$ to converge, $C_P$ (and similarly $C_V$) must go to zero at $T = 0$.
   - This requirement can be compared with the classical (equipartition) heat capacity and with the Einstein heat capacity of [Section 2.6](../chapter-02/section-06.md).

5. **Unattainability of Absolute Zero** (Nernst form)
   - No finite sequence of thermodynamic operations can cool a system to exactly $T = 0$. Each successive cooling step removes less entropy than the previous one; the process converges but never reaches zero.

6. **Absolute (Third-Law) Entropies**
   - With $S(0) = 0$ fixed for ideal crystals, absolute entropies are computed by integrating $C_P/T$ from $0$ to $T$ and adding phase-transition contributions:

   ```{math}
   S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}}.
   ```

   - Tabulated $S^\circ(298.15\ \mathrm{K})$ values ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/), [NIST-JANAF Thermochemical Tables](https://janaf.nist.gov/)) are the basis for computing reaction entropies and Gibbs energies.

---

### Section 5.3: Ammonia Formation

See [Section 5.3](section-03.md).

1. **Thermodynamics vs. Kinetics**
   - $\Delta G$ determines which direction is favored at equilibrium; activation barriers determine *how fast* equilibrium is reached. These are independent questions.
   - For $\mathrm{N_2 + 3\,H_2 \to 2\,NH_3}$ at 298.15 K: $\Delta_r H^\circ \approx -92\ \mathrm{kJ\,mol^{-1}}$, $\Delta_r S^\circ \approx -198\ \mathrm{J\,K^{-1}\,mol^{-1}}$, $\Delta_r G^\circ \approx -33\ \mathrm{kJ\,mol^{-1}}$. Thermodynamically spontaneous but kinetically inert because of the $\mathrm{N\equiv N}$ triple bond.
   - A **catalyst** lowers the activation barrier without changing $\Delta_r G^\circ$ or the equilibrium constant — it accelerates forward and reverse reactions equally.

2. **Temperature Effect for $\Delta_r S^\circ < 0$**
   - From $\Delta_r G^\circ = \Delta_r H^\circ - T\,\Delta_r S^\circ$, a negative $\Delta_r S^\circ$ means $-T\,\Delta_r S^\circ > 0$ grows with $T$, driving $\Delta_r G^\circ$ upward.
   - At high enough $T$, $\Delta_r G^\circ$ changes sign and the reaction becomes non-spontaneous at 1 bar.

3. **Temperature Corrections to $\Delta H$ and $\Delta S$**
   - From constant-$P$ integration:

   ```{math}
   H(T_f) - H(T_i) = \int_{T_i}^{T_f} C_P\,dT,\qquad
   S(T_f) - S(T_i) = \int_{T_i}^{T_f} \frac{C_P}{T}\,dT.
   ```

   - With constant $C_P$: $\Delta H = C_P\,\Delta T$ and $\Delta S = C_P\ln(T_f/T_i)$.

4. **Ideal-Gas $C_P$ from Active Degrees of Freedom**
   - Using $f$ translational + rotational degrees of freedom (rigid-rotor, no vibrations):

   ```{math}
   C_P = \left(\frac{f}{2} + 1\right)R.
   ```

   - Linear $\mathrm{N_2, H_2}$: $f = 5$, $C_P = \tfrac{7}{2}R$.
   - Nonlinear $\mathrm{NH_3}$: $f = 6$, $C_P = 4R$.
   - For the ammonia reaction: $\Delta_r C_P = 2(4R) - (7R/2) - 3(7R/2) = -6R$.

5. **Result at 500 °C**
   - Applying constant-$C_P$ corrections from 298.15 K to 773.15 K gives $\Delta_r H^\circ(773) \approx -115.7\ \mathrm{kJ\,mol^{-1}}$, $\Delta_r S^\circ(773) \approx -245.5\ \mathrm{J\,K^{-1}\,mol^{-1}}$, and

   ```{math}
   \Delta_r G^\circ(773) \approx +74.1\ \mathrm{kJ\,mol^{-1}}.
   ```

   - The reaction is non-spontaneous at 1 bar — the central dilemma of Haber–Bosch.

6. **Pressure Dependence of $G$ (Ideal Gas)**
   - From $dG\big|_T = V\,dP$ with $V = RT/P$:

   ```{math}
   G(T, P) - G^\circ(T) = RT\ln\!\left(\frac{P}{P^\circ}\right).
   ```

   - For a reaction with stoichiometric change $\Delta\nu$ in moles of gas:

   ```{math}
   \Delta_r G(T, P) \approx \Delta_r G^\circ(T) + \Delta\nu\,RT\ln\!\left(\frac{P}{P^\circ}\right).
   ```

   - For $\mathrm{N_2 + 3\,H_2 \to 2\,NH_3}$: $\Delta\nu = 2 - 4 = -2$, so increasing $P$ lowers $\Delta_r G$ (Le Châtelier: high pressure favors the side with fewer gas molecules).

7. **Restoring Spontaneity with Pressure**
   - Setting $\Delta_r G = 0$ and solving for $P$:

   ```{math}
   P = P^\circ\exp\!\left(-\frac{\Delta_r G^\circ(T)}{\Delta\nu\,RT}\right).
   ```

   - At 500 °C: $P \approx 320\ \mathrm{bar}$. At 400 °C: $P \approx 86\ \mathrm{bar}$. The Haber–Bosch operating window (400–500 °C, 150–350 bar) is the engineering compromise predicted by these estimates.

8. **Engineering Knobs Summary**
   - **Catalyst**: speeds up both directions; does *not* change $\Delta_r G^\circ$ or $K$.
   - **Temperature**: raises rate (Arrhenius) but pushes $\Delta_r G^\circ$ upward when $\Delta_r S^\circ < 0$.
   - **Pressure**: compensates for the temperature penalty when $\Delta\nu < 0$.

---

## 2. Checklist of Most Important Equations

Below is a unified list of the major equations from Sections 5.1–5.3.

A. **Fundamental Inequality**

```{math}
dU \le T\,dS - P\,dV.
```

- **Applicability**: closed system with $PV$-only work, in thermal and mechanical equilibrium with its surroundings. Equality when every change in the system is reversible; strict inequality when an irreversible process takes place inside the system.

---

B. **Thermodynamic Potentials and Their Differentials** (fixed $N$, $PV$-only work)

```{math}
\begin{aligned}
dU &= T\,dS - P\,dV,\\
dH &= T\,dS + V\,dP,\\
dA &= -S\,dT - P\,dV,\\
dG &= -S\,dT + V\,dP.
\end{aligned}
```

- Natural variables: $U(S,V)$, $H(S,P)$, $A(T,V)$, $G(T,P)$.

---

C. **Extremum Principles** (equilibrium conditions)

```{math}
dU\big|_{S,V} \le 0,\quad
dH\big|_{S,P} \le 0,\quad
dA\big|_{T,V} \le 0,\quad
dG\big|_{T,P} \le 0.
```

- The potential whose natural variables match the imposed constraints is minimized at equilibrium.

---

D. **Helmholtz Free Energy as a Master Potential**

```{math}
A = -k_{\mathrm B}T\ln Q,\qquad
S = -\left(\frac{\partial A}{\partial T}\right)_V,\qquad
P = -\left(\frac{\partial A}{\partial V}\right)_T.
```

- **Applicability**: canonical ensemble. All equilibrium thermodynamic quantities follow by differentiation of a single function $A(T,V,N)$.

---

E. **Pressure from the Partition Function**

```{math}
P = k_{\mathrm B}T\!\left(\frac{\partial\ln Q}{\partial V}\right)_T.
```

- **Applicability**: canonical ensemble. Recovers $PV = Nk_{\mathrm B}T$ for the monatomic ideal gas.

---

F. **Sackur–Tetrode Entropy** (monatomic ideal gas)

```{math}
S = Nk_{\mathrm B}\!\left[\ln\!\left(\frac{V}{N\Lambda^3}\right) + \frac{5}{2}\right],\qquad
\Lambda = \frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}.
```

- **Applicability**: ideal monatomic gas in the classical (high-$T$, low-density) limit where Stirling's approximation and Boltzmann statistics apply.

---

G. **Gibbs Free Energy and Non-$PV$ Work**

```{math}
dG\big|_{T,P} = \delta w_{\text{non-}PV,\mathrm{rev}}.
```

- **Applicability**: reversible process at constant $T$ and $P$. $\Delta G$ equals the reversible non-$PV$ work.

---

H. **Third Law**

```{math}
S(0\,\mathrm{K}) = 0 \quad\text{(ideal crystal)},\qquad
S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0 \quad\text{(degenerate ground state)}.
```

- Fixes the entropy reference point; residual entropy reflects ground-state multiplicity.

---

I. **Absolute (Third-Law) Entropy**

```{math}
S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}}.
```

- **Applicability**: any substance whose heat capacity and phase-transition enthalpies are known from 0 to $T$; supplies tabulated $S^\circ(298.15\ \mathrm{K})$ values.

---

J. **Temperature Corrections to Reaction Thermodynamics** (constant $C_P$)

```{math}
\Delta_r H^\circ(T_f) = \Delta_r H^\circ(T_i) + \Delta_r C_P(T_f - T_i),
```

```{math}
\Delta_r S^\circ(T_f) = \Delta_r S^\circ(T_i) + \Delta_r C_P\ln\!\left(\frac{T_f}{T_i}\right).
```

- **Applicability**: reaction data from a reference temperature $T_i$ to a target temperature $T_f$, assuming temperature-independent heat capacities over the interval.

---

K. **Ideal-Gas Pressure Dependence of $G$**

```{math}
G(T,P) - G^\circ(T) = RT\ln\!\left(\frac{P}{P^\circ}\right).
```

- **Applicability**: pure ideal gas at constant $T$. For a reaction with stoichiometric gas change $\Delta\nu$:

```{math}
\Delta_r G(T,P) \approx \Delta_r G^\circ(T) + \Delta\nu\,RT\ln\!\left(\frac{P}{P^\circ}\right).
```

- For $\Delta\nu < 0$, higher pressure lowers $\Delta_r G$ and favors products.

## 3. Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. A van der Waals fluid from its Helmholtz free energy

Draws on [Section 5.1](section-01.md) and [Section 1.4](../chapter-01/section-04.md).

A model for $n$ moles ($N$ molecules) of a van der Waals fluid has the Helmholtz free energy

```{math}
A(T,V)=-nRT\left[\ln\left(\frac{V-nb_\text{m}}{N\Lambda^3}\right)+1\right]-\frac{a_\text{m}n^2}{V}+n\,\phi(T),
```

where $\Lambda$ is the thermal wavelength, which depends only on $T$, and $\phi(T)$ collects the contributions of the rotations and vibrations of the molecules, which also depend only on $T$. For carbon dioxide, $a_\text{m}=3.6551\ \mathrm{L^2\,bar\,mol^{-2}}$ and $b_\text{m}=0.042816\ \mathrm{L\,mol^{-1}}$ (the values used in Section 1.4). Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$ and $1\ \mathrm{L\,bar}=100\ \mathrm{J}$.

**(a)** Use $P=-(\partial A/\partial V)_T$ to find the equation of state.

**(b)** The fluid expands isothermally from $V_1$ to $V_2$. Find expressions for $\Delta A$, for $\Delta S$ from $S=-(\partial A/\partial T)_V$, and for $\Delta U$.

**(c)** Evaluate $\Delta A$, $\Delta S$, and $\Delta U$ for one mole of carbon dioxide that expands from $0.300\ \mathrm{L}$ to $1.00\ \mathrm{L}$ at $320\ \mathrm{K}$.

**(d)** At constant temperature, $\Delta A$ is the reversible work ([Section 5.1](section-01.md)). Compare the reversible work of part (c) with that of an ideal gas between the same volumes.

````{admonition} Solution
:class: dropdown

1. **Equation of state**

   Only the first two terms of $A$ depend on $V$:

   ```{math}
   P=-\left(\frac{\partial A}{\partial V}\right)_T=\frac{nRT}{V-nb_\text{m}}-\frac{a_\text{m}n^2}{V^2}.
   ```

   This is the van der Waals equation of state of [Section 1.4](../chapter-01/section-04.md).

2. **Change in $A$ at constant temperature**

   The terms $\ln(N\Lambda^3)$ and $n\,\phi(T)$ do not change, so

   ```{math}
   \Delta A=-nRT\ln\left(\frac{V_2-nb_\text{m}}{V_1-nb_\text{m}}\right)-a_\text{m}n^2\left(\frac{1}{V_2}-\frac{1}{V_1}\right).
   ```

3. **Change in $S$ and in $U$**

   The entropy is $S=-(\partial A/\partial T)_V$. The term $-a_\text{m}n^2/V$ does not depend on $T$. Of the remaining terms, only $-nRT\ln(V-nb_\text{m})$ depends on $V$, and it contributes

   ```{math}
   -\frac{\partial}{\partial T}\Big[-nRT\ln\left(V-nb_\text{m}\right)\Big]=nR\ln\left(V-nb_\text{m}\right).
   ```

   to the entropy. The entropy change of the isothermal expansion is therefore

   ```{math}
   \Delta S=nR\ln\left(\frac{V_2-nb_\text{m}}{V_1-nb_\text{m}}\right).
   ```

   From $A=U-TS$ at constant temperature, $\Delta U=\Delta A+T\,\Delta S$, and the logarithmic terms cancel:

   ```{math}
   \Delta U=-a_\text{m}n^2\left(\frac{1}{V_2}-\frac{1}{V_1}\right).
   ```

4. **Numbers for carbon dioxide**

   With $n=1\ \mathrm{mol}$, $nRT=(8.3145\ \mathrm{J\,K^{-1}})(320\ \mathrm{K})=2660.6\ \mathrm{J}$, and

   ```{math}
   \ln\left(\frac{1.00-0.042816}{0.300-0.042816}\right)=\ln(3.7218)=1.3142,
   ```

   ```{math}
   \begin{aligned}
   -a_\text{m}n^2\left(\frac{1}{V_2}-\frac{1}{V_1}\right)&=(3.6551\ \mathrm{L^2\,bar})\left(\frac{1}{0.300\ \mathrm{L}}-\frac{1}{1.00\ \mathrm{L}}\right)\\
   &=8.529\ \mathrm{L\,bar}=853\ \mathrm{J}.
   \end{aligned}
   ```

   Then

   ```{math}
   \begin{aligned}
   \Delta S&=(8.3145\ \mathrm{J\,K^{-1}})(1.3142)=10.93\ \mathrm{J\,K^{-1}},\\
   \Delta A&=-(2660.6\ \mathrm{J})(1.3142)+853\ \mathrm{J}=-3497+853=-2644\ \mathrm{J},\\
   \Delta U&=+853\ \mathrm{J}.
   \end{aligned}
   ```

5. **Comparison with an ideal gas**

   For an ideal gas, $a_\text{m}=b_\text{m}=0$, and

   ```{math}
   w_{\mathrm{rev}}=\Delta A=-nRT\ln\left(\frac{V_2}{V_1}\right)=-(2660.6\ \mathrm{J})\ln\left(\frac{1.00}{0.300}\right)=-3203\ \mathrm{J}.
   ```

   The van der Waals fluid has $w_{\mathrm{rev}}=\Delta A=-2644\ \mathrm{J}$. The term with $b_\text{m}$ changes the logarithmic part from $-3203\ \mathrm{J}$ to $-3497\ \mathrm{J}$, and the term with $a_\text{m}$ adds $+853\ \mathrm{J}$.

**Result.** $A$ returns the van der Waals equation of state. For the expansion of carbon dioxide, $\Delta S=10.93\ \mathrm{J\,K^{-1}}$, $\Delta U=+853\ \mathrm{J}$, and $\Delta A=w_{\mathrm{rev}}=-2644\ \mathrm{J}$, against $-3203\ \mathrm{J}$ for an ideal gas.
````

### Problem 2. Standard entropy of oxygen, mode by mode

Draws on [Section 5.1](section-01.md), [Section 4.3](../chapter-04/section-03.md), and Sections [2.6](../chapter-02/section-06.md) to [2.8](../chapter-02/section-08.md).

The tables of Sections 2.6 to 2.8 give, for $\mathrm{O_2}$, $\Theta_{\mathrm{rot}}=2.080\ \mathrm{K}$ with $\sigma=2$, $\Theta_{\mathrm{vib}}=2274\ \mathrm{K}$, and a ground electronic level of degeneracy $g_1=3$ whose first excited level is too high to be populated at room temperature. The molar mass is $31.998\ \mathrm{g\,mol^{-1}}$. Treat $\mathrm{O_2}$ as an ideal gas at $T=298.15\ \mathrm{K}$ and $P=1\ \mathrm{bar}$.

Use $h=6.62607\times10^{-34}\ \mathrm{J\,s}$, $k_{\mathrm B}=1.380649\times10^{-23}\ \mathrm{J\,K^{-1}}$, $N_A=6.02214\times10^{23}\ \mathrm{mol^{-1}}$, and $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find the translational entropy per mole from the Sackur–Tetrode equation.

**(b)** For each of the other modes, the molar entropy is $S_m=U_m/T+R\ln q$, with $U_m$ and $q$ of that mode. Find the rotational entropy with $q_{\mathrm{rot}}=T/(\sigma\Theta_{\mathrm{rot}})$ and $U_m=RT$.

**(c)** Find the vibrational entropy and the electronic entropy.

**(d)** Add the four contributions, and compare the total with $S^\circ(298.15\ \mathrm{K})=205.147\ \mathrm{J\,mol^{-1}\,K^{-1}}$ from the [NIST-JANAF tables](https://janaf.nist.gov/).

````{admonition} Solution
:class: dropdown

1. **Thermal wavelength and volume per molecule**

   The mass of one molecule is $m=(0.031998\ \mathrm{kg\,mol^{-1}})/N_A=5.3134\times10^{-26}\ \mathrm{kg}$.

   ```{math}
   \begin{aligned}
   \Lambda&=\frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}=1.7874\times10^{-11}\ \mathrm{m},\\
   \frac{V}{N}&=\frac{k_{\mathrm B}T}{P}=4.1164\times10^{-26}\ \mathrm{m^3}.
   \end{aligned}
   ```

2. **Translational entropy**

   ```{math}
   \frac{V}{N\Lambda^3}=\frac{4.1164\times10^{-26}\ \mathrm{m^3}}{(1.7874\times10^{-11}\ \mathrm{m})^3}=7.209\times10^{6},
   ```

   ```{math}
   \begin{aligned}
   S_{\mathrm{trans}}&=R\left[\ln\left(\frac{V}{N\Lambda^3}\right)+\frac{5}{2}\right]\\
   &=(8.3145)(15.791+2.500)=152.08\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   \end{aligned}
   ```

   The factor $1/N!$ of the partition function of the gas is counted here, so the other modes contribute $U_m/T+R\ln q$ each.

3. **Rotational entropy**

   ```{math}
   q_{\mathrm{rot}}=\frac{298.15}{(2)(2.080)}=71.67,
   ```

   ```{math}
   \begin{aligned}
   S_{\mathrm{rot}}&=\frac{U_m}{T}+R\ln q_{\mathrm{rot}}=R\,(1+\ln 71.67)\\
   &=(8.3145)(5.2721)=43.83\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   \end{aligned}
   ```

4. **Vibrational entropy**

   From [Section 2.6](../chapter-02/section-06.md), with $x=\Theta_{\mathrm{vib}}/T$,

   ```{math}
   \frac{U_m}{RT}=\frac{x}{2}+\frac{x}{e^{x}-1},
   \qquad
   \ln q_{\mathrm{vib}}=-\frac{x}{2}-\ln\left(1-e^{-x}\right).
   ```

   The zero-point terms cancel in the entropy:

   ```{math}
   S_{\mathrm{vib}}=R\left[\frac{x}{e^{x}-1}-\ln\left(1-e^{-x}\right)\right].
   ```

   Here $x=2274/298.15=7.627$ and $e^{x}=2053$, so

   ```{math}
   S_{\mathrm{vib}}=(8.3145)(0.003717+0.000487)=0.035\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   ```

5. **Electronic entropy**

   With $q_{\mathrm{elec}}=g_1=3$ and no electronic energy above the ground level,

   ```{math}
   S_{\mathrm{elec}}=R\ln 3=9.134\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   ```

6. **Total**

   ```{math}
   S^\circ=152.08+43.83+0.035+9.134=205.08\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   ```

   The tabulated value is $205.147\ \mathrm{J\,mol^{-1}\,K^{-1}}$, which is $0.07\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (0.03%) higher.

**Result.** $S^\circ(298.15\ \mathrm{K})=205.08\ \mathrm{J\,mol^{-1}\,K^{-1}}$, against the tabulated $205.147\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

| Mode | $S_m$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) |
|---|---:|
| Translation | 152.08 |
| Rotation | 43.83 |
| Vibration | 0.035 |
| Electronic | 9.134 |
| Total | 205.08 |
````

### Problem 3. Residual entropy from counting orientations

Draws on [Section 5.2](section-02.md) and [Section 4.3](../chapter-04/section-03.md).

Suppose that the molecules of a crystal freeze into their orientations independently of one another.

- In solid nitrous oxide, the linear molecule $\mathrm{NNO}$ can lie in either of two directions ($\mathrm{NNO}$ or $\mathrm{ONN}$).
- In solid perchloryl fluoride, $\mathrm{FClO_3}$, the chlorine atom is surrounded by one fluorine atom and three oxygen atoms at the corners of a tetrahedron, and the fluorine atom can occupy any of the four corners.

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find $\Omega_0$ for $N$ molecules and the molar residual entropy of each crystal if all orientations are equally likely.

**(b)** Suppose instead that the two orientations of a linear molecule have probabilities $p$ and $1-p$. Use the Gibbs entropy to find the molar residual entropy as a function of $p$, and evaluate it for $p=0.50$ and $p=0.75$.

**(c)** [Section 5.2](section-02.md) quotes a measured residual entropy of about $4.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for $\mathrm{CO}$. Find, to two figures, the value of $p$ at which the model of part (b) gives this entropy.

````{admonition} Solution
:class: dropdown

1. **Equally likely orientations**

   Each molecule has $\omega$ orientations, and the molecules are independent, so $\Omega_0=\omega^N$ and

   ```{math}
   S_{\mathrm{res}}=k_{\mathrm B}\ln\Omega_0=Nk_{\mathrm B}\ln\omega,
   \qquad
   S_{\mathrm{res,m}}=R\ln\omega.
   ```

   For $\mathrm{N_2O}$, $\omega=2$ and $\Omega_0=2^N$. For $\mathrm{FClO_3}$, $\omega=4$ and $\Omega_0=4^N$.

   ```{math}
   \begin{aligned}
   S_{\mathrm{res,m}}(\mathrm{N_2O})&=R\ln 2=5.763\ \mathrm{J\,mol^{-1}\,K^{-1}},\\
   S_{\mathrm{res,m}}(\mathrm{FClO_3})&=R\ln 4=11.526\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   \end{aligned}
   ```

2. **Unequal probabilities**

   For one molecule with two orientations, the Gibbs entropy is $-k_{\mathrm B}\left[p\ln p+(1-p)\ln(1-p)\right]$. The entropies of independent molecules add ([Section 4.3](../chapter-04/section-03.md)), so for one mole

   ```{math}
   S_{\mathrm{res,m}}=-R\left[p\ln p+(1-p)\ln(1-p)\right].
   ```

   For $p=0.50$ the bracket is $-\ln 2$, and the result of step 1 returns: $S_{\mathrm{res,m}}=R\ln 2=5.763\ \mathrm{J\,mol^{-1}\,K^{-1}}$. For $p=0.75$,

   ```{math}
   \begin{aligned}
   S_{\mathrm{res,m}}&=-R\left[0.75\ln(0.75)+0.25\ln(0.25)\right]\\
   &=(8.3145)(0.2158+0.3466)=4.68\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   \end{aligned}
   ```

3. **The value of $p$ for $\mathrm{CO}$**

   The equation $-R\left[p\ln p+(1-p)\ln(1-p)\right]=4.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$ has no solution in closed form. The entropy falls as $p$ rises above $0.50$, and step 2 shows that $p$ is slightly above $0.75$:

   | $p$ | $S_{\mathrm{res,m}}$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) |
   |---:|---:|
   | 0.75 | 4.68 |
   | 0.76 | 4.58 |
   | 0.77 | 4.48 |

   The measured value lies between the first two rows, closer to the second, so $p\approx0.76$.

**Result.** $S_{\mathrm{res,m}}=R\ln 2=5.763\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for two orientations and $R\ln 4=11.526\ \mathrm{J\,mol^{-1}\,K^{-1}}$ for four. With unequal probabilities, $S_{\mathrm{res,m}}=-R\left[p\ln p+(1-p)\ln(1-p)\right]$, which is $4.68\ \mathrm{J\,mol^{-1}\,K^{-1}}$ at $p=0.75$ and matches the measured value for $\mathrm{CO}$ at $p\approx0.76$. This model keeps the molecules independent; Section 5.2 describes the same shortfall through correlations between neighboring molecules.
````

### Problem 4. When does steam reforming become spontaneous?

Draws on [Section 5.3](section-03.md).

Most industrial hydrogen is made by steam reforming of methane:

```{math}
\mathrm{CH_4}(g)+\mathrm{H_2O}(g)\;\longrightarrow\;\mathrm{CO}(g)+3\,\mathrm{H_2}(g).
```

The table lists data at $298.15\ \mathrm{K}$ from the [NIST-JANAF tables](https://janaf.nist.gov/).

| Species | $\Delta_f H^\circ$ ($\mathrm{kJ\,mol^{-1}}$) | $S^\circ$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) |
|---|---:|---:|
| $\mathrm{CH_4}(g)$ | −74.873 | 186.251 |
| $\mathrm{H_2O}(g)$ | −241.826 | 188.834 |
| $\mathrm{CO}(g)$ | −110.527 | 197.653 |
| $\mathrm{H_2}(g)$ | 0 | 130.680 |

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find $\Delta_r H^\circ$, $\Delta_r S^\circ$, and $\Delta_r G^\circ$ at $298.15\ \mathrm{K}$.

**(b)** Treat $\Delta_r H^\circ$ and $\Delta_r S^\circ$ as independent of temperature, and find the temperature at which $\Delta_r G^\circ=0$.

**(c)** Find $\Delta_r G^\circ$ at $1200\ \mathrm{K}$ with the constant heat capacities of [Section 5.3](section-03.md): $C_P=\tfrac{7}{2}R$ for the linear molecules and $C_P=4R$ for the nonlinear molecules. Compare the result with $-77.92\ \mathrm{kJ\,mol^{-1}}$ from the NIST-JANAF tables.

**(d)** Reformers operate at about $25\ \mathrm{bar}$. Find $\Delta_r G$ at $1200\ \mathrm{K}$ with every species at $25\ \mathrm{bar}$, starting from your result of part (c).

````{admonition} Solution
:class: dropdown

1. **Reaction quantities at $298.15\ \mathrm{K}$**

   ```{math}
   \begin{aligned}
   \Delta_r H^\circ&=\left[-110.527+3(0)\right]-\left[-74.873-241.826\right]\\
   &=+206.172\ \mathrm{kJ\,mol^{-1}},\\
   \Delta_r S^\circ&=\left[197.653+3(130.680)\right]-\left[186.251+188.834\right]\\
   &=+214.608\ \mathrm{J\,mol^{-1}\,K^{-1}},\\
   \Delta_r G^\circ&=206.172-(298.15)(0.214608)=+142.19\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

2. **Temperature at which $\Delta_r G^\circ=0$**

   ```{math}
   T=\frac{\Delta_r H^\circ}{\Delta_r S^\circ}=\frac{206\,172\ \mathrm{J\,mol^{-1}}}{214.608\ \mathrm{J\,mol^{-1}\,K^{-1}}}=961\ \mathrm{K}.
   ```

   With constant $\Delta_r H^\circ$ and $\Delta_r S^\circ$, $\Delta_r G^\circ$ is positive below this temperature and negative above it.

3. **Heat capacity of reaction**

   $\mathrm{CO}$ and $\mathrm{H_2}$ are linear; $\mathrm{CH_4}$ and $\mathrm{H_2O}$ are nonlinear.

   ```{math}
   \Delta_r C_P=\left[\tfrac{7}{2}R+3\left(\tfrac{7}{2}R\right)\right]-\left[4R+4R\right]=6R=49.887\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   ```

4. **Corrections to $1200\ \mathrm{K}$**

   With $T_i=298.15\ \mathrm{K}$ and $T_f=1200\ \mathrm{K}$,

   ```{math}
   \begin{aligned}
   \Delta_r H^\circ(T_f)&=206.172+(0.049887)(1200-298.15)\\
   &=206.172+44.991=251.16\ \mathrm{kJ\,mol^{-1}},\\
   \Delta_r S^\circ(T_f)&=214.608+(49.887)\ln\left(\frac{1200}{298.15}\right)\\
   &=214.608+69.467=284.07\ \mathrm{J\,mol^{-1}\,K^{-1}},\\
   \Delta_r G^\circ(T_f)&=251.16-(1200)(0.28407)\\
   &=251.16-340.89=-89.73\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

   The NIST-JANAF value is $-77.92\ \mathrm{kJ\,mol^{-1}}$. The estimate has the right sign and is $11.8\ \mathrm{kJ\,mol^{-1}}$ too negative. The model leaves out the vibrations: in the same tables the heat capacity of methane at $1200\ \mathrm{K}$ is $78.8\ \mathrm{J\,mol^{-1}\,K^{-1}}$, more than twice $4R$.

5. **Pressure**

   The reaction turns two moles of gas into four, so $\Delta\nu=(1+3)-(1+1)=+2$. With every species at $P=25\ \mathrm{bar}$,

   ```{math}
   \begin{aligned}
   \Delta_r G&=\Delta_r G^\circ+\Delta\nu\,RT\ln\left(\frac{P}{P^\circ}\right)\\
   &=-89.73+(2)(0.0083145)(1200)\ln(25)\\
   &=-89.73+64.23=-25.5\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

**Result.** At $298.15\ \mathrm{K}$, $\Delta_r H^\circ=+206.17\ \mathrm{kJ\,mol^{-1}}$, $\Delta_r S^\circ=+214.61\ \mathrm{J\,mol^{-1}\,K^{-1}}$, and $\Delta_r G^\circ=+142.19\ \mathrm{kJ\,mol^{-1}}$. $\Delta_r G^\circ$ changes sign near $961\ \mathrm{K}$. At $1200\ \mathrm{K}$ the constant-$C_P$ model gives $\Delta_r G^\circ=-89.7\ \mathrm{kJ\,mol^{-1}}$ (tabulated: $-77.9\ \mathrm{kJ\,mol^{-1}}$), and $\Delta_r G=-25.5\ \mathrm{kJ\,mol^{-1}}$ with every species at $25\ \mathrm{bar}$.
````

## 4. Further Reading

**Textbooks**

- @mcquarrie1997, Chapter 21 (entropy and the third law of thermodynamics) and Chapter 22, Sections 22.1 to 22.7 (Helmholtz and Gibbs energies).
- @dill2011, Chapter 8 (laboratory conditions and free energies) and Chapter 9 (Maxwell's relations and mixtures).
- @hill1986, Chapter 2 (the second and third laws in terms of ensembles).

**Original papers**

- @clayton1932 measure the entropy of carbon monoxide by calorimetry and compare it with the value calculated from spectroscopic data. The difference is the residual entropy of [Section 5.2](section-02.md).
- @pauling1935 counts the arrangements of the hydrogen atoms in ice, the model of Example 2 of [Section 5.2](section-02.md).
- @haber1920 gives his account of the synthesis of ammonia from nitrogen and hydrogen ([Section 5.3](section-03.md)).
