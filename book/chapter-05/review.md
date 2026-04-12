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

### Section 5.1: Free Energy

1. **Fundamental Inequality**
   - Combining the First Law with the Clausius inequality $dS \ge \delta q/T$ and the work bound $\delta w \ge -P\,dV$ gives

   ```{math}
   dU \le T\,dS - P\,dV,
   ```

   with equality for reversible changes. This single inequality encodes the direction of spontaneous evolution for a closed system with $PV$-only work.

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

   - Defining $A$ and $G$ extends the strategy used in Section 3.3 to define $H$.

4. **Conjugate Variables as Derivatives**
   - At equilibrium, intensive variables are read off as slopes of a potential along its natural-variable axes. For example:

   ```{math}
   T = \left(\frac{\partial U}{\partial S}\right)_V,\qquad
   P = -\left(\frac{\partial U}{\partial V}\right)_S,
   ```

   with analogous identities for $H$, $A$, and $G$.

5. **Euler Relation** (fixed $N$)
   - For a simple system with fixed composition, $U(S,V)$ is homogeneous of degree 1 in its extensive arguments. Euler's theorem gives

   ```{math}
   U = TS - PV.
   ```

   - Allowing $N$ to vary adds a $\mu N$ term; this will be developed when the chemical potential is introduced.

6. **Helmholtz Free Energy as a Master Potential**
   - The bridge to statistical mechanics from Section 4.3 is

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
   dG\big|_{T,P} = \delta w_{\text{non-}PV,\text{rev}}.
   ```

   - $\Delta G$ measures the reversible non-$PV$ work available from a process — the "useful" work beyond unavoidable $PV$ work against the atmosphere. For irreversible processes, the useful work output is strictly less than $|\Delta G|$.

---

### Section 5.2: Third Law

1. **Planck's Statement of the Third Law**
   - As $T \to 0\,\mathrm{K}$, the entropy of any pure crystalline substance approaches a constant. For an **ideal crystal** (no defects, impurities, or disorder), that constant is zero:

   ```{math}
   S(0\,\mathrm{K}) = 0 \quad\text{(ideal crystal)}.
   ```

   - This fixes the otherwise-undetermined integration constant in $\Delta S = \int \delta q_{\mathrm{rev}}/T$.

2. **Statistical-Mechanical Justification**
   - From Boltzmann's formula $S = k_{\mathrm B}\ln\Omega$ (Section 4.3), a unique ground state has $\Omega_0 = 1$, so $S(0) = k_{\mathrm B}\ln 1 = 0$. The third law is the statement that a system with a nondegenerate ground state has zero entropy at absolute zero.

3. **Residual Entropy**
   - If the ground state is degenerate or the system has frozen-in disorder ($\Omega_0 > 1$), entropy approaches a nonzero constant:

   ```{math}
   S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0.
   ```

   - This does **not** violate the third law: Planck's statement says $S \to$ a constant, not necessarily zero.
   - Examples:
     - **Solid CO** (two nearly isoenergetic orientations per molecule): $\Omega_0 = 2^N$, giving $S_{\mathrm{res,m}} = R\ln 2 \approx 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (experiment: $\approx 4.6$, indicating partial ordering).
     - **Ice** (Pauling's model with ice rules): $\Omega_0 \approx (3/2)^N$, giving $S_{\mathrm{res,m}} = R\ln(3/2) \approx 3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}$ (experiment: $\approx 3.41$).

4. **Heat Capacities Must Vanish as $T \to 0$**
   - For the integral $S(T) = S(0) + \int_0^T (C_P/T')\,dT'$ to converge, $C_P$ (and similarly $C_V$) must go to zero at $T = 0$.
   - Classical equipartition gives a temperature-independent $C_V$ and therefore **violates** the third law. Quantum statistical mechanics resolves the contradiction: the Einstein and Debye heat capacities (Section 2.6) vanish as $T \to 0$ because vibrational modes freeze out.

5. **Unattainability of Absolute Zero** (Nernst form)
   - No finite sequence of thermodynamic operations can cool a system to exactly $T = 0$. Each successive cooling step removes less entropy than the previous one; the process converges but never reaches zero.

6. **Absolute (Third-Law) Entropies**
   - With $S(0) = 0$ fixed for ideal crystals, absolute entropies are computed by integrating $C_P/T$ from $0$ to $T$ and adding phase-transition contributions:

   ```{math}
   S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}}.
   ```

   - Tabulated $S^\circ(298.15\ \mathrm{K})$ values (NIST WebBook, JANAF) are the basis for computing reaction entropies and Gibbs energies.

---

### Section 5.3: Ammonia Formation

1. **Thermodynamics vs. Kinetics**
   - $\Delta G$ determines which direction is favored at equilibrium; activation barriers determine *how fast* equilibrium is reached. These are independent questions.
   - For $\mathrm{N_2 + 3\,H_2 \to 2\,NH_3}$ at 298.15 K: $\Delta H_r^\circ \approx -92\ \mathrm{kJ\,mol^{-1}}$, $\Delta S_r^\circ \approx -198\ \mathrm{J\,K^{-1}\,mol^{-1}}$, $\Delta G_r^\circ \approx -33\ \mathrm{kJ\,mol^{-1}}$. Thermodynamically spontaneous but kinetically inert because of the $\mathrm{N\equiv N}$ triple bond.
   - A **catalyst** lowers the activation barrier without changing $\Delta G_r^\circ$ or the equilibrium constant — it accelerates forward and reverse reactions equally.

2. **Temperature Effect for $\Delta S_r^\circ < 0$**
   - From $\Delta G_r^\circ = \Delta H_r^\circ - T\,\Delta S_r^\circ$, a negative $\Delta S_r^\circ$ means $-T\,\Delta S_r^\circ > 0$ grows with $T$, driving $\Delta G_r^\circ$ upward.
   - At high enough $T$, $\Delta G_r^\circ$ changes sign and the reaction becomes non-spontaneous at 1 bar.

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
   - For the ammonia reaction: $\Delta C_{P,r} = 2(4R) - (7R/2) - 3(7R/2) = -6R$.

5. **Result at 500 °C**
   - Applying constant-$C_P$ corrections from 298.15 K to 773.15 K gives $\Delta H_r^\circ(773) \approx -115.6\ \mathrm{kJ\,mol^{-1}}$, $\Delta S_r^\circ(773) \approx -245.4\ \mathrm{J\,K^{-1}\,mol^{-1}}$, and

   ```{math}
   \Delta G_r^\circ(773) \approx +74.1\ \mathrm{kJ\,mol^{-1}}.
   ```

   - The reaction is non-spontaneous at 1 bar — the central dilemma of Haber–Bosch.

6. **Pressure Dependence of $G$ (Ideal Gas)**
   - From $dG\big|_T = V\,dP$ with $V = RT/P$:

   ```{math}
   G(T, P) - G^\circ(T) = RT\ln\!\left(\frac{P}{P^\circ}\right).
   ```

   - For a reaction with stoichiometric change $\Delta\nu$ in moles of gas:

   ```{math}
   \Delta G_r(T, P) \approx \Delta G_r^\circ(T) + \Delta\nu\,RT\ln\!\left(\frac{P}{P^\circ}\right).
   ```

   - For $\mathrm{N_2 + 3\,H_2 \to 2\,NH_3}$: $\Delta\nu = 2 - 4 = -2$, so increasing $P$ lowers $\Delta G_r$ (Le Châtelier: high pressure favors the side with fewer gas molecules).

7. **Restoring Spontaneity with Pressure**
   - Setting $\Delta G_r = 0$ and solving for $P$:

   ```{math}
   P = P^\circ\exp\!\left(-\frac{\Delta G_r^\circ(T)}{\Delta\nu\,RT}\right).
   ```

   - At 500 °C: $P \approx 320\ \mathrm{bar}$. At 400 °C: $P \approx 87\ \mathrm{bar}$. The Haber–Bosch operating window (400–500 °C, 150–350 bar) is the engineering compromise predicted by these estimates.

8. **Engineering Knobs Summary**
   - **Catalyst**: speeds up both directions; does *not* change $\Delta G_r^\circ$ or $K$.
   - **Temperature**: raises rate (Arrhenius) but pushes $\Delta G_r^\circ$ upward when $\Delta S_r^\circ < 0$.
   - **Pressure**: compensates for the temperature penalty when $\Delta\nu < 0$.

---

## 2. Checklist of Most Important Equations

Below is a unified list of the major equations from Sections 5.1–5.3.

A. **Fundamental Inequality**

```{math}
dU \le T\,dS - P\,dV.
```

- **Applicability**: closed system with $PV$-only work. Equality for reversible processes; strict inequality for spontaneous irreversible processes.

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
dG\big|_{T,P} = \delta w_{\text{non-}PV,\text{rev}}.
```

- **Applicability**: reversible process at constant $T$ and $P$. $\Delta G$ equals the reversible non-$PV$ work; for irreversible processes, $|\Delta G|$ is an upper bound on useful work output.

---

H. **Third Law**

```{math}
S(0\,\mathrm{K}) = 0 \quad\text{(ideal crystal)},\qquad
S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0 \quad\text{(degenerate ground state)}.
```

- Fixes the entropy reference point; residual entropy reflects ground-state multiplicity, not a failure of the third law.

---

I. **Absolute (Third-Law) Entropy**

```{math}
S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}}.
```

- **Applicability**: any substance whose heat capacity and phase-transition enthalpies are known from 0 to $T$; supplies tabulated $S^\circ(298.15\ \mathrm{K})$ values.

---

J. **Temperature Corrections to Reaction Thermodynamics** (constant $C_P$)

```{math}
\Delta H_r^\circ(T_f) = \Delta H_r^\circ(T_i) + \Delta C_{P,r}(T_f - T_i),
```

```{math}
\Delta S_r^\circ(T_f) = \Delta S_r^\circ(T_i) + \Delta C_{P,r}\ln\!\left(\frac{T_f}{T_i}\right).
```

- **Applicability**: reaction data from a reference temperature $T_i$ to a target temperature $T_f$, assuming temperature-independent heat capacities over the interval.

---

K. **Ideal-Gas Pressure Dependence of $G$**

```{math}
G(T,P) - G^\circ(T) = RT\ln\!\left(\frac{P}{P^\circ}\right).
```

- **Applicability**: pure ideal gas at constant $T$. For a reaction with stoichiometric gas change $\Delta\nu$:

```{math}
\Delta G_r(T,P) \approx \Delta G_r^\circ(T) + \Delta\nu\,RT\ln\!\left(\frac{P}{P^\circ}\right).
```

- For $\Delta\nu < 0$, higher pressure lowers $\Delta G_r$ and favors products.
