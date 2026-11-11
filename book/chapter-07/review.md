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

### Section 7.1: Equilibrium Constant

1. **Extent of Reaction**
   - For a reaction $\nu_A A + \nu_B B \rightleftharpoons \nu_Y Y + \nu_Z Z$, composition changes are tracked by a single scalar $\xi$ (units of moles) via

   ```{math}
   dn_i = \nu_i\,d\xi,
   ```

   with the **signed stoichiometry convention**: $\nu_i < 0$ for reactants and $\nu_i > 0$ for products.
   - $\xi$ collapses the multiple $dn_i$ into one progress variable, enabling reaction thermodynamics to be written as a derivative in a single direction.

2. **Gibbs Free Energy of Reaction**
   - At fixed $T,P$, the Gibbs differential reduces to $dG = \sum_i \mu_i\,dn_i$; substituting $dn_i = \nu_i\,d\xi$ gives

   ```{math}
   \left(\frac{\partial G}{\partial \xi}\right)_{T,P} = \sum_i \nu_i\,\mu_i \equiv \Delta_r G.
   ```

   - $\Delta_r G$ is the slope of $G$ along the reaction coordinate and determines the spontaneous direction at fixed $T,P$.

3. **Chemical Equilibrium Condition**
   - Equilibrium is the stationary point of $G$ at fixed $T,P$:

   ```{math}
   \Delta_r G = \sum_i \nu_i\,\mu_i = 0.
   ```

   - This is the multi-species analog of the phase-equilibrium condition $\mu_\alpha = \mu_\beta$ from § 6.1: instead of moving matter between two phases, matter is redistributed among several chemical species.
   - **Spontaneity**: $\Delta_r G < 0$ drives the reaction forward, $\Delta_r G > 0$ backward, $\Delta_r G = 0$ is equilibrium.

4. **Ideal-Gas Chemical Potential and the Dimensionless Argument of the Logarithm**
   - For a species in an ideal-gas mixture,

   ```{math}
   \mu_i(T,P_i) = \mu_i^\circ(T) + RT\ln\!\left(\frac{P_i}{P^\circ}\right),
   ```

   with $P^\circ = 1\ \mathrm{bar}$.
   - The **ratio** $P_i/P^\circ$ is what makes the logarithm unit-independent; real-gas and solution generalizations preserve this structure by replacing $P_i/P^\circ$ with an activity $a_i$.

5. **Reaction Quotient and Equilibrium Constant**
   - Substituting $\mu_i = \mu_i^\circ + RT\ln(P_i/P^\circ)$ into $\Delta_r G = \sum_i \nu_i\mu_i$ and grouping terms gives

   ```{math}
   \Delta_r G = \Delta_r G^\circ + RT\ln Q_p,
   \qquad
   Q_p \equiv \prod_i\!\left(\frac{P_i}{P^\circ}\right)^{\!\nu_i}.
   ```

   - At equilibrium, $Q_p = K$ and $\Delta_r G = 0$, giving the **bridge between thermochemistry and equilibrium**:

   ```{math}
   \Delta_r G^\circ = -RT\ln K,
   \qquad
   K = \exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right).
   ```

   - Comparing $Q_p$ to $K$ predicts direction: $Q_p < K$ drives products, $Q_p > K$ drives reactants.

6. **Why $G(\xi)$ Has an Interior Minimum**
   - For the NO$_2$ dimerization at $P = P^\circ$, writing $G(\xi)$ explicitly yields

   ```{math}
   G(\xi) - G(0) = \xi\,\Delta_r G^\circ + RT\sum_i n_i(\xi)\ln y_i(\xi).
   ```

   - The first term is linear in $\xi$ (the "bookkeeping" drop in standard chemical potentials); the second is the **ideal-mixing entropy**, which diverges toward $-\infty$ at pure reactants *and* pure products.
   - The tug-of-war between these two contributions is *the* reason equilibrium is generically an interior minimum rather than pure reactants or pure products, even when $|\Delta_r G^\circ| \gg RT$.

7. **Temperature Dependence: the van't Hoff Equation**
   - Starting from the Gibbs–Helmholtz relation and $\Delta_r G^\circ = -RT\ln K$:

   ```{math}
   \left(\frac{\partial \ln K}{\partial T}\right)_P = \frac{\Delta_r H^\circ}{RT^2}.
   ```

   - Assuming $\Delta_r H^\circ$ is approximately constant over $[T_1, T_2]$:

   ```{math}
   \ln\!\left[\frac{K(T_2)}{K(T_1)}\right] = -\frac{\Delta_r H^\circ}{R}\!\left(\frac{1}{T_2} - \frac{1}{T_1}\right).
   ```

   - Exothermic reactions ($\Delta_r H^\circ < 0$) have $K$ that *decreases* with $T$; endothermic reactions have $K$ that *increases* with $T$. This is the quantitative form of Le Châtelier's principle for temperature.
   - The mathematical structure is identical to the integrated Clausius–Clapeyron equation of § 6.2 — a consequence of the shared Gibbs–Helmholtz origin.

8. **Computing $K$ from Tabulated Thermochemistry**
   - The standard workflow for any gas-phase reaction at a specified $T$:
     1. Balance the reaction and fix a consistent standard state ($P^\circ = 1\ \mathrm{bar}$, usually at 298.15 K).
     2. Pull $\Delta_f H_i^\circ$ and $S_i^\circ$ (or $H_i^\circ(T)$ and $S_i^\circ(T)$ if $T \neq 298.15\ \mathrm{K}$) from NIST WebBook, NIST–JANAF, or ATcT.
     3. Form $\Delta_r H^\circ = \sum_i \nu_i H_i^\circ$ and $\Delta_r S^\circ = \sum_i \nu_i S_i^\circ$ using signed stoichiometry.
     4. Combine into $\Delta_r G^\circ(T) = \Delta_r H^\circ(T) - T\,\Delta_r S^\circ(T)$, and finally $K(T) = \exp(-\Delta_r G^\circ(T)/RT)$.
   - Common pitfalls: mixing phases (e.g., $\mathrm{H_2O(\ell)}$ vs. $\mathrm{H_2O(g)}$), mismatched standard states, and the near-universal unit mismatch between $\Delta H^\circ$ (kJ/mol) and $S^\circ$ (J/mol·K).

---

### Section 7.2: Equilibrium Constants from Microscopic Properties

1. **The Right Potential at Constant $T,V$**
   - At fixed $T$ and $V$, the Helmholtz free energy is the master potential (§ 5.1), with differential

   ```{math}
   dA = -S\,dT - P\,dV + \sum_i \mu_i\,dn_i.
   ```

   - Equilibrium is the stationary point of $A$, yielding

   ```{math}
   \sum_i \nu_i\,\mu_i = 0,
   ```

   the same condition as in § 7.1. What changes between fixed $T,P$ and fixed $T,V$ is *which potential is stationary*, not where equilibrium sits — $\mu_i$ is intensive, so the equilibrium composition is the same either way.

2. **Partition Function of an Ideal-Gas Mixture**
   - From factorization over independent subsystems (§ 2.4) and the indistinguishability correction (§ 2.5),

   ```{math}
   Q(T,V,\{N_i\}) = \prod_i \frac{q_i(T,V)^{N_i}}{N_i!},
   ```

   where $q_i$ is the single-molecule partition function for species $i$. The factorization is exact for an ideal gas because molecules of different species do not interact.

3. **Chemical Potential from $A = -k_{\mathrm B}T\ln Q$**
   - Using $A = -k_{\mathrm B}T\ln Q$ (§ 4.3) and differentiating with respect to $N_i$ under Stirling's approximation gives

   ```{math}
   \mu_i = -RT\ln\!\left(\frac{q_i}{N_i}\right) = RT\ln\!\left(\frac{N_i}{q_i}\right),
   ```

   on a per-mole basis. This is the microscopic counterpart of $\mu_i = \mu_i^\circ(T) + RT\ln(P_i/P^\circ)$ from § 7.1: the macroscopic expression writes $\mu_i$ in terms of a tabulated $\mu_i^\circ(T)$ and a partial pressure; the microscopic expression writes the same $\mu_i$ in terms of a counted ratio $q_i/N_i$.

4. **Equilibrium Constant in Terms of Partition Functions**
   - Inserting $\mu_i = -RT\ln(q_i/N_i)$ into $\sum_i \nu_i\mu_i = 0$ and exponentiating:

   ```{math}
   \prod_i N_i^{\nu_i} = \prod_i q_i^{\nu_i}.
   ```

   - Dividing by $V^{\Delta\nu}$ with $\Delta\nu \equiv \sum_i \nu_i$ gives the concentration-based equilibrium constant

   ```{math}
   K_c = \prod_i\!\left(\frac{q_i}{V}\right)^{\!\nu_i},
   ```

   and the pressure-based form follows from $P_i = c_i k_{\mathrm B}T$:

   ```{math}
   K_p = K_c\,\left(\frac{k_{\mathrm B}T}{P^\circ}\right)^{\!\Delta\nu}.
   ```

   - When $\Delta\nu = 0$, the prefactor is unity and $K_p = K_c$; the uncancelled factors of $V$ (and hence pressure) reappear only when $\Delta\nu \neq 0$.

5. **Diatomic Single-Molecule Partition Function**
   - In the rigid-rotor/harmonic-oscillator approximation (§ 2.8), $q = q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}$:

   ```{math}
   q_{\text{diatomic}}
   = \frac{V}{\Lambda^3}\,
     \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}\,
     \frac{e^{-\Theta_{\mathrm{vib}}/(2T)}}{1-e^{-\Theta_{\mathrm{vib}}/T}}\,
     g_1\,e^{D_0/(RT)},
   \qquad
   \Lambda \equiv \sqrt{\frac{h^2}{2\pi m k_{\mathrm B}T}}.
   ```

   - **Energy-zero convention**: the zero of molecular energy is the dissociated atoms at rest, so the bound $v=0$ level sits at $-D_0$ and the bond energy enters the electronic factor as $e^{+D_0/RT}$. An alternative convention absorbs the zero-point energy into the electronic factor and drops the $e^{-\Theta_{\mathrm{vib}}/(2T)}$ prefactor; both conventions yield the same $K$ when used consistently.

6. **$\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ as a $\Delta\nu = 0$ Benchmark**
   - With $\Delta\nu = 0$, volume factors cancel exactly and

   ```{math}
   K = \frac{q_{HI}^2}{q_{H_2}\,q_{I_2}}
   ```

   depends only on intrinsic molecular properties. Expanding with Eq. (F) below separates $K$ into four physically meaningful factors:

   ```{math}
   K =
   \underbrace{\left(\frac{m_{HI}^2}{m_{H_2}\,m_{I_2}}\right)^{3/2}}_{\text{translational}}
   \;
   \underbrace{\frac{\sigma_{H_2}\sigma_{I_2}}{\sigma_{HI}^2}\,
   \frac{\Theta_{\mathrm{rot}}^{H_2}\,\Theta_{\mathrm{rot}}^{I_2}}{\left(\Theta_{\mathrm{rot}}^{HI}\right)^2}}_{\text{rotational}}
   \;
   \underbrace{\frac{(1-e^{-\Theta_{\mathrm{vib}}^{H_2}/T})(1-e^{-\Theta_{\mathrm{vib}}^{I_2}/T})}{(1-e^{-\Theta_{\mathrm{vib}}^{HI}/T})^2}}_{\text{vibrational}}
   \;
   \underbrace{\exp\!\left(\frac{2D_0^{HI}-D_0^{H_2}-D_0^{I_2}}{RT}\right)}_{\text{bond-energy}}.
   ```

   - Zero-point $e^{-\Theta_{\mathrm{vib}}/(2T)}$ factors cancel in pairs under the chosen energy zero, leaving only the "hot-mode" $(1-e^{-\Theta_{\mathrm{vib}}/T})$ factors in the vibrational block.

7. **Physical Reading of the Four Factors at 700 K**
   - Plugging spectroscopic constants for H$_2$, I$_2$, and HI gives $K(700\ \mathrm{K}) \approx 62$, to be compared with the experimental Bodenstein value $\approx 54$ — about 15 % agreement, consistent with the rigid-rotor/harmonic-oscillator approximation and uncertainties in $\Theta_{\mathrm{rot}}$, $\Theta_{\mathrm{vib}}$, and $D_0$.
   - **Translational** ($\approx 181$): dominated by the mass mismatch between H and I, since $m_{HI}^2/(m_{H_2}m_{I_2}) = (m_H+m_I)^2/(4m_H m_I)$ is driven far from unity when $m_I \gg m_H$.
   - **Rotational** ($\approx 0.22$): I$_2$'s tiny $\Theta_{\mathrm{rot}} \approx 0.05\ \mathrm{K}$ makes $q_{\mathrm{rot}}(\mathrm{I_2})$ extremely large, favoring reactants; the symmetry-number prefactor $\sigma_{H_2}\sigma_{I_2}/\sigma_{HI}^2 = 4$ nudges toward products but does not overcome this.
   - **Vibrational** ($\approx 0.36$): HI's stiff H–I stretch remains nearly frozen at 700 K while I$_2$'s low-frequency mode is partially excited, combining unfavorably for products.
   - **Bond-energy** ($\approx 4.3$): $2D_0^{HI} - D_0^{H_2} - D_0^{I_2} \approx +8.5\ \mathrm{kJ/mol}$ — the two HI bonds together are slightly stronger than one H–H plus one I–I bond.
   - Every digit of $K$ traces back to a specific molecular property; the relative sizes of the four factors identify *which* molecular feature most strongly determines the equilibrium.

8. **Consistency with the Macroscopic Route**
   - A student who computes $\Delta_r G^\circ$ from thermochemical tables (§ 7.1.6) and a student who evaluates $K$ from spectroscopic constants (§ 7.2.8) arrive at the same number by very different paths. The macroscopic route answers "what is $K$?" efficiently; the microscopic route answers "*why* is $K$ what it is?" by decomposing it into separable molecular contributions. Together they close the statistical-mechanical bridge built across Chapters 2, 4, 5, and 7.

---

## 2. Checklist of Most Important Equations

Below is a unified list of the major equations from Sections 7.1–7.2.

A. **Extent of Reaction**

```{math}
dn_i = \nu_i\,d\xi,
\qquad
\Delta\nu \equiv \sum_i \nu_i.
```

- **Applicability**: any reaction written with signed stoichiometric coefficients ($\nu_i < 0$ for reactants, $\nu_i > 0$ for products). $\xi$ serves as a single progress variable along the reaction coordinate.

---

B. **Gibbs Free Energy of Reaction**

```{math}
\Delta_r G \equiv \sum_i \nu_i\,\mu_i
= \left(\frac{\partial G}{\partial \xi}\right)_{T,P}.
```

- **Applicability**: any reacting system at constant $T,P$. Signs the spontaneous direction: $\Delta_r G < 0$ forward, $\Delta_r G > 0$ backward.

---

C. **Chemical Equilibrium Condition**

```{math}
\sum_i \nu_i\,\mu_i = 0.
```

- **Applicability**: any reacting system at equilibrium, at constant $T,P$ *or* constant $T,V$. The equilibrium composition is the same in both cases; only the stationary potential ($G$ vs. $A$) differs.

---

D. **Ideal-Gas Chemical Potential and Reaction Quotient**

```{math}
\mu_i(T,P_i) = \mu_i^\circ(T) + RT\ln\!\left(\frac{P_i}{P^\circ}\right),
\qquad
Q_p = \prod_i\!\left(\frac{P_i}{P^\circ}\right)^{\!\nu_i}.
```

- **Applicability**: ideal-gas mixtures with $P^\circ = 1\ \mathrm{bar}$. For non-ideal systems, $P_i/P^\circ$ generalizes to an activity $a_i$, preserving the same structural form.

---

E. **Relation Between $\Delta_r G$, $\Delta_r G^\circ$, and $K$**

```{math}
\Delta_r G = \Delta_r G^\circ + RT\ln Q_p,
\qquad
\Delta_r G^\circ = -RT\ln K,
\qquad
K = \exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right).
```

- **Applicability**: ideal-gas reactions. Bridges tabulated standard thermochemistry to the equilibrium constant and to reaction direction via the comparison $Q_p \lessgtr K$.

---

F. **van't Hoff Equation** (differential and integrated forms)

```{math}
\left(\frac{\partial \ln K}{\partial T}\right)_P = \frac{\Delta_r H^\circ}{RT^2},
\qquad
\ln\!\left[\frac{K(T_2)}{K(T_1)}\right] = -\frac{\Delta_r H^\circ}{R}\!\left(\frac{1}{T_2}-\frac{1}{T_1}\right).
```

- **Applicability**: temperature dependence of $K$; integrated form assumes $\Delta_r H^\circ$ is approximately constant on $[T_1, T_2]$. Structurally identical to the integrated Clausius–Clapeyron equation of § 6.2 (both follow from Gibbs–Helmholtz).

---

G. **Partition Function of an Ideal-Gas Mixture**

```{math}
Q(T,V,\{N_i\}) = \prod_i \frac{q_i(T,V)^{N_i}}{N_i!}.
```

- **Applicability**: ideal gases (non-interacting species); follows from factorization over independent subsystems (§ 2.4) with the indistinguishability correction (§ 2.5).

---

H. **Chemical Potential from the Partition Function**

```{math}
\mu_i = -RT\ln\!\left(\frac{q_i}{N_i}\right)
\qquad (\text{per-mole, via Stirling}).
```

- **Applicability**: ideal-gas species in the thermodynamic limit, where Stirling's approximation is valid. Microscopic counterpart of the macroscopic expression in (D).

---

I. **Equilibrium Constant from Partition Functions**

```{math}
\prod_i N_i^{\nu_i} = \prod_i q_i^{\nu_i},
\qquad
K_c = \prod_i\!\left(\frac{q_i}{V}\right)^{\!\nu_i},
\qquad
K_p = K_c\,\left(\frac{k_{\mathrm B}T}{P^\circ}\right)^{\!\Delta\nu}.
```

- **Applicability**: ideal-gas reactions. When $\Delta\nu = 0$, $K_p = K_c$ and volume (hence pressure) dependence cancels in the ideal-gas model.

---

J. **Diatomic Single-Molecule Partition Function**

```{math}
q_{\text{diatomic}}
= \frac{V}{\Lambda^3}\,
  \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}\,
  \frac{e^{-\Theta_{\mathrm{vib}}/(2T)}}{1-e^{-\Theta_{\mathrm{vib}}/T}}\,
  g_1\,e^{D_0/(RT)}.
```

- **Applicability**: diatomic ideal gas in the rigid-rotor/harmonic-oscillator approximation (§ 2.8), with the energy zero chosen at the dissociated atoms at rest. $\sigma = 2$ for homonuclear diatomics, $\sigma = 1$ for heteronuclear.

---

K. **$K$ for $\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ from Molecular Parameters**

```{math}
K =
\left(\frac{m_{HI}^2}{m_{H_2}\,m_{I_2}}\right)^{\!3/2}
\frac{\sigma_{H_2}\sigma_{I_2}}{\sigma_{HI}^2}\,
\frac{\Theta_{\mathrm{rot}}^{H_2}\,\Theta_{\mathrm{rot}}^{I_2}}{(\Theta_{\mathrm{rot}}^{HI})^2}\,
\frac{(1-e^{-\Theta_{\mathrm{vib}}^{H_2}/T})(1-e^{-\Theta_{\mathrm{vib}}^{I_2}/T})}{(1-e^{-\Theta_{\mathrm{vib}}^{HI}/T})^2}\,
\exp\!\left(\frac{2D_0^{HI}-D_0^{H_2}-D_0^{I_2}}{RT}\right).
```

- **Applicability**: the specific $\Delta\nu = 0$ benchmark reaction in the rigid-rotor/harmonic-oscillator approximation. Separates $K$ into translational (mass), rotational (symmetry and $\Theta_{\mathrm{rot}}$), vibrational (hot-mode), and bond-energy factors — each traceable to identifiable molecular structure.
