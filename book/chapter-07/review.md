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

See [Section 7.1](section-01.md).

1. **Extent of Reaction**
   - For a reaction $a\,\mathrm{A} + b\,\mathrm{B} \rightleftharpoons y\,\mathrm{Y} + z\,\mathrm{Z}$, composition changes are tracked by a single scalar $\xi$ (units of moles) via

   ```{math}
   dn_i = \nu_i\,d\xi,
   ```

   with the **signed stoichiometry convention**: $\nu_i < 0$ for reactants ($\nu_{\mathrm{A}} = -a$, $\nu_{\mathrm{B}} = -b$) and $\nu_i > 0$ for products ($\nu_{\mathrm{Y}} = y$, $\nu_{\mathrm{Z}} = z$).
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

   - This is the multi-species analog of the phase-equilibrium condition $\mu_\alpha = \mu_\beta$ of [Section 6.1](../chapter-06/section-01.md): instead of moving matter between two phases, matter is redistributed among several chemical species.
   - **Spontaneity**: $\Delta_r G < 0$ drives the reaction forward, $\Delta_r G > 0$ backward, $\Delta_r G = 0$ is equilibrium.

4. **Ideal-Gas Chemical Potential and the Standard Pressure**
   - For a species in an ideal-gas mixture,

   ```{math}
   \mu_i(T,P_i) = \mu_i^\circ(T) + RT\ln\!\left(\frac{P_i}{P^\circ}\right),
   ```

   with $P^\circ = 1\ \mathrm{bar}$.
   - The logarithm contains the **ratio** $P_i/P^\circ$, not the pressure itself (Concept Check 2 of [Section 7.1](section-01.md)). For real gases and solutions, the ratio is replaced by an activity $a_i$.

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

   - The first term is linear in $\xi$ (the "bookkeeping" drop in standard chemical potentials); the second is the **ideal-mixing** term, $-T\,\Delta_{\mathrm{mix}}S$, which is zero for pure reactants and for pure products and negative in between.
   - The minimum of $G(\xi)$ lies where the slope $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$ is zero. Concept Check 4 of [Section 7.1](section-01.md) asks why this point lies between pure reactants and pure products even when $|\Delta_r G^\circ| \gg RT$.

7. **Temperature Dependence: the van't Hoff Equation**
   - Starting from the Gibbs–Helmholtz relation, $\left(\partial (G/T)/\partial T\right)_P = -H/T^2$, and $\Delta_r G^\circ = -RT\ln K$:

   ```{math}
   \left(\frac{\partial \ln K}{\partial T}\right)_P = \frac{\Delta_r H^\circ}{RT^2}.
   ```

   - Assuming $\Delta_r H^\circ$ is approximately constant over $[T_1, T_2]$:

   ```{math}
   \ln\!\left[\frac{K(T_2)}{K(T_1)}\right] = -\frac{\Delta_r H^\circ}{R}\!\left(\frac{1}{T_2} - \frac{1}{T_1}\right).
   ```

   - Exothermic reactions ($\Delta_r H^\circ < 0$) have $K$ that *decreases* with $T$; endothermic reactions have $K$ that *increases* with $T$. This is the quantitative form of Le Châtelier's principle for temperature.
   - The mathematical structure is identical to the integrated Clausius–Clapeyron equation of [Section 6.2](../chapter-06/section-02.md), with $K$ in place of the vapor pressure and $\Delta_r H^\circ$ in place of $\Delta \bar h_{\mathrm{vap}}$.

8. **Computing $K$ from Tabulated Thermochemistry**
   - The standard workflow for any gas-phase reaction at a specified $T$:
     1. Balance the reaction and fix a consistent standard state ($P^\circ = 1\ \mathrm{bar}$, usually at 298.15 K).
     2. Pull $\Delta_f H_i^\circ$ and $S_i^\circ$ (or $H_i^\circ(T)$ and $S_i^\circ(T)$ if $T \neq 298.15\ \mathrm{K}$) from NIST WebBook, NIST–JANAF, or ATcT.
     3. Form $\Delta_r H^\circ = \sum_i \nu_i H_i^\circ$ and $\Delta_r S^\circ = \sum_i \nu_i S_i^\circ$ using signed stoichiometry.
     4. Combine into $\Delta_r G^\circ(T) = \Delta_r H^\circ(T) - T\,\Delta_r S^\circ(T)$, and finally $K(T) = \exp(-\Delta_r G^\circ(T)/RT)$.
   - Common pitfalls: mixing phases (e.g., $\mathrm{H_2O(l)}$ vs. $\mathrm{H_2O(g)}$), mismatched standard states, and the near-universal unit mismatch between $\Delta H^\circ$ (kJ/mol) and $S^\circ$ (J/mol·K).

---

### Section 7.2: Equilibrium Constants from Microscopic Properties

See [Section 7.2](section-02.md).

1. **The Right Potential at Constant $T,V$**
   - At fixed $T$ and $V$, the Helmholtz free energy is the master potential ([Section 5.1](../chapter-05/section-01.md)), with differential

   ```{math}
   dA = -S\,dT - P\,dV + \sum_i \mu_i\,dn_i.
   ```

   - Equilibrium is the stationary point of $A$, yielding

   ```{math}
   \sum_i \nu_i\,\mu_i = 0,
   ```

   a condition of the same form as that of [Section 7.1](section-01.md). What changes between fixed $T,P$ and fixed $T,V$ is *which potential is stationary*. The chemical potential is the same quantity in both: $\mu_i = (\partial G/\partial n_i)_{T,P,n_{j\neq i}} = (\partial A/\partial n_i)_{T,V,n_{j\neq i}}$.

2. **Partition Function of an Ideal-Gas Mixture**
   - From factorization over independent subsystems and the indistinguishability correction ([Section 2.4](../chapter-02/section-04.md)),

   ```{math}
   Q(T,V,\{N_i\}) = \prod_i \frac{q_i(T,V)^{N_i}}{N_i!},
   ```

   where $q_i$ is the single-molecule partition function for species $i$. The factorization is exact for an ideal gas because molecules of different species do not interact.

3. **Chemical Potential from $A = -k_{\mathrm B}T\ln Q$**
   - Using $A = -k_{\mathrm B}T\ln Q$ ([Section 4.3](../chapter-04/section-03.md)) and differentiating with respect to $N_i$ under Stirling's approximation gives

   ```{math}
   \mu_i = -RT\ln\!\left(\frac{q_i}{N_i}\right) = RT\ln\!\left(\frac{N_i}{q_i}\right),
   ```

   on a per-mole basis. This is the microscopic counterpart of $\mu_i = \mu_i^\circ(T) + RT\ln(P_i/P^\circ)$ from [Section 7.1](section-01.md): the macroscopic expression writes $\mu_i$ in terms of a tabulated $\mu_i^\circ(T)$ and a partial pressure; the microscopic expression writes the same $\mu_i$ in terms of a counted ratio $q_i/N_i$.

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
   - In the rigid-rotor/harmonic-oscillator approximation ([Section 2.8](../chapter-02/section-08.md)), $q = q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}$:

   ```{math}
   q_{\mathrm{diatomic}}
   = \frac{V}{\Lambda^3}\,
     \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}\,
     \frac{1}{1-e^{-\Theta_{\mathrm{vib}}/T}}\,
     g_1\,e^{D_0/(RT)},
   \qquad
   \Lambda \equiv \sqrt{\frac{h^2}{2\pi m k_{\mathrm B}T}}.
   ```

   - **Energy-zero convention**: the zero of molecular energy is the dissociated atoms at rest, so the bound $v=0$ level sits at $-D_0$ and the bond energy enters the electronic factor as $e^{+D_0/RT}$. Vibrational energies are then measured from $v=0$, which gives $q_{\mathrm{vib}} = 1/(1-e^{-\Theta_{\mathrm{vib}}/T})$ in place of the form of [Section 2.6](../chapter-02/section-06.md). Either zero may be used, provided that $q_{\mathrm{vib}}$ and the bond energy refer to the same one (Concept Check 3 of [Section 7.2](section-02.md)).

6. **$\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ as a $\Delta\nu = 0$ Benchmark**
   - With $\Delta\nu = 0$, volume factors cancel exactly and

   ```{math}
   K = \frac{q_{\mathrm{HI}}^2}{q_{\mathrm{H_2}}\,q_{\mathrm{I_2}}}
   ```

   depends only on intrinsic molecular properties. Expanding with Eq. (J) below separates $K$ into four physically meaningful factors:

   ```{math}
   \begin{aligned}
   K ={}&
   \underbrace{\left(\frac{m_{\mathrm{HI}}^2}{m_{\mathrm{H_2}}\,m_{\mathrm{I_2}}}\right)^{3/2}}_{\text{translational}}
   \;
   \underbrace{\frac{\sigma_{\mathrm{H_2}}\sigma_{\mathrm{I_2}}}{\sigma_{\mathrm{HI}}^2}\,
   \frac{\Theta_{\mathrm{rot}}^{\mathrm{H_2}}\,\Theta_{\mathrm{rot}}^{\mathrm{I_2}}}{\left(\Theta_{\mathrm{rot}}^{\mathrm{HI}}\right)^2}}_{\text{rotational}}
   \\[1ex]
   &\times
   \underbrace{\frac{(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{H_2}}/T})(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{I_2}}/T})}{(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{HI}}/T})^2}}_{\text{vibrational}}
   \;
   \underbrace{\exp\!\left(\frac{2D_0^{\mathrm{HI}}-D_0^{\mathrm{H_2}}-D_0^{\mathrm{I_2}}}{RT}\right)}_{\text{bond-energy}}.
   \end{aligned}
   ```

   - With the energy zero at the separated atoms, the zero-point energies are contained in the $D_0$ values, and the vibrational block holds only the "hot-mode" $(1-e^{-\Theta_{\mathrm{vib}}/T})$ factors.

7. **Physical Reading of the Four Factors at 700 K**
   - Plugging spectroscopic constants for H$_2$, I$_2$, and HI gives $K(700\ \mathrm{K}) \approx 62$, to be compared with 56.8 from the NIST-JANAF tables and with a measured value near 54. Agreement within 10 to 15 % is what the rigid-rotor/harmonic-oscillator approximation can be expected to give; the largest single correction, about 4 %, comes from the high-temperature form of $q_{\mathrm{rot}}$ for $\mathrm{H_2}$.
   - **Translational** ($\approx 181$): dominated by the mass mismatch between H and I, since $m_{\mathrm{HI}}^2/(m_{\mathrm{H_2}}m_{\mathrm{I_2}}) = (m_{\mathrm{H}}+m_{\mathrm{I}})^2/(4m_{\mathrm{H}} m_{\mathrm{I}})$ is driven far from unity when $m_{\mathrm{I}} \gg m_{\mathrm{H}}$.
   - **Rotational** ($\approx 0.22$): I$_2$'s tiny $\Theta_{\mathrm{rot}} \approx 0.05\ \mathrm{K}$ makes $q_{\mathrm{rot}}(\mathrm{I_2})$ extremely large, favoring reactants; the symmetry-number prefactor $\sigma_{\mathrm{H_2}}\sigma_{\mathrm{I_2}}/\sigma_{\mathrm{HI}}^2 = 4$ nudges toward products but does not overcome this.
   - **Vibrational** ($\approx 0.36$): the stiff stretches of H$_2$ and HI remain nearly frozen at 700 K, so their factors are within 1 % of unity; the factor is essentially that of I$_2$, whose low-frequency mode is partially excited.
   - **Bond-energy** ($\approx 4.3$): $2D_0^{\mathrm{HI}} - D_0^{\mathrm{H_2}} - D_0^{\mathrm{I_2}} \approx +8.5\ \mathrm{kJ\,mol^{-1}}$ — the two HI bonds together are slightly stronger than one H–H plus one I–I bond.
   - Every digit of $K$ traces back to a specific molecular property; the relative sizes of the four factors identify *which* molecular feature most strongly determines the equilibrium.

8. **Consistency with the Macroscopic Route**
   - A student who computes $\Delta_r G^\circ$ from thermochemical tables ([Section 7.1](section-01.md)) and a student who evaluates $K$ from spectroscopic constants ([Section 7.2](section-02.md)) arrive at nearly the same number by very different paths. The macroscopic route answers "what is $K$?" efficiently; the microscopic route answers "*why* is $K$ what it is?" by decomposing it into separable molecular contributions. Together they close the statistical-mechanical bridge built across Modules 2, 4, 5, and 7.

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

- **Applicability**: any reacting system at equilibrium, at constant $T,P$ *or* constant $T,V$. The condition has the same form in both cases; the stationary potential ($G$ vs. $A$) differs.

---

D. **Ideal-Gas Chemical Potential and Reaction Quotient**

```{math}
\mu_i(T,P_i) = \mu_i^\circ(T) + RT\ln\!\left(\frac{P_i}{P^\circ}\right),
\qquad
Q_p = \prod_i\!\left(\frac{P_i}{P^\circ}\right)^{\!\nu_i}.
```

- **Applicability**: ideal-gas mixtures with $P^\circ = 1\ \mathrm{bar}$. For non-ideal systems, $P_i/P^\circ$ is replaced by an activity $a_i$.

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

- **Applicability**: temperature dependence of $K$; integrated form assumes $\Delta_r H^\circ$ is approximately constant on $[T_1, T_2]$. Structurally identical to the integrated Clausius–Clapeyron equation of [Section 6.2](../chapter-06/section-02.md).

---

G. **Partition Function of an Ideal-Gas Mixture**

```{math}
Q(T,V,\{N_i\}) = \prod_i \frac{q_i(T,V)^{N_i}}{N_i!}.
```

- **Applicability**: ideal gases (non-interacting species); follows from factorization over independent subsystems with the indistinguishability correction ([Section 2.4](../chapter-02/section-04.md)).

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
q_{\mathrm{diatomic}}
= \frac{V}{\Lambda^3}\,
  \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}\,
  \frac{1}{1-e^{-\Theta_{\mathrm{vib}}/T}}\,
  g_1\,e^{D_0/(RT)}.
```

- **Applicability**: diatomic ideal gas in the rigid-rotor/harmonic-oscillator approximation ([Section 2.8](../chapter-02/section-08.md)), with the energy zero chosen at the dissociated atoms at rest, so that vibrational energies are measured from $v=0$. $\sigma = 2$ for homonuclear diatomics, $\sigma = 1$ for heteronuclear.

---

K. **$K$ for $\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ from Molecular Parameters**

```{math}
\begin{aligned}
K ={}&
\left(\frac{m_{\mathrm{HI}}^2}{m_{\mathrm{H_2}}\,m_{\mathrm{I_2}}}\right)^{\!3/2}
\frac{\sigma_{\mathrm{H_2}}\sigma_{\mathrm{I_2}}}{\sigma_{\mathrm{HI}}^2}\,
\frac{\Theta_{\mathrm{rot}}^{\mathrm{H_2}}\,\Theta_{\mathrm{rot}}^{\mathrm{I_2}}}{(\Theta_{\mathrm{rot}}^{\mathrm{HI}})^2}
\\[1ex]
&\times
\frac{(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{H_2}}/T})(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{I_2}}/T})}{(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{HI}}/T})^2}\,
\exp\!\left(\frac{2D_0^{\mathrm{HI}}-D_0^{\mathrm{H_2}}-D_0^{\mathrm{I_2}}}{RT}\right).
\end{aligned}
```

- **Applicability**: the specific $\Delta\nu = 0$ benchmark reaction in the rigid-rotor/harmonic-oscillator approximation. Separates $K$ into translational (mass), rotational (symmetry and $\Theta_{\mathrm{rot}}$), vibrational (hot-mode), and bond-energy factors — each traceable to identifiable molecular structure.

## 3. Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. Nitric oxide in hot air

Draws on [Section 7.1](section-01.md).

Nitric oxide forms from the nitrogen and oxygen of air in flames and engines:

```{math}
\mathrm{N_2}(g)+\mathrm{O_2}(g)\;\rightleftharpoons\;2\,\mathrm{NO}(g).
```

Model air as $4.00\ \mathrm{mol}$ of $\mathrm{N_2}$ and $1.00\ \mathrm{mol}$ of $\mathrm{O_2}$. The table lists data for $\mathrm{NO}(g)$ from the [NIST-JANAF tables](https://janaf.nist.gov/). Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

| $T$ (K) | $\Delta_f G^\circ$ ($\mathrm{kJ\,mol^{-1}}$) | $\Delta_f H^\circ$ ($\mathrm{kJ\,mol^{-1}}$) |
|---:|---:|---:|
| 1000 | 77.775 | 90.437 |
| 2000 | 65.060 | 90.494 |

**(a)** Find $\Delta_r G^\circ$ and $K$ at $2000\ \mathrm{K}$.

**(b)** Write $K$ in terms of the extent of reaction at equilibrium, $\xi_{\mathrm{eq}}$. Find $\xi_{\mathrm{eq}}$ and the mole fraction of $\mathrm{NO}$ at $2000\ \mathrm{K}$.

**(c)** Estimate $K$ at $2500\ \mathrm{K}$ with the van't Hoff equation and $\Delta_r H^\circ$ at $2000\ \mathrm{K}$. Compare the estimate with $K=3.52\times10^{-3}$ from the same tables.

**(d)** The equilibrium mixture of part (b) is cooled to $1000\ \mathrm{K}$ so quickly that its composition does not change. Find $Q_p$, $K$, and $\Delta_r G$ at $1000\ \mathrm{K}$. In which direction is the reaction spontaneous?

````{admonition} Solution
:class: dropdown

1. **$\Delta_r G^\circ$ and $K$ at $2000\ \mathrm{K}$**

   $\mathrm{N_2}$ and $\mathrm{O_2}$ are elements in their reference states, so $\Delta_f G^\circ=0$ for both and

   ```{math}
   \begin{aligned}
   \Delta_r G^\circ&=2\,\Delta_f G^\circ(\mathrm{NO})=2(65.060)=130.12\ \mathrm{kJ\,mol^{-1}},\\
   K&=\exp\left(-\frac{\Delta_r G^\circ}{RT}\right)=\exp\left(-\frac{130\,120}{(8.3145)(2000)}\right)\\
   &=e^{-7.8249}=4.00\times10^{-4}.
   \end{aligned}
   ```

2. **$K$ in terms of $\xi_{\mathrm{eq}}$**

   With $dn_i=\nu_i\,d\xi$, the amounts are $n_{\mathrm{N_2}}=4-\xi$, $n_{\mathrm{O_2}}=1-\xi$, and $n_{\mathrm{NO}}=2\xi$, and the total is $5\ \mathrm{mol}$ at every $\xi$. With $P_i=y_iP$,

   ```{math}
   K=\frac{(P_{\mathrm{NO}}/P^\circ)^2}{(P_{\mathrm{N_2}}/P^\circ)(P_{\mathrm{O_2}}/P^\circ)}
   =\frac{y_{\mathrm{NO}}^2}{y_{\mathrm{N_2}}\,y_{\mathrm{O_2}}}
   =\frac{4\,\xi_{\mathrm{eq}}^2}{(4-\xi_{\mathrm{eq}})(1-\xi_{\mathrm{eq}})}.
   ```

   The total pressure cancels because $\Delta\nu=0$.

3. **Equilibrium composition**

   Clearing the fraction gives a quadratic equation,

   ```{math}
   (4-K)\,\xi_{\mathrm{eq}}^2+5K\,\xi_{\mathrm{eq}}-4K=0 ,
   ```

   whose positive root is

   ```{math}
   \begin{aligned}
   \xi_{\mathrm{eq}}&=\frac{-5K+\sqrt{25K^2+16K(4-K)}}{2(4-K)}\\
   &=\frac{-0.001998+0.159937}{7.9992}=0.0197\ \mathrm{mol}.
   \end{aligned}
   ```

   The mole fraction of nitric oxide is $y_{\mathrm{NO}}=2\xi_{\mathrm{eq}}/5=0.0079$. Because $\xi_{\mathrm{eq}}\ll1$, the estimate $\xi_{\mathrm{eq}}\approx\sqrt{K}=0.0200$ is already close.

4. **van't Hoff estimate at $2500\ \mathrm{K}$**

   $\Delta_r H^\circ=2\,\Delta_f H^\circ(\mathrm{NO})=2(90.494)=180.99\ \mathrm{kJ\,mol^{-1}}$ at $2000\ \mathrm{K}$. With $T_1=2000\ \mathrm{K}$ and $T_2=2500\ \mathrm{K}$,

   ```{math}
   \begin{aligned}
   \ln\left[\frac{K(T_2)}{K(T_1)}\right]&=-\frac{\Delta_r H^\circ}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right)\\
   &=-(21\,768\ \mathrm{K})\left(-1.000\times10^{-4}\ \mathrm{K^{-1}}\right)=2.1768,\\
   K(T_2)&=(4.00\times10^{-4})\,e^{2.1768}=3.52\times10^{-3}.
   \end{aligned}
   ```

   The estimate agrees with the tabulated $3.52\times10^{-3}$ to the digits shown. The reaction is endothermic, and $K$ rises by a factor of 8.8 over $500\ \mathrm{K}$.

5. **The cooled mixture**

   The composition is that of step 3, so $Q_p$ has the value that $K$ had at $2000\ \mathrm{K}$: $Q_p=4.00\times10^{-4}$. At $1000\ \mathrm{K}$, $\Delta_r G^\circ=2(77.775)=155.55\ \mathrm{kJ\,mol^{-1}}$ and

   ```{math}
   \begin{aligned}
   K&=\exp\left(-\frac{155\,550}{(8.3145)(1000)}\right)=e^{-18.708}=7.5\times10^{-9},\\
   \Delta_r G&=\Delta_r G^\circ+RT\ln Q_p\\
   &=155.55+(8.3145\times10^{-3})(1000)\ln(4.00\times10^{-4})\\
   &=155.55-65.06=+90.49\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

   $Q_p>K$ and $\Delta_r G>0$, so the reverse reaction is spontaneous: nitric oxide should decompose into $\mathrm{N_2}$ and $\mathrm{O_2}$. This result gives the direction of the change and says nothing about its rate.

**Result.** $K=4.00\times10^{-4}$ at $2000\ \mathrm{K}$, with $\xi_{\mathrm{eq}}=0.0197\ \mathrm{mol}$ and $y_{\mathrm{NO}}=0.0079$. The van't Hoff equation gives $K=3.52\times10^{-3}$ at $2500\ \mathrm{K}$. At $1000\ \mathrm{K}$ the frozen mixture has $Q_p=4.00\times10^{-4}>K=7.5\times10^{-9}$ and $\Delta_r G=+90.49\ \mathrm{kJ\,mol^{-1}}$, so the reverse reaction is spontaneous.
````

### Problem 2. Pressure and the extent of a dissociation

Draws on [Section 7.1](section-01.md).

Phosgene dissociates when it is heated:

```{math}
\mathrm{COCl_2}(g)\;\rightleftharpoons\;\mathrm{CO}(g)+\mathrm{Cl_2}(g).
```

At $700\ \mathrm{K}$ the [NIST-JANAF tables](https://janaf.nist.gov/) give $\Delta_f G^\circ=-187.054\ \mathrm{kJ\,mol^{-1}}$ for $\mathrm{COCl_2}(g)$ and $-173.518\ \mathrm{kJ\,mol^{-1}}$ for $\mathrm{CO}(g)$. A vessel initially holds $1.00\ \mathrm{mol}$ of phosgene and no products. Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$ and $P^\circ=1\ \mathrm{bar}$.

**(a)** Find $\Delta_r G^\circ$ and $K$ at $700\ \mathrm{K}$.

**(b)** Write $K$ in terms of $\xi_{\mathrm{eq}}$ and the total pressure $P$, and solve the expression for $\xi_{\mathrm{eq}}$.

**(c)** Find $\xi_{\mathrm{eq}}$ and the three partial pressures at $P=1.00\ \mathrm{bar}$ and at $P=0.100\ \mathrm{bar}$.

````{admonition} Solution
:class: dropdown

1. **$\Delta_r G^\circ$ and $K$**

   $\mathrm{Cl_2}(g)$ is an element in its reference state, so its $\Delta_f G^\circ$ is zero.

   ```{math}
   \begin{aligned}
   \Delta_r G^\circ&=\left[-173.518+0\right]-\left[-187.054\right]=+13.536\ \mathrm{kJ\,mol^{-1}},\\
   K&=\exp\left(-\frac{13\,536}{(8.3145)(700)}\right)=e^{-2.3257}=0.0977.
   \end{aligned}
   ```

2. **$K$ in terms of $\xi_{\mathrm{eq}}$ and $P$**

   The amounts are $n_{\mathrm{COCl_2}}=1-\xi$ and $n_{\mathrm{CO}}=n_{\mathrm{Cl_2}}=\xi$, with $n_{\mathrm{tot}}=1+\xi$. The mole fractions are

   ```{math}
   y_{\mathrm{COCl_2}}=\frac{1-\xi}{1+\xi},\qquad y_{\mathrm{CO}}=y_{\mathrm{Cl_2}}=\frac{\xi}{1+\xi}.
   ```

   With $P_i=y_iP$,

   ```{math}
   K=\frac{(P_{\mathrm{CO}}/P^\circ)(P_{\mathrm{Cl_2}}/P^\circ)}{P_{\mathrm{COCl_2}}/P^\circ}
   =\frac{\xi_{\mathrm{eq}}^2}{(1-\xi_{\mathrm{eq}})(1+\xi_{\mathrm{eq}})}\,\frac{P}{P^\circ}
   =\frac{\xi_{\mathrm{eq}}^2}{1-\xi_{\mathrm{eq}}^2}\,\frac{P}{P^\circ}.
   ```

   Solving for the extent of reaction gives

   ```{math}
   \xi_{\mathrm{eq}}=\sqrt{\frac{K}{K+P/P^\circ}}.
   ```

3. **Two pressures**

   ```{math}
   \begin{aligned}
   P=1.00\ \mathrm{bar}:\quad \xi_{\mathrm{eq}}&=\sqrt{\frac{0.0977}{0.0977+1.00}}=0.298\ \mathrm{mol},\\
   P=0.100\ \mathrm{bar}:\quad \xi_{\mathrm{eq}}&=\sqrt{\frac{0.0977}{0.0977+0.100}}=0.703\ \mathrm{mol}.
   \end{aligned}
   ```

   The partial pressures follow from $P_i=y_iP$.

   | $P$ (bar) | $\xi_{\mathrm{eq}}$ (mol) | $P_{\mathrm{COCl_2}}$ (bar) | $P_{\mathrm{CO}}=P_{\mathrm{Cl_2}}$ (bar) |
   |---:|---:|---:|---:|
   | 1.00 | 0.298 | 0.540 | 0.230 |
   | 0.100 | 0.703 | 0.0174 | 0.0413 |

   As a check, $(0.230)^2/0.540=0.098$ and $(0.0413)^2/0.0174=0.098$: both mixtures have $Q_p=K$. The dissociation has $\Delta\nu=+1$, and a tenfold drop in pressure raises $\xi_{\mathrm{eq}}$ from $0.298$ to $0.703\ \mathrm{mol}$ while $K$ stays the same.

**Result.** $\Delta_r G^\circ=+13.536\ \mathrm{kJ\,mol^{-1}}$ and $K=0.0977$ at $700\ \mathrm{K}$. $\xi_{\mathrm{eq}}=\sqrt{K/(K+P/P^\circ)}$, which is $0.298\ \mathrm{mol}$ at $1.00\ \mathrm{bar}$ and $0.703\ \mathrm{mol}$ at $0.100\ \mathrm{bar}$.
````

### Problem 3. The Gibbs energy minimum of hot carbon dioxide

Draws on [Section 7.1](section-01.md).

At $3000\ \mathrm{K}$ carbon dioxide dissociates appreciably:

```{math}
\mathrm{CO_2}(g)\;\rightleftharpoons\;\mathrm{CO}(g)+\tfrac{1}{2}\,\mathrm{O_2}(g).
```

At this temperature the [NIST-JANAF tables](https://janaf.nist.gov/) give $\Delta_f G^\circ=-395.461\ \mathrm{kJ\,mol^{-1}}$ for $\mathrm{CO_2}(g)$ and $-367.826\ \mathrm{kJ\,mol^{-1}}$ for $\mathrm{CO}(g)$. The system starts as $1.00\ \mathrm{mol}$ of $\mathrm{CO_2}$ and stays at $3000\ \mathrm{K}$ and $P=P^\circ=1\ \mathrm{bar}$. Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find $\Delta_r G^\circ$ and $K$.

**(b)** Write the amounts $n_i(\xi)$ and the mole fractions $y_i(\xi)$. Then write $G(\xi)-G(0)$ in the form used for $\mathrm{NO_2}$ in [Section 7.1](section-01.md).

**(c)** Write $Q_p$ as a function of $\xi$, and confirm that $\xi_{\mathrm{eq}}=0.438\ \mathrm{mol}$ satisfies $Q_p=K$.

**(d)** Evaluate $G(\xi)-G(0)$ at $\xi_{\mathrm{eq}}$ and at $\xi=1\ \mathrm{mol}$. Evaluate the slope $\Delta_r G$ at $\xi=0.200\ \mathrm{mol}$ and at $\xi=0.700\ \mathrm{mol}$.

````{admonition} Solution
:class: dropdown

1. **$\Delta_r G^\circ$ and $K$**

   ```{math}
   \begin{aligned}
   \Delta_r G^\circ&=\left[-367.826+\tfrac{1}{2}(0)\right]-\left[-395.461\right]=+27.635\ \mathrm{kJ\,mol^{-1}},\\
   K&=\exp\left(-\frac{27\,635}{(8.3145)(3000)}\right)=e^{-1.1079}=0.3303.
   \end{aligned}
   ```

2. **Amounts, mole fractions, and $G(\xi)$**

   With $\nu_{\mathrm{CO_2}}=-1$, $\nu_{\mathrm{CO}}=+1$, and $\nu_{\mathrm{O_2}}=+\tfrac{1}{2}$,

   | Species | $n_i$ (mol) | $y_i$ |
   |---|---|---|
   | $\mathrm{CO_2}$ | $1-\xi$ | $(1-\xi)/(1+\xi/2)$ |
   | $\mathrm{CO}$ | $\xi$ | $\xi/(1+\xi/2)$ |
   | $\mathrm{O_2}$ | $\xi/2$ | $(\xi/2)/(1+\xi/2)$ |

   and $n_{\mathrm{tot}}=1+\xi/2$. At $P=P^\circ$,

   ```{math}
   \begin{aligned}
   G(\xi)-G(0)={}&\xi\,\Delta_r G^\circ\\
   &+RT\left[(1-\xi)\ln y_{\mathrm{CO_2}}+\xi\ln y_{\mathrm{CO}}+\frac{\xi}{2}\ln y_{\mathrm{O_2}}\right].
   \end{aligned}
   ```

   The first term is the standard term, and the second is the ideal-mixing term.

3. **The equilibrium condition**

   At $P=P^\circ$,

   ```{math}
   Q_p=\frac{y_{\mathrm{CO}}\,y_{\mathrm{O_2}}^{1/2}}{y_{\mathrm{CO_2}}}
   =\frac{\xi}{1-\xi}\left(\frac{\xi/2}{1+\xi/2}\right)^{1/2}.
   ```

   At $\xi=0.438\ \mathrm{mol}$ the mole fractions are $y_{\mathrm{CO_2}}=0.562/1.219=0.4610$, $y_{\mathrm{CO}}=0.438/1.219=0.3593$, and $y_{\mathrm{O_2}}=0.219/1.219=0.1797$, so

   ```{math}
   Q_p=\frac{(0.3593)(0.1797)^{1/2}}{0.4610}=0.3304 .
   ```

   This equals $K=0.3303$ to three figures, so $\xi_{\mathrm{eq}}=0.438\ \mathrm{mol}$.

4. **$G$ at the minimum and at complete reaction**

   Here $RT=(8.3145)(3000)=24.94\ \mathrm{kJ\,mol^{-1}}$. At $\xi_{\mathrm{eq}}$,

   ```{math}
   \begin{aligned}
   G-G(0)={}&(0.438)(27.635)\\
   &+(24.94)\left[0.562\ln0.4610+0.438\ln0.3593+0.219\ln0.1797\right]\\
   ={}&12.10+(24.94)(-1.2594)=12.10-31.41=-19.31\ \mathrm{kJ}.
   \end{aligned}
   ```

   At $\xi=1\ \mathrm{mol}$ the mixture is $1\ \mathrm{mol}$ of $\mathrm{CO}$ and $\tfrac{1}{2}\ \mathrm{mol}$ of $\mathrm{O_2}$, with $y_{\mathrm{CO}}=\tfrac{2}{3}$ and $y_{\mathrm{O_2}}=\tfrac{1}{3}$:

   ```{math}
   \begin{aligned}
   G(1)-G(0)&=27.635+(24.94)\left[\ln\tfrac{2}{3}+\tfrac{1}{2}\ln\tfrac{1}{3}\right]\\
   &=27.635+(24.94)(-0.9548)=27.63-23.82=+3.82\ \mathrm{kJ}.
   \end{aligned}
   ```

5. **Slopes on the two sides of the minimum**

   The slope is $\Delta_r G=\Delta_r G^\circ+RT\ln Q_p$, with $Q_p$ from step 3.

   | $\xi$ (mol) | $Q_p$ | $RT\ln Q_p$ ($\mathrm{kJ\,mol^{-1}}$) | $\Delta_r G$ ($\mathrm{kJ\,mol^{-1}}$) |
   |---:|---:|---:|---:|
   | 0.200 | 0.0754 | −64.48 | −36.85 |
   | 0.438 | 0.3303 | −27.63 | 0.00 |
   | 0.700 | 1.188 | +4.30 | +31.93 |

   The slope is negative at $\xi=0.200\ \mathrm{mol}$, zero at $\xi_{\mathrm{eq}}$, and positive at $\xi=0.700\ \mathrm{mol}$.

**Result.** $K=0.3303$ and $\xi_{\mathrm{eq}}=0.438\ \mathrm{mol}$. Relative to pure $\mathrm{CO_2}$, $G$ is $-19.31\ \mathrm{kJ}$ at $\xi_{\mathrm{eq}}$ and $+3.82\ \mathrm{kJ}$ at $\xi=1\ \mathrm{mol}$. $\Delta_r G$ is $-36.85\ \mathrm{kJ\,mol^{-1}}$ at $\xi=0.200\ \mathrm{mol}$ and $+31.93\ \mathrm{kJ\,mol^{-1}}$ at $\xi=0.700\ \mathrm{mol}$.
````

### Problem 4. Dissociation of iodine from molecular constants

Draws on [Section 7.2](section-02.md).

Above about $1000\ \mathrm{K}$, iodine vapor contains a measurable fraction of atoms:

```{math}
\mathrm{I_2}(g)\;\rightleftharpoons\;2\,\mathrm{I}(g).
```

The table lists the constants of $\mathrm{I_2}$ from [Section 7.2](section-02.md) and those of the atom, whose ground level has the degeneracy $g_1=4$.

| Species | $M$ ($\mathrm{g\,mol^{-1}}$) | $\Theta_{\mathrm{rot}}$ (K) | $\Theta_{\mathrm{vib}}$ (K) | $D_0$ ($\mathrm{kJ\,mol^{-1}}$) | $\sigma$ | $g_1$ |
| :-- | --: | --: | --: | --: | --: | --: |
| $\mathrm{I_2}$ | 253.808 | 0.05377 | 309 | 148.81 | 2 | 1 |
| $\mathrm{I}$ | 126.904 | | | | | 4 |

Take the zero of energy at the separated atoms, as in [Section 7.2](section-02.md), so that $q_{\mathrm{I}}=(V/\Lambda_{\mathrm{I}}^3)\,g_1$. Use $h=6.6261\times10^{-34}\ \mathrm{J\,s}$, $k_{\mathrm B}=1.3806\times10^{-23}\ \mathrm{J\,K^{-1}}$, $N_{\mathrm A}=6.0221\times10^{23}\ \mathrm{mol^{-1}}$, $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$, and $P^\circ=1\ \mathrm{bar}$.

**(a)** Find $q/V$ for the atom and for the molecule at $1000\ \mathrm{K}$.

**(b)** Find $K_c$ and $K_p$ at $1000\ \mathrm{K}$.

**(c)** Find $\Delta_r G^\circ$ from $K_p$. Compare it with the [NIST-JANAF](https://janaf.nist.gov/) value $\Delta_f G^\circ(\mathrm{I},g)=24.039\ \mathrm{kJ\,mol^{-1}}$ at $1000\ \mathrm{K}$, for which the reference state of iodine is $\mathrm{I_2}(g)$.

**(d)** One mole of $\mathrm{I_2}$ comes to equilibrium at $1000\ \mathrm{K}$ and a total pressure of $1.00\ \mathrm{bar}$. Find $\xi_{\mathrm{eq}}$ and the two partial pressures.

````{admonition} Solution
:class: dropdown

1. **The atom**

   The mass of one atom is $m_{\mathrm{I}}=M/N_{\mathrm A}=2.107\times10^{-25}\ \mathrm{kg}$, so

   ```{math}
   \begin{aligned}
   \Lambda_{\mathrm{I}}&=\frac{h}{\sqrt{2\pi m_{\mathrm{I}}k_{\mathrm B}T}}\\
   &=\frac{6.6261\times10^{-34}}{\sqrt{2\pi(2.107\times10^{-25})(1.3806\times10^{-23})(1000)}}\\
   &=4.901\times10^{-12}\ \mathrm{m},\\
   \frac{q_{\mathrm{I}}}{V}&=\frac{g_1}{\Lambda_{\mathrm{I}}^3}=\frac{4}{(4.901\times10^{-12}\ \mathrm{m})^3}=3.398\times10^{34}\ \mathrm{m^{-3}}.
   \end{aligned}
   ```

2. **The molecule**

   The molecule has twice the mass, so $\Lambda_{\mathrm{I_2}}=\Lambda_{\mathrm{I}}/\sqrt{2}=3.465\times10^{-12}\ \mathrm{m}$ and $1/\Lambda_{\mathrm{I_2}}^3=2.403\times10^{34}\ \mathrm{m^{-3}}$. The other three factors of the diatomic partition function are

   ```{math}
   \begin{aligned}
   \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}&=\frac{1000}{2(0.05377)}=9299,\\
   \frac{1}{1-e^{-\Theta_{\mathrm{vib}}/T}}&=\frac{1}{1-e^{-0.309}}=3.762,\\
   g_1\,e^{D_0/(RT)}&=\exp\left(\frac{148\,810}{(8.3145)(1000)}\right)=e^{17.8976}=5.927\times10^{7}.
   \end{aligned}
   ```

   Their product with the translational factor is

   ```{math}
   \begin{aligned}
   \frac{q_{\mathrm{I_2}}}{V}&=(2.403\times10^{34})(9299)(3.762)(5.927\times10^{7})\ \mathrm{m^{-3}}\\
   &=4.982\times10^{46}\ \mathrm{m^{-3}}.
   \end{aligned}
   ```

3. **$K_c$ and $K_p$**

   For this reaction $\nu_{\mathrm{I}}=+2$, $\nu_{\mathrm{I_2}}=-1$, and $\Delta\nu=+1$.

   ```{math}
   \begin{aligned}
   K_c&=\frac{(q_{\mathrm{I}}/V)^2}{q_{\mathrm{I_2}}/V}=\frac{(3.398\times10^{34})^2}{4.982\times10^{46}}=2.318\times10^{22}\ \mathrm{m^{-3}},\\
   K_p&=K_c\,\frac{k_{\mathrm B}T}{P^\circ}=(2.318\times10^{22}\ \mathrm{m^{-3}})\,\frac{(1.3806\times10^{-23})(1000)\ \mathrm{J}}{10^{5}\ \mathrm{Pa}}\\
   &=3.20\times10^{-3}.
   \end{aligned}
   ```

4. **$\Delta_r G^\circ$ and the tabulated value**

   ```{math}
   \begin{aligned}
   \Delta_r G^\circ=-RT\ln K_p&=-(8.3145)(1000)\ln(3.20\times10^{-3})\ \mathrm{J\,mol^{-1}}\\
   &=+47.76\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

   The tables give $\Delta_r G^\circ=2(24.039)=48.08\ \mathrm{kJ\,mol^{-1}}$, which corresponds to $K=3.08\times10^{-3}$. The model value of $\Delta_r G^\circ$ is about $0.3\ \mathrm{kJ\,mol^{-1}}$ lower, and its $K_p$ is $4\ \%$ higher.

5. **Composition at $1.00\ \mathrm{bar}$**

   The amounts are $n_{\mathrm{I_2}}=1-\xi$ and $n_{\mathrm{I}}=2\xi$, with $n_{\mathrm{tot}}=1+\xi$. Then

   ```{math}
   K_p=\frac{(P_{\mathrm{I}}/P^\circ)^2}{P_{\mathrm{I_2}}/P^\circ}
   =\frac{4\,\xi_{\mathrm{eq}}^2}{1-\xi_{\mathrm{eq}}^2}\,\frac{P}{P^\circ}
   \qquad\Longrightarrow\qquad
   \xi_{\mathrm{eq}}=\sqrt{\frac{K_p}{K_p+4P/P^\circ}}.
   ```

   With $K_p=3.20\times10^{-3}$ and $P=P^\circ$,

   ```{math}
   \begin{aligned}
   \xi_{\mathrm{eq}}&=\sqrt{\frac{3.20\times10^{-3}}{4.0032}}=0.0283\ \mathrm{mol},\\
   P_{\mathrm{I}}&=\frac{2\xi_{\mathrm{eq}}}{1+\xi_{\mathrm{eq}}}\,P=0.0550\ \mathrm{bar},\qquad
   P_{\mathrm{I_2}}=\frac{1-\xi_{\mathrm{eq}}}{1+\xi_{\mathrm{eq}}}\,P=0.945\ \mathrm{bar}.
   \end{aligned}
   ```

   About $2.8\ \%$ of the molecules are dissociated.

**Result.** At $1000\ \mathrm{K}$, $q_{\mathrm{I}}/V=3.398\times10^{34}\ \mathrm{m^{-3}}$ and $q_{\mathrm{I_2}}/V=4.982\times10^{46}\ \mathrm{m^{-3}}$, which give $K_c=2.318\times10^{22}\ \mathrm{m^{-3}}$, $K_p=3.20\times10^{-3}$, and $\Delta_r G^\circ=+47.76\ \mathrm{kJ\,mol^{-1}}$ (NIST-JANAF: $48.08\ \mathrm{kJ\,mol^{-1}}$). At $1.00\ \mathrm{bar}$, $\xi_{\mathrm{eq}}=0.0283\ \mathrm{mol}$, $P_{\mathrm{I}}=0.0550\ \mathrm{bar}$, and $P_{\mathrm{I_2}}=0.945\ \mathrm{bar}$.
````
