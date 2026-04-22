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


# 7.2. Equilibrium Constants from Microscopic Properties

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Section 7.1 derived the equilibrium constant from macroscopic thermochemistry: given $\Delta_r H^\circ$ and $\Delta_r S^\circ$ (e.g., from NIST/JANAF/ATcT tables), one computes $\Delta_r G^\circ$ and then $K = \exp(-\Delta_r G^\circ/RT)$.

In this section we derive the same $K$ from **microscopic** information — molecular masses, rotational and vibrational constants, electronic degeneracies, and bond dissociation energies. This completes the statistical-mechanical bridge built in Chapters 2, 4, and 5: Chapter 4 connected the canonical partition function to the Helmholtz free energy, $A = -k_{\mathrm B}T\ln Q$; Chapter 5 showed that $A$ is the master potential at constant $T,V$, generating $S$, $U$, $P$, and — as we now develop — the chemical potentials $\mu_i$. Once we have $\mu_i$ in terms of the single-molecule partition functions $q_i$, the equilibrium condition $\sum_i\nu_i\mu_i = 0$ gives $K$ directly in terms of molecular properties.

The payoff: for any gas-phase reaction, we can predict $K(T)$ from spectroscopic data alone, and we can *see* which molecular features (masses, moments of inertia, bond strengths) drive the equilibrium.

```{admonition} Notation note
:class: note

In this section, partition functions are $Q$ (many-particle, canonical) and $q$ (single-particle), consistent with Chapters 2, 4, and 5 and the [course-wide conventions](../notation.md). The reaction quotient $Q_p$ of § 7.1 does not appear in this section, so no notational collision arises. All partition functions here are of the single- or many-particle kind.
```

Learning objectives:

- Write the equilibrium condition at fixed $T,V$ as $\sum_i \nu_i\mu_i=0$ from minimizing the Helmholtz free energy.
- Use $Q = \prod_i q_i^{N_i}/N_i!$ for an ideal-gas mixture to derive $\mu_i$ from $A = -k_{\mathrm B}T\ln Q$.
- Express equilibrium constants in terms of single-molecule partition functions, and relate $K_p$ and $K_c$ through $\Delta\nu$.
- Identify how translational, rotational, vibrational, and electronic factors contribute to $K$ for a specific reaction, and compute $K(T)$ numerically from molecular constants.

## Core Ideas and Derivations

### 7.2.1 The right potential for equilibrium at constant $T,V$

Consider a general gas-phase reaction

```{math}
\nu_A A(g) + \nu_B B(g) \rightleftharpoons \nu_Y Y(g) + \nu_Z Z(g),
```

under the **ideal-gas** assumption used throughout this chapter.

Section 5.1 established that at constant $T$ and $V$, the relevant thermodynamic potential is the **Helmholtz free energy** $A$, with differential

```{math}
dA = -S\,dT - P\,dV + \sum_i \mu_i\,dn_i.
```

At constant $T$ and $V$,

```{math}
dA = \sum_i \mu_i\,dn_i = \left(\sum_i \nu_i\,\mu_i\right)d\xi,
```

where the second equality uses $dn_i = \nu_i\,d\xi$ from § 7.1.1. Equilibrium at fixed $T,V$ occurs at a minimum of $A$, so

```{math}
\left(\frac{\partial A}{\partial \xi}\right)_{T,V} = 0
\qquad\Longrightarrow\qquad
\boxed{\sum_i \nu_i\,\mu_i = 0.}
```

This is the **chemical equilibrium condition** written in terms of chemical potentials. It is the same condition as in § 7.1 — the equilibrium position of a reaction does not depend on whether we hold $P$ or $V$ fixed, because $\mu_i$ is intensive. What changes between the two cases is only *which potential is stationary*: $G$ at fixed $T,P$, $A$ at fixed $T,V$.

---

### 7.2.2 Partition function of an ideal-gas mixture

From § 2.4 (factorization over independent subsystems) and § 2.5 (indistinguishable particles), the canonical partition function of a mixture of ideal gases factorizes:

```{math}
:label: eq:Q-mixture
\boxed{
Q(T,V,\{N_i\}) = \prod_i Q_i(T,V,N_i),
\qquad
Q_i(T,V,N_i) = \frac{q_i(T,V)^{N_i}}{N_i!},
}
```

where $q_i$ is the **single-molecule** partition function for species $i$. Each species is treated independently because ideal-gas molecules don't interact.

---

### 7.2.3 Chemical potentials from the partition function

Section 4.3 established

```{math}
:label: eq:A-from-Q
A = -k_{\mathrm B}T\ln Q,
```

and § 5.1 developed $A$ as the "master potential" at fixed $T,V$, with $S$, $U$, $P$ all following from its derivatives. The chemical potential $\mu_i$ is the natural extension of that derivative machinery to composition:

```{math}
\mu_i^{\text{(per molecule)}}
\equiv \left(\frac{\partial A}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}
= -k_{\mathrm B}T\left(\frac{\partial \ln Q}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}.
```

Using Eq. {eq}`eq:Q-mixture`,

```{math}
\ln Q = \sum_i \left(N_i\ln q_i - \ln N_i!\right),
```

and applying Stirling's approximation in the thermodynamic limit ($\ln N_i! \approx N_i\ln N_i - N_i$, as introduced in § 2.5),

```{math}
\left(\frac{\partial \ln Q}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}
= \ln q_i - \ln N_i = \ln\!\left(\frac{q_i}{N_i}\right).
```

So the per-molecule chemical potential is

```{math}
\mu_i^{\text{(per molecule)}} = -k_{\mathrm B}T\ln\!\left(\frac{q_i}{N_i}\right).
```

Multiplying by Avogadro's number $N_{\mathrm A}$ to express $\mu_i$ on a **per-mole** basis (as used throughout § 7.1),

```{math}
:label: eq:mu-from-q
\boxed{
\mu_i = -RT\ln\!\left(\frac{q_i}{N_i}\right) = RT\ln\!\left(\frac{N_i}{q_i}\right),
}
```

where $R = N_{\mathrm A}k_{\mathrm B}$. Equation {eq}`eq:mu-from-q` is the microscopic counterpart of the thermodynamic expression $\mu_i = \mu_i^\circ(T) + RT\ln(P_i/P^\circ)$ from § 7.1.3: both express the chemical potential of a gas-phase species, one in terms of partial pressure and a reference $\mu_i^\circ(T)$, the other in terms of microscopic counting through $q_i/N_i$.

---

### 7.2.4 Equilibrium constants from partition functions

Inserting Eq. {eq}`eq:mu-from-q` into the equilibrium condition $\sum_i \nu_i\mu_i=0$:

```{math}
\sum_i \nu_i\ln\!\left(\frac{q_i}{N_i}\right)=0
\qquad\Longrightarrow\qquad
\prod_i\left(\frac{q_i}{N_i}\right)^{\nu_i}=1.
```

Equivalently,

```{math}
:label: eq:K-partition-numbers
\boxed{
\prod_i N_i^{\nu_i} = \prod_i q_i^{\nu_i}.
}
```

#### Relating to $K_c$ and $K_p$

Define a number concentration $c_i = N_i/V$. Dividing Eq. {eq}`eq:K-partition-numbers` by $V^{\Delta\nu}$ (with $\Delta\nu \equiv \sum_i \nu_i$) gives

```{math}
K_c \equiv \prod_i c_i^{\nu_i} = \prod_i\!\left(\frac{q_i}{V}\right)^{\nu_i}.
```

For an ideal gas, $P_i = c_i\, k_{\mathrm B}T$, so the pressure-based equilibrium constant is

```{math}
:label: eq:Kp-Kc
\boxed{
K_p \equiv \prod_i\left(\frac{P_i}{P^\circ}\right)^{\nu_i} = K_c\,\left(\frac{k_{\mathrm B}T}{P^\circ}\right)^{\Delta\nu}.
}
```

The prefactor $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ converts between "molecules per unit volume" and "pressure ratios relative to $P^\circ$." When $\Delta\nu = 0$, this prefactor is unity and $K_p = K_c$ is dimensionless directly.

So, **once you can compute** the single-molecule partition functions $q_i$, you can compute $K_c$ (and then $K_p$).

---

### 7.2.5 Example: $\mathrm{HI(g)}$ formation

Consider

```{math}
\mathrm{H_2(g) + I_2(g) \rightleftharpoons 2HI(g)}.
```

Here,

```{math}
\Delta\nu = 2 - 1 - 1 = 0,
```

so $K_p = K_c$ and the equilibrium constant is

```{math}
:label: eq:K-HI-ratio
\boxed{
K = \frac{q_{HI}^2}{q_{H_2}\,q_{I_2}}.
}
```

The $\Delta\nu = 0$ choice is pedagogically convenient: it makes volume factors cancel exactly (see the Worked Example below), so $K$ depends only on intrinsic molecular properties (masses, rotational and vibrational constants, bond strengths). With $\Delta\nu \neq 0$ reactions, one would need the $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ factor from Eq. {eq}`eq:Kp-Kc` to convert between $K_c$ and $K_p$.

---

### 7.2.6 Single-molecule partition function for a diatomic ideal gas

Within the rigid-rotor/harmonic-oscillator approximations developed in Chapter 2 (see § 2.8 for the summary table), the single-molecule partition function for a diatomic gas factors as

```{math}
q = q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}.
```

```{admonition} Energy-zero convention
:class: note

A partition function depends on *where* one places the zero of energy. In this section we choose the zero of molecular energy to lie at the **dissociated atoms at rest** — i.e., the ground vibrational level of a bound diatomic sits at energy $-D_0$, where $D_0$ is the bond dissociation energy measured from $v=0$. This choice puts the bond energy directly in the electronic factor as $e^{+D_0/k_{\mathrm B}T}$ (per molecule) or $e^{+D_0/RT}$ (per mole), and keeps the vibrational factor in its "ground-state-at-zero" form $e^{-\Theta_{\mathrm{vib}}/(2T)}/(1-e^{-\Theta_{\mathrm{vib}}/T})$. Some texts instead absorb the $v=0$ zero-point energy into the electronic factor and use $q_{\mathrm{vib}} = 1/(1-e^{-\Theta_{\mathrm{vib}}/T})$; both conventions give the *same* $K$ provided they are used consistently. We use the first convention throughout.
```

Putting the four factors together (cf. § 2.8):

```{math}
:label: eq:q-diatomic
\boxed{
q_{\text{diatomic}}
=
\underbrace{\frac{V}{\Lambda^3}}_{q_{\mathrm{trans}}}
\;\underbrace{\frac{T}{\sigma\,\Theta_{\mathrm{rot}}}}_{q_{\mathrm{rot}}}
\;\underbrace{\frac{e^{-\Theta_{\mathrm{vib}}/(2T)}}{1-e^{-\Theta_{\mathrm{vib}}/T}}}_{q_{\mathrm{vib}}}
\;\underbrace{g_1\,e^{D_0/(RT)}}_{q_{\mathrm{elec}}}.
}
```

where the thermal de Broglie wavelength is

```{math}
\Lambda \equiv \sqrt{\frac{h^2}{2\pi m k_{\mathrm B}T}}.
```

#### Notes on parameters

- $m$: molecular mass.
- $\sigma$: rotational symmetry number — $\sigma=2$ for homonuclear diatomics (H$_2$, I$_2$), $\sigma=1$ for heteronuclear (HI). See § 2.8.
- $\Theta_{\mathrm{rot}} \equiv \hbar^2/(2 I k_{\mathrm B})$: characteristic rotational temperature (from § 2.7).
- $\Theta_{\mathrm{vib}} \equiv h\nu/k_{\mathrm B}$: characteristic vibrational temperature (from § 2.6).
- $g_1$: ground-state electronic degeneracy. For the species in the HI reaction, all three ground states are $^1\Sigma^+$ singlets, so $g_1 = 1$ for each.
- $D_0$: bond dissociation energy referenced from the $v=0$ level (per-mole units when combined with $RT$, per-molecule when combined with $k_{\mathrm B}T$).

---

### 7.2.7 $K$ for $\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ in terms of molecular parameters

Substituting Eq. {eq}`eq:q-diatomic` into $K = q_{HI}^2/(q_{H_2}\,q_{I_2})$, and using the fact that volume cancels when $\Delta\nu=0$ (see Worked Example), yields

```{math}
:label: eq:K-HI-full
\boxed{
K
=
\left(\frac{m_{HI}^2}{m_{H_2}\,m_{I_2}}\right)^{3/2}
\;
\frac{\sigma_{H_2}\sigma_{I_2}}{\sigma_{HI}^2}\,
\frac{\Theta_{\mathrm{rot}}^{H_2}\,\Theta_{\mathrm{rot}}^{I_2}}{\left(\Theta_{\mathrm{rot}}^{HI}\right)^2}
\;
\frac{\bigl(1-e^{-\Theta_{\mathrm{vib}}^{H_2}/T}\bigr)\bigl(1-e^{-\Theta_{\mathrm{vib}}^{I_2}/T}\bigr)}{\bigl(1-e^{-\Theta_{\mathrm{vib}}^{HI}/T}\bigr)^2}
\;
\exp\!\left(\frac{2D_0^{HI}-D_0^{H_2}-D_0^{I_2}}{RT}\right).
}
```

With $\sigma_{H_2}=\sigma_{I_2}=2$ and $\sigma_{HI}=1$, the symmetry-number prefactor is $\sigma_{H_2}\sigma_{I_2}/\sigma_{HI}^2 = 4$; vibrational zero-point factors have cancelled in pairs under the chosen energy zero, leaving only the $(1-e^{-\Theta_{\mathrm{vib}}/T})$ "hot-mode" factors.

Equation {eq}`eq:K-HI-full` makes the microscopic physics very explicit:

- **Translation** contributes the $m^{3/2}$ mass dependence.
- **Rotation** contributes the symmetry-number and $\Theta_{\mathrm{rot}}$ factors.
- **Vibration** contributes the $(1-e^{-\Theta_{\mathrm{vib}}/T})$ factors that encode how many vibrational quanta are thermally accessible.
- **Bond strengths** enter through the dissociation energies $D_0$ in the exponential; the numerator $2D_0^{HI}-D_0^{H_2}-D_0^{I_2}$ is essentially the negative of the reaction internal energy at 0 K.

---

### 7.2.8 Numerical evaluation: $K(700\,\mathrm{K})$ from spectroscopic data

To see the microscopic-to-macroscopic bridge at work, we evaluate Eq. {eq}`eq:K-HI-full` at $T = 700\ \mathrm{K}$ using standard spectroscopic constants (e.g., from Herzberg or the NIST WebBook):

| species | $M$ (g/mol) | $\Theta_{\mathrm{rot}}$ (K) | $\Theta_{\mathrm{vib}}$ (K) | $D_0$ (kJ/mol) | $\sigma$ | $g_1$ |
| :-- | --: | --: | --: | --: | --: | --: |
| $\mathrm{H_2}$ | 2.016 | 87.6 | 6332 | 432.07 | 2 | 1 |
| $\mathrm{I_2}$ | 253.808 | 0.0537 | 308 | 148.81 | 2 | 1 |
| $\mathrm{HI}$ | 127.912 | 9.25 | 3266 | 294.67 | 1 | 1 |

```{code-cell} ipython3
import numpy as np

R = 8.314462618      # J mol^-1 K^-1

# Spectroscopic / thermochemical constants
mol = {
    "H2": dict(M=2.016e-3,   Theta_rot=87.6,   Theta_vib=6332.0, D0=432.07e3, sigma=2),
    "I2": dict(M=253.808e-3, Theta_rot=0.0537, Theta_vib=308.0,  D0=148.81e3, sigma=2),
    "HI": dict(M=127.912e-3, Theta_rot=9.25,   Theta_vib=3266.0, D0=294.67e3, sigma=1),
}
T = 700.0  # K

# Mass factor (translational)
mass_factor = (mol["HI"]["M"]**2 / (mol["H2"]["M"] * mol["I2"]["M"]))**1.5

# Rotational factor
rot_factor = (mol["H2"]["sigma"] * mol["I2"]["sigma"] / mol["HI"]["sigma"]**2) \
             * (mol["H2"]["Theta_rot"] * mol["I2"]["Theta_rot"]
                / mol["HI"]["Theta_rot"]**2)

# Vibrational factor: "hot-mode" (1 - exp(-Th_vib/T)) factors
def vfac(Tvib): return 1 - np.exp(-Tvib/T)
vib_factor = (vfac(mol["H2"]["Theta_vib"]) * vfac(mol["I2"]["Theta_vib"])
              / vfac(mol["HI"]["Theta_vib"])**2)

# Bond-energy factor: exp((2 D0_HI - D0_H2 - D0_I2)/RT)
dD0 = 2*mol["HI"]["D0"] - mol["H2"]["D0"] - mol["I2"]["D0"]  # J/mol
bond_factor = np.exp(dD0 / (R*T))

K = mass_factor * rot_factor * vib_factor * bond_factor

print(f"T = {T:.0f} K")
print(f"  Mass (translational) factor:       {mass_factor:10.3f}")
print(f"  Rotational factor:                 {rot_factor:10.4f}")
print(f"  Vibrational factor:                {vib_factor:10.4f}")
print(f"  Bond-energy factor exp(dD0/RT):    {bond_factor:10.4f}")
print(f"  2 D0(HI) - D0(H2) - D0(I2):        {dD0/1e3:+.2f} kJ/mol")
print()
print(f"  K(700 K) from molecular constants: {K:.1f}")
print(f"  Experimental K (Bodenstein, ~700 K): ~54")
```

**Discussion.** The predicted $K \approx 62$ compares to an experimental value near 54 at 700 K. Agreement at the ~15 % level is about what one should expect from the rigid-rotor/harmonic-oscillator approximation combined with uncertainties in $\Theta_{\mathrm{rot}}$, $\Theta_{\mathrm{vib}}$, and $D_0$ — corrections for anharmonicity, vibration–rotation coupling, and centrifugal distortion (none of which are in our model) would each shift $K$ by a few percent.

More importantly, the decomposition makes the physics transparent:

- The **translational mass factor** ($\approx 181$) is the largest single contribution. The ratio $m_{HI}^2/(m_{H_2}m_{I_2}) = (m_H+m_I)^2/(4 m_H m_I)$ is driven far from unity when the two atomic masses differ by a lot — for H and I, $m_I/m_H \approx 127$, so $m_{HI}^2 \gg 4\,m_H m_I$. If instead the atoms had similar masses (e.g., $\mathrm{D_2+X_2}\rightleftharpoons 2\mathrm{DX}$ with $m_D\sim m_X$), this factor would be close to 1.
- The **rotational factor** ($\approx 0.22$) is well below unity because I$_2$ is unusually heavy and has an extraordinarily small $\Theta_{\mathrm{rot}} \approx 0.05$ K; the corresponding $q_{\mathrm{rot}}(\mathrm{I_2})$ is very large, which favors keeping mass in I$_2$ rather than redistributing it into HI. The factor-of-4 symmetry-number prefactor (from $\sigma_{H_2}\sigma_{I_2}/\sigma_{HI}^2$) nudges equilibrium toward products but doesn't overcome this.
- The **vibrational factor** ($\approx 0.36$) is also below unity. HI has a stiff H–I stretch that is still nearly frozen out at 700 K (numerator small for HI), while I$_2$'s low-frequency mode is partially excited (numerator moderate for I$_2$), and these combine unfavorably for products.
- The **bond-energy factor** ($\approx 4.3$) contains the "chemistry": $2D_0^{HI} - D_0^{H_2} - D_0^{I_2} \approx +8.5\,\mathrm{kJ/mol}$ says the two HI bonds together are slightly stronger than one H–H plus one I–I bond, a modest exothermicity that at 700 K gives just over a factor of 4 in $K$.

The lesson: every digit of $K$ can be traced back to a specific molecular property, and the relative sizes of the four factors tell us which part of molecular structure matters most for this particular equilibrium. A student who computed $\Delta_r G^\circ$ from thermochemical tables (§ 7.1.6) and a student who computed $K$ from spectroscopic constants (§ 7.2.8) would arrive at the same number by very different routes — and each route teaches something different about *why* the equilibrium sits where it does.

---

## Worked Example

### Why volume cancels when $\Delta\nu=0$: $\mathrm{H_2+I_2\rightleftharpoons 2HI}$

Writing each single-molecule partition function as $q_i = q_{\mathrm{trans},i}\,q_{\mathrm{int},i}$ with $q_{\mathrm{trans},i}=V/\Lambda_i^3$ (and $q_{\mathrm{int},i}$ collecting the rotational, vibrational, and electronic factors, all of which are $V$-independent):

```{math}
K = \frac{q_{HI}^2}{q_{H_2}\,q_{I_2}}
= \frac{(V/\Lambda_{HI}^3)^2\,q_{\mathrm{int},HI}^2}{(V/\Lambda_{H_2}^3)(V/\Lambda_{I_2}^3)\,q_{\mathrm{int},H_2}\,q_{\mathrm{int},I_2}}.
```

Because $\Delta\nu = 2-1-1 = 0$, the volume factors cancel:

```{math}
\frac{V^2}{V\cdot V}=1,
```

leaving

```{math}
K=\left(\frac{\Lambda_{H_2}^3\Lambda_{I_2}^3}{\Lambda_{HI}^6}\right)\,
\frac{q_{\mathrm{int},HI}^2}{q_{\mathrm{int},H_2}\,q_{\mathrm{int},I_2}}.
```

Since $\Lambda\propto m^{-1/2}$, the translational contribution collapses to the mass factor $\left(m_{HI}^2/(m_{H_2}m_{I_2})\right)^{3/2}$ seen in Eq. {eq}`eq:K-HI-full`.

**Result.** For reactions with $\Delta\nu=0$, $K$ is independent of volume (and pressure) in the ideal-gas model; microscopic physics enters through masses and internal partition functions only. For reactions with $\Delta\nu\ne 0$, the uncancelled factors of $V$ reappear as the $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ conversion in Eq. {eq}`eq:Kp-Kc`.

## Concept Checks

1. Why is the Helmholtz free energy the relevant potential for equilibrium at constant $T,V$, and why does the resulting condition $\sum_i\nu_i\mu_i=0$ match the condition derived in § 7.1 at constant $T,P$?
2. Where does Stirling's approximation enter the derivation of $\mu_i = -RT\ln(q_i/N_i)$, and what physical limit makes it reasonable?
3. In Eq. {eq}`eq:q-diatomic`, how would the vibrational factor change if we chose the energy zero at the bottom of the potential well rather than at $v=0$? Show that $K$ is unchanged by this rechoice (as it must be).
4. Vibrational contributions to $K$ often matter most at high temperatures. Why? Reason from the temperature dependence of $(1-e^{-\Theta_{\mathrm{vib}}/T})$.
5. Suppose all three species in a $\Delta\nu=0$ reaction had identical masses, moments of inertia, vibrational frequencies, and $g_1$ values, but different $D_0$ values. What would $K$ reduce to, and what does this limit teach you about when "bond counting" arguments are quantitatively trustworthy?

## Key Takeaways

- Microscopic partition functions provide a route to equilibrium constants by linking $A = -k_{\mathrm B}T\ln Q$ (§ 4.3) to chemical potentials via $\mu_i = (\partial A/\partial N_i)_{T,V}$.
- For ideal-gas mixtures, $\mu_i = -RT\ln(q_i/N_i)$ (per-mole), and the equilibrium condition $\sum_i\nu_i\mu_i=0$ delivers $K$ as a ratio of single-molecule partition functions.
- Translational, rotational, vibrational, and electronic factors contribute in separable ways under the common rigid-rotor/harmonic-oscillator approximations of Chapter 2.
- For $\Delta\nu=0$ gas reactions, volume (and hence pressure) dependence cancels in the ideal-gas equilibrium constant; for $\Delta\nu\ne 0$, the $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ factor converts between $K_c$ and $K_p$.
- For the $\mathrm{H_2+I_2\rightleftharpoons 2HI}$ benchmark, the four factors in Eq. {eq}`eq:K-HI-full` predict $K(700\,\mathrm{K})\approx 62$, in ~15 % agreement with the experimental value of ~54 — a concrete illustration of how molecular structure alone determines the position of equilibrium.
