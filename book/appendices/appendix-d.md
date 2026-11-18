---
numbering:
  enumerator: D.%s
---

# Appendix D. Data Sources

The measured and tabulated numbers in these notes come from the sources below. For each source, this appendix states what it tabulates, where the notes use it, which conventions to check before using a value, and how to cite it. NIST is the National Institute of Standards and Technology, and JANAF stands for Joint Army-Navy-Air Force.

## Which Source to Use

```{list-table} Sources for the quantities used in these notes
:header-rows: 1
:name: data-sources

* - Quantity
  - Source
* - $\Delta_f H^\circ$, $S^\circ$, and $C_P^\circ$ of a species at $298.15\ \mathrm{K}$
  - [NIST Chemistry WebBook](#section:webbook) or [NIST-JANAF Thermochemical Tables](#section:janaf)
* - $H^\circ(T) - H^\circ(298.15\ \mathrm{K})$, $S^\circ(T)$, and $\Delta_f G^\circ(T)$ at other temperatures
  - [NIST-JANAF Thermochemical Tables](#section:janaf)
* - $\Delta_f H^\circ$ with a stated uncertainty
  - [Active Thermochemical Tables](#section:atct)
* - Vibrational wavenumbers, rotational constants, and bond lengths of molecules
  - [NIST Chemistry WebBook](#section:webbook) (diatomic molecules) or [CCCBDB](#section:cccbdb) (polyatomic molecules)
* - Energy levels of atoms
  - [NIST Atomic Spectra Database](#section:asd)
* - Density, saturation pressure, and enthalpy of vaporization of a fluid
  - [NIST Chemistry WebBook](#section:webbook), fluid properties
* - Physical constants
  - [CODATA recommended values](#section:constants)
```

(section:webbook)=
## NIST Chemistry WebBook

**Address.** [https://webbook.nist.gov/chemistry/](https://webbook.nist.gov/chemistry/); the fluid properties are at [https://webbook.nist.gov/chemistry/fluid/](https://webbook.nist.gov/chemistry/fluid/).

**What it tabulates.** For each species, the WebBook collects thermochemical data ($\Delta_f H^\circ$, $S^\circ$, and the coefficients of the Shomate equation for $C_P^\circ(T)$), phase-change data (boiling and triple points, enthalpies of vaporization and sublimation, vapor-pressure equations), and the spectroscopic constants of diatomic molecules. The fluid-properties pages give $P$, $\bar v$, $\bar h$, $\bar s$, and the heat capacities of more than 70 fluids along isotherms, isobars, and the saturation curve.

**Where the notes use it.** The isotherms and isobars of Sections [1.3](../chapter-01/section-03.md) and [1.4](../chapter-01/section-04.md), the molecular constants of Sections [2.6](../chapter-02/section-06.md) to [2.8](../chapter-02/section-08.md) and [7.2](../chapter-07/section-02.md), the heat capacities of [Section 5.3](../chapter-05/section-03.md), and the saturation properties of Sections [6.1](../chapter-06/section-01.md) and [6.2](../chapter-06/section-02.md).

**Check before using a value.**

- The phase. The gas and the condensed phases of a species have separate tables, and $\Delta_f H^\circ$ and $S^\circ$ differ between them.
- The range of a Shomate equation. Each set of coefficients holds between two stated temperatures, and the equation is written in the variable $t = T/(1000\ \mathrm{K})$.
- The reference state of the fluid properties. On the fluid pages, $\bar h$ and $\bar s$ are measured from a reference state that is chosen separately for each fluid, and the page offers four conventions for it. Only a difference between two states of the same fluid has a meaning. These values cannot be combined with enthalpies of formation or with third-law entropies, and they cannot be compared between two fluids.
- The origin of the fluid properties. They are calculated from a reference equation of state fitted to many measurements (for water, the IAPWS-95 formulation of @wagner2002), so they vary smoothly and carry no scatter.

**How to cite it.** Cite the database [@nist_webbook], and for fluid properties the chapter that supplies them [@lemmon_fluids]. Give the date on which the value was retrieved.

(section:janaf)=
## NIST-JANAF Thermochemical Tables

**Address.** [https://janaf.nist.gov/](https://janaf.nist.gov/)

**What it tabulates.** One table for each species and phase, with rows at $0$, $100$, $200$, $250$, $298.15$, and $300\ \mathrm{K}$ and at regular intervals above. The columns are $C_P^\circ$, $S^\circ$, $-[G^\circ - H^\circ(T_r)]/T$, $H^\circ - H^\circ(T_r)$, $\Delta_f H^\circ$, $\Delta_f G^\circ$, and $\log K_f$. The online tables reproduce the fourth edition of the printed tables [@chase1998].

**Where the notes use it.** The reaction enthalpies of [Section 3.3](../chapter-03/section-03.md), the free energies and absolute entropies of Sections [5.1](../chapter-05/section-01.md) and [5.2](../chapter-05/section-02.md), the heat capacities of water in [Section 6.1](../chapter-06/section-01.md), and the equilibrium constants of Sections [7.1](../chapter-07/section-01.md) and [7.2](../chapter-07/section-02.md).

**Check before using a value.**

- The standard-state pressure is $P^\circ = 0.1\ \mathrm{MPa} = 1\ \mathrm{bar}$, and the reference temperature is $T_r = 298.15\ \mathrm{K}$. Both are printed at the head of each table.
- The units differ between columns: $C_P^\circ$, $S^\circ$, and $-[G^\circ - H^\circ(T_r)]/T$ are in $\mathrm{J\,K^{-1}\,mol^{-1}}$, and the enthalpies and Gibbs energies are in $\mathrm{kJ\,mol^{-1}}$.
- $\log K_f$ is the base-10 logarithm of the equilibrium constant of the formation reaction of the species from its elements.
- $\Delta_f H^\circ$ and $\Delta_f G^\circ$ of an element in its reference state are zero at every temperature.

**How to cite it.** Cite the online database [@nist_janaf] or the printed tables [@chase1998].

(section:atct)=
## Active Thermochemical Tables

**Address.** [https://atct.anl.gov/](https://atct.anl.gov/)

**What it tabulates.** Enthalpies of formation at $0\ \mathrm{K}$ and $298.15\ \mathrm{K}$, each with an uncertainty. The values are obtained by solving one network of measured and calculated reaction enthalpies for all species at once, and not species by species [@ruscic2004].

**Where the notes use it.** Sections [3.3](../chapter-03/section-03.md), [5.3](../chapter-05/section-03.md), and [7.1](../chapter-07/section-01.md) name it as a source of $\Delta_f H^\circ$.

**Check before using a value.**

- The version. The values are revised whenever the network grows, so the same species can have slightly different values in two versions.
- The quantity. The tables give enthalpies of formation. Take entropies and heat capacities from the NIST-JANAF tables or the NIST Chemistry WebBook.

**How to cite it.** Cite the version of the Thermochemical Network printed on the page that gave the value, together with @ruscic2004.

(section:cccbdb)=
## NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB)

**Address.** [https://cccbdb.nist.gov/](https://cccbdb.nist.gov/)

**What it tabulates.** Measured and calculated properties of gas-phase atoms and small molecules: geometries, rotational constants, vibrational wavenumbers, and enthalpies of formation.

**Where the notes use it.** The rotational constants and vibrational wavenumbers of polyatomic molecules in [Section 2.8](../chapter-02/section-08.md) and in the [Module 2 review](../chapter-02/review.md).

**Check before using a value.**

- Measured or calculated. Each property has pages of experimental values and pages of values calculated with many quantum-chemical methods. The notes use the experimental values.
- The units. Rotational constants and vibrational frequencies are given as wavenumbers in $\mathrm{cm^{-1}}$. Multiply by $hc/k_{\mathrm{B}} = 1.43878\ \mathrm{cm\,K}$ to obtain the characteristic temperature, as in [Section 2.7](../chapter-02/section-07.md).

**How to cite it.** Cite the database with its release number [@johnson_cccbdb].

(section:asd)=
## NIST Atomic Spectra Database

**Address.** [https://www.nist.gov/pml/atomic-spectra-database](https://www.nist.gov/pml/atomic-spectra-database)

**What it tabulates.** The energy levels of atoms and atomic ions, each with its configuration, term, and total angular momentum quantum number $J$, and the wavelengths and probabilities of the transitions between them.

**Where the notes use it.** The electronic partition functions of [Section 2.8](../chapter-02/section-08.md) and the practice problems of the [Module 2](../chapter-02/review.md) and [Module 4](../chapter-04/review.md) reviews.

**Check before using a value.**

- The units. Level energies are wavenumbers in $\mathrm{cm^{-1}}$ above the lowest level.
- The degeneracy. A level with quantum number $J$ has the degeneracy $g = 2J + 1$.

**How to cite it.** Cite the database with its version number [@kramida_asd].

(section:constants)=
## Physical Constants and Other Sources

The links at the physical constants in the notes lead to the NIST pages for the CODATA recommended values [@mohr2025]. Since the 2019 revision of the SI, $h$, $k_{\mathrm{B}}$, $N_{\mathrm{A}}$, and the elementary charge have exact defined values, so $R = N_{\mathrm{A}} k_{\mathrm{B}}$ is exact as well.

The phase diagram of water in [Section 6.1](../chapter-06/section-01.md) follows the equations of the International Association for the Properties of Water and Steam (IAPWS), whose formulation for the fluid phases is IAPWS-95 [@wagner2002]. The vapor pressures of [](#fig:cc-studio) in [Section 6.2](../chapter-06/section-02.md) were calculated with the CoolProp library [@bell2014].
