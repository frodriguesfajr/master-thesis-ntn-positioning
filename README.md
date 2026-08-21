# NTN Positioning Performance Based on the Cramér–Rao Bound

MATLAB simulation framework developed for the master's dissertation:

**“Análise do Desempenho de Posicionamento em Redes Não Terrestres Baseada no Limite de Cramér–Rao”**

This repository contains the simulation code used to investigate positioning performance in heterogeneous non-terrestrial networks (NTNs) composed of **HAPS, LEO, MEO, and GEO transmitters**.

The analysis considers two main factors:

* transmitter geometry, represented by the **Position Dilution of Precision (PDOP)**;
* statistical quality of pseudorange measurements, represented by the measurement covariance matrix $\mathbf{R}$.

These effects are jointly incorporated through the Fisher Information Matrix:

$$\mathbf{J}=\mathbf{H}^{T}\mathbf{R}^{-1}\mathbf{H}$$

The positional Cramér–Rao bound is obtained from:

$$
\mathcal{B}_{\mathrm{CRB}} =
\sqrt{\mathrm{tr}\left([\mathbf{J}^{-1}]_{1:3,1:3}\right)}
$$

---

## Repository Structure

```text
.
├── scenario1.m
├── scenario2.m
├── scenario3.m
├── scenario3a.m
│
├── results_scenario_01/
├── results_scenario_02/
├── results_scenario_03/
└── results_scenario_04/
```

Each MATLAB script is self-contained and includes the local auxiliary functions required for the corresponding simulation.

Generated results may include:

* `.csv` tables;
* `.png` figures;
* `.mat` files containing simulation variables and selected geometries.

---

# Simulation Scenarios

## Scenario 1 — Isolated NTN Architectures

File:

```text
scenario1.m
```

HAPS, LEO, MEO, and GEO architectures are evaluated separately.

The scenario analyzes:

* PDOP;
* positional Cramér–Rao bound;
* WLS positioning error;
* Monte Carlo RMSE;
* sensitivity to the reference $C/N_0$;
* empirical error distributions;
* positioning-error percentiles.

For each architecture, a controlled synthetic geometry is generated within predefined elevation, propagation-distance, and PDOP ranges.

The effective carrier-to-noise-density ratio is modeled as:

$$(C/N_0)_{\mathrm{eff}}=(C/N_0)_{\mathrm{ref}}+\Delta(C/N_0)$$

where architecture-dependent offsets are used to represent different measurement-quality conditions.

The nominal pseudorange standard deviation is modeled as:

$$\sigma_{\rho}=\frac{c}{2\pi\beta\sqrt{(C/N_0)*{\mathrm{linear}}T*{\mathrm{coh}}}}$$

Monte Carlo simulations are used to compare the WLS positioning error with the theoretical positional bound.

The main purpose of this scenario is to show that positioning performance depends not only on transmitter geometry but also on the statistical quality of the pseudorange measurements.

---

## Scenario 2 — Greedy Transmitter Selection

File:

```text
scenario2.m
```

A heterogeneous candidate set composed of:

$$4,\mathrm{HAPS}+8,\mathrm{LEO}+8,\mathrm{MEO}+6,\mathrm{GEO}=26$$

transmitters is considered.

The objective is to construct a low-cardinality subset satisfying:

$$\mathcal{B}_{\mathrm{CRB}}\leq\epsilon$$

with:

$$\epsilon = 1 \mathrm{m}$$

At least four transmitters are required because the positioning state contains three spatial coordinates and one receiver-clock term.

### Initial subset

The initialization performs an exhaustive search over every four-transmitter combination:

$$\binom{26}{4}=14,950$$

For each combination, the positional Cramér–Rao bound is evaluated.

The initial subset is therefore obtained from:

```math
\mathcal{S}_4
=
\arg\min_{\mathcal{A}\subseteq\mathcal{C},\;|\mathcal{A}|=4}
\mathcal{B}_{\mathrm{CRB}}(\mathcal{A})
```


### Greedy expansion

After initialization, one transmitter is added at each iteration.

The selected transmitter is the candidate that produces the lowest positional bound when added to the current subset:

$$ m^\star=
\underset{
m\in\mathcal{C}\setminus\mathcal{S}
}{
\arg\min
}
;
\mathcal{B}_{\mathrm{CRB}}
\left(
\mathcal{S}\cup{m}
\right)
$$

The procedure continues until:

$$
\mathcal{B}_{\mathrm{CRB}}
\leq
\epsilon
$$

or until no additional candidates remain.

The greedy procedure substantially reduces computational cost compared with evaluating every possible transmitter subset.

However, the method does **not** guarantee the globally optimal subset. It should be interpreted as a computationally efficient heuristic for constructing a configuration that satisfies the positioning requirement.

---

## Scenario 3 — Sensitivity to Measurement Quality

File:

```text
scenario3.m
```

The candidate geometry is kept fixed while the reference $C/N_0$ is varied from:

$$
30
\leq
(C/N_0)_{\mathrm{ref}}
\leq
50
\quad
\mathrm{dB!-!Hz}
$$

For each value of reference $C/N_0$:

1. the pseudorange uncertainties are recalculated;
2. the covariance matrix $\mathbf{R}$ changes;
3. the greedy transmitter-selection procedure is executed;
4. the final number and composition of selected transmitters are recorded.

Because the geometry is kept fixed throughout the experiment, this scenario isolates the effect of measurement quality on transmitter selection.

The relative offsets among the architectures remain constant:

$$
\Delta(C/N_0)_{\mathrm{HAPS}} = +8,\mathrm{dB}
$$

$$
\Delta(C/N_0)_{\mathrm{LEO}} = +4,\mathrm{dB}
$$

$$
\Delta(C/N_0)_{\mathrm{MEO}} = 0,\mathrm{dB}
$$

$$
\Delta(C/N_0)_{\mathrm{GEO}} = -4,\mathrm{dB}
$$

Therefore, improving the reference $C/N_0$ improves the overall measurement quality while preserving the relative measurement-quality differences between architectures.

The main purpose of this experiment is to evaluate how improved measurement conditions reduce the number of transmitters required to satisfy:

$$
\mathcal{B}_{\mathrm{CRB}}
\leq
1,\mathrm{m}
$$

---

## Scenario 3a — HAPS Availability and LEO Compensation

File:

```text
scenario3a.m
```

This experiment investigates how transmitter availability affects the configuration selected by the greedy algorithm.

The analysis is divided into two parts.

---

### Part A — HAPS Restriction

The maximum number of HAPS allowed in the selected subset is varied as:

$$
N_{\mathrm{HAPS,max}}
\in
{4,2,0}
$$

The remaining candidate availability is kept at:

$$
8,\mathrm{LEO}
+
8,\mathrm{MEO}
+
6,\mathrm{GEO}
$$

The imposed restriction is:

$$
N_{\mathrm{HAPS}}(\mathcal{S})
\leq
N_{\mathrm{HAPS,max}}
$$

This means that the algorithm can select at most the specified number of HAPS from the available HAPS candidates.

The experiment evaluates the effects of this restriction on:

* selected architecture composition;
* total number of selected transmitters;
* PDOP;
* positional Cramér–Rao bound;
* propagation distance;
* one-way propagation delay.

---

### Part B — LEO Compensation

HAPS are excluded from the selected configuration:

$$
N_{\mathrm{HAPS,max}}=0
$$

The number of available LEO candidates is progressively increased:

$$
N_{\mathrm{LEO}}
\in
{8,12,16,20}
$$

while MEO and GEO availability remains fixed.

The purpose is to investigate whether additional LEO candidates can compensate for the absence of HAPS and recover the positioning requirement:

$$
\mathcal{B}_{\mathrm{CRB}}
\leq
1,\mathrm{m}
$$

The LEO candidate sets are nested. Therefore:

$$
\mathcal{C}*{8}
\subset
\mathcal{C}*{12}
\subset
\mathcal{C}*{16}
\subset
\mathcal{C}*{20}
$$

This allows the experiment to represent a progressive increase in LEO availability rather than independent candidate sets for each case.

---

## Defense Presentation

The presentation used for the master's thesis defense is available in:

`docs/defense/master_thesis_defense.pdf`

The slides summarize the motivation, methodology, simulation scenarios,
greedy transmitter-selection strategy, main results, and conclusions.

---

# Greedy Selection Strategy

The greedy strategy used in this repository is motivated by the need to reduce the computational cost associated with exhaustive combinatorial search.

The approach is conceptually related to greedy antenna-selection techniques found in the signal-processing literature, including:

**M. O. K. Mendonça, P. S. R. Diniz, T. N. Ferreira, and L. Lovisolo,
“Antenna Selection in Massive MIMO Based on Greedy Algorithms,”
IEEE Transactions on Wireless Communications, 2020.**

The algorithm implemented in this repository is **not a direct implementation of Matching Pursuit**.

In Matching Pursuit, a target vector is approximated using a sparse combination of dictionary elements. At each iteration, the selected codeword is typically the one most correlated with the current approximation residue.

In the positioning problem considered here, the criterion is instead based directly on the positional Cramér–Rao bound.

The common element between the approaches is the **iterative greedy construction of a subset**, while the selection metric is adapted to the positioning problem.

---

# Positioning Model

The positioning state contains three spatial coordinates and one receiver-clock term:

$$\boldsymbol{\theta}=\begin{bmatrix}x & y & z & b\end{bmatrix}^{T}$$

A pseudorange measurement associated with transmitter $i$ can be represented as:

$$\rho_i=|\mathbf{p}_i-\mathbf{p}|+b+\varepsilon_i$$

where:

* $\mathbf{p}$ is the receiver position;
* $\mathbf{p}_i$ is the position of transmitter $i$;
* $b$ is the receiver-clock term expressed in distance units;
* $\varepsilon_i$ represents the pseudorange measurement error.

After linearization, the geometry matrix is denoted by $\mathbf{H}$.

The pseudorange-error covariance matrix is:

$$\mathbf{R}=\rm{diag}\left(\sigma_{\rho_1}^{2},\sigma_{\rho_2}^{2},\ldots,\sigma_{\rho_M}^{2}\right)$$

The Fisher Information Matrix is:

$$\mathbf{J}=\mathbf{H}^{T}\mathbf{R}^{-1}\mathbf{H}$$

Therefore, the metric used by the greedy algorithm jointly depends on:

$$
\boxed{
\text{transmitter geometry}
}
$$

and:

$$
\boxed{
\text{pseudorange measurement quality}
}
$$

---

# PDOP

The Position Dilution of Precision is calculated from:

$$\mathbf{Q}=(\mathbf{H}^{T}\mathbf{H})^{-1}$$

and:

$$\mathrm{PDOP}=\sqrt{Q_{xx}+Q_{yy}+Q_{zz}}$$

or equivalently:

$$\mathrm{PDOP}=\sqrt{\rm{tr}\left([\mathbf{Q}]_{1:3,1:3}\right)}$$

PDOP characterizes the influence of geometry alone.

In contrast, the positional Cramér–Rao bound also incorporates the statistical quality of the measurements through $\mathbf{R}$.

---

# Reproducibility

The synthetic transmitter geometries are generated using MATLAB pseudorandom-number generators.

For reproducible experiments, a fixed random seed should be defined before geometry generation.

For example:

```matlab
rng(2,'twister');
```

Using a fixed seed ensures that the same pseudorandom sequence and, consequently, the same controlled geometry realization can be reproduced.

The results of the greedy-selection scenarios should therefore be interpreted as analyses of a **controlled and reproducible realization of the candidate geometry**.

They are not intended to represent a statistical characterization of every possible NTN geometry.

A Monte Carlo analysis over independent geometry realizations is a natural extension of this work.

Such an analysis could evaluate:

* probability of satisfying the positioning requirement;
* distribution of the number of selected transmitters;
* frequency of selection of each architecture;
* distribution of final PDOP;
* distribution of the positional Cramér–Rao bound;
* robustness of selected configurations to geometric variations.

For a large number of geometry realizations, the exhaustive initialization with four transmitters can become computationally demanding.

For example, with 26 candidates, each realization initially requires evaluation of:

$$\binom{26}{4}=14,950$$

four-transmitter subsets.

For $10,000$ independent geometry realizations, this initial stage alone would require approximately:

$$10,000\times14,950=149,500,000$$

subset evaluations.

Such an analysis would therefore benefit from optimized code and parallel computation.

---

# Important Modeling Assumptions

## Synthetic Geometries

The transmitter geometries used in the simulations are synthetic and controlled.

The architecture labels HAPS, LEO, MEO, and GEO define parameters such as:

* elevation masks;
* transmitter–receiver distance ranges;
* target PDOP ranges;
* relative $C/N_0$ offsets.

The scripts do **not** propagate operational satellite constellations or implement complete orbital dynamics.

Therefore, the simulations should be interpreted as controlled positioning scenarios rather than instantaneous representations of specific operational constellations.

---

## Measurement Quality

The pseudorange uncertainty is modeled from the effective $C/N_0$.

Propagation distance does not directly determine $C/N_0$ through a complete communication link-budget model.

Instead, architecture-dependent relative $C/N_0$ offsets are introduced to represent different measurement-quality conditions.

---

## Propagation Delay

The propagation delay reported in Scenario 3a is calculated as:

$$\tau=\frac{d}{c}$$

where:

* $d$ is the geometric transmitter–receiver distance;
* $c$ is the speed of light.

This quantity represents only the **one-way geometric propagation delay**.

It does not include:

* processing delay;
* queueing;
* routing;
* retransmissions;
* medium-access delay;
* other network latency components.

The propagation delay is currently used as an output metric and does not participate in the greedy transmitter-selection criterion.

---

# Running the Simulations

Open MATLAB in the repository directory.

To execute Scenario 1:

```matlab
scenario1
```

To execute Scenario 2:

```matlab
scenario2
```

To execute Scenario 3:

```matlab
scenario3
```

To execute Scenario 3a:

```matlab
scenario3a
```

Each script automatically creates its corresponding result directory.

---

# Outputs

Typical output files include:

```text
*.csv
*.png
*.mat
```

### CSV files

Contain:

* simulation parameters;
* selected configurations;
* PDOP values;
* positional Cramér–Rao bounds;
* greedy-algorithm histories;
* empirical error statistics.

### PNG files

Contain the figures generated by the simulation scripts.

### MAT files

Preserve simulation variables such as:

* generated geometries;
* transmitter positions;
* architecture labels;
* selected subsets;
* greedy-selection histories;
* numerical results.

The `.mat` files can therefore be useful for reproducing the exact geometry realization associated with a given result.

---

# Software

The simulations were developed in **MATLAB**.

No orbital propagation software is required for the synthetic geometry generation implemented in the current version of the repository.

---

# Scope and Future Extensions

The purpose of this repository is to provide a reproducible simulation framework for studying the interaction between:

$$
\boxed{
\text{transmitter geometry}
}
$$

and:

$$
\boxed{
\text{pseudorange measurement quality}
}
$$

in heterogeneous non-terrestrial positioning networks.

Possible extensions include:

* Monte Carlo analysis over geometry realizations;
* time-varying transmitter geometries;
* orbital propagation models;
* more detailed pseudorange-error models;
* complete propagation and link-budget models;
* availability constraints;
* latency constraints;
* computational-cost constraints;
* integration between NTN architectures and conventional GNSS;
* comparison of greedy selection with other optimization strategies.

---

# References

The theoretical background and simulation methodology used in this repository are based on established references in GNSS positioning, statistical estimation, non-terrestrial networks, satellite systems, and greedy subset-selection methods.

## GNSS and Positioning

1. **P. Misra and P. Enge**,
   *Global Positioning System: Signals, Measurements, and Performance*,
   2nd ed., Ganga-Jamuna Press, 2006.

2. **E. D. Kaplan and C. J. Hegarty**,
   *Understanding GPS/GNSS: Principles and Applications*,
   3rd ed., Artech House, 2017.

3. **P. D. Groves**,
   *Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems*,
   2nd ed., Artech House, 2013.

4. **J. B.-Y. Tsui**,
   *Fundamentals of Global Positioning System Receivers: A Software Approach*,
   2nd ed., John Wiley & Sons, 2005.
   DOI: `10.1002/0471712582`

5. **P. J. G. Teunissen and O. Montenbruck**, eds.,
   *Springer Handbook of Global Navigation Satellite Systems*,
   Springer, 2017.
   DOI: `10.1007/978-3-319-42928-1`

---

## Statistical Estimation and Cramér–Rao Bound

6. **S. M. Kay**,
   *Fundamentals of Statistical Signal Processing, Volume I: Estimation Theory*,
   Prentice Hall, 1993.

This reference provides the main theoretical foundation for estimation theory, the Fisher Information Matrix, and the Cramér–Rao lower bound used in this work.

7. **H. V. Poor**,
   *An Introduction to Signal Detection and Estimation*,
   2nd ed., Springer, 1994.

8. **H. L. Van Trees**,
   *Detection, Estimation, and Modulation Theory, Part I*,
   John Wiley & Sons, 2001.

---

## Non-Terrestrial Networks and Satellite Systems

9. **3rd Generation Partnership Project (3GPP)**,
   *Study on New Radio (NR) to Support Non-Terrestrial Networks*,
   Technical Report TR 38.811, Release 15.

10. **M. M. Azari et al.**,
    “Evolution of Non-Terrestrial Networks From 5G to 6G: A Survey,”
    *IEEE Communications Surveys & Tutorials*,
    vol. 24, no. 4, pp. 2633–2672, 2022.
    DOI: `10.1109/COMST.2022.3199901`

11. **O. Kodheli et al.**,
    “Satellite Communications in the New Space Era: A Survey and Future Challenges,”
    *IEEE Communications Surveys & Tutorials*,
    vol. 23, no. 1, pp. 70–109, 2021.
    DOI: `10.1109/COMST.2020.3028247`

12. **D. A. Vallado**,
    *Fundamentals of Astrodynamics and Applications*,
    4th ed., Microcosm Press, 2013.

13. **J. R. Wertz, D. F. Everett, and J. J. Puschell**, eds.,
    *Space Mission Engineering: The New SMAD*,
    Microcosm Press, 2011.

14. **G. Maral and M. Bousquet**,
    *Satellite Communications Systems: Systems, Techniques and Technology*,
    5th ed., John Wiley & Sons, 2009.

---

## Greedy Transmitter Selection

15. **M. O. K. Mendonça, P. S. R. Diniz, T. N. Ferreira, and L. Lovisolo**,
    “Antenna Selection in Massive MIMO Based on Greedy Algorithms,”
    *IEEE Transactions on Wireless Communications*,
    vol. 19, no. 3, pp. 1868–1881, 2020.
    DOI: `10.1109/TWC.2019.2959317`

This work provides methodological motivation for using greedy subset-selection strategies to reduce the computational complexity of combinatorial selection problems.

The algorithm implemented in this repository is not a direct implementation of the Matching Pursuit algorithms proposed by Mendonça et al. Instead, the greedy subset-construction principle is adapted to transmitter selection for positioning, using the positional Cramér–Rao bound as the selection criterion.

16. **M. O. K. de Mendonça**,
    *Greedy Algorithms and Machine Learning for Communications*,
    Ph.D. dissertation, COPPE/UFRJ, Electrical Engineering Program,
    Rio de Janeiro, Brazil, 2022.

---

## Additional Positioning and NTN References

17. **Y. Yang, W. Gao, S. Guo, Y. Mao, and Y. Yang**,
    “Introduction to BeiDou-3 Navigation Satellite System,”
    *NAVIGATION*, vol. 66, no. 1, pp. 7–18, 2019.
    DOI: `10.1002/navi.291`

18. **C. Pinell, F. S. Prol, M. Z. H. Bhuiyan, and J. Praks**,
    “Receiver Architectures for Positioning with Low Earth Orbit Satellite Signals: A Survey,”
    *EURASIP Journal on Advances in Signal Processing*,
    vol. 2023, article 60, 2023.
    DOI: `10.1186/s13634-023-01022-1`

---

For the complete bibliography and detailed theoretical discussion, please refer to the master's dissertation.

---

# Author

**Filipe Augusto Jesus Rodrigues**

Master's research developed at the Federal University of Rio de Janeiro (UFRJ/COPPE), Electrical Engineering Program.
