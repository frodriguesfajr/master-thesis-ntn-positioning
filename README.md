# NTN Positioning Performance Based on the Cramér–Rao Bound

MATLAB simulation framework developed for the master's dissertation:

**“Análise do Desempenho de Posicionamento em Redes Não Terrestres Baseada no Limite de Cramér–Rao”**

This repository contains the simulation code used to investigate positioning performance in heterogeneous non-terrestrial networks (NTNs) composed of **HAPS, LEO, MEO, and GEO transmitters**.

The analysis considers two main factors:

* transmitter geometry, represented by the **Position Dilution of Precision (PDOP)**;
* statistical quality of pseudorange measurements, represented by the measurement covariance matrix (\mathbf{R}).

These effects are jointly incorporated through the Fisher Information Matrix

[
\mathbf{J} =
\mathbf{H}^{T}\mathbf{R}^{-1}\mathbf{H},
]

from which the positional Cramér–Rao lower bound is obtained as

[
\mathcal{B}_{\mathrm{CRB}}
==========================

\sqrt{
\operatorname{tr}
\left(
[\mathbf{J}^{-1}]_{1:3,1:3}
\right)
}.
]

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

Generated results include:

* `.csv` tables;
* `.png` figures;
* `.mat` files containing the simulation variables and selected geometries.

---

## Simulation Scenarios

### Scenario 1 — Isolated NTN Architectures

`scenario1.m`

HAPS, LEO, MEO, and GEO architectures are evaluated separately.

The scenario analyzes:

* PDOP;
* positional Cramér–Rao bound;
* WLS positioning error;
* Monte Carlo RMSE;
* sensitivity to the reference (C/N_0);
* empirical error distributions and percentiles.

For each architecture, a controlled synthetic geometry is generated within predefined elevation, propagation-distance, and PDOP ranges.

The effective carrier-to-noise-density ratio is modeled as

[
(C/N_0)_{\mathrm{eff}}
======================

(C/N_0)_{\mathrm{ref}}
+
\Delta(C/N_0),
]

where architecture-dependent offsets are used to represent different measurement-quality conditions.

The nominal pseudorange standard deviation is modeled as

[
\sigma_{\rho}
=============

\frac{c}
{2\pi\beta
\sqrt{(C/N_0)*{\mathrm{linear}}T*{\mathrm{coh}}}}.
]

Monte Carlo simulations are used to compare the WLS positioning error with the theoretical positional bound.

---

### Scenario 2 — Greedy Transmitter Selection

`scenario2.m`

A heterogeneous candidate set composed of

[
4\text{ HAPS}
+
8\text{ LEO}
+
8\text{ MEO}
+
6\text{ GEO}
============

26
]

transmitters is considered.

The objective is to construct a low-cardinality subset satisfying

[
\mathcal{B}_{\mathrm{CRB}}
\leq
\epsilon,
]

with

[
\epsilon = 1,\mathrm{m}.
]

At least four transmitters are required because the positioning state contains three spatial coordinates and one receiver-clock term.

The initialization performs an exhaustive search over all

[
\binom{26}{4}=14,950
]

four-transmitter combinations and selects the subset with the lowest positional bound.

After initialization, a greedy procedure adds one transmitter at a time. At each iteration, the candidate producing the lowest value of (\mathcal{B}_{\mathrm{CRB}}) is selected.

The method significantly reduces the search complexity compared with evaluating every possible subset.

The greedy procedure does **not** guarantee the globally optimal subset. It provides a computationally efficient heuristic for constructing configurations satisfying the positioning requirement.

---

### Scenario 3 — Sensitivity to Measurement Quality

`scenario3.m`

The candidate geometry is kept fixed while the reference (C/N_0) is varied from

[
30\text{ to }50\ \mathrm{dB!-!Hz}.
]

For each value, the pseudorange uncertainties are recalculated and the greedy transmitter-selection procedure is executed.

Because the geometry is fixed, this experiment isolates the effect of measurement quality on the number of transmitters required to satisfy

[
\mathcal{B}_{\mathrm{CRB}}\leq1,\mathrm{m}.
]

The relative (C/N_0) offsets among HAPS, LEO, MEO, and GEO remain fixed throughout the experiment.

Consequently, the experiment primarily investigates how improving the overall measurement quality allows the positioning requirement to be reached with fewer transmitters.

---

### Scenario 3a / Scenario 4 — HAPS Availability and LEO Compensation

`scenario3a.m`

This experiment evaluates the sensitivity of the selected configuration to transmitter availability.

#### Part A — HAPS Restriction

The maximum number of HAPS allowed in the selected configuration is varied as

[
N_{\mathrm{HAPS,max}}
\in
{4,2,0}.
]

The remaining candidate availability is kept at:

[
8\text{ LEO},
\quad
8\text{ MEO},
\quad
6\text{ GEO}.
]

The experiment evaluates how limiting the use of HAPS changes:

* the number of selected transmitters;
* the architecture composition;
* PDOP;
* positional Cramér–Rao bound;
* propagation distance;
* one-way propagation delay.

The HAPS restriction is implemented as a **maximum allowed number of HAPS in the selected subset**.

#### Part B — LEO Compensation

HAPS are excluded from the selected configuration and the number of available LEO candidates is progressively increased:

[
N_{\mathrm{LEO}}
\in
{8,12,16,20}.
]

The purpose is to investigate whether increasing LEO availability can compensate for the absence of HAPS and recover the positioning requirement.

The candidate sets are nested: the configurations with additional LEO transmitters extend the previously available candidate set rather than generating an independent geometry for each case.

---

## Greedy Selection Strategy

The greedy approach used in this repository is motivated by the need to avoid exhaustive evaluation of a combinatorial number of transmitter subsets.

The strategy is conceptually related to greedy antenna-selection approaches found in the signal-processing literature, including:

**M. O. K. Mendonça, P. S. R. Diniz, T. N. Ferreira, and L. Lovisolo,
“Antenna Selection in Massive MIMO Based on Greedy Algorithms,”
IEEE Transactions on Wireless Communications, 2020.**

The algorithm implemented here is **not a direct implementation of Matching Pursuit**.

In Matching Pursuit, candidates are selected according to their correlation with an approximation residue.

In this work, the candidate-selection criterion is instead the positional Cramér–Rao bound:

[
m^\star
=======

\underset{m\in\mathcal C\setminus\mathcal S}
{\arg\min}
;
\mathcal B_{\mathrm{CRB}}
(\mathcal S\cup{m}).
]

Thus, the common element is the **greedy construction of a subset**, while the optimization criterion is specifically adapted to the positioning problem.

---

## Reproducibility

The transmitter geometries are generated using MATLAB pseudorandom-number generators.

For reproducible experiments, a fixed random seed should be used before geometry generation, for example:

```matlab
rng(2,'twister');
```

Using a fixed seed ensures that the same synthetic geometry is generated in repeated executions.

The results in Scenarios 2, 3, and 3a should therefore be interpreted as analyses of a **controlled and reproducible realization of the candidate geometry**, rather than as statistical characterizations of all possible NTN geometries.

A Monte Carlo analysis over multiple independent geometry realizations is a natural extension of this work.

Such an analysis would allow statistical evaluation of quantities such as:

* probability of satisfying the positioning requirement;
* distribution of the number of selected transmitters;
* frequency with which each architecture is selected;
* robustness of the greedy solution to geometric variations.

For a large number of geometry realizations, the exhaustive four-transmitter initialization becomes computationally demanding and would benefit from optimized and parallel implementations.

---

## Important Modeling Assumptions

The simulations use **synthetic controlled transmitter geometries**.

The HAPS, LEO, MEO, and GEO labels define scenario characteristics such as:

* elevation masks;
* propagation-distance ranges;
* target PDOP ranges;
* relative (C/N_0) conditions.

The scripts do **not** propagate operational satellite constellations or model complete orbital dynamics.

Propagation distance does not directly determine (C/N_0) through a complete link-budget model. Instead, measurement quality is modeled separately using architecture-dependent (C/N_0) offsets.

The calculated propagation delay is

[
\tau=\frac{d}{c},
]

and represents only the **one-way geometric propagation delay**. It does not include processing, routing, queueing, retransmission, or other network latency components.

---

## Running the Simulations

Open MATLAB in the repository directory and execute the desired script:

```matlab
scenario1
```

or

```matlab
scenario2
```

```matlab
scenario3
```

```matlab
scenario3a
```

Each script automatically creates its corresponding output directory and saves tables, figures, and MATLAB workspace variables.

For long Monte Carlo experiments, verify the number of realizations configured in the script before execution.

---

## Outputs

Typical output files include:

```text
*.csv    Simulation tables and algorithm histories
*.png    Generated figures
*.mat    Saved geometries, selections, and simulation variables
```

The `.mat` files are particularly useful for preserving the exact geometry realization associated with a given simulation result.

---

## Scope

The purpose of this repository is to provide a reproducible research implementation for studying the interaction between:

[
\boxed{\text{transmitter geometry}}
]

and

[
\boxed{\text{pseudorange measurement quality}}
]

in heterogeneous non-terrestrial positioning networks.

The positional Cramér–Rao bound provides the common theoretical metric used to compare configurations and guide transmitter selection.
