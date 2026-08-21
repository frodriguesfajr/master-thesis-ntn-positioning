# Reproducibility

This directory contains information required to reproduce the MATLAB simulations associated with the master's dissertation:

**“Análise do Desempenho de Posicionamento em Redes Não Terrestres Baseada no Limite de Cramér–Rao”**

The simulations investigate positioning performance using synthetic and controlled geometries composed of HAPS, LEO, MEO, and GEO transmitters.

The main scripts are:

```text
scenario1.m
scenario2.m
scenario3.m
scenario3a.m
```

---

# 1. Reproducibility Philosophy

Two different concepts of reproducibility are relevant for these simulations:

1. **Exact reproduction of a previously generated result**
2. **Generation of a new but reproducible random realization**

These two cases should not be confused.

---

## 1.1 Exact Reproduction

Some scenarios generate transmitter geometries using MATLAB pseudorandom-number functions such as:

```matlab
rand
randn
randperm
```

If a simulation was originally executed without a fixed random seed, defining a seed afterward does **not** reproduce the original realization.

For example, adding:

```matlab
rng(2,'twister');
```

to a script after the original simulation was performed does not guarantee reproduction of the geometry previously used in the dissertation.

For exact reproduction of an already generated result, the preferred approach is to use the corresponding saved `.mat` file containing the original geometry and simulation variables.

Examples include:

```text
results_scenario_01/scenario_01_results.mat
results_scenario_02/scenario_02_results.mat
results_scenario_03/scenario_03_results.mat
results_scenario_04/scenario_04_results.mat
```

These files may contain variables such as:

```text
SatPool
SatByArch
SatLEO_extra
archLabel
selected
allSelected
allHistories
```

and therefore preserve the geometry realization associated with the corresponding simulation output.

---

## 1.2 Reproducible New Realizations

For new simulations, a fixed random seed should be defined before any random geometry or measurement generation.

For example:

```matlab
rng(2,'twister');
```

This ensures that repeated executions using the same MATLAB random-number generator and the same code generate the same pseudorandom sequence.

A recommended structure is:

```matlab
close all;
clear;
clc;

rng(2,'twister');
```

The seed should be defined **before** any call to:

```matlab
rand
randn
randperm
```

---

# 2. MATLAB Environment

The simulations were developed in MATLAB.

The scripts use standard MATLAB functionality including:

```text
table
writetable
save
nchoosek
inputParser
exportgraphics
```

No external orbital propagation software is required for the current synthetic geometry model.

The scripts can be executed using either:

* MATLAB Desktop;
* MATLAB Online.

For reproducibility purposes, the MATLAB release used to generate final published results should ideally be recorded.

It is recommended to save the following information together with the final results:

```matlab
version
computer
```

For example:

```matlab
fprintf('MATLAB version:\n');
disp(version);

fprintf('Computer architecture:\n');
disp(computer);
```

---

# 3. Random Geometry Generation

The NTN transmitter positions used in the simulations are synthetic.

For each architecture, random:

* azimuths;
* elevation angles;
* transmitter–receiver distances

are generated subject to predefined constraints.

The geometry generator searches for a realization with PDOP inside a specified interval.

Conceptually:

```text
random geometry
      |
      v
elevation-mask constraint
      |
      v
distance-range constraint
      |
      v
PDOP evaluation
      |
      v
accept geometry if target interval is satisfied
```

The resulting geometry therefore represents a **controlled synthetic realization**, rather than an instantaneous state of a propagated operational constellation.

---

# 4. Scenario 1 — Isolated Architectures

Script:

```text
scenario1.m
```

Scenario 1 evaluates HAPS, LEO, MEO, and GEO separately.

The experiment includes Monte Carlo simulations of pseudorange measurements and WLS positioning.

The number of Monte Carlo experiments is controlled by:

```matlab
Nexpe
```

For code testing, a smaller value may be used:

```matlab
Nexpe = 500;
```

or:

```matlab
Nexpe = 1000;
```

For the final high-statistics simulation, use the value adopted for the dissertation results, for example:

```matlab
Nexpe = 50000;
```

The value of `Nexpe` should always be verified before generating final figures or tables.

For a reproducible Monte Carlo run, define a seed before geometry and noise generation, for example:

```matlab
rng(1,'twister');
```

The seed used for a new reproducible run should be documented together with the generated results.

---

# 5. Scenario 2 — Greedy Transmitter Selection

Script:

```text
scenario2.m
```

The base candidate set contains:

```math
4\,\mathrm{HAPS}
+
8\,\mathrm{LEO}
+
8\,\mathrm{MEO}
+
6\,\mathrm{GEO}
=
26
```

candidates.

The algorithm starts by evaluating all four-transmitter combinations:

```math
\binom{26}{4}
=
14\,950
```

The best four-transmitter subset according to the positional Cramér–Rao bound is used as the initialization of the greedy procedure.

After initialization, one transmitter is added at each iteration until the target:

```math
\mathcal{B}_{\mathrm{CRB}}(\mathcal{S})
\leq
\epsilon
```

is satisfied or all candidates are exhausted.

Because the candidate geometry is randomly generated, a different geometry realization may produce:

* different transmitter IDs;
* a different initial subset;
* a different greedy sequence;
* a different architecture composition;
* a different final number of selected transmitters.

For new reproducible experiments, use:

```matlab
rng(2,'twister');
```

before candidate-geometry generation.

However, if the original dissertation result was generated without a fixed seed, the saved:

```text
scenario_02_results.mat
```

should be used for exact reproduction of that original geometry.

---

# 6. Scenario 3 — Sensitivity to Reference C/N0

Script:

```text
scenario3.m
```

Scenario 3 generates one candidate geometry and keeps it fixed while varying the reference carrier-to-noise-density ratio.

The reference range is:

```math
30
\leq
(C/N_0)_{\mathrm{ref}}
\leq
50
\quad
\mathrm{dB\!-\!Hz}
```

Keeping the geometry fixed is essential because the objective is to isolate the effect of measurement quality.

Therefore, within one execution:

```text
same candidate geometry
        +
different C/N0 values
        =
measurement-quality sensitivity analysis
```

A new execution without a fixed seed may generate a different candidate geometry and therefore a different curve.

For new reproducible experiments:

```matlab
rng(2,'twister');
```

should be defined before geometry generation.

For exact reproduction of a previously generated experiment, use the saved:

```text
scenario_03_results.mat
```

when available.

---

# 7. Scenario 3a / Scenario 4 — HAPS Availability and LEO Compensation

Script:

```text
scenario3a.m
```

The script currently stores its results under:

```text
results_scenario_04/
```

and saves:

```text
scenario_04_results.mat
```

The analysis is divided into two parts.

---

## Part A — HAPS Restriction

The maximum number of HAPS allowed in the selected subset is:

```math
N_{\mathrm{HAPS,max}}
\in
\{4,2,0\}
```

The experiment does not necessarily remove specific HAPS from the candidate pool.

Instead, the selection constraint is:

```math
N_{\mathrm{HAPS}}(\mathcal{S})
\leq
N_{\mathrm{HAPS,max}}
```

Thus, for example, when:

```math
N_{\mathrm{HAPS,max}}=2
```

the algorithm can select the best two HAPS among the available HAPS candidates.

---

## Part B — LEO Compensation

HAPS are excluded from the selected configuration and LEO availability is increased according to:

```math
N_{\mathrm{LEO}}
\in
\{8,12,16,20\}
```

The candidate LEO sets are nested.

Therefore:

```math
\mathcal{C}_{8}
\subset
\mathcal{C}_{12}
\subset
\mathcal{C}_{16}
\subset
\mathcal{C}_{20}
```

This ensures that the analysis represents progressively increasing LEO availability rather than generating a completely independent LEO configuration for every case.

Because both the base geometries and the additional LEO geometry are random, a fixed seed is required for reproducible new experiments.

Recommended:

```matlab
rng(2,'twister');
```

For exact reproduction of a previously generated experiment, use:

```text
results_scenario_04/scenario_04_results.mat
```

when available.

---

# 8. Measurement Model

The nominal pseudorange uncertainty is modeled from the effective carrier-to-noise-density ratio.

For an architecture-dependent offset:

```math
(C/N_0)_{\mathrm{eff}}
=
(C/N_0)_{\mathrm{ref}}
+
\Delta(C/N_0)
```

the pseudorange standard deviation is:

```math
\sigma_{\rho}
=
\frac{c}
{2\pi\beta
\sqrt{
(C/N_0)_{\mathrm{linear}}
T_{\mathrm{coh}}
}}
```

The measurement covariance matrix is:

```math
\mathbf{R}
=
\mathrm{diag}
\left(
\sigma_{\rho_1}^{2},
\ldots,
\sigma_{\rho_M}^{2}
\right)
```

and the Fisher Information Matrix is:

```math
\mathbf{J}
=
\mathbf{H}^{T}
\mathbf{R}^{-1}
\mathbf{H}
```

The positional bound is calculated as:

```math
\mathcal{B}_{\mathrm{CRB}}
=
\sqrt{
\mathrm{tr}
\left(
[\mathbf{J}^{-1}]_{1:3,1:3}
\right)
}
```

---

# 9. Important Parameters to Record

When generating final results, the following parameters should be recorded together with the output:

```text
random seed
MATLAB version
simulation script version / Git commit
number of Monte Carlo experiments
C/N0 range
C/N0 architecture offsets
coherent integration time
code bandwidth
positioning target epsilon
number of candidates per architecture
PDOP target ranges
elevation masks
propagation-distance ranges
```

This allows a simulation result to be associated with a specific numerical configuration.

---

# 10. Git Version Control

For research reproducibility, each set of final simulation results should ideally be associated with a Git commit.

Before generating final results:

```bash
git status
```

Commit the exact version of the scripts:

```bash
git add .
git commit -m "Simulation configuration for final results"
```

The corresponding commit can be obtained using:

```bash
git rev-parse HEAD
```

The commit hash provides an unambiguous reference to the exact source-code version used in the experiment.

---

# 11. Suggested Result Metadata

A simple text file can be stored with each final simulation run.

Example:

```text
Simulation: Scenario 02
Date: YYYY-MM-DD
MATLAB version: R20XXx
Git commit: <commit-hash>
Random seed: 2
RNG generator: twister
C/N0 reference: 50 dB-Hz
Target epsilon: 1 m
Candidates: 4 HAPS + 8 LEO + 8 MEO + 6 GEO
```

For Scenario 1, also record:

```text
Monte Carlo realizations: 50000
```

This makes it much easier to reproduce published tables and figures later.

---

# 12. Exact Results vs. Statistical Generalization

The simulations in Scenarios 2, 3, and 3a should be interpreted primarily as controlled experiments on specific geometry realizations.

A fixed random seed provides:

```text
repeatability
+
controlled comparison
+
debugging capability
```

but it does **not** provide statistical characterization over all possible geometries.

A statistical robustness analysis would require repeating the complete geometry-generation and selection process over many independent realizations.

For example:

```math
N_{\mathrm{MC}} = 10\,000
```

geometry realizations.

For Scenario 2, the four-transmitter initialization alone would then require approximately:

```math
10\,000
\times
14\,950
=
149\,500\,000
```

subset evaluations.

This does not include the subsequent greedy iterations.

A large-scale geometry Monte Carlo analysis would therefore benefit from:

* code optimization;
* parallel processing;
* vectorization;
* reuse of intermediate calculations;
* more efficient subset-search strategies.

A GPU would only provide benefits if the implementation were explicitly redesigned to use GPU-compatible operations. The current MATLAB implementation is primarily limited by combinatorial search and CPU-side iterative computation.

---

# 13. Recommended Interpretation

The current reproducible simulations should be described as:

> **Controlled and reproducible realizations of synthetic NTN geometries used to isolate the effects of transmitter geometry, measurement quality, and availability constraints.**

They should not be described as:

> **A statistical representation of every possible HAPS, LEO, MEO, or GEO configuration.**

A Monte Carlo analysis over the geometry itself is considered a natural extension of the present work.

---

# 14. Recommended Workflow

For a new experiment:

```text
1. Pull the desired Git commit
2. Verify simulation parameters
3. Define the random seed
4. Execute the MATLAB script
5. Save generated .mat results
6. Save CSV tables and figures
7. Record MATLAB version
8. Record Git commit
9. Record RNG seed
```

For an exact reconstruction of a previously generated result:

```text
1. Checkout the corresponding Git commit
2. Load the saved .mat geometry/results
3. Avoid regenerating the random geometry
4. Recreate tables and figures from the stored variables
```

---

# 15. Scope

The reproducibility procedure described here is intended to preserve the numerical experiments associated with this research.

The most important distinction is:

```text
fixed RNG seed
        |
        v
reproducible new realization
```

versus:

```text
saved original geometry
        |
        v
exact reproduction of an existing result
```

Both approaches are useful, but they serve different purposes.
