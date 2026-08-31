SCENARIO 4 - REPRODUCIBILITY FILES
==================================

This folder contains the MATLAB scripts and TLE files required to
reproduce the numerical results of Scenario 4.

Scenario:
Hybrid NTN positioning geometry composed of:
- 4 synthetic quasi-stationary HAPS
- 8 real Starlink LEO satellites
- 7 real O3b/O3b mPOWER MEO satellites
- 7 real GEO satellites

Total candidate set: 26 transmitters.

Evaluation epoch:
01-Aug-2026 12:00:00 UTC

User position:
Latitude  = -22.8596582 deg
Longitude = -43.2303236 deg
Height    = 10 m

MATLAB REQUIREMENTS
-------------------
MATLAB with Satellite Communications Toolbox, required for:
- satelliteScenario
- satellite
- SGP4 propagation
- ECEF satellite states

No STK installation is required to reproduce the numerical results.
STK was used only for visual inspection of the generated geometry.

FILES
-----
scenario4_build_candidate_pool.m
    Builds the final 26-transmitter candidate set at the common epoch.
    Propagates the LEO, MEO and GEO TLEs using SGP4 and generates the
    four synthetic HAPS positions.

scenario4a.m
    Applies the pseudorange-quality model and performs transmitter
    selection using exhaustive initialization with four transmitters
    followed by forward greedy selection.

scenario4b.m
    Repeats the selection considering HAPS availability constraints:
    maximum 4 HAPS, maximum 2 HAPS, and no HAPS.

TLE files:
scenario4_LEO8_eval_2026-08-01T120000Z_combined.tle
scenario4_MEO7_eval_2026-08-01T120000Z_combined.tle
scenario4_GEO7_v2_eval_2026-08-01T120000Z_combined.tle

EXECUTION
---------
1. Set this folder as the MATLAB Current Folder.

2. Run:
   scenario4_build_candidate_pool.m

3. Run:
   scenario4a.m

4. Run:
   scenario4b.m

The scripts automatically create the output directories containing
the numerical results.

MAIN SIMULATION PARAMETERS
--------------------------
C/N0 reference = 50 dB-Hz

Architecture offsets:
HAPS = +8 dB
LEO  = +4 dB
MEO  =  0 dB
GEO  = -4 dB

beta = 1.023 MHz
Tcoh = 20 ms

Target positional bound = 1 m

Elevation masks:
HAPS = 15 deg
LEO  = 5 deg
MEO  = 5 deg
GEO  = 5 deg

NOTE
----
The C/N0 values are controlled simulation parameters representing
relative measurement-quality differences between architectures.
They do not represent a complete link-budget calculation.

The supplied TLE files correspond to the historical orbital elements
selected for the common evaluation epoch.