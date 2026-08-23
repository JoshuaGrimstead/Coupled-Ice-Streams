# Variability in a Model of Divergent and Convergent Coupled Ice Streams
Code used to support the MSc dissertation of Joshua Grimstead and subsequently
Grimstead et al. (2026). Ice stream model code originally developed by Alex Robel and initially adapted by Kolja Kypke for a divergent coupling.


This repository contains a box model of an ice stream with divergent and convergent
topologies, together with the analysis used to characterise its dynamics:
time series, Poincaré sections, a high-resolution parameter sweep in
(`T_s2`, `L1`), maximal Lyapunov exponents by Rosenstein's method, and
bifurcation diagrams in `L2` and `T_s2`. A convergent-topology variant of
the model is included for comparison.

>

---

## Requirements

- MATLAB [R2025b] or later
- Statistics and Machine Learning Toolbox (`pdist2`, `prctile`)
- Parallel Computing Toolbox (`parfor`, `parpool`) — for the parameter
  sweep and the Rosenstein grid

The parameter sweep and the Lyapunov grid were run on an HPC resource using 48 cores. Both call `parpool(48)`; reduce this to suit your
machine. Everything else runs on a desktop.

---

## Repository layout

### Model

| File | Purpose |
|---|---|
| `Ice_Stream_Box_Model_Divergent.m` | The divergent-topology model: state derivatives and event function |
| `Ice_Stream_Box_Model_Convergent.m` | The convergent-topology model: state derivatives and event function |
| `p_base_values.mat` | Shared physical parameters for divergent model |
| `p_base_values_conv.mat` | Shared physical parameters for convergent model|
| `peak_Vice_Divergent.m` | Find peaks in ice stream volume for the divegent model |
| `peak_Vice_Convergent.m` | Find peaks in ice stream volume for the convergent model |

### Single-run diagnostics

| File | Purpose |
|---|---|
| `Ice_Stream_Box_Model_Divergent_Timeseries.m` | Integrates one parameter set and plots the state through time |
| `Ice_Stream_Box_Model_Divergent_Poincare.m` | Poincaré section for one parameter set |
| `Ice_Stream_Box_Model_Divergent_Transient_Test.m` | Sensitivity of the solution to `L2` and `L3`, used to set the transient discarded before analysis |
| `Ice_Stream_Box_Model_Convergent_Timeseries.m` | Time series for the convergent model |


### Parameter sweep

| File | Purpose |
|---|---|
| `Ice_Stream_Box_Model_Divergent_Parameter_Sweep.m` | High-resolution sweep over (`T_s2`, `L1`); writes Poincaré crossings and crossing times for every cell |

### Lyapunov exponents

| File | Purpose |
|---|---|
| `Ice_Stream_Box_Model_Divergent_Rosenstein_Cell.m` | Maximal Lyapunov exponent for a single Poincaré section |
| `Ice_Stream_Box_Model_Divergent_Rosenstein_Grid.m` | Applies the above across the sweep grid in parallel |
| `Ice_Stream_Box_Model_Divergent_Rosenstein_Run.m` | Loads the sweep output, runs the grid, saves the result |

### Bifurcation diagrams

| File | Purpose |
|---|---|
| `Ice_Stream_Box_Model_Divergent_Length_Bif_Hysteresis.m` | Bifurcation diagram varying `L2` demonstrating hysteresis |
| `Ice_Stream_Box_Model_Divergent_Temp_Bif.m` | Bifurcation diagram varying `T_s2` |
| `Ice_Stream_Box_Model_Convergent_Length_Bif_Hysteresis.m` | As above, convergent model |
| `Ice_Stream_Box_Model_Convergent_Temp_Bif.m` | As above, convergent model |

### Plotting

| File | Purpose |
|---|---|
| `Ice_Stream_Box_Model_Divergent_Frequency_Plot.m` | Plots from the parameter sweep comparing the frequency of box 2 and box 3 void ratios |
| `Ice_Stream_Box_Model_Divergent_Overlay_Switches.m` | Sweep and Lyapunov maps with the dynamical switches overlaid |

---

## Reproducing the results

The sweep is the expensive step; everything downstream reads its output,
so the stages can be run independently once it exists.

1. **Parameter sweep** — run
   `Ice_Stream_Box_Model_Divergent_Parameter_Sweep.m`. Produces
   `raw_simulation2_output3_HPC.mat`, containing `PC_results_e2_all` and
   `crossing_times_e2_mat` (one Poincaré section and crossing-time vector
   per grid cell) along with `T_s2_vals`, `L1_vals` and
   `ThresholdCrossingIndiciesMat`.

2. **Lyapunov exponents** — run
   `Ice_Stream_Box_Model_Divergent_Rosenstein_Run.m`. Produces
   `rosenstein_grid_output3.mat` with `lambda`, `lambda_phys`,
   and `cloud_rel` .


4. **Complexity and winding** — `Ice_Stream_Box_Model_Divergent_Frequency_Plot.m` produces plots and
   `arnold_tongue_complexity_winding_frequencyratio.mat`.

5. **Plots** — run
   `Ice_Stream_Box_Model_Divergent_Frequency_Plot.m` and
   `Ice_Stream_Box_Model_Divergent_Overlay_Switches.m` for
   parameter sweep plots. All other plots are contained within their respetive code files.

---

## Reading the Lyapunov output

Three matrices, each `[num_T x num_L]`, aligned with `T_s2_vals` and
`L1_vals`:

- **`lambda`** — the exponent **per crossing**: the slope of the mean
  log-separation curve, i.e. the exponent of the Poincaré map itself.
  Because return times vary across the grid, cells are not directly
  comparable in physical time on this map.
- **`lambda_phys`** — the same exponent per kyr, obtained by dividing by
  the cell's mean return time. `NaN` where crossing times were missing or
  inconsistent, which does not invalidate `lambda`.
- **`cloud_rel`** — the normalised size of the attractor, the largest
  per-column spread relative to that column's own magnitude.


### Method notes

The exponent follows Rosenstein et al. (1993). Two choices are specific
to this application and worth stating:

- **The fit window is set by the attractor's own size.** A Poincaré
  section has a finite diameter, so nearest-neighbour separations cannot
  grow past it and the log-separation curve must saturate whatever the
  dynamics. The fit therefore stops where the curve reaches `SatFrac` of
  that diameter rather than at a visually chosen scaling region.
- **The fit starts at k = 1, not k = 0.** The first step reflects the
  choice of neighbour rather than the dynamics: the initial separation is
  a nearest-neighbour distance and so lies far below typical separation,
  and pairs revert to typical separation within one step whatever the
  cell is doing. Left in, it inflates the exponent most where the cloud
  is clustered — a period-2 orbit has a small initial separation within a
  cluster and typical cluster spread thereafter — so low-period cells
  near a boundary would be reported as chaotic.

State columns are selected by `Dims`, which defaults to
`[1 5 6 7 8 9 10]` as 2 can become arbitrarily large 
and all other columns cannot show chaotic behaviour.

---

## Data

Given their size, the .mat outputs must be regenerated. Base parameter values are provided, see **Model**.

