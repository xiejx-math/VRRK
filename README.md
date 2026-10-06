# VRRK-AS and AVRRK-AS

MATLAB implementations accompanying the article **“Efficient variance
reduction for randomized Kaczmarz method”** by Jiaxin Xie, Hou-Duo Qi, and
Deren Han.

The repository contains the row-action variance-reduced randomized Kaczmarz
method (VRRK-AS), its accelerated variant (AVRRK-AS), and the REK and RCD
baselines used in the numerical experiments. All reported test systems are
inconsistent.

## Requirements

- MATLAB R2022a or newer.
- Internet access is needed only when the SuiteSparse matrices are not already
  present under `data/suitesparse/`.
- No add-on MATLAB toolbox is required.

## Quick start

Clone or download the repository, open MATLAB in the repository root, and run:

```matlab
setup
run_smoke_test
run_regression_tests
```

The regression tests cover plateau handling, arithmetic medians, percentiles,
terminal states, stopping, stability failure, theoretical parameters and symmetric
matrix input. The smoke test uses a small synthetic problem and checks that VRRK-AS,
AVRRK-AS, REK, and RCD all return finite trajectories.

## Numerical experiment demos

Run `setup` once in each new MATLAB session, followed by one of these scripts:

```matlab
Demo_comparison_VRRK_AVRRK_quantiles
Demo_comparison_full_rank_quantiles
Demo_comparison_rank_deficient_quantiles
Demo_comparison_SuiteSparse_full_rank_quantiles
```

- `Demo_comparison_VRRK_AVRRK_quantiles.m` compares VRRK-AS and AVRRK-AS on
  one manually configurable synthetic problem.
- `Demo_comparison_full_rank_quantiles.m` compares VRRK-AS, AVRRK-AS, REK,
  and RCD on full-column-rank synthetic systems.
- `Demo_comparison_rank_deficient_quantiles.m` compares VRRK-AS, AVRRK-AS,
  and REK on rank-deficient synthetic systems. RCD is intentionally omitted.
- `Demo_comparison_SuiteSparse_full_rank_quantiles.m` compares all four
  methods on `ash958`, `well1033`, and `illc1033`.

Each script has a clearly marked parameter section near the top. Publication
settings use 20 independent trials and may take substantial time. For a short
test, set `runTime=2`, reduce `maxEpoch` (or `maxIter`), and use shorter tuning
budgets such as `[2 5 10]`.

Check `autoTune` before running: the current two-method demo preserves a
manual-parameter configuration (`autoTune=false`). Manual parameters are not
automatically stability-validated and may not transfer to other matrices.

The synthetic REK/RCD comparisons use explicitly named budget multipliers to
let slower baselines continue beyond the primary row-method budget. These are
stopping safeguards only; the plotted horizontal coordinates always use the
method-dependent epoch definitions below.

## Parameter tuning

VRRK-AS uses the dimensionless multiplier

```text
alpha = cAlpha / L_D,
L_D = max_i ||a_i||_2^2 / p_i.
```

When `autoTune=true`, the code tunes VRRK-AS and AVRRK-AS on an independently
generated problem instance. It uses staged successive halving for
`cAlpha`, `eta`, `cTheta`, `cGamma`, and `theta2`, followed by a longer
two-seed, fixed-budget stability check for both methods. A failed check raises
an error and saves the evaluation log, instead of returning an unvalidated
configuration. Selected parameters and the complete tuning log are saved
with the experiment results. Set `autoTune=false` to use the manual values in
the corresponding script.

If all surviving candidates fail, the tuner revisits the complete stage grid
at its longest budget, then tries step multipliers down to 1/256 of the
original values. Every accepted candidate must reduce terminal RSE on both
streams. Recovery is bounded and logged; these checks do not prove convergence
or identify globally optimal parameters.

The default accelerated parameters are **empirically tuned**, not the parameter
schedule in Theorem 4.1. Run `Demo_theoretical_parameters` for a small example
using `struct('mode',"theory")`. That mode computes a dense SVD, supports
singleton sampling, and is deliberately limited to small validation matrices.
AVRRK histories record its updated `y` sequence.

## Figures and statistics

The epoch definitions are

```text
VRRK-AS and AVRRK-AS: k/m
RCD:                   k/n
REK:                   k/m + k/n
```

In epoch plots, the bold line is the pointwise median RSE over the independent
trials. The dark band is the 25th--75th percentile range, and the light band
is the minimum--maximum range. Time plots use exactly the same epoch-aligned
RSE statistics, with the median elapsed time at each checkpoint as the x
coordinate. Thus both plots have identical RSE endpoints and bands. Time and
RSE use the same previous recorded checkpoint on the shared epoch grid.
Each method ends at the earliest final epoch among its trials; missing tails
are neither zero-padded nor extrapolated. Bands describe RSE variation, not
time variation. The separate
`time_curve_statistics` helper remains available for first-recorded-hit analyses,
but is not used by the comparison plots.
Default curves use only the range covered by all trials of each method;
their endpoints need not equal the terminal-summary medians.

All comparison demos draw samples inside the timed loops using cached CDFs.
The elapsed time includes sampling and updates but excludes initialization,
problem generation, parameter tuning, and RSE/normal-equation diagnostics.
For deterministic testing, the two row solvers still accept an explicit index
vector; such runs are marked `samplingIncluded=false` and must not be mixed
with the default timings.

The work safeguard `cfg.maxPasses` counts accessed nonzeros divided by
`nnz(A)`; it is distinct from the plotted epochs for irregular sparse
matrices. Optional `cfg.maxIterations` supplies a separate iteration cap.
Stopping is checked at diagnostic checkpoints and at the final state.
Histories contain `xFinal`, `stopReason`, and the actual final iteration.

Generated figures, CSV tables, tuning logs, and MAT files are written below
`results/`, which is ignored by Git. Each run uses a unique subdirectory and
saves its initial configuration as `run_settings.mat`.

To redraw saved trials without rerunning tuning or algorithms, run:

```matlab
Replot_saved_results
```

Choose an experiment folder or one case subfolder. The program saves new PNG,
FIG, and plotting-statistics MAT files under a unique `replot_double_median_*`
folder, preserving the original data. Identical epoch grids retain all recorded
checkpoints; otherwise a shared 501-point grid uses previous recorded values.

## SuiteSparse data

The SuiteSparse demo downloads Matrix Market archives directly from the
[SuiteSparse Matrix Collection](https://sparse.tamu.edu/) when necessary.
Downloaded data are stored below `data/suitesparse/` and are not committed to
this repository.

## Repository layout

```text
.
|-- Demo_*.m             numerical experiment scripts
|-- run_smoke_test.m     fast installation check
|-- setup.m              adds the source directory to the MATLAB path
|-- src/                 algorithms, tuning, plotting, and data utilities
|-- data/README.md       SuiteSparse download information
`-- results/             generated locally and ignored by Git
```

See [CHANGELOG.md](CHANGELOG.md) for corrections affecting old results.
Results from older timing implementations must be rerun before timing comparisons.

## Preparing a GitHub release

```matlab
run_smoke_test
run_regression_tests
package_release
```

Upload the contents of `VRRK-codes-GitHub-release.zip` as the repository root.
The archive contains only allowlisted source, documentation, and test fixtures;
it excludes local results and downloaded matrices. Existing demo settings are
editable example configurations, not a guarantee of reproducing every figure
in the manuscript. Plot-only updates can reuse saved histories; changes to
algorithms, tuning, or timing require rerunning the relevant experiments.

## License

This code is released under the MIT License. See `LICENSE`.
