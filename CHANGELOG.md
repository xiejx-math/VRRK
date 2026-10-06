# Release audit, v1.2.0 (2026-10-06)

- Validate plotting coordinates and use their actual shared interval without
  extrapolation. Identical trial grids retain every recorded checkpoint.
- Reject invalid synthetic dimensions, nonpositive inconsistency levels,
  unsupported solution types, and full-row-rank inconsistent constructions.
  Validate problem vector dimensions and finite values.
- Add source-only release packaging and document saved-result replotting and
  bounded tuning recovery. Preserve user-selected demo configurations.
- These changes do not alter solver updates. Replot saved histories when
  adopting the checkpoint-preserving plot behavior.

# Epoch-aligned double-median plots (2026-09-15)

- Supersedes the fixed-time plot below. All comparison time plots now use
  median elapsed time versus median RSE at shared epoch checkpoints.
  Epoch and time plots reuse identical RSE statistics and bands, including
  the endpoint. No missing-tail zero padding or extrapolation is performed.
- Algorithm updates, tuning, timing, and stopping rules are unchanged.

# Pointwise elapsed-time plots (2026-09-14)

- All four comparison demos now plot RSE statistics at common elapsed times:
  median lines, minimum--maximum light bands, and quartile dark bands.
  Each method uses only the time interval covered by all its trials, with
  previous-value alignment. All lines are drawn above the shaded bands.
- The first-attainment-time helper remains available but is no longer used
  by these plots. Saved histories can be replotted without rerunning solvers.

# Tuning recovery correction (2026-09-14)

- Candidates must reduce terminal RSE on both tuning streams, not merely
  remain below the divergence threshold. If a stage loses all candidates,
  revisit its full original grid at the longest tuning budget, then retry
  with steps scaled by 1/2 through 1/256 if needed. Stop at the first scale
  with a valid candidate and select by the existing history score.
- Recovery evaluations are logged with a `_recovery` stage and actual
  parameters. Final independent-stream stability validation remains required
  for both methods, with nine step scales and at least the longest tuning budget.
- Tuning ignores demo iteration caps. Failure logs include the training
  problem for reproduction. Exhausted recovery still raises an error.
- This is empirical validation, not a convergence guarantee or global
  parameter optimum. Re-run tuning and comparisons together for new results.

# Scientific reproducibility corrections (2026-09-13)

This revision changes benchmark timings and plotted statistics. Re-run the
experiments before replacing manuscript figures; previous CSVs and timings
must not be presented as results from this version.

- All four comparison demos draw samples inside solver timers. VRRK and
  AVRRK use the same reset RNG seed, while REK/RCD have independent seeds.
  CDF construction is outside timers; cached-CDF binary-search sampling is
  inside timers for all methods.
- Iteration time excludes initialization, problem construction, tuning and
  recorded RSE/normal-equation diagnostic calculations. It includes updates,
  sampling, loop/control overhead, and history-call overhead outside the
  explicitly excluded diagnostic block. It is elapsed time, not CPU time.
- All solvers save the actual terminal state, xFinal, stopReason and tolerance.
  Tolerance checks occur at recording checkpoints and the final state;
  tol=0 disables tolerance stopping. Optional cfg.maxIterations provides an
  iteration cap in addition to the accessed-nonzero work safeguard.
- The default AVRRK mode remains empirical. A separate small-matrix
  Demo_theoretical_parameters demonstrates the singleton formulas of
  Theorem 4.1. Its SVD and parameter construction are outside iteration time.
  AVRRK histories report the updated y sequence (stored as xFinal).
- Time plots use the first RECORDED target hit with no interpolation. Their
  resolution is limited by the recording interval. Optional unreached
  targets return Inf and reachedCount.
- Epoch statistics use the original RSE scale, including the arithmetic
  median for even trial counts. Percentiles use midpoint plotting positions,
  matching the traditional MATLAB prctile convention rather than type 7.
- Default plots use only the interval shared by all runs of each method.
  A curve endpoint is not necessarily the median terminal error. Failed
  trajectories must be investigated, not silently dropped from plots.
- Both methods undergo two-seed, fixed-budget stability validation.
  Every candidate evaluation logs its parameters, stage, seed and outcome.
  No valid/stable candidate raises VRRK:TuningFailed and saves a MAT log.
  Empirical stability checks do not establish a mathematical guarantee.
- Each demo creates a unique results subdirectory with case information,
  and saves run_settings.mat before running. Previously generated results
  are not overwritten.
- Matrix Market symmetric off-diagonal entries are mirrored correctly;
  unsupported header formats are explicitly rejected.

Validation on MATLAB R2022a (2026-09-13): run_regression_tests and
run_smoke_test both passed. The standalone Code Analyzer also reported no
diagnostics. These are small regression checks, not a rerun of the 20-trial
publication experiments or a validation of all manuscript conclusions.
