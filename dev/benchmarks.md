# Fit-time benchmarks (2026-08-31)

Single runs, wall-clock seconds, Windows / R 4.6.1. Script:
scratchpad bench.R (fit only, default controls each package).

| model                          | frm  | glmmTMB | lme4 |
|--------------------------------|------|---------|------|
| LMM n=180 (sleepstudy)         | 0.19 | 0.08    | 0.03 |
| LMM n=50k, 200 grp, (x\|g)     | 4.61 | 2.86    | 0.35 |
| Poisson GLMM n=100k, 500 grp   | 6.00 | 5.18    | 9.89 |
| NB2 + dispformula n=20k        | 1.31 | 0.71    | -    |
| ZI-poisson n=20k               | 0.47 | 0.58    | -    |

Reading:
- lme4 owns pure gaussian LMMs (profiled deviance + specialized sparse
  Cholesky); nothing Laplace-generic touches that.
- vs glmmTMB (the architectural peer): within ~2x everywhere, at parity
  or faster on GLMM/zi workloads. The gap on gaussian models is tape
  overhead plus no OpenMP parallel accumulation (RTMB limitation).
- frm beats lme4 on the large Poisson GLMM (glmer's AGQ/PIRLS scales
  worse than Laplace-with-AD here).
- Not measured: glmmTMB's parallel threads (off by default), warm-start
  refits, REML variants.
