<!-- BEGIN GENERATED: dev/rel066-counts.R -->

### Ungated suite, all eight packages

From `dev/release/suite.log`. SUITE ran 327 of 327

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 193 | 13494 | 0 | 0 | 163 | 0 |
| frmtmb.eam | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.sample | 43 | 2213 | 0 | 0 | 4 | 0 |
| frmtmb.coupling | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.learn | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.latent | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.spline | 15 | 553 | 0 | 0 | 1 | 0 |
| **total** | **327** | **19880** | **0** | **0** | **192** | **0** |

Files with a skip: frmtmb `test-bcm-bart.R` (1), frmtmb
`test-bcm-binomial.R` (6), frmtmb `test-bcm-data-analysis.R` (6),
frmtmb `test-bcm-esp.R` (2), frmtmb `test-bcm-gaussian.R` (3), frmtmb
`test-bcm-gcm.R` (2), frmtmb `test-bcm-latent-mixtures.R` (5), frmtmb
`test-bcm-model-selection.R` (3), frmtmb `test-bcm-mpt.R` (2), frmtmb
`test-bcm-psychophysics.R` (2), frmtmb `test-bcm-retention.R` (3),
frmtmb `test-bcm-signal-detection.R` (3), frmtmb `test-bcm-simple.R`
(1), frmtmb `test-brms-agreement.R` (2), frmtmb
`test-brms-likelihood.R` (48), frmtmb `test-brms-methods.R` (46),
frmtmb `test-brms-priors.R` (12), frmtmb `test-drmtmb-agreement.R`
(13), frmtmb `test-fuzz.R` (1), frmtmb `test-rl-example.R` (2),
frmtmb.eam `test-scale.R` (3), frmtmb.sample `test-loo.R` (2),
frmtmb.sample `test-sampling-ported.R` (1), frmtmb.sample
`test-scale.R` (1), frmtmb.coupling `test-scale.R` (5), frmtmb.learn
`test-scale.R` (2), frmtmb.learn `test-stan-identity.R` (11),
frmtmb.latent `test-scale.R` (2), frmtmb.ode `test-scale.R` (1),
frmtmb.spline `test-scale.R` (1).

### Gated tier, with the brms, drmTMB and fuzz gates set

From `dev/release/gated.log`. GATED ran 39 of 39

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 31 | 2810 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 1 | 71 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 7 | 525 | 0 | 0 | 0 | 0 |
| **total** | **39** | **3406** | **0** | **0** | **0** | **0** |

### Scale tier (FRMTMB_SCALE_TESTS set)

From `dev/release/scale.log`. SCALE ran 7 of 7

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb.eam | 1 | 18 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 1 | 2 | 0 | 0 | 0 | 1 |
| frmtmb.latent | 1 | 6 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 1 | 6 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 1 | 3 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 1 | 14 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 1 | 2 | 0 | 0 | 0 | 0 |
| **total** | **7** | **51** | **0** | **0** | **0** | **1** |

### frmtmb.ode's scale file, rerun after the wrap

From `dev/release/scale-ode-rerun.log`. SCALE-ODE-RERUN ran 1 of 1

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb.ode | 1 | 2 | 0 | 0 | 0 | 0 |
| **total** | **1** | **2** | **0** | **0** | **0** | **0** |

### frmtmb.ode's ungated suite, rerun after the new helper

From `dev/release/ode-recheck.log`. ODE-RECHECK ran 11 of 11

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb.ode | 11 | 547 | 0 | 0 | 1 | 0 |
| **total** | **11** | **547** | **0** | **0** | **1** | **0** |

Files with a skip: frmtmb.ode `test-scale.R` (1).

### Ported brms suite, verdicts asserted

From `dev/rel066-log/tier.txt`. TIER ran 15 of 15 files;
pass 634, fail 0, error 0, skip 0.

### The ledger (`dev/brmsport-ledger.tsv`)

| outcome | lane defects alone | release |
|---|---|---|
| pass | 306 | 368 |
| cannot transfer | 126 | 67 |
| defect | 21 | 14 |
| divergence | 41 | 45 |
| **total** | **494** | **494** |

Bin 1 passes: 306 of 494 before, 368 of 494 after.
Rows whose outcome or class moved: 79.

### Escaped warnings (`dev/release/warnscan-066/`)

- ungated: 327 files scanned, 0 with any line
- gated: 39 files scanned, 0 with any line
- escaped warnings: 0

### R CMD check --as-cran

From `dev/release/check-066/<pkg>/check.log`.

- frmtmb: Status: 1 NOTE (* checking HTML version of manual ... [19s] NOTE)
- frmtmb.eam: Status: 1 NOTE (* checking HTML version of manual ... NOTE)
- frmtmb.learn: Status: OK
- frmtmb.ode: Status: 1 NOTE (* checking HTML version of manual ... NOTE)
- frmtmb.sample: Status: OK
<!-- END GENERATED -->
