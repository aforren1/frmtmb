<!-- BEGIN GENERATED: dev/rel067-counts.R -->

### Ungated suite, all eight packages

From `dev/release/suite.log`. SUITE ran 338 of 338

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 202 | 14081 | 0 | 0 | 169 | 0 |
| frmtmb.eam | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.sample | 45 | 2299 | 0 | 0 | 4 | 0 |
| frmtmb.coupling | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.learn | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.latent | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.spline | 15 | 553 | 0 | 0 | 1 | 0 |
| **total** | **338** | **20553** | **0** | **0** | **198** | **0** |

Files with a skip: frmtmb `test-bcm-bart.R` (1), frmtmb
`test-bcm-binomial.R` (6), frmtmb `test-bcm-data-analysis.R` (6),
frmtmb `test-bcm-esp.R` (2), frmtmb `test-bcm-gaussian.R` (3), frmtmb
`test-bcm-gcm.R` (2), frmtmb `test-bcm-latent-mixtures.R` (5), frmtmb
`test-bcm-model-selection.R` (3), frmtmb `test-bcm-mpt.R` (2), frmtmb
`test-bcm-psychophysics.R` (2), frmtmb `test-bcm-retention.R` (3),
frmtmb `test-bcm-signal-detection.R` (3), frmtmb `test-bcm-simple.R`
(1), frmtmb `test-brms-agreement.R` (2), frmtmb
`test-brms-likelihood.R` (52), frmtmb `test-brms-methods.R` (47),
frmtmb `test-brms-priors.R` (12), frmtmb `test-drmtmb-agreement.R`
(13), frmtmb `test-fuzz.R` (1), frmtmb `test-offset-grid.R` (1), frmtmb
`test-rl-example.R` (2), frmtmb.eam `test-scale.R` (3), frmtmb.sample
`test-loo.R` (2), frmtmb.sample `test-sampling-ported.R` (1),
frmtmb.sample `test-scale.R` (1), frmtmb.coupling `test-scale.R` (5),
frmtmb.learn `test-scale.R` (2), frmtmb.learn `test-stan-identity.R`
(11), frmtmb.latent `test-scale.R` (2), frmtmb.ode `test-scale.R` (1),
frmtmb.spline `test-scale.R` (1).

### Gated tier, with the brms, drmTMB and fuzz gates set

From `dev/release/gated.log`. GATED ran 39 of 39

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 31 | 2891 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 7 | 525 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 1 | 71 | 0 | 0 | 0 | 0 |
| **total** | **39** | **3487** | **0** | **0** | **0** | **0** |

### Scale tier (FRMTMB_SCALE_TESTS set)

From `dev/release/scale.log`. SCALE ran 7 of 7

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb.eam | 1 | 18 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 1 | 3 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 1 | 6 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 1 | 14 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 1 | 6 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 1 | 2 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 1 | 2 | 0 | 0 | 0 | 0 |
| **total** | **7** | **51** | **0** | **0** | **0** | **0** |

### Ported brms suite, verdicts asserted

From `dev/rel067-log/tier.txt`. TIER ran 15 of 15 files;
pass 634, fail 0, error 0, skip 0.

### The ledger (`dev/brmsport-ledger.tsv`)

| outcome | 0.66.0 | 0.67.0 |
|---|---|---|
| pass | 368 | 387 |
| cannot transfer | 67 | 53 |
| defect | 14 | 2 |
| divergence | 45 | 52 |
| **total** | **494** | **494** |

Bin 1 passes: 368 of 494 before, 387 of 494 after.
Rows whose outcome or class moved: 30.

### Escaped warnings (`dev/release/warnscan-067/`)

- ungated: 338 files scanned; 0 with an escaped warning;
  0 the scan could not run
- gated: 39 files scanned; 0 with an escaped warning;
  0 the scan could not run
- escaped warnings: 0

### R CMD check --as-cran

From `dev/release/check-067/<pkg>/check.log`.

- frmtmb: Status: 1 NOTE (* checking HTML version of manual ... [20s] NOTE)
- frmtmb.sample: Status: OK
<!-- END GENERATED -->
