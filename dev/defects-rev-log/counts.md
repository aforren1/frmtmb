| run | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb, every file, gates open (lane) | SUITE lane frmtmb ran 183 of 183 files | 14459 | 0 | 0 | 14 | 0 |
| frmtmb.sample, every file, gates open (lane) | SUITE lane frmtmb.sample ran 36 of 36 files | 2168 | 0 | 0 | 1 | 0 |
| frmtmb.latent (lane core) | SUITE lane frmtmb.latent ran 10 of 10 files | 359 | 0 | 0 | 2 | 0 |
| frmtmb.coupling (lane core) | SUITE lane frmtmb.coupling ran 11 of 11 files | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam (lane core) | SUITE lane frmtmb.eam ran 29 of 29 files | 1733 | 0 | 0 | 3 | 0 |
| frmtmb.learn (lane core) | SUITE lane frmtmb.learn ran 15 of 15 files | 500 | 0 | 0 | 2 | 4 |
| frmtmb.ode (lane core) | SUITE lane frmtmb.ode ran 11 of 11 files | 547 | 0 | 0 | 1 | 0 |
| frmtmb.spline (lane core) | SUITE lane frmtmb.spline ran 15 of 15 files | 553 | 0 | 0 | 1 | 0 |
| changed frmtmb files (base, rellib-r3) | SUITE base frmtmb ran 5 of 5 files | 999 | 15 | 17 | 0 | 1 |
| changed frmtmb.sample file (base, rellib-r3) | SUITE base frmtmb.sample ran 1 of 1 files | 148 | 0 | 1 | 0 | 0 |

Per file, base and lane, for the files the lane changed or added:

| file | base pass/fail/err/skip | lane pass/fail/err/skip |
|---|---|---|
| test-brms-parity-defects.R | 2/9/15/0 | 78/0/0/0 |
| test-arg-refusal.R | 117/0/0/0 | 117/0/0/0 |
| test-compat.R | 634/0/1/0 | 640/0/0/0 |
| test-mv-gaps.R | 55/0/1/0 | 58/0/0/0 |
| test-prior-compat.R | 191/6/0/0 | 197/0/0/0 |
| frmtmb.sample test-draws-methods.R | 148/0/1/0 | 156/0/0/0 |

Ported tier, verdicts asserted (inside the runs above): 15 files, pass 632, fail 0, error 0, skip 0, warn 0

Mutants (dev/defects-rev-mutant.R):

- MUTANT M1 ord_linear_per_threshold test-brms-parity-defects.R pass= 76 fail= 2 err= 0; caught by: fitted(scale = "linear") of a cs() fit is one layer per threshold 
- MUTANT M1 ord_linear_per_threshold test-brms-suite-methods.R pass= 162 fail= 0 err= 0 
- MUTANT M2 prior_coef_rows test-brms-parity-defects.R pass= 76 fail= 2 err= 0; caught by: the prior table spells a smooth's column as brms does 
- MUTANT M2 prior_coef_rows test-brms-suite-priors.R pass= 34 fail= 2 err= 0; caught by: default_prior returns correct fixed effect names for GAMMs 
- MUTANT M3 predict_new_level_draw test-brms-parity-defects.R pass= 77 fail= 1 err= 0; caught by: predict() takes brms's sample_new_levels = "old_levels" 
- MUTANT M3 predict_new_level_draw test-brms-suite-methods.R pass= 162 fail= 0 err= 0 
- MUTANT M4 check_newdata_frame test-brms-parity-defects.R pass= 77 fail= 1 err= 0; caught by: numeric codes of a fitted factor in newdata are its levels 
- MUTANT M4 check_newdata_frame test-brms-suite-data-helpers.R pass= 5 fail= 1 err= 0; caught by: validate_newdata handles factors correctly 
- MUTANT M5 residuals_newdata test-brms-parity-defects.R pass= 76 fail= 2 err= 0; caught by: residuals() takes newdata and answers a multivariate fit 
- MUTANT M6 fitted_point test-brms-parity-defects.R pass= 73 fail= 1 err= 1; caught by: fitted(scale = "linear") of a cs() fit is one layer per threshold 
- MUTANT M6 fitted_point test-brms-suite-methods.R pass= 161 fail= 1 err= 0; caught by: fitted has reasonable outputs 
- MUTANT M7 predict_new_level_spec test-brms-parity-defects.R pass= 78 fail= 0 err= 0 
- MUTANT M7 predict_new_level_spec test-brms-suite-methods.R pass= 162 fail= 0 err= 0 
- MUTANT M8 summary_gp_frame test-brms-parity-defects.R pass= 77 fail= 1 err= 0; caught by: summary() reports a gp() term's hyperparameters, as brms does 
