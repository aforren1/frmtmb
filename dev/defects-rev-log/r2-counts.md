| run | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb, every file, gates open (lane) | SUITE lane frmtmb ran 183 of 183 files | 14482 | 0 | 0 | 14 | 0 |
| frmtmb.sample, every file, gates open (lane) | SUITE lane frmtmb.sample ran 36 of 36 files | 2173 | 0 | 0 | 1 | 0 |
| frmtmb.latent (lane core) | SUITE lane frmtmb.latent ran 10 of 10 files | 359 | 0 | 0 | 2 | 0 |
| frmtmb.coupling (lane core) | SUITE lane frmtmb.coupling ran 11 of 11 files | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam (lane core) | SUITE lane frmtmb.eam ran 29 of 29 files | 1733 | 0 | 0 | 3 | 0 |
| frmtmb.learn (lane core) | SUITE lane frmtmb.learn ran 15 of 15 files | 500 | 0 | 0 | 2 | 4 |
| frmtmb.ode (lane core) | SUITE lane frmtmb.ode ran 11 of 11 files | 547 | 0 | 0 | 1 | 0 |
| frmtmb.spline (lane core) | SUITE lane frmtmb.spline ran 15 of 15 files | 553 | 0 | 0 | 1 | 0 |
| changed frmtmb files (base, rellib-r3) | SUITE base frmtmb ran 1 of 1 files | 9 | 14 | 19 | 0 | 1 |
| changed frmtmb.sample file (base, rellib-r3) | SUITE base frmtmb.sample ran 1 of 1 files | 149 | 2 | 1 | 0 | 0 |

Per file, base and lane, for the files the lane changed or added:

| file | base pass/fail/err/skip | lane pass/fail/err/skip |
|---|---|---|
| test-brms-parity-defects.R | 9/14/19/0 | 101/0/0/0 |
| test-arg-refusal.R | not run | 117/0/0/0 |
| test-compat.R | not run | 640/0/0/0 |
| test-mv-gaps.R | not run | 58/0/0/0 |
| test-prior-compat.R | not run | 197/0/0/0 |
| frmtmb.sample test-draws-methods.R | 149/2/1/0 | 159/0/0/0 |

Ported tier, verdicts asserted (inside the runs above): 15 files, pass 634, fail 0, error 0, skip 0, warn 0

Mutants (dev/defects-rev-mutant.R):

