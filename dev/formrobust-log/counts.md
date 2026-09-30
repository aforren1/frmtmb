### Final ungated suite (NOT_CRAN true) (`dev/formrobust-log/suite2-files`)

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 199 | 199 | 13646 | 0 | 0 | 163 | 0 |
| frmtmb.coupling | 11 | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 10 | 359 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 15 | 429 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 11 | 547 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 44 | 44 | 2223 | 0 | 0 | 4 | 0 |
| frmtmb.spline | 15 | 15 | 553 | 0 | 0 | 1 | 0 |

Libraries the files loaded (package path):

- 199 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.coupling
- 29 x C:/Users/adf44/source/r/rellib-r4/frmtmb.eam
- 10 x C:/Users/adf44/source/r/rellib-r4/frmtmb.latent
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.learn
- 11 x C:/Users/adf44/source/r/rellib-r4/frmtmb.ode
- 44 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample
- 15 x C:/Users/adf44/source/r/rellib-r4/frmtmb.spline

### Gated tier (FRMTMB_BRMS_FIT_TESTS true) (`dev/formrobust-log/gated1-files`)

| package | files | with RESULT | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | 31 | 31 | 2802 | 8 | 0 | 0 | 0 |
| frmtmb.sample | 5 | 5 | 201 | 0 | 0 | 0 | 0 |

Files with a failure, an error, a warning or no RESULT line:

- frmtmb test-brms-suite-brm.R: pass=22 fail=1 err=0 skip=0 warn=0
- frmtmb test-brms-suite-methods.R: pass=159 fail=3 err=0 skip=0 warn=0
- frmtmb test-brms-suite-standata.R: pass=83 fail=4 err=0 skip=0 warn=0

Libraries the files loaded (package path):

- 31 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb
- 5 x C:/Users/adf44/source/r/wt-formrobust-lib/frmtmb.sample

### The new test files, before (rellib-r4) and after (lane library)

| file | before | after |
|---|---|---|
| frmtmb--test-arma-na-newdata | pass=0 fail=0 err=4 skip=0 warn=0 | pass=13 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-aterm-expr | pass=0 fail=0 err=8 skip=0 warn=0 | pass=49 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-bernoulli-coding | pass=0 fail=0 err=5 skip=0 warn=0 | pass=37 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-brms-api-formrobust | pass=1 fail=0 err=6 skip=0 warn=0 | pass=26 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-offset-grid | pass=1 fail=1 err=3 skip=0 warn=0 | pass=15 fail=0 err=0 skip=0 warn=0 |
| frmtmb--test-update-pool | pass=1 fail=1 err=2 skip=0 warn=0 | pass=11 fail=0 err=0 skip=0 warn=0 |
| frmtmb.sample--test-formrobust-draws | pass=0 fail=0 err=2 skip=0 warn=0 | pass=10 fail=0 err=0 skip=0 warn=0 |

