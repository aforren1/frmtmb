| package | section (as filed) | source | record | bullet |
|---|---|---|---|---|
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | `mo()` codes an integer variable from its minimum, as brms |
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | `mo()` takes negative integers, as brms does; 0.75.0 refused |
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | New data for an integer `mo()` variable are coded from the f |
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | `conditional_effects()` holds a numeric `mo()` variable at |
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | `emmeans` refuses a reference grid that holds an integer `mo |
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | The standard error of a prediction includes the uncertainty  |
| frmtmb | Breaking changes | mocode | `dev/mocode-findings.md` | A fit saved by frmtmb 0.75.0 or earlier with an integer `mo( |
| frmtmb | Breaking changes | acprior | `dev/acprior-findings.md` | One residual-correlation parameter has one prior slot, as in |
| frmtmb | Breaking changes | acprior | `dev/acprior-findings.md` | The prior table bounds the residual-correlation classes as b |
| frmtmb | Breaking changes | acprior | `dev/acprior-findings.md` | A prior on `sigma` under a first-order `cov = TRUE` `ar()`,  |
| frmtmb | New features | acprior | `dev/acprior-findings.md` | Three functions join the interface for a sampling extension |
| frmtmb | Bug fixes | acprior | `dev/acprior-findings.md` | A bound at an end of a residual-correlation parameter's inte |
| frmtmb | Bug fixes | acprior | `dev/acprior-findings.md` | `?frmtmb-autocor` said that `set_prior()` cannot target the |
| frmtmb | New features | loglik | `dev/loglik-findings.md` | New in the sampling API: `autocor_row_loglik(fit, dpv, resp) |
| frmtmb | New features | loglik | `dev/loglik-findings.md` | `frmtmb_register_autocor()` takes `jacobian =`, the Jacobian |
| frmtmb.sample | Other changes | mocode | `dev/mocode-findings.md` | With the frmtmb that codes an integer `mo()` variable from |
| frmtmb.sample | Breaking changes | acprior | `dev/acprior-findings.md` | The default priors of a residual-correlation term are brms's |
| frmtmb.sample | Breaking changes | acprior | `dev/acprior-findings.md` | `sigma` under a first-order `cov = TRUE` `ar()`, `ma()` or ` |
| frmtmb.sample | Breaking changes | acprior | `dev/acprior-findings.md` | A prior of your own keeps the bound brms declares for its sl |
| frmtmb.sample | Breaking changes | acprior | `dev/acprior-findings.md` | A prior on `class = "theta", coef = "thetaac_1"` is refused  |
| frmtmb.sample | Bug fixes | acprior | `dev/acprior-findings.md` | A bounded residual-correlation parameter is sampled through  |
| frmtmb.sample | Bug fixes | acprior | `dev/acprior-findings.md` | When the ML estimate of `cosy` is below brms's bound at 0, t |
| frmtmb.sample | New features | loglik | `dev/loglik-findings.md` | `log_lik()` of draws no longer refuses a residual covariance |
| frmtmb.spatial | Bug fixes | acprior | `dev/acprior-findings.md` | The `default_prior()` table of a `sar()` model passes back a |
| frmtmb.spatial | New features | loglik | `dev/loglik-findings.md` | `sar()` gives its whitening map's Jacobian, `(I - rho W) / s |
