# Fix round 1: the NEWS of core, frmtmb.spatial and frmtmb.sample. Part
# A's edits go to the worktree and the stage-A copy, part B's to the
# worktree only. Rscript <this> <worktree> <stageA>
args <- commandArgs(TRUE)
wt <- args[1]
st <- args[2]
edit <- function(root, p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
core_a <- list(
  c("  in brms; a vector is its diagonal. The log-likelihood equals brms's
  Stan program at frmtmb's estimates to 7e-15 (gaussian) and 8e-7
  (student, at `nu` near 5e9), and the brms gradient there is below
  3e-6. A brms fit of 64 rows agrees with the maximum likelihood fit
  to within 0.33 standard errors in every parameter. `fitted()`,
  `residuals()` (pearson divides by `sigma * sqrt(diag(M))`),
  `simulate()`, `predict()` and `autocor_matrix()` work on the fitted
  rows. `newdata` is refused there, because `M` covers the fitted rows
  only; brms refuses it too unless `newdata2` holds a matrix over the
  new rows. `frm_linpred()` takes `newdata`. frmtmb refused `fcor()`
  by name before.",
    "  in brms; a vector is its diagonal. The log-likelihood equals brms's
  Stan program at frmtmb's estimates to 7e-15 (gaussian, with a brms
  gradient there below 3e-6), and at set parameter values with
  `nu = 5` to 1.4e-14 (student). A brms fit of 64 rows agrees with the
  maximum likelihood fit to within 0.33 standard errors in every
  parameter. `fitted()`, `residuals()` (pearson divides by
  `sigma * sqrt(diag(M))`), `simulate()`, `predict()` and
  `autocor_matrix()` work on the fitted rows. They also take `newdata`
  with as many rows as the fit, read as the fitted rows in their order,
  as brms reads it, so a counterfactual (`x + 1` everywhere) works;
  `fitted()` there equals brms's `posterior_epred()` at the same
  parameter values. Other `newdata` is refused, because `M` covers the
  fitted rows only; brms refuses it too unless `newdata2` holds a
  matrix over the new rows. `frm_linpred()` takes any `newdata`. frmtmb
  refused `fcor()` by name before."),
  c("- brms 2.23.0 checks that `M` is positive definite with
  `min(eigen(M)$values <= 0)`, which passes a matrix with a negative
  eigenvalue. frmtmb refuses every `M` that is not positive definite.",
    "- brms 2.23.0 checks that `M` is positive definite with
  `min(eigen(M)$values <= 0)`, which passes a matrix with a negative
  eigenvalue; its fit of such an `M` then fails in Stan, where the
  Cholesky factor of `M` is taken. frmtmb refuses every `M` that is not
  positive definite, before the fit."),
  c("- brms's `sar()` is now supported in frmtmb.spatial. Without that
  package, `sar()` stops with a message that names it, where frmtmb
  refused it as unsupported before.",
    "- brms's `sar()` is now supported in frmtmb.spatial. Without that
  package, `sar()` in the formula or in `bf(autocor = ~ sar(M))` stops
  with a message that names it, where frmtmb refused it as unsupported
  before. `set_prior()` takes `class = \"lagsar\"` and
  `class = \"errorsar\"`, as brms does."),
  c("  cannot express such a term: it adds a latent vector to the linear
  predictor, and a lag SAR changes the mean of the response itself.
  `frmtmb_local_registries()` saves and restores the new registry.",
    "  cannot express such a term: it adds a latent vector to the linear
  predictor, and a lag SAR changes the mean of the response itself.
  `frmtmb_local_registries()` saves and restores the new registry.
- The map gets one `sigma` when the term sets `sigma_varies = FALSE`.
  `unwhiten()` takes a matrix of innovations, one column each, so the
  pearson residuals and `autocor_matrix()` solve the term's system once
  for every column: on a 50 x 50 lattice under `sar()` (2500 rows) they
  take 0.48 s and 4.3 s, where one solve per column took 16 s each.
  When the model is built, frmtmb checks the map once: `unwhiten()` must
  give one value per row and one column per column of innovations,
  `whiten()` one innovation per row and one finite log-Jacobian, and
  `whiten(unwhiten(u))` must give back `u` to a relative 1e-6. A map
  that fails is refused with the term's name.
- A term's `bounds` slot states the interval of a parameter that
  `natural` maps by `lower + (upper - lower) plogis(theta)`.
  `set_prior()` then takes the parameter's name as its class (`lb` and
  `ub` narrow the interval), and frmtmb.sample's default prior on it is
  flat on the interval, as brms's on `lagsar`. A class prior on a
  parameter without `bounds` is refused, with the internal parameter to
  write it on instead. New export `frmtmb_autocor_pars()` lists the
  parameters of a fit's registered terms, with their intervals.
- A slot of `frmtmb_register_autocor()` or
  `frmtmb_register_covstruct()` can be a function of `...`."))
core_b <- list(
  c("- New export `frmtmb_covstruct_info()`: what the registry knows about
  one covariance structure (its package, whether it is a field, whether
  its levels are coupled, whether it is full rank). frmtmb.sample reads
  it instead of keeping its own list of field structures.
- `frmtmb_register_covstruct(full_rank =)` states whether a structure's
  covariance is positive definite at some parameter value. Without it,
  frmtmb reads the property off the structure's covariance at the
  start: a bar structure with uncoupled levels and a positive definite
  covariance there is full rank.
- autoscale leaves the random slopes of a structure another package
  registers as they are, because nothing tells frmtmb that its
  parameters scale with the column, and now says so in a message that
  names the structure and the column. It said nothing before.",
    "- New export `frmtmb_covstruct_info()`: what the registry knows about
  one covariance structure (its package, whether it is a field, whether
  its levels are coupled, whether it is full rank). frmtmb.sample reads
  it instead of keeping its own list of field structures. A name that
  core lists as an extension's (`\"car\"`, `\"spde\"`) loads that
  extension.
- `frmtmb_register_covstruct(full_rank =)` states whether an effect a
  row owns can reach any value with a finite density: a covariance that
  is positive definite at some parameter value and a density that is
  positive everywhere. A structure that states nothing leaves the
  weights check open for it, because its covariance does not tell: a
  latent with bounded support keeps a weighted model's maximum. A
  package whose structure is full rank should state
  `full_rank = TRUE`.
- autoscale leaves the random slopes of every structure but `us()` and
  `diag()` as they are, frmtmb's own `homdiag()` or `cs()` as much as a
  structure another package registers, because nothing tells frmtmb
  that their parameters scale with the column. It now says so in a
  message that names the structure and the column. It said nothing
  before."),
  c("  structure's `full_rank`. A Gamma model with weights 1.5 and a
  registered iid block on one level per row is refused now, as the
  same model with `(1 | obs)` is; before, it fitted to a logLik of 849
  with a shape of 5e14 and two warnings. A registered structure that is
  not full rank is not refused.",
    "  structure's `full_rank`. A Gamma model with weights 1.5 and a
  registered iid block on one level per row that states
  `full_rank = TRUE` is refused now, as the same model with `(1 | obs)`
  is; before, it fitted to a logLik of 849 with a shape of 5e14 and two
  warnings. A registered structure that is not full rank, or states
  nothing, is not refused.
- The Student-t latents of `gr(dist = \"student\")` were not full rank
  for the weights check, so it failed open for them. A Gamma model with
  weights 1.5 and `(1 | gr(obs, dist = \"student\"))` fitted to a logLik
  of 4.7e7 with a shape of 8e19 and two warnings; it is refused now, as
  the same model with `(1 | obs)` is, and the refusal names the term as
  written."))
for (root in c(wt, st)) edit(root, "NEWS.md", core_a)
edit(wt, "NEWS.md", core_b)

spatial <- list(
  c("  `-N/2 log(2 pi)` (`-N/2 log(pi)` for student) that brms leaves out,
  to 1.4e-14 (gaussian) and 4.4e-6 (student, at `nu` near 5e9). brms
  fits of 64 rows agree with the maximum likelihood fits to within
  0.44 standard errors in every parameter. `fitted()` reports brms's
  `posterior_epred()`, `(I - rho W)^-1 mu` for the lag form.
  `summary()` and `confint_varcorr()` report `lagsar` or `errorsar`.",
    "  `-N/2 log(2 pi)` (`-N/2 log(pi)` for student) that brms leaves out,
  to 1.4e-14 (gaussian, at the estimates) and 4e-15 (student, at set
  parameter values with `nu = 5`). brms fits of 64 rows agree with the
  maximum likelihood fits to within 0.44 standard errors in every
  parameter. `fitted()` reports brms's `posterior_epred()`,
  `(I - rho W)^-1 mu` for the lag form. `summary()` and
  `confint_varcorr()` report `lagsar` or `errorsar`, and `set_prior()`
  takes them as classes: a density on `rho`, carried to the internal
  parameter.
- `fitted()`, `predict()` and `simulate()` take `newdata` with as many
  rows as the fit, read as the fitted locations in their order, as
  brms reads it. A counterfactual such as `x + 1` everywhere then gives
  the lag model's spillover, `(I - rho W)^-1 mu_new`, which equals
  brms's `posterior_epred()` on the same `newdata` at the same
  parameter values. Other `newdata` is refused, as brms refuses it
  without `newdata2`."),
  c("- `rho` is estimated on `(1 / min eigenvalue, 1 / max eigenvalue)` of
  `M`, where `I - rho M` has a positive determinant. brms bounds it by
  `(min eigenvalue, max eigenvalue)`, which for a row-standardized `M`
  cuts off negative values below the smallest eigenvalue, and for a
  matrix that is not row-standardized reaches past
  `1 / max eigenvalue`, where brms's log-determinant is the log of a
  negative number.",
    "- `rho` is estimated on `(1 / min eigenvalue, 1 / max eigenvalue)` of
  `M`, where `I - rho M` has a positive determinant. brms bounds it by
  `(min eigenvalue, max eigenvalue)`. For a row-standardized `M` that
  cuts off negative values below the smallest eigenvalue, so frmtmb can
  estimate `rho` beyond brms's bound, and the two are then different
  models. On a row-standardized 8 x 8 queen lattice (eigenvalues -0.5006
  to 1) with data drawn at `rho = -0.8`, frmtmb's estimate is -1.103,
  with a logLik 3.74 above the best inside brms's interval; brms's
  posterior median is -0.443, with a fifth of its draws within 0.02 of
  its bound. For a matrix that is not row-standardized, brms's interval
  reaches past `1 / max eigenvalue`, where its log-determinant is not a
  number and Stan rejects the proposal, so brms's support there is
  frmtmb's interval, at the cost of rejected proposals.
- A matrix whose smallest or largest eigenvalue repeats for all rows
  but one, such as a complete graph, gives a warning: at that end of
  the interval `I - rho M` keeps one direction only, and with an
  intercept the likelihood has no maximum there, so `rho` runs to the
  end. The model is not refused, because without an intercept it has
  its maximum."))
for (root in c(wt, st)) edit(root, "extensions/frmtmb.spatial/NEWS.md", spatial)

sample <- list(
  c("- The draws of a model with frmtmb.spatial's `sar()` carry brms's
  `lagsar` or `errorsar` column, derived from the sampled internal
  parameter, as `ar[1]` and `cosy` are.",
    "- The draws of a model with frmtmb.spatial's `sar()` carry brms's
  `lagsar` or `errorsar` column, derived from the sampled internal
  parameter, as `ar[1]` and `cosy` are. The names come from core's
  registry (`frmtmb::frmtmb_autocor_pars()`), so a term another package
  registers gets its own column the same way.
- The default prior on `lagsar` and `errorsar` is brms's: flat on
  `rho` over its interval, which `frm_sample()` announces as
  `uniform(lower, upper)`. The same holds for a bounded parameter of a
  term another package registers. `set_prior(class = \"lagsar\")` and
  `class = \"errorsar\"` replace it."))
for (root in c(wt, st)) edit(root, "extensions/frmtmb.sample/NEWS.md", sample)
cat("ok\n")
