# The sentences over 25 words in the 0.75.0 NEWS sections
# (dev/rel075-news-long.R, dev/rel075-log/news-long-1.txt: core 29,
# frmtmb.latent 1, frmtmb.sample 2, frmtmb.survival 2, frmtmb.spatial
# 8), split into shorter ones, as ASD-STE100 asks, and the edits of fact
# the merge needs (a release a lane could not number, a code span longer
# than a line, and the residual terms the merged refusal covers).
# dev/rel075-news-merge.R applies each pair (lane text, new text) to the
# top section of the file named, before it rewraps the section; a run of
# whitespace in the lane text matches any whitespace, so a pair does not
# depend on where the lane wrapped its lines. The numbers and the facts
# are the lanes'; only the sentence boundaries and the order of clauses
# change, except where a comment says otherwise. Adapted from
# dev/rel074-news-split.R.
list(
  "NEWS.md" = list(
    c(r"---[Only the model can tell a declared simplex
      (`frmtmb_family(simplex_pars =)`) from a distributional parameter, so
      `set_prior()` now accepts the row for every class, and the refusal
      repeats the spelling as written.]---",
      r"---[Only the model can tell a declared simplex
      (`frmtmb_family(simplex_pars =)`) from a distributional parameter.
      So `set_prior()` now accepts the row for every class, and the
      refusal repeats the spelling as written.]---"),
    c(r"---[A brms prior object with such a row is refused the
      same way, where it was refused before as "a brms prior table has 1 row
      with no faithful frmtmb spelling".]---",
      r"---[A brms prior object with such a row is refused the same way.
      Before, it was refused as "a brms prior table has 1 row with no
      faithful frmtmb spelling".]---"),
    c(r"---[To keep `exp(x)` from being dropped as a covariance structure,
      the parser renames the call to `.frm_exp(x)`, and since 0.15.0 that
      alias was the name of the design column.]---",
      r"---[To keep `exp(x)` from being dropped as a covariance structure,
      the parser renames the call to `.frm_exp(x)`. Since 0.15.0 that
      alias was the name of the design column.]---"),
    c(r"---[So the coefficient was `.frm_expx` in `fixef()`, `vcov()`,
      `summary()`, `coef()`, `ranef()` and `VarCorr()`, `b_.frm_expx` in
      `variables()` and in the draws of `frmtmb.sample::frm_sample()`, and
      `.frm_exp(x)` in `confint()` and `default_prior()`.]---",
      r"---[So the coefficient was `.frm_expx` in `fixef()`, `vcov()`,
      `summary()`, `coef()`, `ranef()` and `VarCorr()`. It was
      `b_.frm_expx` in `variables()` and in the draws of
      `frmtmb.sample::frm_sample()`, and `.frm_exp(x)` in `confint()` and
      `default_prior()`.]---"),
    c(r"---[The `coef` column held the design column, such as `log(xp)`,
      `I(x^2)` or `poly(x, 2)1`, where brms 2.23.0 lists `logxp`, `IxE2`
      and `polyx21`, and `set_prior()` refused brms's name where the rename
      is more than dropped parentheses, such as `IxE2`, `polyx21` or
      `expxD2`.]---",
      r"---[The `coef` column held the design column, such as `log(xp)`,
      `I(x^2)` or `poly(x, 2)1`, where brms 2.23.0 lists `logxp`, `IxE2`
      and `polyx21`. And `set_prior()` refused brms's name where the
      rename is more than dropped parentheses, such as `IxE2`, `polyx21`
      or `expxD2`.]---"),
    c(r"---[The design column name still works in `set_prior(coef =)`, for
      class `"b"`, class `"sd"` and a `cs()` term alike, and the two
      spellings name one slot, so writing both is refused as a
      duplicate.]---",
      r"---[The design column name still works in `set_prior(coef =)`, for
      class `"b"`, class `"sd"` and a `cs()` term alike. The two spellings
      name one slot, so writing both is refused as a duplicate.]---"),
    c(r"---[`default_prior()` orders them as brms's does, by response,
      distributional parameter, nonlinear parameter, class, group and coef,
      so `b` comes before `cor`, `Intercept` and `sd`; `validate_prior()`
      and `prior_summary()` order them by class, group, response,
      distributional parameter, nonlinear parameter and coef, as brms's
      do.]---",
      r"---[`default_prior()` orders them as brms's does, by response,
      distributional parameter, nonlinear parameter, class, group and
      coef, so `b` comes before `cor`, `Intercept` and `sd`.
      `validate_prior()` and `prior_summary()` order them as brms's do, by
      class, group, response, distributional parameter, nonlinear
      parameter and coef.]---"),
    c(r"---[registers an addition term that carries no data and sets an
      argument of the family instead, as brms's `bhaz(df = 6)` sets the
      baseline basis of `cox()`.]---",
      r"---[registers an addition term that carries no data and sets an
      argument of the family instead. brms's `bhaz(df = 6)`, for example,
      sets the baseline basis of `cox()`.]---"),
    c(r"---[now say what each kind of term does: a registered data term
      carries a column to the family's density, and an option term adds
      nothing to the likelihood and sets an argument of the family of its
      own response.]---",
      r"---[now say what each kind of term does. A registered data term
      carries a column to the family's density. An option term adds
      nothing to the likelihood and sets an argument of the family of its
      own response.]---"),
    c(r"---[The residuals of a `gaussian()` or `student()` response are one
      draw with covariance `sigma^2 M`, or `D M D` with `D = diag(sigma)`
      under a predicted `sigma ~ ...`, as in brms 2.23.0.]---",
      r"---[As in brms 2.23.0, the residuals of a `gaussian()` or
      `student()` response are one draw with covariance `sigma^2 M`.
      Under a predicted `sigma ~ ...` the covariance is `D M D` with
      `D = diag(sigma)`.]---"),
    c(r"---[The log-likelihood equals brms's Stan program at frmtmb's
      estimates to 7e-15 (gaussian, with a brms gradient there below
      3e-6), and at set parameter values with `nu = 5` to 1.4e-14
      (student).]---",
      r"---[The log-likelihood equals brms's Stan program to 7e-15 at
      frmtmb's estimates (gaussian, with a brms gradient there below
      3e-6). It equals it to 1.4e-14 at set parameter values with
      `nu = 5` (student).]---"),
    c(r"---[They also take `newdata` with one row per row of `M`, read in
      their order, as brms reads it, so a counterfactual (`x + 1`
      everywhere) works; when the fit dropped rows for missing values,
      `newdata` covers them too, and a row with a missing predictor gives
      `NA`, as in brms; `fitted()` there equals brms's `posterior_epred()`
      at the same parameter values.]---",
      r"---[They also take `newdata` with one row per row of `M`, read in
      their order, as brms reads it. So a counterfactual (`x + 1`
      everywhere) works. When the fit dropped rows for missing values,
      `newdata` covers them too, and a row with a missing predictor gives
      `NA`, as in brms. `fitted()` there equals brms's
      `posterior_epred()` at the same parameter values.]---"),
    c(r"---[which passes a matrix with a negative eigenvalue; its fit of
      such an `M` then fails in Stan,]---",
      r"---[which passes a matrix with a negative eigenvalue. Its fit of
      such an `M` then fails in Stan,]---"),
    c(r"---[A package registers a term that correlates the residuals of a
      gaussian or student response, as one affine map from the response to
      independent innovations and its inverse.]---",
      r"---[A package registers a term that correlates the residuals of a
      gaussian or student response. The term is one affine map from the
      response to independent innovations, with its inverse.]---"),
    c(r"---[for every column: on a 50 x 50 lattice under `sar()` (2500 rows)
      they take 0.48 s and 4.3 s, where one solve per column took 16 s
      each.]---",
      r"---[for every column. On a 50 x 50 lattice under `sar()` (2500
      rows), they take 0.48 s and 4.3 s. One solve per column took 16 s
      each.]---"),
    c(r"---[frmtmb checks the map once: `unwhiten()` must give one value
      per row and one column per column of innovations, `whiten()` one
      innovation per row and one finite log-Jacobian, and
      `whiten(unwhiten(u))` must give back `u` to a relative 1e-6.]---",
      r"---[frmtmb checks the map once. `unwhiten()` must give one value
      per row and one column per column of innovations. `whiten()` must
      give one innovation per row and one finite log-Jacobian.
      `whiten(unwhiten(u))` must give back `u` to a relative 1e-6.]---"),
    c(r"---[`set_prior()` then takes the parameter's name as its class (`lb`
      and `ub` narrow the interval), and frmtmb.sample's default prior on
      it is flat on the interval, as brms's on `lagsar`.]---",
      r"---[`set_prior()` then takes the parameter's name as its class, and
      `lb` and `ub` narrow the interval. frmtmb.sample's default prior on
      it is flat on the interval, as brms's on `lagsar` is.]---"),
    c(r"---[autoscale leaves the random slopes of every structure but `us()`
      and `diag()` as they are, frmtmb's own `homdiag()` or `cs()` as much
      as a structure another package registers, because nothing tells
      frmtmb that their parameters scale with the column.]---",
      r"---[autoscale leaves the random slopes of every structure but
      `us()` and `diag()` as they are. That holds for frmtmb's own
      `homdiag()` or `cs()` as much as for a structure another package
      registers. Nothing tells frmtmb that their parameters scale with the
      column.]---"),
    # the merged refusal covers a residual correlation term too
    # (dev/rel075-resolve.R, item 1)
    c(r"---[Every method of a fit that reads the model's families or
      structures checks them all first, so]---",
      r"---[Every method of a fit that reads the model's families,
      structures or residual correlation terms checks them all first. So]---"),
    c(r"---[The refusal for a `cox()` fit an older frmtmb built names the
      other missing packages too, and where its new home, frmtmb.survival,
      is installed but cannot be loaded, it says to reinstall it.]---",
      r"---[The refusal for a `cox()` fit an older frmtmb built names the
      other missing packages too. Where its new home, frmtmb.survival, is
      installed but cannot be loaded, the refusal says to reinstall
      it.]---"),
    c(r"---[Every handler that reads an error as "no answer" now lets this
      refusal through, so `predict()` no longer warns first that "the
      covariance of the estimates could not be recovered from the Hessian",
      and `cooks.distance()` no longer says that "all deletion refits
      failed".]---",
      r"---[Every handler that reads an error as "no answer" now lets this
      refusal through. So `predict()` no longer warns first that "the
      covariance of the estimates could not be recovered from the
      Hessian". And `cooks.distance()` no longer says that "all deletion
      refits failed".]---"),
    c(r"---[for every model, and a `mixture_mvn()` fit is none of the
      three.]---",
      r"---[for every model. A `mixture_mvn()` fit is none of the
      three.]---"),
    c(r"---[The row is now `exp(x)`, with `Df 1` and AIC 625.69, the AIC of
      `y ~ z` refitted, and `step()` drops `z`, which lowers the AIC to
      605.84, where it kept both terms.]---",
      r"---[The row is now `exp(x)`, with `Df 1` and AIC 625.69, the AIC of
      `y ~ z` refitted. And `step()` drops `z`, which lowers the AIC to
      605.84; before, it kept both terms.]---"),
    c(r"---[name the same coefficient, but in 0.74.0 the pair was taken as
      two specifications,]---",
      r"---[name the same coefficient. But in 0.74.0 the pair was taken as
      two specifications,]---"),
    # a code span longer than a line, in two
    c(r"---[for example `mvbf(bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1, nl =
      TRUE), bf(y ~ c + d * x, c ~ exp(z), d ~ 1, nl = TRUE))`,]---",
      r"---[for example `mvbf(f1, f2)` with
      `f1 = bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE)` and
      `f2 = bf(y ~ c + d * x, c ~ exp(z), d ~ 1, nl = TRUE)`,]---"),
    c(r"---[The weights check, which refuses a weighted likelihood with no
      maximum, named its full-rank structures in a list, so it failed open
      for every structure another package registers.]---",
      r"---[The weights check refuses a weighted likelihood with no maximum.
      It named its full-rank structures in a list, so it failed open for
      every structure another package registers.]---"),
    c(r"---[A Gamma model with weights 1.5 and a registered iid block on one
      level per row that states `full_rank = TRUE` is refused now, as the
      same model with `(1 | obs)` is; before, it fitted to a logLik of 849
      with a shape of 5e14 and two warnings.]---",
      r"---[Take a Gamma model with weights 1.5 and a registered iid block
      on one level per row that states `full_rank = TRUE`. It is refused
      now, as the same model with `(1 | obs)` is. Before, it fitted to a
      logLik of 849 with a shape of 5e14 and two warnings.]---"),
    c(r"---[fitted to a logLik of 4.7e7 with a shape of 8e19 and two
      warnings; it is refused now,]---",
      r"---[fitted to a logLik of 4.7e7 with a shape of 8e19 and two
      warnings. It is refused now,]---"),
    c(r"---[frmtmb's R CMD check (four platforms) and its coverage job fail
      when a test skips because a suggested package is installed but
      cannot be loaded, or is not installed, as the extension jobs have
      since 0.74.0, and when testthat's output has no summary line.]---",
      r"---[frmtmb's R CMD check (four platforms) and its coverage job fail
      when a test skips because a suggested package is missing or cannot
      be loaded. The extension jobs have done so since 0.74.0. Every job
      with this gate also fails when testthat's output has no summary
      line.]---"),
    c(r"---[The GEV density test reads extraDistr, a suggested package, in
      place of evd, which frmtmb does not declare, so it runs on every
      platform; through 0.74.0 it skipped on all four R CMD check jobs and
      in coverage.]---",
      r"---[The GEV density test reads extraDistr, a suggested package, in
      place of evd, which frmtmb does not declare. So it runs on every
      platform. Through 0.74.0 it skipped on all four R CMD check jobs and
      in coverage.]---"),
    # release fixes (dev/rel075-release-fix.R, items 2 and 3): what the
    # merged code does for a bounded parameter's internal name and for
    # a residual term's name in frmtmb_covstruct_info()
    c(r"---[A class prior on a parameter without `bounds` is refused, with
      the internal parameter to write it on instead.]---",
      r"---[A class prior on a parameter without `bounds` is refused, with
      the internal parameter to write it on instead. A parameter with
      `bounds` has one prior slot, its class: its internal row is not in
      `default_prior()`, and a prior on it is refused.]---"),
    c(r"---[A name that core lists as an extension's (`"car"`, `"spde"`)
      loads that extension.]---",
      r"---[A name that core lists as an extension's (`"car"`, `"spde"`)
      loads that extension. A residual correlation term's name, such as
      `"sar"`, is refused with a pointer to `frmtmb_autocor_pars()`.]---"),
    # the test now reads a saved sar() fit too (dev/rel075-resolve.R,
    # item 6)
    c(r"---[The pkgdown job, the one that installs every extension, also
      runs the test that reads a saved `cox()` and `car()` fit where
      neither extension can be loaded, and fails if any part of it
      skips.]---",
      r"---[The pkgdown job installs every extension. It also runs the
      saved-fit test of `cox()`, `car()` and `sar()` with neither
      extension loadable, and fails if any part of it skips.]---")
  ),
  "extensions/frmtmb.latent/NEWS.md" = list(
    c(r"---[The `mixture_mvn()` tests compare each result as a ratio to a
      quantity the run measures (a log-likelihood's size, a standard error,
      the unit roundoff), not to an absolute tolerance, and each still
      fails on a constructed wrong answer.]---",
      r"---[The `mixture_mvn()` tests compare each result as a ratio to a
      quantity the run measures, not to an absolute tolerance. The
      quantity is a log-likelihood's size, a standard error or the unit
      roundoff. Each test still fails on a constructed wrong answer.]---")
  ),
  "extensions/frmtmb.sample/NEWS.md" = list(
    # a lane cannot write the release number
    c(r"---[Needs the frmtmb release that exports `frmtmb_covstruct_info()`.]---",
      r"---[Needs frmtmb 0.75.0, which exports `frmtmb_covstruct_info()`.]---"),
    c(r"---[The change is in frmtmb, which builds the names; the new test of
      it needs that frmtmb.]---",
      r"---[The change is in frmtmb, which builds the names, and the new test
      of it needs frmtmb 0.75.0.]---"),
    c(r"---[`sigma_Intercept`, `x`, `sigma_expz`; the draws' `fixef()`
      listed]---",
      r"---[`sigma_Intercept`, `x`, `sigma_expz`. The draws' `fixef()`
      listed]---"),
    c(r"---[was mapped onto that block's names, failed the map's own check
      and sent every draw through the slow path:]---",
      r"---[was mapped onto that block's names. It failed the map's own
      check and sent every draw through the slow path:]---")
  ),
  "extensions/frmtmb.survival/NEWS.md" = list(
    c(r"---[`bhaz(gr = g)`, brms's stratified Cox model with one baseline per
      level of `g`, stops with a message that says `cox()` has one
      baseline, and that tells you to fit each stratum as a model of its
      own or to add `g` to the formula.]---",
      r"---[`bhaz(gr = g)`, brms's stratified Cox model with one baseline
      per level of `g`, stops. The message says that `cox()` has one
      baseline. It tells you to fit each stratum as a model of its own, or
      to add `g` to the formula.]---"),
    c(r"---[62 rows lie below the lower knot and 1 above; at equal weights
      brms's cumulative baseline hazard reaches -1.08 there, so the
      survivor function exceeds one, and the hazard is negative for some
      weights.]---",
      r"---[62 rows lie below the lower knot and 1 above. At equal weights,
      brms's cumulative baseline hazard reaches -1.08 there, so the
      survivor function exceeds one. The hazard is negative for some
      weights.]---")
  ),
  "extensions/frmtmb.spatial/NEWS.md" = list(
    c(r"---[`sar(M, type = )`, brms's spatial simultaneous autoregressive
      term, in its lag and error forms, for `gaussian()` and `student()`
      with a constant `sigma` and `nu`, as in brms 2.23.0.]---",
      r"---[`sar(M, type = )`, brms's spatial simultaneous autoregressive
      term, in its lag and error forms. As in brms 2.23.0, it is for
      `gaussian()` and `student()` with a constant `sigma` and
      `nu`.]---"),
    c(r"---[The log-likelihood equals brms's Stan program at frmtmb's
      estimates, less the `-N/2 log(2 pi)` (`-N/2 log(pi)` for student)
      that brms leaves out, to 1.4e-14 (gaussian, at the estimates) and
      4e-15 (student, at set parameter values with `nu = 5`).]---",
      r"---[brms leaves out `-N/2 log(2 pi)` (`-N/2 log(pi)` for student).
      Less that constant, the log-likelihood equals brms's Stan program to
      1.4e-14 (gaussian, at frmtmb's estimates) and 4e-15 (student, at set
      parameter values with `nu = 5`).]---"),
    c(r"---[read as the locations in their order, as brms reads it; when the
      fit dropped rows for missing values,]---",
      r"---[read as the locations in their order, as brms reads it. When the
      fit dropped rows for missing values,]---"),
    c(r"---[gives the lag model's spillover, `(I - rho W)^-1 mu_new`, which
      equals brms's]---",
      r"---[gives the lag model's spillover, `(I - rho W)^-1 mu_new`. That
      equals brms's]---"),
    c(r"---[cuts off negative values below the smallest eigenvalue, so frmtmb
      can estimate]---",
      r"---[cuts off negative values below the smallest eigenvalue. So
      frmtmb can estimate]---"),
    c(r"---[On a row-standardized 8 x 8 queen lattice (eigenvalues -0.5006 to
      1) with data drawn at `rho = -0.8`, frmtmb's estimate is -1.103, with
      a logLik 3.74 above the best inside brms's interval; brms's posterior
      median is -0.443, with a fifth of its draws within 0.02 of its
      bound.]---",
      r"---[An example is a row-standardized 8 x 8 queen lattice
      (eigenvalues -0.5006 to 1) with data drawn at `rho = -0.8`.
      frmtmb's estimate is -1.103, with a logLik 3.74 above the best inside
      brms's interval. brms's posterior median is -0.443, with a fifth of
      its draws within 0.02 of its bound.]---"),
    c(r"---[where its log-determinant is not a number and Stan rejects the
      proposal, so brms's support there]---",
      r"---[where its log-determinant is not a number and Stan rejects the
      proposal. So brms's support there]---"),
    c(r"---[gives a warning: at that end of the interval `I - rho M` keeps
      one direction only, and with an intercept the likelihood has no
      maximum there, so `rho` runs to the end.]---",
      r"---[gives a warning. At that end of the interval `I - rho M` keeps
      one direction only. With an intercept the likelihood has no maximum
      there, so `rho` runs to the end.]---")
  )
)
