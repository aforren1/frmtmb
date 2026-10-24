# The sentences over 25 words in the 0.74.0 NEWS sections
# (dev/rel074-news-long.R, dev/rel074-log/news-long-1.txt: core 13,
# frmtmb.latent 3, frmtmb.survival 1, frmtmb.spatial 2), split into
# shorter ones, as ASD-STE100 asks, and four edits of fact the merge
# needs. dev/rel074-news-merge.R applies each pair (lane text, new text)
# to the top section of the file named, before it rewraps the section;
# a run of whitespace in the lane text matches any whitespace, so a pair
# does not depend on where the lane wrapped its lines. The numbers and
# the facts are the lanes'; only the sentence boundaries and the order
# of clauses change, except where a comment says otherwise.
list(
  "NEWS.md" = list(
    c(r"---[A fit is the same fit: on eight `cox()` models (censored, a
      frailty, `df = 6` with the identity link, truncated, a `sbhaz`
      prior, and the three spellings) the log likelihood, the estimates
      and the standard errors are `identical()` to frmtmb 0.73.0's, and
      `get_prior()` gives the same table.]---",
      r"---[A fit is the same fit. The tests have eight `cox()` models:
      censored, a frailty, `df = 6` with the identity link, truncated, a
      `sbhaz` prior, and the three spellings. On each, the log
      likelihood, the estimates and the standard errors are
      `identical()` to frmtmb 0.73.0's, and `get_prior()` gives the same
      table.]---"),
    # the splitsurv review, finding 8: validate_prior() refuses it too
    c(r"---[`frm()` refuses the row instead, and names the simplex
      classes, as `set_prior()` did.]---",
      r"---[`frm()` refuses the row instead, and names the simplex
      classes, as `set_prior()` did. `validate_prior()` refuses it
      too.]---"),
    # the install line is a code span longer than a line; the README
    # has it
    c(r"---[To fit a multivariate normal mixture, install frmtmb.latent
      (`remotes::install_github("aforren1/frmtmb", subdir =
      "extensions/frmtmb.latent")`) and run `library(frmtmb.latent)`.]---",
      r"---[To fit a multivariate normal mixture, install frmtmb.latent
      and run `library(frmtmb.latent)`.]---"),
    c(r"---[The family, its arguments and its fits did not change: on the
      59 fits of its tests, the log-likelihood, the estimates, the
      standard errors, `predict()`, `fitted()`, `mixture_probs()` and
      `simulate()` are identical to the last bit on the same
      machine.]---",
      r"---[The family, its arguments and its fits did not change. On the
      59 fits of its tests, these are identical to the last bit on the
      same machine: the log-likelihood, the estimates, the standard
      errors, `predict()`, `fitted()`, `mixture_probs()` and
      `simulate()`.]---"),
    c(r"---[The fields are the same: on the models of the moved tests the
      log-likelihood, every estimate and standard error, the fitted
      values with and without an unseen location, the `ranef()`
      conditional standard deviations, `confint_varcorr()` and
      `simulate()` are identical to frmtmb 0.73.0's, and
      `frmtmb.sample::frm_sample()` gives the same draws from the same
      seed.]---",
      r"---[The fields are the same. On the models of the moved tests,
      these are identical to frmtmb 0.73.0's: the log-likelihood, every
      estimate and standard error, the fitted values with and without an
      unseen location, the `ranef()` conditional standard deviations,
      `confint_varcorr()` and `simulate()`.
      `frmtmb.sample::frm_sample()` gives the same draws from the same
      seed.]---"),
    c(r"---[A declared simplex is read in the chart of a `mo()` simplex,
      searched off a pole and off a face the likelihood rises from, drawn
      by frmtmb.sample in the softmax chart under brms's `dirichlet(1)`,
      named as brms names its draws (`sbhaz[k]`) and given a prior class
      of its own name in `get_prior()` and `set_prior()`.]---",
      r"---[A declared simplex is read in the chart of a `mo()` simplex,
      and searched off a pole and off a face the likelihood rises from.
      frmtmb.sample draws it in the softmax chart under brms's
      `dirichlet(1)`, and names its draws as brms does (`sbhaz[k]`). It
      has a prior class of its own name in `get_prior()` and
      `set_prior()`.]---"),
    c(r"---[`frm()` refuses a declaration that names no extra parameter
      the family returns, or one with no coordinate, and names it;
      before, such a typo ended in "not enough degrees of freedom to
      define contrasts".]---",
      r"---[`frm()` refuses, by name, a declaration that names no extra
      parameter the family returns, or one with no coordinate. Before,
      such a typo ended in "not enough degrees of freedom to define
      contrasts".]---"),
    c(r"---[A saved fit whose family was built by an extension package
      that this session cannot load now stops, at the first method that
      needs the family's own functions (`summary()`, `vcov()`,
      `confint()`, `fitted()`, `residuals()`, `simulate()`, `update()`,
      `frm_curve()` and the rest), with "The cox family of this model is
      from the extension package frmtmb.survival, which is not installed
      ...".]---",
      r"---[Some saved fits have a family that an extension package built.
      When this session cannot load that package, such a fit now stops at
      the first method that needs the family's own functions. These
      include `summary()`, `vcov()`, `confint()`, `fitted()`,
      `residuals()`, `simulate()`, `update()` and `frm_curve()`. The
      message names the
      package: "The cox family of this model is from the extension
      package frmtmb.survival, which is not installed ...".]---"),
    c(r"---[`us_sigma()` builds the covariance matrix that the parameter
      vector of a `us()` block encodes, and it is AD-safe, so a family in
      another package can carry a free covariance matrix as an extra
      parameter; `us_theta()` gives the parameter vector of a covariance
      matrix, for a start.]---",
      r"---[`us_sigma()` builds the covariance matrix that the parameter
      vector of a `us()` block encodes. It is AD-safe, so a family in
      another package can carry a free covariance matrix as an extra
      parameter. `us_theta()` gives the parameter vector of a covariance
      matrix, for a start.]---"),
    c(r"---[a bar structure written `name(terms | g)`, or a field written
      as a predictor special such as `car(M, gr = g)`, whose precision or
      covariance may depend on parameters and on fixed data the special
      names.]---",
      r"---[a bar structure written `name(terms | g)`, or a field written
      as a predictor special such as `car(M, gr = g)`. Its precision or
      covariance may depend on parameters and on fixed data that the
      special names.]---"),
    c(r"---[One entry carries everything the objective, the standard
      errors and their checks, `VarCorr()`, `confint_varcorr()`,
      `ranef()`, `predict()` with and without new levels, `simulate()`,
      `conditional_effects()`, `vcov_cluster()` and frmtmb.sample read,
      and the built-in structures are entries of the same shape.]---",
      r"---[One entry carries everything that frmtmb and frmtmb.sample
      read: the objective, the standard errors and their checks,
      `VarCorr()`, `confint_varcorr()`, `ranef()`, `predict()` with and
      without new levels, `simulate()`, `conditional_effects()` and
      `vcov_cluster()`. The built-in structures are entries of the same
      shape.]---"),
    c(r"---[A name the rest of frmtmb, brms or another extension owns, or a
      function formulas call as a term, is refused (`poly`, `ns`, `log`,
      `gaussian`, `dec`, ...), and a name another package registered is
      replaced only with `replace = TRUE`.]---",
      r"---[The call refuses a name that the rest of frmtmb, brms or
      another extension owns, such as `gaussian` or `dec`. It also
      refuses a function that formulas call as a term, such as `poly`,
      `ns` or `log`. A name that another package registered is replaced
      only with `replace = TRUE`.]---"),
    c(r"---[`frm_compat()` resolves the rules that packages contribute in
      the order of the packages' names, so its answer no longer depends
      on the order of the `library()` calls.]---",
      # the release review, finding 7: the order of the pairs still
      # moves with the load order; their status and note do not
      r"---[`frm_compat()` resolves the rules that packages contribute in
      the order of the packages' names. So the status and note it gives
      a pair no longer depend on the order of the `library()`
      calls.]---"),
    # cox() is frmtmb.survival's after the merge
    c(r"---[`bhaz()` (the basis is set in `cox(df =, degree =,
      intercept =)`, and]---",
      r"---[`bhaz()` (the basis is set in `cox(df =, degree =,
      intercept =)` of frmtmb.survival, and]---"),
    c(r"---[Every family, addition term and formula special of brms
      2.23.0, in its plainest spelling and in the others brms documents,
      now resolves in core, resolves in a named extension, or has such a
      refusal; before, `custom`, `bhaz()`, `dec()` and bare `gr(g)` did
      not.]---",
      r"---[Every family, addition term and formula special of brms
      2.23.0 now resolves in core, resolves in a named extension, or has
      such a refusal. This holds for its plainest spelling and for the
      others that brms documents. Before, `custom`, `bhaz()`, `dec()` and
      bare `gr(g)` did not.]---"),
    c(r"---[A gaussian model whose field has one row per level, such as
      `car()` over areas with one observation each, no longer warns that
      its variance is confounded with the residual sd.]---",
      r"---[A gaussian model can have a field with one row per level, such
      as `car()` over areas with one observation each. Such a model no
      longer warns that its variance is confounded with the residual
      sd.]---"),
    # the release review, finding 8: 86 was the count with
    # frmtmb.spline loaded on 0.73.0; the claim is per pair
    c(r"---[`frm_compat()` reports the same 86 `cox` pairs, status and note
      alike.]---",
      r"---[Every `cox` pair of `frm_compat()` keeps its status and
      note.]---"),
    # the release review, finding 3: a cox() fit an older frmtmb saved
    # is refused by name now (old_core_family_home(), R/extension-home.R)
    c(r"---[and `get_prior()` gives the same table.]---",
      r"---[and `get_prior()` gives the same table. A `cox()` fit saved by
      frmtmb 0.73.0 or earlier must be refitted: frmtmb refuses it by
      name and says to refit it with frmtmb.survival loaded.]---")
  ),
  "extensions/frmtmb.latent/NEWS.md" = list(
    c(r"---[It fits mixtures of multivariate normal distributions to an n
      x D matrix response, with mclust's eight covariance models, class
      means that take covariates, and covariate-dependent mixing
      weights.]---",
      r"---[It fits mixtures of multivariate normal distributions to an n
      x D matrix response. It has mclust's eight covariance models, class
      means that take covariates, and covariate-dependent mixing
      weights.]---"),
    c(r"---[Nothing about the family changed in the move: on the 59 fits
      of its tests, the log-likelihood, the estimates, the standard
      errors, `predict()`, `fitted()`, `mixture_probs()` and
      `simulate()` are identical to the last bit to the ones frmtmb
      gave, on the same machine.]---",
      r"---[Nothing about the family changed in the move. On the 59 fits
      of its tests, these are identical to the last bit to the ones
      frmtmb gave, on the same machine: the log-likelihood, the
      estimates, the standard errors, `predict()`, `fitted()`,
      `mixture_probs()` and `simulate()`.]---"),
    # the release number the lane could not write
    c(r"---[`hmm()` reads core's one list of field terms (`car()` and
      `spde()`, now parsed through core's covariance-structure registry)
      when it decides whether a distributional parameter is constant, so
      it needs the frmtmb release that carries
      `frmtmb_register_covstruct()`.]---",
      r"---[`hmm()` reads core's one list of field terms when it decides
      whether a distributional parameter is constant. That list holds
      `car()` and `spde()`, which core now parses through its
      covariance-structure registry. So the package needs frmtmb 0.74.0,
      which has `frmtmb_register_covstruct()`.]---")
  ),
  "extensions/frmtmb.survival/NEWS.md" = list(
    c(r"---[The curve functions that read a fit, `frm_curve()`,
      `frm_curve_deriv()` and `frm_curve_feature()`, are in frmtmb now,
      and still refuse a fit whose cumulative hazard is not monotone,
      through the family's `post$curve_check` hook.]---",
      r"---[The curve functions that read a fit, `frm_curve()`,
      `frm_curve_deriv()` and `frm_curve_feature()`, are in frmtmb now.
      They still refuse a fit whose cumulative hazard is not monotone,
      through the family's `post$curve_check` hook.]---")
  ),
  "extensions/frmtmb.spatial/NEWS.md" = list(
    c(r"---[New package. `car(M, gr = g, type = )`, brms's conditional
      autoregressive field in its four types (`"escar"`, `"esicar"`,
      `"icar"`, `"bym2"`), and `spde(fem, gr = node)`, a Matern field
      over a finite-element mesh, moved here from frmtmb.]---",
      r"---[New package. Two fields moved here from frmtmb.
      `car(M, gr = g, type = )` is brms's conditional autoregressive
      field in its four types (`"escar"`, `"esicar"`, `"icar"`,
      `"bym2"`). `spde(fem, gr = node)` is a Matern field over a
      finite-element mesh.]---"),
    c(r"---[Fits are unchanged: on the moved tests' models the
      log-likelihood, every estimate and standard error, `fitted()` with
      and without an unseen location, `ranef()` conditional standard
      deviations, `confint_varcorr()` and `simulate()` are identical to
      frmtmb 0.73.0's.]---",
      r"---[Fits are unchanged. On the moved tests' models, these are
      identical to frmtmb 0.73.0's: the log-likelihood, every estimate
      and standard error, `fitted()` with and without an unseen location,
      `ranef()` conditional standard deviations, `confint_varcorr()` and
      `simulate()`.]---")
  )
)
