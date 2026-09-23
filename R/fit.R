#' Fit a model
#'
#' Fits a model specified with [bf()] by maximum likelihood, using the
#' Laplace approximation for random effects through RTMB.
#'
#' @param formula A `frmtmb_formula` from [bf()] (with a family attached
#'   via `+`), or a plain formula combined with the `family` argument.
#'   With neither, the family is `gaussian()`.
#' @param data A data frame. A `tibble`, a `data.table`, or a plain
#'   named list of equal-length columns is accepted as well, since each
#'   reaches [stats::model.frame()] unchanged. A matrix column is a
#'   supported model variable and enters the design as its own block of
#'   columns (this is how a functional predictor or a `cbind()` term is
#'   written). A list column is not a model variable and is refused by
#'   `model.frame()` if the formula names one, though it may sit unused
#'   in `data`.
#' @param data2 A named list of objects that are not columns of `data`:
#'   the adjacency matrix of `car()`, the mesh triple of `spde()`, and
#'   the matrices of `gr(prec = )`, `gr(cov = )` and `equalto()`. This
#'   is brms's `data2` argument, with one deliberate extension: brms
#'   accepts a bare name only, while frmtmb also evaluates compound
#'   expressions with `data2` in front of the data mask, so
#'   `gr(g, cov = solve(Q))` finds `Q` there.
#'   Each structural expression resolves from `data2` first, then
#'   `data`, then the formula environment; that last step is what a
#'   model written before `data2` relied on, so old code keeps working.
#'   Prefer `data2`: its objects are stored on the fit, so `saveRDS()`
#'   and a later `refit()`, `influence()` or `update()` in a fresh
#'   session do not need the calling environment to still exist.
#' @param family A family: a `frmtmb_family`, a [stats::family] object,
#'   a family constructor with or without its parentheses (for example
#'   `gaussian()`, `poisson`, `cumulative`), or a family name as a
#'   string. It overrides a family already attached to `formula`; in a
#'   multivariate model it fills only the responses that have none.
#'   The default, `NULL`, means `gaussian()` - the brms, `lme4` and
#'   `glmmTMB` convention - so `frm(y ~ x, data = d)` is a linear
#'   model.
#' @param REML If `TRUE`, integrate the `mu` fixed effects out of the
#'   likelihood along with the random effects (restricted maximum
#'   likelihood). For a gaussian model the Laplace approximation is
#'   exact for that integral, so this IS classical
#'   Patterson-Thompson REML. For any other family no error-contrast
#'   derivation exists; what is computed is the Laplace-approximated
#'   integrated likelihood for the variance parameters under a flat
#'   prior on the fixed effects, which agrees to first order with the
#'   Cox-Reid adjusted profile likelihood (the general-model REML
#'   analogue). glmmTMB computes the same quantity. It usually
#'   reduces the finite-sample downward bias of variance components,
#'   which is REML's purpose, but it is an approximation stacked on
#'   an approximation, not an exactness result. The usual restriction
#'   carries over: [anova()] compares REML fits only when their
#'   fixed-effect designs agree.
#'
#'   Distributional-parameter coefficients (a `sigma ~ z` model's, for
#'   example) are deliberately NOT integrated: they belong to the
#'   variance-parameter set, exactly as a `varExp()` coefficient does
#'   in `nlme::gls(method = "REML")`, and the two agree to optimizer
#'   precision for gaussian models. Random effects appearing in a
#'   distributional parameter's own formula are integrated like any
#'   other latent. `mgcv`'s `method = "REML"` for a location-scale
#'   family such as `gaulss()` integrates the coefficients of every
#'   linear predictor, so for a smooth in `sigma` the two REML criteria
#'   are different and give different curves and smoothing parameters.
#' @param start Optional named list of starting values, one entry per
#'   parameter component (`beta`, `betad`, `theta`, ...).
#'   [par_template()] returns that list already filled with the
#'   defaults, for a fitted model or for a formula and data alone; edit
#'   it and pass it back.
#'
#'   A component is read by NAME when its vector has names, so
#'   `start = list(beta = c("(Intercept)" = 5000))` moves that one
#'   coefficient and leaves the rest at their defaults; parentheses may
#'   be dropped, as everywhere a parameter is named. A vector with no
#'   names is positional and must be full length. Naming some entries
#'   of one vector and not others is refused.
#' @param control A list from [frmtmb_control()].
#' @param se If `TRUE`, run [RTMB::sdreport()] at fit time. The default
#'   (`FALSE`) defers it until standard errors are first needed
#'   (`summary`, `vcov`, `confint`, `frm_linpred(se.fit = TRUE)`), which cuts
#'   roughly a quarter off fit time in fit-and-predict or bootstrap
#'   loops. The deferred report is cached, so nothing is computed twice.
#' @param na.action How to handle missing values, as in [stats::lm()]
#'   (default [stats::na.omit]). Rows dropped for missingness are
#'   reported in a message; wrap the call in `suppressMessages()` to
#'   silence it.
#' @param prior Optional [set_prior()] specification. This makes the
#'   fit MAP / regularized ML (glmmTMB's `priors=` in spirit): useful
#'   for stabilizing singular variance components or separating
#'   binomials. The reported logLik/AIC then include the prior terms
#'   and are penalized quantities, and `anova()` comparisons across
#'   different priors are meaningless. A `brmsprior` object built by
#'   brms's own `prior()` is translated row by row, so priors copied
#'   out of a brms script work whichever package's `prior()` was in
#'   scope. The argument takes brms's spelling, `prior`; the `priors`
#'   of releases before 0.43 is gone rather than aliased, and a call
#'   still using it fails as an unused argument.
#'
#'   This is also where HARD BOUNDS are written, as `lb`/`ub` on a
#'   specification that may carry no distribution at all:
#'   `set_prior("", nlpar = "guess", lb = 0, ub = 1)` is a box
#'   constraint and nothing else. The `lower`/`upper` arguments of
#'   releases before 0.49 are gone rather than aliased, and a call still
#'   using them fails as an unused argument; every outer parameter they
#'   could reach has a `set_prior()` class, down to a single internal
#'   covariance parameter (`class = "theta", coef = "thetaac_1"`). See
#'   [set_prior()]'s Hard bounds section.
#'
#'   *The spelling is brms's; the SEMANTICS are not.* A prior here is a
#'   penalty on the likelihood and the answer is one mode. It is not a
#'   posterior, and no interval this fit reports is a credible
#'   interval. The two land close where the data dominates: a kidney
#'   frailty model given brms's own priors returns `sd(patient)` 0.38
#'   against brms's posterior mean of 0.40. That is a measurement on
#'   one model, though, not a property of the translation. Use
#'   `frmtmb.sample::frm_sample()` when the posterior is the answer.
#'
#'   *A prior with a location places a nonlinear start.* brms uses its
#'   priors to place the sampler, and `frm()` now reads them the same
#'   way, for nonlinear parameters only: where `start` does not set it,
#'   a nonlinear coefficient begins at the location of a `normal()`,
#'   `student_t()` or `cauchy()` prior on it, and a message names what
#'   was placed. This is what lets the insurance-loss growth curve fit
#'   from its brms priors with no `start` at all. Every other parameter
#'   keeps the start it always had.
#'
#'   Without such a prior a nonlinear model still needs `start`.
#'   `frm()` optimizes rather than samples, so it evaluates the
#'   objective AT the starting values, and a nonlinear body is rarely
#'   defined at zero. [par_template()] shows the names and the values
#'   in force.
#' @param quadrature If `TRUE`, marginalize each scalar random effect by
#'   adaptive Gauss-Kronrod quadrature instead of the Laplace
#'   approximation (the `glmer(nAGQ = k)` analogue; matches it in
#'   tests). Worth it for Bernoulli responses with small clusters,
#'   where Laplace biases variance components. Scalar random-intercept
#'   models only, and not with `mi()`, `trunc()`, `REML = TRUE`, or
#'   `frmtmb_control(profile = TRUE)`. A plain Laplace fit runs first
#'   and the quadrature tape is built at its optimum: the
#'   Gauss-Kronrod rescaling is fixed when the tape is built, so the
#'   starting point decides whether the marginalized objective is
#'   finite at all. That fit also supplies the conditional modes, which
#'   the marginalized objective no longer carries, so `ranef()`,
#'   `fitted()` and `predict()` work as usual.
#' @param importance Number of importance draws per group, or `0`
#'   (default) for none. A positive count replaces the Laplace
#'   approximation with an importance-sampling correction of it
#'   (Skaug and Fournier 2006, the ADMB `-is` option): each group's
#'   marginal likelihood is re-estimated by reweighting draws from the
#'   Laplace Gaussian, which is unbiased for the exact integral rather
#'   than accurate to `O(n^-1)`. Like `quadrature = TRUE` it is worth
#'   it for binary responses in small clusters, and unlike it, the
#'   blocks may have any dimension: this is the correction for a
#'   correlated random slope, where quadrature cannot go.
#'
#'   The cost is one tape of `N` stacked copies of the likelihood, so
#'   time and memory grow linearly in the draw count. An odd count is
#'   rounded up to the next even one, because the draws are taken in
#'   antithetic pairs (measured on this package's own probe design,
#'   the pairing is worth two to three times the draws it costs).
#'
#'   Scope, with everything else refused by name: any number of
#'   random-effect blocks, of any dimension, provided they share ONE
#'   grouping factor and one set of levels, in any covariance
#'   structure that is Gaussian within a level and independent
#'   between levels; one response; a family with a rowwise density.
#'   Not with `quadrature`, `REML = TRUE`,
#'   `frmtmb_control(profile = TRUE)`, `mi()`, a residual correlation
#'   term, `rescor`, a nonlinear predictor, or `cs()`. `trunc()` and
#'   `cens()` are supported: the draws sit near the conditional mode,
#'   which is why the truncation normalizer that underflows at the
#'   Gauss-Kronrod nodes does not underflow here.
#'
#'   Distributional regression writes several blocks by construction:
#'   `(1 | g)` in `mu` and `(1 | g)` in `sigma` are two blocks over
#'   `g`, and a level's coefficients from both are drawn together
#'   from one joint proposal. Blocks over DIFFERENT grouping factors,
#'   crossed or nested, are refused.
#'
#'   The fit records the draw count, the seed, the rounds taken, the
#'   per-group effective sample sizes and the Monte Carlo standard
#'   error of the corrected log-likelihood in `fit$importance`;
#'   `summary()` prints them. `logLik()`, `AIC()`, `vcov()` and
#'   `confint()` all come from the corrected objective, with two
#'   caveats worth knowing. The reported estimate is the optimum of the
#'   previous round's proposal, so the covariance is a Hessian at a
#'   point that is stationary only to `fit$importance$grad`. And
#'   `confint(method = "profile")` walks the FROZEN proposal instead of
#'   refitting, so it warns when a bound lands where the weights no
#'   longer cover the integrand. The seed, the round cap and the
#'   effective-sample-size threshold are [frmtmb_control()] settings.
#' @param dry_run `"spec"` returns the parsed intermediate representation
#'   without touching `data`; `"frame"` returns the assembled design
#'   matrices and parameter template without fitting; `"objective"`
#'   additionally tapes the objective and returns an UNFITTED object
#'   carrying it, which is what `frmtmb.sample::frm_sample()` samples
#'   when it is given a formula rather than a fit. Methods that report a
#'   maximum-likelihood quantity refuse on that object.
#' @param verbose Report fit progress; a shortcut for
#'   `control = frmtmb_control(verbose =)`, whose value wins when both
#'   are given. See [frmtmb_control()] for the levels and the output.
#' @return An object of class `frmtmb_fit`. It is a list, and two of its
#'   elements are read directly often enough to name here. `fit$data2`
#'   is the `data2` list of known matrices. `fit$data` is the MODEL
#'   FRAME, the object [model.frame()] returns: one column per term the
#'   formula names, under the term's own spelling, so a model with
#'   `offset(Age)` carries a literal `offset(Age)` column and a data
#'   column no term uses is absent. brms keeps its VALIDATED RAW DATA
#'   under the same name, so `names(fit$data)` and `ncol(fit$data)`
#'   differ from brms on any model with a transformed term, and
#'   `newdata = fit$data` is not the brms idiom it looks like. The
#'   element exists because brms-shaped code reads `fit$data` and used
#'   to reach the `data2` list through `$`'s partial matching; the model
#'   frame is the object this package already had. It is the same object
#'   in memory that `fit$frame[["data_frame"]]` holds, but R's
#'   serializer does not deduplicate a shared value, so a SAVED fit
#'   carries the frame twice: 5.3% more raw bytes and 14.0% more gzipped
#'   on a 20,000-row fit, 0.16% and 1.49% on a 240-row one
#'   (`dev/adefects-log/p1-dollar.txt`). Everything else the fit carries
#'   has an accessor, and the accessor is the supported route.
#'
#' @section What a nonlinear body sees:
#' The body of a nonlinear formula - the response formula under
#' `nl = TRUE`, or an [nlf()] body - is evaluated as a small language
#' for the AD tape, not as ordinary R. A function it CALLS is looked up
#' in RTMB first, so a bare `pnorm()`, `qgamma()`, `plogis()`,
#' `besselK()` or `matrix()` in a body is RTMB's tape-capable version.
#' The `RTMB::` prefix is no longer needed, and writing it changes
#' nothing: both spellings tape to the same function and give the same
#' fit.
#'
#' RTMB WINS. It wins over the search path and over a function of the
#' same name that you defined yourself, because the body is a language
#' with its own vocabulary rather than a piece of your session. To reach
#' a different `pnorm` inside a body, qualify it: `stats::pnorm()` is
#' the escape hatch, and so is `myPkg::pnorm()`.
#'
#' Nothing else changes. Only the names RTMB replaces are shadowed, and
#' only where the body calls them, so a helper function of your own, an
#' object the body reads, and every name that is not in RTMB still
#' resolve in the formula's environment exactly as before. A name the
#' body READS is never shadowed either, so a data frame called `df` or a
#' vector called `order` keeps its meaning.
#'
#' The rule stops at the end-user surface. A [frmtmb_family()]
#' log-density and a [frmtmb_structure()] log-likelihood are ordinary R
#' functions written by an extension author, and they keep resolving
#' lexically: qualify with `RTMB::` there, and see [frmtmb_ad_overload()]
#' for the rest of that contract.
#'
#' @section Monotonic effects:
#' `mo(x)` fits an ordinal predictor without assuming its categories are
#' equally spaced. `x` is an ordered factor or a non-negative integer
#' vector with at least three categories, coded `0..D`. The term
#' contributes `b * D * sum(zeta[1:k])` at category `k`, where `zeta` is
#' a simplex of `D` non-negative steps summing to one. `b` is therefore
#' the AVERAGE step, on the coefficient scale of any other predictor,
#' and `zeta` says how that total is distributed over the categories.
#' `b` appears in `fixef()` under the term's label (`mox`, or `mox:z`
#' for an interaction); the simplex is held as its `D - 1` free softmax
#' coordinates in a `zeta<j>` component of [par_template()], and
#' `summary()` prints those coordinates rather than the simplex.
#'
#' Every monotonic TERM gets its own simplex. `y ~ mo(x) * z` fits two:
#' one shape for the main effect and one for the interaction, because
#' the two terms describe different shapes and there is no reason for
#' the interaction to bend the way the main effect does. The simplexes
#' are numbered in the order [brms::brm()] enumerates its special terms,
#' which is `stats::terms()` order and so lists every main effect before
#' any interaction. The j-th monotonic term of a linear predictor is
#' therefore the one brms calls `simo_<j>` and names `<label><1>` in the
#' `simo` rows of its `get_prior()`.
#'
#' The `zeta<j>` NUMBER is not that j in general. Simplexes continue the
#' numbering of whatever parameters the family contributed first, so
#' `zeta<j>` is brms's `simo_<j>` only for a family that declares none:
#' an ordinal fit spends slot 1 on its thresholds and a `cox()` fit on
#' its baseline hazard, and in both the first monotonic simplex is
#' `zeta2`. Read the terms in order rather than parsing the number.
#'
#' frmtmb has no prior class for a simplex; brms's flat `dirichlet(1)`
#' is the one prior it cannot be told to drop, and it contributes only a
#' constant, so the two log-densities agree up to `lgamma(D)` per
#' simplex.
#'
#' `mo()` may be crossed with a single numeric multiplier and no more:
#' `mo(x) * z` and `mo(x):z` are fitted, `mo(x) * f` for a factor `f` is
#' refused (the simplex carries one coefficient and a contrast expansion
#' has no column to go in), and `mo(x):mo(w)` is refused outright.
#'
#' @section The Laplace approximation, and how to check it:
#' Random effects are integrated out by the Laplace approximation,
#' which assumes the integrand is close to Gaussian around the
#' conditional mode; the Wald intervals `confint()` reports assume the
#' log-likelihood is close to quadratic at the optimum. Both degrade in
#' the same places: variance components estimated from few groups, and
#' binary data in small clusters.
#'
#' Three remedies are in this package. `confint(method = "profile")`
#' replaces the quadratic assumption with a profile likelihood,
#' [frm_bootstrap()] with a resampling distribution, and
#' `quadrature = TRUE` replaces the Laplace approximation itself with
#' adaptive quadrature (the test suite checks that fit against
#' `lme4::glmer(nAGQ = 25)` and GLMMadaptive in exactly the regime where
#' the Laplace fit is biased).
#'
#' To MEASURE the violation rather than route around it, install the
#' `frmtmb.sample` package and call `frmtmb.sample::check_laplace()`. It
#' runs NUTS on this very objective and reports how far the posterior
#' mean sits from the maximum-likelihood estimate in posterior standard
#' deviations, and the ratio of the posterior standard deviation to the
#' Wald standard error. `vignette("diagnostics")` works through it.
#'
#' @srrstats {RE1.4} The assumptions the fit rests on are documented,
#'   and the consequences of violating them are both documented and
#'   testable. The section above names the two approximations, the
#'   regimes where each degrades, the three remedies inside this
#'   package (`confint(method = "profile")`, [frm_bootstrap()] and
#'   `quadrature = TRUE`), and the direct measurement in the companion
#'   sampling package, which runs NUTS on the same objective and
#'   reports the shift of the posterior mean from the
#'   maximum-likelihood estimate in posterior standard deviations and
#'   the ratio of the posterior standard deviation to the Wald standard
#'   error. The quadrature fit is checked against
#'   `lme4::glmer(nAGQ = 25)` and GLMMadaptive in the regime where the
#'   Laplace fit is biased; `vignette("diagnostics")` works through the
#'   whole question.
#' @srrstats {RE1.0} Models are specified through a formula interface:
#'   [bf()] builds a `frmtmb_formula` from one or more R formulas and a
#'   family attaches with `+`. A plain formula plus `family =` is also
#'   accepted. There is no matrix-only entry point.
#' @srrstats {G2.4,G2.4a,G2.4b,G2.4c,G2.4e} Type conversion during frame
#'   assembly is explicit, never implicit. The response is converted with
#'   `as.numeric()`, or `storage.mode(y) <- "double"` for a matrix
#'   response; an ordinal factor response is converted from factor with
#'   `as.numeric()` and a two-level binomial factor with
#'   `as.numeric(y) - 1`; `mo()` category codes use `as.integer()`;
#'   grouping factors use `as.integer()` for the level index, and a
#'   grouping value is matched to its level through `as.character()`, so
#'   an integer, a factor and a character grouping column index alike;
#'   addition terms other than `cens()` use `as.numeric()`.
#' @srrstats {G2.5} Where a factor input is expected, the expected kind is
#'   checked and documented. `mo()` requires an ordered factor and errors
#'   otherwise ("mo(): factor variables must be ordered factors"). An
#'   ordinal family takes level order as category order, and refuses an
#'   unordered factor response, naming its level order: that order is
#'   alphabetical unless the user set it, so the model could differ from
#'   the one intended. brms refuses it too. A factor response for a
#'   non-ordinal, non-categorical, non-binomial family
#'   is refused. The compatibility registry states the requirement in
#'   prose and it is rendered in `vignette("compatibility")`.
#' @srrstats {G2.6} One-dimensional responses are pre-processed to a plain
#'   numeric vector regardless of the class they arrive in: a factor, a
#'   one-column matrix, a `scale()`d matrix carrying attributes, or a bare
#'   vector all reach the objective as `as.numeric(as.vector(y))`.
#' @srrstats {G2.8} Pre-processing funnels every input into two internal
#'   classes before any analytic code runs: `parse_spec()` produces a
#'   data-free `frmtmb_spec`, and `assemble_frame()` produces a
#'   `frmtmb_frame` holding the design matrices and the parameter
#'   template. Every sub-function downstream of assembly sees only those
#'   two classes. Both are reachable for inspection through `dry_run`.
#' @srrstats {G2.13} Missing data is checked during pre-processing, before
#'   anything is passed to the optimizer. Rows are removed by `na.action`
#'   inside `stats::model.frame()`; assembly then errors if any missing
#'   value remains in a model variable other than one declared `mi()` or
#'   modelled by a family that declares `keep_na`, where the `NA` is the
#'   estimand and is carried deliberately, and errors if no complete
#'   observation is left.
#' @srrstats {G2.14,G2.14a} `na.action` lets the user choose how missing
#'   data is handled, following [stats::lm()]. `stats::na.fail` errors on
#'   missing data, `stats::na.omit` (the default) and `stats::na.exclude`
#'   drop the affected rows, and `na.exclude` pads `fitted()`,
#'   `residuals()`, `predict()`, and `simulate()` back to the input
#'   length with `NA` in the original positions. A column declared
#'   `mi()`, and the response of a family that declares `keep_na`,
#'   bypass `na.action` by design: their `NA`s are what the model
#'   estimates, so assembly keeps those rows and applies the requested
#'   action to the rest.
#' @srrstats {G2.14b} Rows dropped for missingness are reported, not
#'   dropped silently: frame assembly emits one `message()` per fit
#'   giving the number of rows removed. It is a message, not a warning,
#'   so `suppressMessages()` silences it for callers that ask for
#'   `na.omit` deliberately, and `na.action()` on the fit still names the
#'   rows.
#' @srrstats {G2.15} No function assumes non-missingness by inheriting a
#'   default. The invariant is established once at the boundary: assembly
#'   errors unless every model variable is free of missing values, so
#'   downstream arithmetic operates on complete data by construction
#'   rather than by defensive `na.rm` flags that would silently change
#'   the estimand. The exceptions are explicit and modelled, not
#'   inherited: an `mi()` column and a `keep_na` family's response keep
#'   their `NA`s, and the code that reads one writes `na.rm = TRUE` at the
#'   call site, where the reader can see which quantity it changes.
#' @srrstats {G2.16} Undefined values are handled separately from missing
#'   ones. The response check is explicitly written as
#'   `any(!is.finite(y) & !is.na(y))`, so `Inf` and `-Inf` are rejected
#'   with their own message, which no `na.action` setting can silence,
#'   while anything `is.na()` calls missing is left to `na.action`.
#'   `NaN` counts as missing to `is.na()`, so it takes the `NA` route
#'   and is dropped (or refused by `na.fail`) rather than reaching the
#'   finiteness check. That is `stats::lm()`'s own division, and it is
#'   the reason the two are described together here: the split is
#'   between values that are undefined and values that are absent, not
#'   between the `Inf` and `NaN` spellings.
#' @srrstats {RE2.1} The processing of missing values is controlled by an
#'   explicit parameter (`na.action`), and `NA`/`NaN` are distinguished
#'   from `Inf` as described under G2.16.
#' @srrstats {RE2.4,RE2.4a} Perfect collinearity among predictors is
#'   detected during pre-processing. Each parametric fixed-effect design
#'   is factorized with `qr()`; when the rank is short of the column
#'   count the aliased columns are named in a `message()` and dropped,
#'   and the null space of the design is stored on the fit so that
#'   `predict()` can refuse rows of `newdata` that are not estimable. A
#'   cheap sparse singular-value screen gates the dense check so that the
#'   sparse and dense backends drop the same columns.
#' @srrstats {RE3.0} Models that fail to converge raise warnings: a
#'   nonzero optimizer status (with the optimizer's own message and, for
#'   a nonlinear model, a hint that `start` was not set), a maximum
#'   absolute gradient above `grad_tol`, and, once the standard-error
#'   machinery has run (at fit time under `se = TRUE`, otherwise on the
#'   first `vcov()`, `summary()` or [diagnose()] call), a Hessian that is
#'   not positive definite or, failing that, non-finite standard errors.
#'   The last two are exclusive: a Hessian that is not positive definite
#'   explains the standard errors, so only the first of the pair is
#'   raised.
#' @srrstats {RE3.1} Those diagnostics are `warning()` conditions, so
#'   `suppressWarnings()` silences them, and the returned object still
#'   carries enough to identify the failure: `fit$opt$convergence` and
#'   `fit$opt$message` hold the optimizer verdict, and [diagnose()]
#'   recomputes the gradient, the Hessian verdict, separation, singular
#'   variance components, and predictor scaling on demand.
#' @srrstats {RE4.0} The return value is a model object of class
#'   `frmtmb_fit`, with the standard `stats` accessor methods defined
#'   for it.
#' @srrstats {RE4.1} `dry_run` generates the model object without fitting
#'   it: `"spec"` returns the parsed representation without touching
#'   `data`, and `"frame"` returns the assembled design matrices and the
#'   parameter template. Both are useful for batch setup and for
#'   inspecting what a formula compiled to.
#' @srrstats {RE4.4} The model specification is recoverable as a formula
#'   through `formula()`, and the original call through `fit$call`.
#' @srrstats {RE4.5} The number of observations used is returned by
#'   `nobs()`, and the rows dropped by `na.action` by `na.action()`.
#' @srrstats {RE4.8} The response and its metadata are retained on the
#'   fit: `fit$frame$y` holds the response per response name,
#'   `fit$frame$y_levels` the original factor levels for ordinal and
#'   categorical responses, and `fit$frame$aterm_values` the addition
#'   terms (`weights`, `trials`, `cens`, `trunc_lb` and `trunc_ub`,
#'   `se`, `vint` and `vreal`, `mi`). An offset belongs to one linear
#'   predictor rather than to the response, so it is kept there, as
#'   `fit$frame$linpreds[[k]]$offset`.
#'   `model.frame()` returns the stored model frame with its row names.
#' @srrstats {RE4.17} `print()` on a `frmtmb_fit` summarizes the model
#'   input (family, formula, grouping structure) and the fitted
#'   coefficients.
#' @srrstats {RE4.18} A distinct `summary()` method is implemented, and it
#'   is computationally non-trivial: with the default `se = FALSE` it
#'   triggers the deferred `RTMB::sdreport()` that produces the standard
#'   errors, then forms z and p values, variance components, group
#'   counts, information criteria, residual correlations, and smooth
#'   effective degrees of freedom. The report is cached on the fit, so a
#'   second `summary()` is cheap.
#'
#' @examples
#' \dontrun{
#' data(sleepstudy, package = "lme4")
#' fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
#'               data = sleepstudy)
#' summary(fit)
#' }
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' summary(fit)
#' fixef(fit)
#' VarCorr(fit)
#'
#' # distributional regression: model sigma too
#' fit2 <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
#' anova(fit, fit2)
#'
#' # a nonlinear model: discover the names, edit, fit. A nonlinear body
#' # is undefined at the zero start, so `start` is not optional here.
#' nd <- data.frame(t = rep(seq(0, 10, length.out = 25), 4))
#' nd$y <- 8 * (1 - exp(-nd$t / 3)) + rnorm(nrow(nd), 0, 0.3)
#' nf <- bf(y ~ asym * (1 - exp(-t / lrc)), asym ~ 1, lrc ~ 1,
#'          nl = TRUE) + gaussian()
#' st <- par_template(nf, data = nd)
#' st
#' st$beta[["asym_(Intercept)"]] <- 5
#' st$beta[["lrc_(Intercept)"]] <- 1
#' fixef(frm(nf, data = nd, start = st))
#'
#' # or let a located prior place the same starts, brms-style
#' fixef(frm(nf, data = nd,
#'           prior = prior(normal(5, 5), nlpar = "asym") +
#'             prior(normal(1, 5), nlpar = "lrc")))
#' @export
frm <- function(formula, data, family = NULL, REML = FALSE, start = NULL,
                control = frmtmb_control(), se = FALSE,
                na.action = stats::na.omit,
                prior = NULL, quadrature = FALSE, importance = 0L,
                data2 = list(), dry_run = NULL, verbose = FALSE) {
  cl <- match.call()
  # a brms prior object is translated at the boundary, so nothing
  # downstream sees anything but a frmtmb_priorlist
  prior <- as_priorlist(prior)
  check_prior_slots(prior)
  # Every one of these used to be read through isTRUE() or an empty
  # names() loop, so a flag set by mistake either fitted a DIFFERENT
  # model in silence or died with an obscure downstream error (REML
  # reached as.logical() inside fit_assembled).
  # `verbose` is deliberately not among them: it takes integer
  # levels as well as TRUE/FALSE and documents that it ignores anything
  # else (see verbose_level()).
  check_flag(REML, "REML")
  check_flag(se, "se")
  check_flag(quadrature, "quadrature")
  check_count(importance, "importance", min = 0L)
  if (!is.null(dry_run)) {
    check_string_choice(dry_run, "dry_run", c("spec", "frame", "objective"))
  }
  if (!is.null(start)) {
    check_named_list(start, "start", "start = list(beta = c(0, 1))")
  }
  # a list is not enough: control = list() passes is.list() and then
  # reaches seq_len(NULL) inside the optimizer loop, where the message
  # is "argument must be coercible to non-negative integer"
  if (!is.list(control)) {
    frm_stop("`control` must be a list from frmtmb_control(), not ",
             arg_desc(control), call. = FALSE)
  }
  ctl_need <- c("optimizer", "optCtrl", "restarts", "grad_tol")
  if (!all(ctl_need %in% names(control))) {
    frm_stop("`control` must come from frmtmb_control(); this list is missing ",
             paste(setdiff(ctl_need, names(control)), collapse = ", "),
             call. = FALSE)
  }
  if (!is.null(prior) && !is.list(prior)) {
    frm_stop("`prior` must be a set_prior() specification or a named list ",
             "of prior objects, not ", arg_desc(prior), call. = FALSE)
  }
  if (!is.function(na.action) && !is.character(na.action)) {
    frm_stop("`na.action` must be a function such as stats::na.omit, or its ",
             "name as a string, not ", arg_desc(na.action), call. = FALSE)
  }
  data2 <- validate_data2(data2)
  # frmtmb_control() leaves verbose unset (NULL), so an explicit control
  # value always wins over the frm() shortcut
  control$verbose <- control$verbose %||% verbose
  vb <- verbose_level(control)
  bform <- as_bform(formula, family)
  # a family whose parameter vocabulary is a property of the data
  # (categorical()'s one predictor per observed category) becomes
  # concrete here, before the grammar that names those dpars is parsed
  bform <- resolve_deferred_families(bform, data)

  if (vb) t0 <- vb_now()
  spec <- parse_spec(bform)
  if (vb) vb_stage("parse", t0)
  if (identical(dry_run, "spec")) return(spec)

  if (vb) t0 <- vb_now()
  frame <- assemble_frame(spec, data, na.action = na.action,
                          sparse_x = isTRUE(control$sparse_x),
                          data2 = data2)
  if (vb) vb_stage("frame", t0, vb_frame_detail(frame))
  # A family with a `family_finalize` slot derives itself from the
  # response during assembly, and the fit stores the spec separately
  # from the frame, so family(fit) would otherwise report the family as
  # written rather than the one that was taped.
  spec <- carry_finalized_responses(spec, frame)
  check_re_structure(spec, frame, control)
  suggest_bernoulli(spec, frame)
  if (identical(dry_run, "frame")) return(frame)

  fit_assembled(spec, frame, bform, cl, REML = REML, start = start,
                control = control, se = se, lower = NULL, upper = NULL,
                prior = prior, quadrature = quadrature,
                importance = importance, data2 = data2,
                objective_only = identical(dry_run, "objective"),
                # the user's own call is the one place a placed start is
                # news; refit(), influence and the pre-fit are not
                announce_start = TRUE)
}

# --- verbose progress reporting --------------------------------------
# Stage lines go through message(): suppressMessages() silences them and
# they can never contaminate stdout results. Every call site is guarded
# by an `if (vb)` so a default fit does no clock reads and builds no
# strings. Nothing here reaches inside the tape - per-evaluation
# printing would dominate the fit it is meant to measure.

#' Resolved level: 0 = silent, 1 = stage progress, 2 = also the
#' optimizer's own iteration trace.
#'
#' @noRd
verbose_level <- function(control) {
  v <- control$verbose
  if (is.null(v) || isFALSE(v)) return(0L)
  if (isTRUE(v)) return(1L)
  v <- suppressWarnings(as.integer(v)[1L])
  if (is.na(v) || v < 0L) 0L else v
}

#' The clock every verbose timing reads: elapsed seconds as one number.
#'
#' @noRd
vb_now <- function() proc.time()[["elapsed"]]

#' One progress line with the package prefix, through `message()`.
#'
#' @noRd
vb_say <- function(...) frm_message("frmtmb: ", ...)

#' One timed stage line, shaped `"frmtmb: <stage> [1.23s]: <detail>"`.
#'
#' @noRd
vb_stage <- function(stage, t0, detail = NULL) {
  frm_message("frmtmb: ", stage, " [", sprintf("%.2f", vb_now() - t0), "s]",
              if (is.null(detail)) "" else paste0(": ", detail))
}

#' A count and its noun, with the plural `s` only when the count needs
#' one.
#'
#' @noRd
vb_plural <- function(n, what) {
  paste0(n, " ", what, if (n != 1L) "s")
}

#' The size of an assembled frame in one clause: observations, linear
#' predictors, and random-effect blocks. The detail of the "frame" stage
#' line.
#'
#' @noRd
vb_frame_detail <- function(frame) {
  paste0(frame[["n_obs"]], " obs, ",
         vb_plural(length(frame[["linpreds"]]), "linear predictor"), ", ",
         vb_plural(length(frame[["re_blocks"]]), "random-effect block"))
}

#' First line of a fit: the family and every mode that changes what the
#' optimizer is solving, so a slow log says which problem it is timing.
#'
#' @noRd
vb_fit_detail <- function(spec, REML, control, quadrature, prior,
                          importance = 0L) {
  fams <- vapply(spec$responses, function(r) r$family[["family"]] %||% "?", "")
  opt <- control$optimizer %||% "nlminb"
  if (is.function(opt)) opt <- "custom"
  flags <- c(if (isTRUE(REML)) "REML" else "ML",
             if (isTRUE(control$profile)) "profile",
             if (isTRUE(quadrature)) "quadrature",
             if (importance > 0L) paste0("importance ", importance),
             if (isTRUE(control$autoscale)) "autoscale",
             if (!is.null(prior)) "prior")
  paste0(paste(unique(fams), collapse = " + "), ", ",
         paste(flags, collapse = ", "), ", ", opt)
}

#' An optimizer result in one clause: the objective, plus the status code
#' when the optimizer did not report success and the count of non-finite
#' trials when there were any.
#'
#' @noRd
vb_opt_detail <- function(opt) {
  paste0("objective ", format(opt$objective, digits = 8),
         if (opt$convergence != 0) {
           paste0(", convergence ", opt$convergence)
         },
         # the NA/NaN warning these trials used to raise is gone (see
         # nlminb_trial_fn()), so the trace is where a user watching the
         # optimizer still sees them
         if (isTRUE(opt$nonfinite_trials > 0L)) {
           paste0(", ", vb_plural(opt$nonfinite_trials,
                                  "non-finite trial"))
         })
}

#' verbose >= 2 turns on the optimizer's own iteration trace, unless the
#' user already asked for one. That trace is printed by nlminb/optim
#' themselves and so goes to stdout, not through message(); a custom
#' optimizer receives optCtrl untouched.
#'
#' @noRd
vb_trace_ctrl <- function(optCtrl, optimizer) {
  optCtrl <- optCtrl %||% list()
  if (!is.null(optCtrl[["trace"]])) return(optCtrl)
  if (identical(optimizer, "nlminb")) {
    optCtrl$trace <- 1L
  } else if (identical(optimizer, "optim")) {
    optCtrl$trace <- 1L
    if (is.null(optCtrl[["REPORT"]])) optCtrl$REPORT <- 1L
  }
  optCtrl
}

#' Fitting core shared by frm() and refit(): objective build through the
#' convergence check. A non-NULL `template` bypasses make_start (warm
#' starts when refitting to a new response).
#'
#' `objective_only = TRUE` stops one step before the optimizer and hands
#' back an UNFITTED object carrying the taped objective, the starting
#' template and the resolved bounds. `frm_sample()` on a formula needs
#' exactly that: NUTS wants the density, not the mode, and every guard
#' above the optimizer (REML, quadrature and the structured families'
#' own refusals) still has to run. The autoscale pre-fit is skipped, because
#' it is an optimization and there is nothing here to warm-start.
#'
#' @noRd
fit_assembled <- function(spec, frame, bform, cl, REML, start, control,
                          se, lower, upper, prior, quadrature,
                          importance = 0L,
                          template = NULL, data2 = list(),
                          objective_only = FALSE,
                          announce_start = FALSE) {
  lower_arg <- lower
  upper_arg <- upper
  vb <- verbose_level(control)
  n_imp <- as.integer(importance %||% 0L)
  if (vb) {
    t_fit <- vb_now()
    vb_say("fit: ", vb_fit_detail(spec, REML, control, quadrature,
                                  prior, n_imp))
  }
  ascale <- if (isTRUE(control$autoscale) && !objective_only) {
    autoscale_plan(frame)
  }
  if (!is.null(ascale) && is.null(template)) {
    # two-stage warm start: fit the standardized frame, back-transform
    # the optimum, and continue below as the ordinary unscaled fit
    # (see R/autoscale.R). Doubles the (cheap) optimization. A caller
    # template (refit and friends) skips the pre-fit but keeps the
    # plan, so the optimizer and sdreport still run in natural units.
    if (vb) t0 <- vb_now()
    template <- autoscale_prefit(spec, frame, bform, cl, REML = REML,
                                 start = start, control = control,
                                 lower = lower, upper = upper,
                                 prior = prior,
                                 quadrature = quadrature, plan = ascale)
    if (vb) vb_stage("autoscale pre-fit", t0)
  }
  if (vb) t0 <- vb_now()
  nll <- build_objective(frame)
  # also read by make_start(): a prior's location places a nonlinear
  # parameter's starting value
  prior_entries <- NULL
  if (!is.null(prior)) {
    # MAP / regularized ML: the optimized objective includes the prior
    # terms, so logLik/AIC are penalized quantities - documented
    ri <- resolve_prior_input(list(frame = frame, spec = spec), prior)
    prior_entries <- ri$entries
    if (length(ri$entries)) {
      nll0 <- nll
      nlp <- neg_log_prior_fn(ri$entries)
      nll <- function(pars) nll0(pars) + nlp(pars)
    }
    # The prior is the ONLY source of bounds now that frm() has no
    # lower/upper of its own. `lower`/`upper` still arrive from the
    # refit paths (allfit, anova, influence, simulate), which forward
    # the box a fit was built with alongside its prior; those two agree
    # by construction, so the merge is idempotent there and exists to
    # keep a caller that forwards only the box from losing it.
    if (length(ri$lower)) {
      lower <- utils::modifyList(as.list(ri$lower),
                                 as.list(lower %||% c()))
      lower <- unlist(lower)
    }
    if (length(ri$upper)) {
      upper <- utils::modifyList(as.list(ri$upper),
                                 as.list(upper %||% c()))
      upper <- unlist(upper)
    }
  }
  if (is.null(template)) {
    template <- make_start(frame, start, prior_entries,
                           announce = announce_start)
  }

  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- c(if (!is.null(template[["b"]])) "b",
              if (!is.null(template[["miss"]])) "miss")
  if (REML) random <- c(random, "beta")
  if (!length(random)) random <- NULL

  # A family whose likelihood does not factorize over rows refuses the
  # fitting options its own structure says it cannot answer, in its own
  # words: each of REML, quadrature and profiling integrates something
  # out with a Laplace approximation about a single inner mode, which a
  # structured likelihood either has no definition for or has several
  # of. The structure also gets to look at where the fit STARTS, which
  # is where a multimodal likelihood's traps are.
  check_structure_fit(spec, frame, template, REML, quadrature, control)
  # ps() blocks refuse the same three options for a reason of their
  # own: each integrates out something a penalized block has already
  # put inside the Laplace approximation.
  check_ps_fit(frame, REML, quadrature, control)

  # The importance correction reweights the SAME Laplace integral, so
  # every guard the Laplace fit passes still applies and only the
  # correction's own restrictions are added here, before any tape.
  if (n_imp > 0L) {
    check_importance_scope(spec, frame, template, REML, quadrature, control)
  }

  integrate <- NULL
  if (isTRUE(quadrature)) {
    if (!is.null(template[["miss"]])) {
      frm_stop("quadrature = TRUE cannot be combined with mi()",
               call. = FALSE)
    }
    # The Gauss-Kronrod rule integrates whatever density the tape
    # produces, so a scalar Student-t latent is in scope. There the
    # result is EXACT rather than approximate, which makes
    # quadrature = TRUE the recommended check on a t-block Laplace fit
    # (see ?frmtmb-student-re and dev/tre-feasibility.md section 4).
    scalar_iid <- vapply(frame[["re_blocks"]], function(bk) {
      bk[["dim"]] == 1L &&
        bk[["covstruct"]] %in% c("us", "diag", "homdiag", "us_t", "diag_t")
    }, TRUE)
    if (!length(scalar_iid) || !all(scalar_iid)) {
      frm_stop("quadrature = TRUE currently supports scalar random ",
               "intercepts only (every block must be a dim-1 us/diag term)",
               call. = FALSE)
    }
    if (REML) {
      frm_stop("quadrature = TRUE cannot be combined with REML = TRUE",
               call. = FALSE)
    }
    if (length(frame[["autocor"]] %||% list())) {
      # the Gauss-Kronrod rule integrates one scalar random effect
      # against a PRODUCT of per-row densities; an R-side residual is a
      # joint density over each group, so no per-row integrand exists
      frm_stop("quadrature = TRUE cannot be combined with the residual ",
               "correlation term ", frame[["autocor"]][[1L]]$label,
               ": the rule integrates a random effect against ",
               "per-observation densities, and this residual is a joint ",
               "density over each group. Use quadrature = FALSE (Laplace) ",
               "or REML = TRUE", call. = FALSE)
    }
    # The truncation normalizer is log(F(ub) - F(lb)) over plain CDFs.
    # The Gauss-Kronrod nodes reach random-effect values where that
    # difference underflows to exactly zero while the density itself is
    # still representable, so the integrand is +Inf there and the
    # marginalized objective comes back -Inf - at the Laplace optimum
    # as well as at the starting values. Laplace never leaves the
    # neighborhood of the mode and is unaffected. Refusing beats
    # reporting logLik = +Inf as a converged fit.
    trunc_resp <- names(which(vapply(
      frame[["aterm_values"]],
      function(a) !is.null(a$trunc_lb) || !is.null(a$trunc_ub), TRUE)))
    if (length(trunc_resp)) {
      frm_stop("quadrature = TRUE cannot be combined with trunc() (",
               paste(trunc_resp, collapse = ", "), "): the truncation ",
               "normalizer underflows at the Gauss-Kronrod nodes and the ",
               "marginalized objective is unbounded. Use quadrature = ",
               "FALSE (Laplace), REML = TRUE, or ",
               "frmtmb_control(profile = TRUE)", call. = FALSE)
    }
    # adaptive Gauss-Kronrod marginalization per scalar random effect
    # (TMB's experimental `integrate`; the nAGQ analogue - matches
    # glmer(nAGQ = 25) in tests). Spec replicated from TMB's GK().
    integrate <- list(b = structure(
      list(dim = 1, adaptive = FALSE, debug = FALSE,
           method = "marginal_gk"),
      class = "GK"
    ))
  }
  profile_arg <- NULL
  if (isTRUE(control$profile)) {
    if (REML) {
      frm_stop("frmtmb_control(profile = TRUE) cannot be combined with ",
               "REML = TRUE (beta is already integrated)", call. = FALSE)
    }
    if (isTRUE(quadrature)) {
      frm_stop("frmtmb_control(profile = TRUE) cannot be combined with ",
               "quadrature = TRUE", call. = FALSE)
    }
    profile_arg <- "beta"
  }
  # Under integrate= the objective is built in two passes: the plain
  # Laplace tape first, then the marginalized one calibrated at its
  # optimum. See quad_fit() for why. The Laplace objective is kept
  # afterwards because it is the only source of the conditional modes.
  # The importance correction needs the same two-pass shape for the
  # same reason: its proposal IS the Laplace Gaussian, so the Laplace
  # tape has to exist to supply the conditional modes and their
  # Hessian, both while fitting and afterwards for ranef().
  lap_obj <- if (is.null(integrate) && n_imp == 0L) NULL else {
    RTMB::MakeADFun(nll, template, random = random, map = frame[["map"]],
                    silent = TRUE)
  }
  obj <- if (is.null(lap_obj)) {
    RTMB::MakeADFun(nll, template, random = random, map = frame[["map"]],
                    profile = profile_arg, silent = TRUE)
  } else lap_obj
  if (vb) {
    vb_stage("tape", t0,
             paste0(length(obj$par), " outer, ",
                    length(obj$env$random), " inner parameters"))
  }
  if (objective_only) {
    if (isTRUE(quadrature)) {
      frm_stop("quadrature = TRUE has no unfitted form: the Gauss-Kronrod ",
               "tape is calibrated at a Laplace OPTIMUM, and stopping ",
               "before the optimizer leaves nothing to calibrate it at. ",
               "Sample the Laplace objective instead (quadrature = FALSE)",
               call. = FALSE)
    }
    if (n_imp > 0L) {
      frm_stop("`importance` has no unfitted form: the proposal is the ",
               "Laplace Gaussian at a conditional mode, and stopping ",
               "before the optimizer leaves no mode to centre it on. ",
               "Sample the Laplace objective instead (importance = 0)",
               call. = FALSE)
    }
    return(unfitted_object(spec, frame, obj, template, bform, cl, REML,
                           prior, data2, control, lower_arg, upper_arg))
  }
  # control must ride along: outer_par_names drops beta under
  # profile = TRUE, and a shim without it misaligns every bound
  bounds <- resolve_bounds(list(frame = frame, REML = REML,
                                control = control), lower, upper)
  # a badly scaled coefficient has a badly scaled gradient too: judge
  # (and steer) the optimizer in per-parameter natural units
  par_units <- if (!is.null(ascale)) {
    autoscale_units(frame, ascale, names(obj$par))
  }
  # the optimizer trace rides on the optimizer's own control list, so
  # keep it out of the control stored on the fit (refit and friends
  # reuse that list and must not inherit a trace)
  ctl_opt <- control
  if (vb >= 2L) {
    ctl_opt$optCtrl <- vb_trace_ctrl(control$optCtrl, control$optimizer)
  }
  imp <- NULL
  # shared by the first attempt, every recovery start and every
  # stationary-point escape, so the trials a failed attempt mapped are
  # still in the count the fit reports
  tally <- new.env(parent = emptyenv())
  tally$n <- 0L
  if (n_imp > 0L) {
    imp <- importance_fit(nll, template, random, frame[["map"]], lap_obj,
                          ctl_opt, bounds, par_units, frame,
                          imp_n_draw(n_imp), vb)
    obj <- imp$obj
    opt <- imp$opt
  } else if (is.null(integrate)) {
    opt <- fit_error_context(
      spec, start, REML, control, quadrature, prior,
      tryCatch(optimize_obj(obj, ctl_opt, bounds, par_units, verbose = vb,
                            tally = tally),
               error = function(e) {
                 rs <- fit_recovery_starts(obj, nll, template, random,
                                           frame[["map"]], frame, start,
                                           ctl_opt, bounds, par_units,
                                           prior_entries)
                 for (lbl in names(rs)) {
                   if (vb) {
                     vb_say("optimizer failed (", conditionMessage(e),
                            "); restarting from ", lbl)
                   }
                   op <- tryCatch(optimize_obj(obj, ctl_opt, bounds,
                                               par_units, verbose = vb,
                                               start_par = rs[[lbl]],
                                               tally = tally),
                                  error = function(e2) NULL)
                   if (!is.null(op)) return(op)
                 }
                 stop(e)
               }))
  } else {
    qf <- quad_fit(nll, template, random, frame[["map"]], integrate, lap_obj,
                   ctl_opt, bounds, par_units, frame, vb)
    obj <- qf$obj
    opt <- qf$opt
  }
  # A family that declares a stationary point gets looked at from both
  # sides of it. Not under importance=: its diagnostics (the effective
  # sample size, the proposal) belong to the optimum importance_fit()
  # returned, and replacing that optimum would leave them describing a
  # point the fit no longer reports.
  if (is.null(imp)) {
    opt <- escape_stationary(obj, opt, frame, ctl_opt, bounds, par_units,
                             verbose = vb, tally = tally)
  }

  # Estimates come cheaply from the parameter list at the optimum;
  # sdreport (a quarter of typical fit time) is computed on demand
  # through sdr_of() unless se = TRUE asked for it now.
  # The corrected tape carries no random effects at all, so its
  # parList() has no slot for them; the Laplace inner solve at the
  # corrected optimum is both the source of the conditional modes and
  # a complete parameter list to start from.
  est <- if (is.null(imp)) obj$env$parList(opt$par) else {
    solved_par_list(lap_obj, opt$par)
  }
  for (nm in names(frame[["par_template"]])) {
    names(est[[nm]]) <- names(frame[["par_template"]][[nm]])
  }
  # under integrate= (quadrature), parList leaves outer components NA;
  # the optimizer vector is authoritative for them either way
  pn <- names(opt$par)
  for (cp in setdiff(unique(pn), random %||% character(0))) {
    pos <- seq_along(frame[["par_template"]][[cp]])
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      pos <- setdiff(pos, frame[["betad_fixed_idx"]])
    }
    est[[cp]][pos] <- unname(opt$par[pn == cp])
  }
  if (!is.null(integrate)) {
    # integrate= removes the random effects from the tape, so this
    # objective has no conditional modes to report: parList() leaves
    # them NA and slides the outer values into their slots. ranef(),
    # fitted() and predict(newdata =) all read them, so recover them
    # from the inner Newton solve of the Laplace objective at the
    # quadrature optimum.
    inner <- solved_par_list(lap_obj, opt$par)
    for (cp in random) est[[cp]][] <- inner[[cp]]
  }
  fit <- structure(
    list(spec = spec, frame = frame, obj = obj, opt = opt, sdr = NULL,
         REML = REML, estimates = est, prior = prior,
         bform = bform, call = cl, data = frame[["data_frame"]],
         data2 = data2,
         control = control, quadrature = isTRUE(quadrature),
         importance = imp_record(imp),
         lower = lower_arg, upper = upper_arg, par_units = par_units,
         cache = new.env(parent = emptyenv())),
    class = "frmtmb_fit"
  )
  if (!is.null(imp)) {
    imp_ess_warning(imp$ess$ess, imp$lay, imp$plan[["n_draw"]],
                    control$importance_ess %||% imp_ess_floor)
  }
  if (se) {
    if (vb) t0 <- vb_now()
    fit$cache$sdr <- autoscale_sdreport(fit)
    if (vb) vb_stage("sdreport", t0)
  }
  chk <- check_convergence(fit, control)
  # Fit-end checks: a family that wants to inspect its own fit, and the
  # knot-span coverage of any ps() block. Both need the finished object,
  # which is why neither can live in family_finalize() or in the frame.
  fit_end_checks(fit)
  if (vb) {
    vb_stage("done", t_fit,
             paste0("objective ", format(fit$opt$objective, digits = 8),
                    ", max|grad| ", format(chk$grad, digits = 3), ", ",
                    vb_plural(length(chk$warnings), "warning"),
                    # the per-stage lines give each run's own count;
                    # this is the total the fit carries
                    if (isTRUE(fit$opt$nonfinite_trials > 0L)) {
                      paste0(", ", vb_plural(fit$opt$nonfinite_trials,
                                             "non-finite trial"))
                    }))
  }
  fit
}

#' The unfitted counterpart of a `frmtmb_fit`: the same spec, frame and
#' taped objective, with the STARTING values in the `estimates` slot and
#' no `opt`, no sdreport and no convergence check.
#'
#' It carries `frmtmb_fit` as a second class on purpose. Everything that
#' reads only `spec`, `frame` and `estimates` - `eval_dpars()`,
#' `predict()`, the whole draws surface through `draws_fit_at()` - is
#' correct on it once a real parameter vector has been written in, which
#' is exactly what a posterior draw does. What is NOT correct is reading
#' the slot as an estimate, so `require_fitted()` guards the accessors
#' that would.
#'
#' @noRd
unfitted_object <- function(spec, frame, obj, template, bform, cl, REML,
                            prior, data2, control, lower, upper) {
  est <- template
  for (nm in names(frame[["par_template"]])) {
    names(est[[nm]]) <- names(frame[["par_template"]][[nm]])
  }
  structure(
    list(spec = spec, frame = frame, obj = obj, opt = NULL, sdr = NULL,
         REML = REML, estimates = est, prior = prior,
         bform = bform, call = cl, data = frame[["data_frame"]],
         data2 = data2,
         control = control, quadrature = FALSE, importance = NULL,
         lower = lower, upper = upper, par_units = NULL,
         cache = new.env(parent = emptyenv())),
    class = c("frmtmb_unfitted", "frmtmb_fit")
  )
}

#' The object really is a fit, before anything reaches into its slots.
#'
#' Without this the first reach lands inside `sdreport()` and the reader
#' gets a message from three packages down that names neither the
#' argument nor the fault.
#'
#' @noRd
require_frmtmb_fit <- function(fit, what) {
  if (inherits(fit, "frmtmb_fit")) return(invisible(NULL))
  frm_stop(what, " needs a model fitted by frm(), not ", arg_desc(fit),
           call. = FALSE)
}

#' Refuse a method that reports a maximum-likelihood quantity on an
#' object that has none. `what` names the method, so the one message
#' still tells the user which call to change.
#'
#' @noRd
require_fitted <- function(fit, what) {
  if (!inherits(fit, "frmtmb_unfitted")) return(invisible(NULL))
  frm_stop(what, " needs a fitted model. This object was assembled by ",
           "frm_sample() from a formula, for sampling only: it holds the ",
           "objective and the starting values, but no optimizer result, no ",
           "mode and no sdreport, so there is no estimate to report. Use ",
           "the draws (summary(), fixef(), VarCorr(), hypothesis(), ",
           "posterior_*()), or fit the model with frm() first",
           call. = FALSE)
}

#' parList() at an outer parameter vector with the inner problem solved
#' there first. fn() runs the inner Newton solve and leaves the
#' conditional modes in last.par, which is parList's default `par`.
#'
#' @noRd
solved_par_list <- function(obj, par) {
  obj$fn(par)
  obj$env$parList(par)
}

#' Which of two quadrature candidates quad_fit() keeps. A stationary one
#' wins outright (the loop stops there); among non-stationary ones the
#' lowest objective is the best point reached, not the first that
#' happened to tape and optimize without breaking.
#'
#' @noRd
quad_keep_best <- function(best, a) {
  if (is.null(a) || is.null(a$obj)) return(best)
  if (isTRUE(a$stationary)) return(a)
  if (is.null(best) || a$opt$objective < best$opt$objective) a else best
}

#' Build and optimize the Gauss-Kronrod (integrate=) objective.
#'
#' TMBad's marginal_gk transform rescales each integrand ONCE: it finds
#' the mode and curvature of the log-integrand by finite differences and
#' bakes that (mu, sigma) pair into the tape as constants, at whichever
#' parameter values `template` happens to hold. Everything downstream
#' depends on that one calibration, so this function has to do two
#' things the transform does not do for itself.
#'
#' 1. Tape at a sensible point. From the cold start the frozen rescaling
#'    sits far from the real conditional mode, and for every family
#'    whose inverse link exponentiates the linear predictor the rescaled
#'    integrand then overflows: obj$fn() is NaN before the optimizer
#'    takes a step (poisson, Gamma and Beta over nested scalar blocks,
#'    Beta over a single one). Gaussian responses survive it only
#'    because their integrand is quadratic wherever it is sampled. So
#'    fit the plain Laplace objective first and tape the marginalized
#'    one at that optimum: the two optima maximize the same marginal
#'    likelihood, one exactly and one to O(n^-1).
#'
#' 2. Recalibrate when the tape expires. A frozen rescaling is only
#'    trustworthy near the point it was made at, so it can run out in
#'    two ways. The optimizer can walk far enough that the rescaled
#'    integrand breaks (RTMB then raises "NA/NaN gradient evaluation"
#'    from inside nlminb), or it can stop somewhere the tape's own
#'    gradient does not vanish - a mixture whose Laplace fit collapses a
#'    mixing weight to exp(-35) does that, and the reported objective is
#'    then a value no neighborhood shares. Either way the answer is to
#'    tape again at the best point reached and carry on, and to keep the
#'    cold template as a last anchor when the Laplace optimum is the bad
#'    one. Each candidate costs a tape, so they are tried in order and
#'    the first stationary result wins; if none is stationary the
#'    candidate with the lowest objective does.
#'
#' 3. Displace the widths when neither anchor holds. What the transform
#'    freezes is a width per random effect, so a calibration that expires
#'    within a step or two is one whose widths are wrong, and the
#'    parameter that sets them is theta. Half a log-SD either way, then a
#'    whole one, is enough where displacement helps at all - it recovers
#'    single-block fits whose every un-displaced calibration dies inside
#'    the optimizer. It cannot recover a nested block, where the outer
#'    integrand is the inner rescaling's output and the objective is NaN
#'    before the optimizer takes a step.
#'
#' @noRd
quad_fit <- function(nll, template, random, map, integrate, lap_obj,
                     control, bounds, par_units, frame = NULL, vb = 0L,
                     rounds = 3L) {
  if (vb) t0 <- vb_now()
  lap_opt <- optimize_obj(lap_obj, control, bounds, par_units)
  if (vb) vb_stage("quadrature warm start", t0, vb_opt_detail(lap_opt))

  # A template holding the outer values `par` plus the conditional
  # modes there, which is what MakeADFun needs to calibrate the tape.
  anchor <- function(par) solved_par_list(lap_obj, par)[names(template)]

  # One build-and-optimize pass. Returns the fit, or - when the tape
  # broke mid-optimization - the best point it reached, so the caller
  # can recalibrate there.
  attempt <- function(tpl) {
    if (vb) t0 <- vb_now()
    o <- RTMB::MakeADFun(nll, tpl, random = random, map = map,
                         integrate = integrate, silent = TRUE)
    f0 <- try(o$fn(o$par), silent = TRUE)
    if (inherits(f0, "try-error") || !is.finite(f0)) return(NULL)
    if (vb) vb_stage("quadrature tape", t0)
    op <- tryCatch(optimize_obj(o, control, bounds, par_units,
                                verbose = vb),
                   error = function(e) NULL)
    if (is.null(op) || !is.finite(op$objective)) {
      pb <- o$env$last.par.best
      if (is.null(pb) || length(pb) != length(o$par) || anyNA(pb)) {
        return(NULL)
      }
      return(list(retry = stats::setNames(as.numeric(pb),
                                          names(o$par))))
    }
    g <- try(max(abs(o$gr(op$par) * (par_units %||% 1))), silent = TRUE)
    if (inherits(g, "try-error")) g <- NA_real_
    list(obj = o, opt = op, stationary = isTRUE(g < control$grad_tol))
  }

  # Calibration points, in the order they are tried. The Laplace optimum
  # first; then the untouched template, because the anchor itself can be
  # the problem (the Laplace optimum can sit on a singular variance
  # component); then the anchor with the integrand widths displaced.
  lap_tpl <- anchor(lap_opt$par)
  seeds <- list(lap_tpl, template)
  if (length(lap_tpl$theta)) {
    seeds <- c(seeds, lapply(c(-0.5, 0.5, -1), function(s) {
      tpl <- lap_tpl
      tpl[["theta"]] <- tpl[["theta"]] + s
      tpl
    }))
  }

  best <- NULL
  for (seed in seeds) {
    tpl <- seed
    for (i in seq_len(rounds)) {
      a <- attempt(tpl)
      best <- quad_keep_best(best, a)
      if (is.null(a) || isTRUE(a$stationary) || is.null(a$retry)) break
      # the frozen rescaling expired where the optimizer walked to:
      # tape again at the best point it managed
      tpl <- anchor(a$retry)
    }
    if (!is.null(best) && isTRUE(best$stationary)) break
  }
  if (is.null(best)) {
    frm_stop(quad_breakdown_message(frame, lap_tpl$theta), call. = FALSE)
  }
  best[c("obj", "opt")]
}

#' The refusal message when every calibration point breaks down. A bare
#' "the objective was non-finite" tells a user nothing they can act on,
#' and the two shapes this failure takes want different answers, so name
#' whichever applies: an iterated integral over nested blocks (the outer
#' integrand is the inner rescaling's output, and no calibration of the
#' outer one can repair that), and a variance component so narrow that
#' the transform's unit-step finite differences cannot measure it.
#'
#' @noRd
quad_breakdown_message <- function(frame, theta) {
  blocks <- frame[["re_blocks"]] %||% list()
  desc <- function(b) paste0("'", b$term_label, "' (", b$n_levels,
                             " levels)")
  sds <- vapply(blocks, function(b) {
    v <- tryCatch(covstruct_registry[[b$covstruct]]$vcov(
      theta[b$theta_idx], b), error = function(e) NULL)
    if (is.null(v)) NA_real_ else sqrt(v[1L, 1L])
  }, 1)
  why <- character(0)
  if (length(blocks) > 1L) {
    why <- c(why, paste0("the model asks for an iterated integral over ",
                         length(blocks), " nested blocks (",
                         paste(vapply(blocks, desc, ""), collapse = ", "),
                         ") on ", frame[["n_obs"]], " observations, and the ",
                         "outer integrand is itself the frozen ",
                         "rescaling of the inner one"))
  }
  sing <- which(is.finite(sds) & sds < 1e-3)
  if (length(sing)) {
    why <- c(why, paste0("the Laplace optimum leaves ",
                         paste(vapply(blocks[sing], desc, ""),
                               collapse = ", "),
                         " on a singular variance component (sd ",
                         paste(format(sds[sing], digits = 2),
                               collapse = ", "),
                         "), which the transform measures by finite ",
                         "differences with a step of 1"))
  }
  if (!length(why)) {
    why <- paste0("the calibration expired on ",
                  paste(vapply(blocks, desc, ""), collapse = ", "),
                  " at every point tried")
  }
  paste0("quadrature = TRUE could not marginalize this model: the ",
         "Gauss-Kronrod objective broke down (a non-finite objective ",
         "or gradient) at every calibration point - the Laplace ",
         "optimum, the best point the optimizer reached, the cold ",
         "start, and the optimum with the integrand widths displaced. ",
         "Here ", paste(why, collapse = "; "),
         ". Refit with quadrature = FALSE: the Laplace approximation ",
         "fits this model (it is what the warm start above already did)")
}

#' The joint precision is the covariance source for parameters outside
#' cov.fixed: REML (beta random) and control profile = TRUE (beta inner).
#'
#' @noRd
needs_jp <- function(fit) {
  fit$REML || isTRUE(fit$control$profile)
}

#' Memoized sdreport: the standard-error machinery (summary, vcov,
#' confint, predict se.fit, diagnose) triggers it on first use.
#'
#' @noRd
sdr_of <- function(fit) {
  require_fitted(fit, paste("The standard-error machinery (summary(),",
                            "vcov(), confint(), frm_linpred(se.fit =))"))
  cache <- fit$cache
  if (is.null(cache$sdr)) {
    cache$sdr <- autoscale_sdreport(fit)
  }
  cache$sdr
}

#' Control parameters for frmtmb fits
#'
#' @param optimizer `"nlminb"` (default), `"optim"` (L-BFGS-B), or a
#'   function with signature `(par, fn, gr, lower, upper, control)`
#'   returning a list with elements `par`, `objective`, `convergence`
#'   (0 = success), and optionally `message` - the hook for optimx,
#'   nloptr, and friends without frmtmb depending on them. Example:
#'   ```
#'   nlopt <- function(par, fn, gr, lower, upper, control) {
#'     r <- nloptr::nloptr(par, fn, gr, lb = lower, ub = upper,
#'                         opts = list(algorithm = "NLOPT_LD_LBFGS",
#'                                     xtol_rel = 1e-10, maxeval = 2000))
#'     list(par = r$solution, objective = r$objective,
#'          convergence = as.integer(r$status < 0), message = r$message)
#'   }
#'   frm(..., control = frmtmb_control(optimizer = nlopt))
#'   ```
#' @param optCtrl Control list passed to the built-in optimizers
#'   ([stats::nlminb()] / [stats::optim()]).
#' @param restarts Number of times to restart the optimizer from the
#'   current optimum while the gradient remains above `grad_tol`.
#' @param grad_tol Warn (and restart) if the maximum absolute gradient at
#'   the optimum exceeds this value.
#' @param profile Experimental: move the primary (`beta`) coefficients
#'   into the inner (Laplace) problem, TMB's `profile` argument - the
#'   analog of `glmmTMBControl(profile = TRUE)` and `glmer(nAGQ = 0)`.
#'   Speeds up models with many fixed effects, and like those it is an
#'   approximation: estimates differ slightly from the exact fit.
#'   Coefficient covariance comes from the joint precision. Not
#'   compatible with `REML = TRUE` or `quadrature = TRUE`;
#'   profile/uniroot `confint()` and `hypothesis(method = "profile")`
#'   need a non-profiled fit.
#' @param sparse_x Build the parametric fixed-effect design matrices as
#'   sparse [Matrix::sparse.model.matrix()] objects, the analog of
#'   `glmmTMB(sparseX =)`. Worth it when a many-level fixed factor makes
#'   the dense design dominate memory; estimates are identical either
#'   way. `model.matrix()` on the fit then returns a sparse matrix.
#' @param autoscale Standardize badly scaled continuous predictors
#'   internally (the lme4 >= 1.1.37 feature): fit a copy of the model
#'   with each qualifying fixed-effect column centered and scaled
#'   (scaled only, in a linear predictor without an intercept), map
#'   that optimum back to the original parameterization exactly, and
#'   warm-start the ordinary fit there. Reported results are always on
#'   the original scale, so every downstream method works unchanged;
#'   the cost is a second (cheap) optimization. Columns qualify when
#'   they are parametric and take more than two distinct values;
#'   intercepts, factor contrasts, smooth bases, and mo()/mi() columns
#'   are never touched, and the whole step is a silent no-op when
#'   nothing qualifies. Compatible with `profile = TRUE`. Under
#'   `prior` or bounds, the first stage applies them to the scaled
#'   coefficients; the second stage is the fit that is reported.
#' @param check_nlev_1 What to do about a scalar random-effect term
#'   whose grouping factor has a single level: `"warning"` (default),
#'   `"ignore"`, or `"stop"`, following lme4's `lmerControl()` check
#'   vocabulary. Such a term has no variance to estimate - the single
#'   level is absorbed by the intercept - and its standard deviation
#'   collapses to zero. Structured blocks over several terms per level
#'   (`ar1()`, `us()`, the spatial covariance structures) are never
#'   flagged: one grouping level there is a single realization of a
#'   field, which is the normal way to write them.
#' @param check_olre What to do about an observation-level random
#'   effect - one level per row - on a gaussian, student or lognormal
#'   response: `"warning"` (default), `"ignore"`, or `"stop"`. Its
#'   variance is confounded with the residual standard deviation, so
#'   only their sum is identified and the split between them is
#'   arbitrary. The check is skipped when `sigma` is not free to absorb
#'   it - a `se()` response or a constant `sigma` - which is the
#'   random-effects meta-analysis, and for discrete families, where an
#'   observation-level term is the usual overdispersion model. It is
#'   also skipped for the known-structure blocks `gr(cov = )`,
#'   `gr(prec = )` and `equalto()`: their level covariance is not
#'   proportional to the identity, so the block and the residual are
#'   separately identified. That is the animal model with one
#'   measurement per individual.
#' @param importance_seed Seed for the standard normal draws of
#'   `frm(importance =)`. The draws are taken from a private random
#'   stream, so the fit neither reads nor disturbs the session's random
#'   state: the same seed gives the same answer to the last bit,
#'   whatever else has drawn random numbers.
#' @param importance_rounds Cap on the number of times
#'   `frm(importance =)` refreezes its proposal. Each round rebuilds
#'   the proposal at the current estimate and optimizes again, and the
#'   loop stops as soon as either the estimates move less than `1e-3`
#'   or the freshly anchored objective is stationary at them to
#'   `grad_tol`. Measured on this package's own designs that takes two
#'   to four rounds, with the move falling about tenfold each time, so
#'   the default leaves one in hand; an unused round costs nothing, and
#'   each round actually taken costs one tape. A fit that uses every
#'   round warns, and the warning distinguishes two cases: moves that
#'   are still shrinking want a larger cap, while moves that are all
#'   the SAME size are a stalled iteration, whose reported shift is
#'   that step times the round count and grows with the cap instead of
#'   settling.
#' @param importance_ess Effective sample size, as a fraction of the
#'   draw count, below which `frm(importance =)` warns and names the
#'   groups. The default `0.25` separates two measured regimes rather
#'   than clearing every design. A proposal that is working sits near
#'   the top of the range: the probe design behind this correction (60
#'   groups of 8 Bernoulli rows, correlated slope) holds a worst group
#'   of `0.95` at 1000 draws, and a gaussian response, where the
#'   correction has nothing to correct, holds exactly `1` at any draw
#'   count. Displacing that probe's two log standard deviations by half
#'   a unit takes its worst group to `0.03`. The designs the warning is
#'   FOR sit between the two, and draws are what move them: a correlated
#'   slope in `mu` plus a `sigma` block over 12 groups of 10 rows holds
#'   a worst group of `0.17` at 500 draws and `0.40` at 2000.
#' @param verbose Report fit progress through [message()], one terse
#'   line per stage with its elapsed seconds, so a slow fit shows where
#'   the time went. `FALSE` (default) is silent and costs nothing.
#'   `TRUE` (or `1`) reports parsing, frame assembly, the autoscale
#'   pre-fit, tape construction, each optimizer run and restart,
#'   `sdreport()` when `se = TRUE`, and a closing line with the
#'   objective, the maximum absolute gradient, and the number of
#'   convergence warnings. The fit itself opens with a line naming the
#'   family, the mode (ML or REML, plus profile, quadrature, autoscale,
#'   prior), and the optimizer. `2` adds the optimizer's own trace, by
#'   setting `optCtrl$trace` unless you set it yourself. That trace is
#'   printed by [stats::nlminb()] / [stats::optim()] to standard
#'   output, not through `message()`, and a custom optimizer function
#'   receives `optCtrl` unchanged, so `verbose` does not reach it.
#'
#'   Use `suppressMessages()` to silence a verbose fit. The refit loops
#'   in [frm_bootstrap()], [influence()], and [frm_allfit()] force
#'   `verbose` off, so a verbose fit does not make them print hundreds
#'   of lines; [refit()] and [frm_multiple()] report normally.
#' @return A list of control settings.
#'
#' @srrstats {RE2.0} Data transformations are documented and can be turned
#'   off. `autoscale` is the only transformation of predictor values, it
#'   is `FALSE` by default, and the documentation states exactly which
#'   columns qualify, which are never touched (intercepts, factor
#'   contrasts, smooth bases, `mo()`/`mi()` columns), and that results are
#'   always mapped back and reported on the original scale. `sparse_x`
#'   changes only the storage of the design and is documented as leaving
#'   estimates identical. Factor handling is delegated to
#'   `stats::model.matrix()` with the contrasts frozen at fit time and
#'   reapplied to `newdata`.
#' @srrstats {RE2.3} `autoscale = TRUE` centers and scales qualifying
#'   continuous fixed-effect columns, that is, converts them to z-scores.
#'   Centering is applied only where an intercept exists to absorb the
#'   shift, and the effect is exercised in `tests/testthat/test-autoscale.R`,
#'   which checks that the scaled and unscaled fits agree.
#' @srrstats {RE3.2} Convergence thresholds have documented defaults:
#'   `grad_tol = 1e-3` on the maximum absolute gradient at the optimum,
#'   and `optCtrl = list(iter.max = 1000, eval.max = 1000)` for the
#'   built-in optimizers. Both appear in the usage section of the manual
#'   page with an `@param` describing them.
#' @srrstats {RE3.3} Those thresholds can be set explicitly:
#'   `grad_tol` for the gradient criterion, `optCtrl` for the optimizer's
#'   own tolerances and iteration caps, `restarts` for how many times to
#'   restart from the current optimum while the gradient is still above
#'   `grad_tol`, and `optimizer` to substitute another optimizer
#'   altogether.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 20)))
#' dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # another optimizer, with its own control list
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'            control = frmtmb_control(optimizer = "optim",
#'                                     optCtrl = list(maxit = 500)))
#' fit$opt$convergence
#'
#' # a tighter gradient criterion, with restarts from the current
#' # optimum until it is met
#' frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'     control = frmtmb_control(grad_tol = 1e-4, restarts = 3))
#'
#' # badly scaled predictors: fit an internally standardized copy first,
#' # then warm-start the reported fit from it
#' dd$xbig <- dd$x * 1e5
#' frm(bf(y ~ xbig + (1 | g)) + gaussian(), data = dd,
#'     control = frmtmb_control(autoscale = TRUE))
#'
#' # the object is a plain list, so it can be built once and reused
#' ctrl <- frmtmb_control(check_nlev_1 = "ignore")
#' ctrl$optimizer
#' @export
frmtmb_control <- function(optimizer = "nlminb",
                           optCtrl = list(iter.max = 1000, eval.max = 1000),
                           restarts = 1, grad_tol = 1e-3,
                           profile = FALSE, sparse_x = FALSE,
                           autoscale = FALSE,
                           check_nlev_1 = c("warning", "ignore", "stop"),
                           check_olre = c("warning", "ignore", "stop"),
                           importance_seed = 1L, importance_rounds = 5L,
                           importance_ess = 0.25,
                           verbose = NULL) {
  # The three flags below reach isTRUE() one line down, which reads a
  # string or a length-2 vector as FALSE; checking here refuses the
  # mistake instead of quietly turning the option off.
  check_flag(profile, "profile")
  check_flag(sparse_x, "sparse_x")
  check_flag(autoscale, "autoscale")
  check_count(restarts, "restarts", min = 0L)
  check_positive(grad_tol, "grad_tol")
  check_count(importance_seed, "importance_seed", min = 0L)
  check_count(importance_rounds, "importance_rounds", min = 1L)
  check_positive(importance_ess, "importance_ess")
  if (importance_ess > 1) {
    frm_stop("`importance_ess` is an effective sample size as a FRACTION ",
             "of the draw count, so it lies in (0, 1]; ",
             format(importance_ess), " asks for more effective draws than ",
             "there are draws", call. = FALSE)
  }
  # Only the SHAPE of `optimizer` is checked here. Which names are
  # known is settled at fit time, on purpose: that refusal carries the
  # family and mode of the fit it was raised from ("raised while
  # fitting: gaussian, ML, nope"), which a check in this constructor
  # could not know. A length-2 value has no such excuse - it used to
  # reach switch() and report "EXPR must be a length 1 vector" once per
  # element.
  if (!is.function(optimizer) &&
        !(is.character(optimizer) && length(optimizer) == 1L &&
            !is.na(optimizer))) {
    frm_stop("`optimizer` must be a single optimizer name or a function, ",
             "not ", arg_desc(optimizer), call. = FALSE)
  }
  if (!is.list(optCtrl)) {
    frm_stop("`optCtrl` must be a list of options for the optimizer, e.g. ",
             "optCtrl = list(iter.max = 1000), not ", arg_desc(optCtrl),
             call. = FALSE)
  }
  # verbose stays NULL when unset, which is how frm(verbose =) knows an
  # explicit control value must win over its own shortcut
  list(optimizer = optimizer, optCtrl = optCtrl, restarts = restarts,
       grad_tol = grad_tol, profile = isTRUE(profile),
       sparse_x = isTRUE(sparse_x), autoscale = isTRUE(autoscale),
       check_nlev_1 = frm_match_arg(check_nlev_1),
       check_olre = frm_match_arg(check_olre),
       importance_seed = as.integer(importance_seed),
       importance_rounds = as.integer(importance_rounds),
       importance_ess = importance_ess,
       verbose = verbose)
}

#' lme4's lmerControl runs a battery of structural checks before the fit
#' and gives each one an ignore/warning/stop setting; these are the two
#' that change what a frmtmb fit MEANS rather than how fast it runs.
#' Both currently fit silently to an answer the user did not ask for.
#' `[lme4 lmerControl checks]`
#'
#' @noRd
re_check_act <- function(what, msg) {
  switch(what %||% "warning",
         ignore = invisible(NULL),
         stop = frm_stop(msg, call. = FALSE),
         frm_warning(msg, call. = FALSE))
}

#' Runs those two checks over the assembled random-effect blocks, before
#' any fitting starts: a scalar term whose grouping factor has one level,
#' and an observation-level term on a response whose residual sd already
#' holds that variance. Each check reports through `re_check_act()`, so
#' the control setting decides between silence, a warning, and an error.
#'
#' @noRd
check_re_structure <- function(spec, frame, control) {
  gaussian_like <- c("gaussian", "student", "lognormal")
  for (bk in frame[["re_blocks"]]) {
    # smooth / gp / hsgp blocks carry a synthetic n_levels of 1 and no
    # grouping levels at all; only real grouping factors are checked
    if (is.null(bk[["levels"]])) next
    # A structured block over several terms per level (ar1, us, the
    # spatial covstructs) is a single realization of a field, and one
    # group level is the normal way to write it; only a SCALAR term
    # loses its variance to a single level.
    if (bk[["n_levels"]] == 1L && bk[["dim"]] == 1L) {
      re_check_act(
        control$check_nlev_1,
        paste0("Grouping factor '", bk[["group_name"]], "' in `",
               bk[["term_label"]], "` has a single level, so its variance is ",
               "not identified and collapses to zero. Drop the term (it ",
               "is absorbed by the intercept), or set ",
               "frmtmb_control(check_nlev_1 = \"ignore\")"))
      next
    }
    lp <- frame[["linpreds"]][[bk[["components"]][[1L]]$lp_key]]
    resp <- spec$responses[[lp[["resp"]]]]
    fam <- resp$family
    if (is.null(fam)) next
    # se() supplies the residual sd row by row and a constant dpar pins
    # it outright; either way sigma is no longer free to absorb the
    # observation-level variance, so the two are identified. That is
    # exactly the random-effects meta-analysis, where the
    # observation-level term IS the between-study variance.
    sigma_free <- is.null(resp$aterms[["se"]]) &&
      is.null(frame[["linpreds"]][[linpred_key(lp[["resp"]],
                                  "sigma")]]$constant)
    # A known-structure block is not an OLRE even with one row per
    # level. Its levels are correlated through the fixed relationship
    # matrix (gr(cov = A), gr(prec = Q)) or its covariance is fixed
    # outright (equalto), so the block covariance is no longer
    # proportional to the identity and the residual sd is not a
    # reparameterization of it. That is the animal model: the additive
    # genetic variance and the residual variance are separately
    # identified precisely BECAUSE A is not the identity.
    structured <- bk[["covstruct"]] %in% c("gr_cov", "gr_prec", "equalto")
    if (sigma_free && !structured &&
        bk[["n_levels"]] == frame[["n_obs"]] && bk[["dim"]] == 1L &&
        fam[["family"]] %in% gaussian_like &&
        bk[["dpar"]] %in% (fam[["primary_dpars"]] %||% "mu")) {
      re_check_act(
        control$check_olre,
        paste0("`", bk[["term_label"]], "` gives every observation its own ",
               "random effect, and for a ", fam[["family"]], " response that ",
               "variance is confounded with the residual sd: only their ",
               "sum is identified, so the split between them is ",
               "arbitrary. Observation-level random effects are ",
               "meaningful for discrete families (overdispersion), not ",
               "here. Set frmtmb_control(check_olre = \"ignore\") to ",
               "keep it"))
    }
  }
  invisible(NULL)
}

#' The objective as nlminb sees it: a NaN at a TRIAL point becomes +Inf.
#'
#' A link whose inverse is undefined on part of the line (`1/mu^2` below
#' zero, `inverse` for a positive mean, `identity` or `log` for a
#' probability) makes the objective NaN wherever a line search overshoots
#' into that part. PORT treats NaN and +Inf alike, backtracking from
#' either, but R's nlminb warns "NA/NaN function evaluation" for each
#' NaN, so a correct `inverse.gaussian()` fit on its default link warned
#' every time. Measured on the same tape, the two values give the same
#' iterates to the last bit and differ only in the warning
#' (dev/famlink-1b-nlminb-inf-log.txt).
#'
#' The START is passed through unchanged: an objective that is NaN
#' where the optimizer begins is a broken model, not an overshoot, and
#' keeps its warning. Every mapped trial is counted into
#' `nonfinite_trials` on the result, so the event stays measurable.
#'
#' `tally`, an environment holding `n`, is counted into as each trial is
#' mapped rather than when nlminb returns: optimize_obj() keeps one total
#' across its restarts, and a run that raises (and is retried from its
#' best point) never returns a result that could carry its count.
#'
#' @noRd
nlminb_trial_fn <- function(fn, tally = NULL) {
  n_eval <- 0L
  n_mapped <- 0L
  list(
    fn = function(par) {
      v <- fn(par)
      n_eval <<- n_eval + 1L
      if (n_eval > 1L && length(v) == 1L && is.nan(v)) {
        n_mapped <<- n_mapped + 1L
        if (!is.null(tally)) tally$n <- tally$n + 1L
        return(Inf)
      }
      v
    },
    count = function() n_mapped
  )
}

#' One optimizer invocation, normalized to nlminb's result shape.
#' par_units (autoscale) carries per-parameter magnitudes into nlminb's
#' scaling hook; the custom-optimizer contract is unchanged.
#'
#' @noRd
run_optimizer <- function(optimizer, par, fn, gr, lower, upper, control,
                          par_units = NULL, tally = NULL) {
  if (is.function(optimizer)) {
    res <- optimizer(par, fn, gr, lower, upper, control)
    need <- c("par", "objective", "convergence")
    if (!all(need %in% names(res))) {
      frm_stop("A custom optimizer must return par, objective, and ",
               "convergence", call. = FALSE)
    }
    res$message <- res$message %||% ""
    return(res)
  }
  switch(optimizer,
    nlminb = {
      fnw <- nlminb_trial_fn(fn, tally)
      res <- stats::nlminb(par, fnw$fn, gr, control = control,
                           # PORT iterates in scale * par units
                           scale = if (is.null(par_units)) 1 else
                             1 / par_units,
                           lower = lower, upper = upper)
      res$nonfinite_trials <- fnw$count()
      res
    },
    optim = {
      ctl <- control[names(control) %in%
                       c("maxit", "factr", "pgtol", "trace", "REPORT")]
      if (is.null(ctl$maxit)) ctl$maxit <- 1000
      # L-BFGS-B ignores parscale (see ?optim); par_units still govern
      # the gradient-based convergence checks
      r <- stats::optim(par, fn, gr, method = "L-BFGS-B",
                        lower = lower, upper = upper, control = ctl)
      list(par = r$par, objective = r$value,
           convergence = r$convergence, message = r$message %||% "")
    },
    frm_stop("Unknown optimizer '", optimizer,
             "' (use \"nlminb\", \"optim\", or a function)", call. = FALSE)
  )
}

#' A parameter list flattened into the outer vector the optimizer
#' iterates on: the components of names(obj$par), in that order, with
#' the dpars that se() maps to constants removed. Same loop
#' autoscale_units() walks, and the inverse of the one that fills
#' `estimates` from opt$par.
#'
#' @noRd
outer_from_template <- function(tpl, obj, frame) {
  pn <- names(obj$par)
  out <- numeric(0)
  for (cp in unique(pn)) {
    v <- tpl[[cp]]
    if (is.null(v)) return(NULL)
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      v <- v[-frame[["betad_fixed_idx"]]]
    }
    out <- c(out, unname(v))
  }
  if (length(out) != length(pn)) return(NULL)
  stats::setNames(out, pn)
}

#' Starting points to try when the optimizer could not get through at
#' all. Both only ever run after a failure, so a healthy fit pays
#' nothing.
#'
#' The start itself can be unusable. The autoscale pre-fit hands back a
#' warm start, and when the standardized fit ran a correlated block to a
#' perfect correlation the mapped-back template has an infinite
#' objective and a NaN gradient - nlminb dies before its first step,
#' while the cold start make_start() would have built fits the same
#' model. So offer the cold start whenever the fit did not begin there.
#'
#' Or the path from a usable start crosses a hole. Under profile = TRUE
#' beta moves into the inner problem, and the outer objective left
#' behind can be undefined where the optimizer must pass (a Gamma shape
#' intercept on its way to exp(34) does it). The plain Laplace objective
#' over the same model - the one every other mode optimizes - has no
#' such barrier, so its optimum is a starting point on the far side of
#' the hole. Same recipe quad_fit() uses to calibrate the Gauss-Kronrod
#' tape at the Laplace optimum.
#'
#' @noRd
fit_recovery_starts <- function(obj, nll, template, random, map, frame,
                                start, control, bounds, par_units,
                                prior_entries = NULL) {
  out <- list()
  differs <- function(p) {
    !is.null(p) && length(p) == length(obj$par) && all(is.finite(p)) &&
      !isTRUE(all.equal(unname(as.numeric(p)), unname(as.numeric(obj$par))))
  }
  cold <- outer_from_template(make_start(frame, start, prior_entries),
                              obj, frame)
  if (differs(cold)) out[["the cold starting values"]] <- cold
  if (isTRUE(control$profile)) {
    plain <- tryCatch({
      o <- RTMB::MakeADFun(nll, template, random = random, map = map,
                           silent = TRUE)
      # bounds and units are resolved for the profiled parameterization,
      # which has no beta in it; this fit only supplies a starting point
      # and the real optimization below applies both
      op <- optimize_obj(o, control)
      outer_from_template(o$env$parList(op$par), obj, frame)
    }, error = function(e) NULL)
    if (differs(plain)) out[["the Laplace optimum"]] <- plain
  }
  # a start outside the model's bounds is rejected by the optimizer
  lapply(out, function(p) pmin(pmax(p, bounds$lower), bounds$upper))
}

#' Second chance after the optimizer aborts the whole fit.
#'
#' nlminb hands RTMB a step that lands outside the region where the
#' likelihood is defined, RTMB raises "NA/NaN gradient evaluation" from
#' inside the optimizer, and the fit dies - on models that fit perfectly
#' well from a different starting point (an ar1 block under
#' profile = TRUE, a wide smooth under autoscale). The tape is not
#' broken: it still holds the best point the line search reached, and
#' restarting there steps around the hole. quad_fit() already does
#' exactly this when a frozen Gauss-Kronrod rescaling expires the same
#' way. One retry only, and only from a point the optimizer preferred to
#' where it started: a second failure is an objective that is genuinely
#' undefined nearby, not an unlucky step.
#'
#' @noRd
optimizer_from_best <- function(obj, par, e, optimizer, bounds, control,
                                par_units, verbose = 0L, tally = NULL) {
  pb <- obj$env$last.par.best
  # last.par.best spans the joint vector when the model has random
  # effects (profiled parameters included); lfixed() selects the outer
  # block the optimizer actually iterates on
  p0 <- if (is.null(pb)) NULL else pb[obj$env$lfixed()]
  if (is.null(p0) || length(p0) != length(par) || anyNA(p0) ||
      isTRUE(all.equal(unname(as.numeric(p0)), unname(as.numeric(par))))) {
    stop(e)
  }
  if (verbose) {
    vb_say("optimizer failed (", conditionMessage(e),
           "); restarting from the best point it reached")
  }
  run_optimizer(optimizer, stats::setNames(as.numeric(p0), names(par)),
                obj$fn, obj$gr, bounds$lower, bounds$upper,
                control$optCtrl, par_units, tally)
}

#' The single entry point to the optimizer for every fit mode: it runs
#' the objective to an optimum, restarts from there while the gradient
#' stays above `grad_tol`, and keeps the better result. Bounds, autoscale
#' units, and the restart-from-best recovery are applied here, so no
#' caller has to repeat them.
#'
#' The result's `nonfinite_trials` is the total over EVERY optimizer run
#' here, not the count of the run that was kept: a restart that replaces
#' `opt` maps few or no trials, so its own count would erase the ones
#' that got the fit there. `tally` lets a caller that retries a failed
#' call (fit_assembled()'s recovery starts) keep counting across calls.
#'
#' @noRd
optimize_obj <- function(obj, control,
                         bounds = list(lower = -Inf, upper = Inf),
                         par_units = NULL, verbose = 0L,
                         start_par = obj$par, tally = NULL) {
  optimizer <- control$optimizer %||% "nlminb"
  # A model with no free outer parameters - every dpar pinned by a
  # constant and every design zero-column, e.g. y | trials(n) ~ 0 - is
  # already at its optimum. nlminb's PORT front end rejects the empty
  # start vector with "'d' must be a nonempty numeric (double) vector",
  # which names nothing the user wrote, so evaluate the template
  # instead and report a converged degenerate fit. [glmmTMB#1325, #1317]
  if (!length(obj$par)) {
    if (verbose) vb_say("no free parameters; evaluating the template")
    return(list(par = obj$par, objective = as.numeric(obj$fn(obj$par)),
                convergence = 0L,
                message = "no free parameters (degenerate model)"))
  }
  if (is.null(tally)) {
    tally <- new.env(parent = emptyenv())
    tally$n <- 0L
  }
  run <- function(par) {
    tryCatch(run_optimizer(optimizer, par, obj$fn, obj$gr,
                           bounds$lower, bounds$upper, control$optCtrl,
                           par_units, tally),
             error = function(e) optimizer_from_best(obj, par, e, optimizer,
                                                     bounds, control,
                                                     par_units, verbose,
                                                     tally))
  }
  if (verbose) t0 <- vb_now()
  opt <- run(start_par)
  if (verbose) vb_stage("optimize", t0, vb_opt_detail(opt))
  for (i in seq_len(control$restarts)) {
    g <- max(abs(obj$gr(opt$par) * (par_units %||% 1)))
    if (is.finite(g) && g < control$grad_tol) break
    if (verbose) t0 <- vb_now()
    opt2 <- run(opt$par)
    if (verbose) {
      vb_stage(paste0("restart ", i), t0,
               paste0(vb_opt_detail(opt2), ", from max|grad| ",
                      format(g, digits = 3)))
    }
    if (opt2$objective <= opt$objective) opt <- opt2
  }
  # only the nlminb path maps trials; optim and a custom optimizer have
  # no count to report, and NULL says so where 0 would claim a measurement
  if (identical(optimizer, "nlminb")) opt$nonfinite_trials <- tally$n
  opt
}

#' Where a linear predictor's coefficients sit in the outer parameter
#' vector the optimizer iterates on, one index per column of its design,
#' or NULL when they are not there at all (`beta` under
#' `profile = TRUE`, or a component whose length does not line up).
#'
#' Same bookkeeping outer_from_template() and autoscale_units() walk:
#' the components of `names(obj$par)` in order, with the `betad` entries
#' that se() maps to constants removed.
#'
#' @noRd
outer_index_of <- function(lp, par_names, frame) {
  cp <- lp[["par"]]
  off <- which(par_names == cp)
  if (!length(off)) return(NULL)
  keep <- seq_along(frame[["par_template"]][[cp]])
  if (identical(cp, "betad") && length(fx <- frame[["betad_fixed_idx"]])) {
    keep <- setdiff(keep, fx)
  }
  if (length(off) != length(keep)) return(NULL)
  pos <- match(lp[["idx"]], keep)
  if (anyNA(pos)) return(NULL)
  off[pos]
}

#' Coefficients that put a linear predictor at a constant `value`
#' everywhere, or NULL when the design cannot be solved.
#'
#' An intercept takes the value and the other columns take zero, which
#' is the exact answer and keeps such a design bit-for-bit where it
#' was. WITHOUT an intercept there is no single coefficient to carry
#' the value, and skipping the design entirely is what let
#' `alpha ~ 0 + grp` start on the stationary point and never leave it:
#' cell-means syntax spans the same columns as `alpha ~ grp` and lost
#' tens of log-likelihood units to it. The least-squares solution of
#' `X b = value` is the generalization: on `0 + grp` it sets every cell
#' to `value`, and on a design with an intercept it returns exactly the
#' intercept answer.
#'
#' A design whose column space cannot represent a constant gets the
#' closest predictor it can, which the caller checks before using.
#'
#' @noRd
predictor_at <- function(X, value) {
  # a sparse() design is a Matrix and not a base matrix, so the shape
  # is judged by dim(): an is.matrix() gate here silently dropped the
  # initializer of every sparse model
  if (!design_ok(X) || length(value) != 1L || !is.finite(value)) {
    return(NULL)
  }
  icpt <- match("(Intercept)", colnames(X))
  if (!is.na(icpt)) {
    b <- rep(0, ncol(X))
    b[icpt] <- value
    return(b)
  }
  Xd <- tryCatch(as.matrix(X), error = function(e) NULL)
  if (!is.matrix(Xd)) return(NULL)
  b <- tryCatch(base::qr.coef(base::qr(Xd), rep(value, nrow(Xd))),
                error = function(e) NULL)
  if (is.null(b) || length(b) != ncol(X)) return(NULL)
  # an aliased column gets NA from qr.coef; it contributes nothing
  b[is.na(b)] <- 0
  if (!all(is.finite(b))) return(NULL)
  as.numeric(b)
}

#' Is this a usable two-dimensional design? Base matrix or Matrix
#' alike, because `sparse()` predictors are the latter.
#'
#' @noRd
design_ok <- function(X) {
  d <- dim(X)
  !is.null(X) && length(d) == 2L && d[1L] > 0L && d[2L] > 0L
}

#' A linear predictor as a plain numeric vector, for a base matrix or a
#' Matrix. `as.vector()` alone does not flatten a Matrix product.
#'
#' @noRd
predictor_eta <- function(X, b) {
  tryCatch(as.numeric(as.matrix(X %*% b)), error = function(e) NULL)
}

#' Coefficients that move a linear predictor AWAY from zero by `value`
#' in root mean square, for a design whose columns cannot represent a
#' constant at all (`alpha ~ 0 + x` with a centred x is the case).
#'
#' There `predictor_at()` returns the least-squares fit of a constant,
#' which on such a design is the zero vector, so the escape would start
#' on the very point it is escaping. Moving along the design's largest
#' column instead gives a start that is genuinely elsewhere, scaled so
#' the predictor's spread is the size the family asked for.
#'
#' @noRd
predictor_away <- function(X, value) {
  if (!design_ok(X) || length(value) != 1L || !is.finite(value) ||
        value == 0) {
    return(NULL)
  }
  nrm <- tryCatch(as.numeric(sqrt(Matrix::colSums(X^2))),
                  error = function(e) NULL)
  if (is.null(nrm) || length(nrm) != ncol(X) ||
        !any(is.finite(nrm) & nrm > 0)) {
    return(NULL)
  }
  j <- which.max(ifelse(is.finite(nrm), nrm, -Inf))
  s <- nrm[j] / sqrt(nrow(X))
  if (!is.finite(s) || s <= 0) return(NULL)
  b <- rep(0, ncol(X))
  b[j] <- value / s
  if (!all(is.finite(b))) return(NULL)
  b
}

#' The starting points that a family's declared stationary points ask
#' for, given where the optimizer actually stopped.
#'
#' A `post$stationary` entry says the likelihood has a stationary point
#' at a known place, at every sample. The skew normal's alpha = 0 is the
#' case this exists for: the information is singular there, so the
#' gradient AND the curvature vanish and a converged fit is
#' indistinguishable from a maximum by any test of the optimum itself.
#' The only way to tell is to look from somewhere else.
#'
#' Empty unless the fit landed on a declared point, so a fit that did
#' not pays one comparison per declaration and nothing else.
#'
#' @noRd
stationary_escapes <- function(obj, opt, frame) {
  out <- list()
  par_names <- names(obj$par)
  for (lp in frame[["linpreds"]]) {
    if (!is.null(lp[["constant"]])) next
    resp <- frame[["spec"]]$responses[[lp[["resp"]]]]
    # `[[` on an ATOMIC declaration raises "subscript out of bounds",
    # so the shape of the whole field is checked before it is indexed:
    # post$stationary = "alpha" used to kill the fit outright
    st_all <- resp$family[["post"]][["stationary"]]
    if (!is.list(st_all)) next
    st <- st_all[[lp[["dpar"]]]]
    # a malformed declaration from a third-party family must not reach
    # the parameter vector: a character `from` would coerce the whole
    # start, and a non-finite `tol` would fire on every fit
    if (!is.list(st)) next
    from <- st[["from"]]
    at <- st[["at"]] %||% 0
    tol <- st[["tol"]] %||% 0.05
    if (!is.numeric(from)) next
    # A non-finite restart is dropped HERE rather than failing later
    # inside predictor_at(): `from = c(Inf, -2)` then yields one start
    # by decision and not by accident, which is what the Rd says.
    from <- from[is.finite(from)]
    if (!length(from) ||
          !is.numeric(at) || length(at) != 1L || !is.finite(at) ||
          !is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
          tol <= 0) next
    idx <- outer_index_of(lp, par_names, frame)
    if (is.null(idx)) next
    X <- lp[["X"]]
    if (!design_ok(X) || ncol(X) != length(idx)) next
    # The stationary point is a property of the PREDICTOR, not of the
    # coefficients: `alpha ~ 0 + grp` sits on it with every cell at
    # zero and no coefficient named "(Intercept)" to read. Judging the
    # fitted predictor covers both spellings, and on an intercept-only
    # design it is the same comparison as before.
    eta <- predictor_eta(X, opt$par[idx])
    if (is.null(eta) || !all(is.finite(eta)) ||
          max(abs(eta - at)) >= tol) next
    # A start counts as an escape when it moves the predictor a fair
    # share of what the family asked for. The yardstick is the
    # requested move and NOT `tol`: a family that declares a wide
    # window would otherwise reject every start it asked for.
    moved <- function(b, v) {
      if (is.null(b)) return(FALSE)
      e <- predictor_eta(X, b)
      if (is.null(e) || !all(is.finite(e))) return(FALSE)
      max(abs(e - at)) >= 0.5 * abs(v - at)
    }
    for (v in from) {
      b <- predictor_at(X, v)
      # a design that cannot represent a constant would "escape" to the
      # point it is already on, which is the one start that never
      # leaves; move along the design's own largest column instead
      if (!moved(b, v)) b <- predictor_away(X, v - at)
      if (!moved(b, v)) next
      p <- opt$par
      p[idx] <- b
      out[[length(out) + 1L]] <- p
    }
  }
  out
}

#' Refit from each escape start and keep the best optimum.
#'
#' Silent. The optimizer's other recoveries (restarts, the
#' restart-from-best after a failure) say nothing either, and this one
#' would speak on roughly half of all skew-normal fits of symmetric
#' data, where alpha = 0 is reached legitimately. `verbose = TRUE`
#' reports it with the other optimizer stages, and the result carries
#' `stationary_escape` so the event stays measurable.
#'
#' @noRd
escape_stationary <- function(obj, opt, frame, control, bounds,
                              par_units = NULL, verbose = 0L,
                              tally = NULL) {
  starts <- stationary_escapes(obj, opt, frame)
  if (!length(starts)) return(opt)
  if (verbose) t0 <- vb_now()
  # what the tally already held before this escape: quad_fit() counts
  # into no tally, so `tally$n` alone is the escape's own trials and
  # would REPLACE a quadrature fit's total with a smaller number
  n_before <- if (is.null(tally)) 0L else tally$n
  best <- opt
  for (p in starts) {
    p <- pmin(pmax(p, bounds$lower), bounds$upper)
    o <- tryCatch(optimize_obj(obj, control, bounds, par_units,
                               start_par = p, tally = tally),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$objective) &&
          o$objective < best$objective) {
      best <- o
    }
  }
  # recorded whenever the escape RAN, gain or no gain: a fit that paid
  # for two restarts and kept its own optimum is the measurement that
  # says what the fallback costs
  gain <- opt$objective - best$objective
  best$stationary_escape <- c(starts = length(starts), gain = gain)
  # the count a fit reports is the total over every run of its
  # objective, so the escape ADDS its own trials to whatever the run
  # being kept already carried rather than overwriting the total
  if (!is.null(tally) && !is.null(best$nonfinite_trials)) {
    added <- max(tally$n - n_before, 0L)
    base_n <- opt$nonfinite_trials %||% 0L
    best$nonfinite_trials <- base_n + added
  }
  if (verbose) {
    vb_stage("stationary escape", t0,
             paste0(length(starts), " restart",
                    if (length(starts) != 1L) "s", ", gain ",
                    format(gain, digits = 3)))
  }
  best
}

#' Does an optimizer failure message read like an undefined objective?
#'
#' The diagnosis is only attached to errors that actually look numerical.
#' Everything raised inside the wrapped expression lands here, including
#' a misspelled optimizer, a custom optimizer that broke its return
#' contract, and errors from user code in a custom family - none of which
#' says anything about the likelihood. Reporting those as an undefined
#' objective sends the user to the wrong remedies and buries the real
#' message under advice.
#'
#' @noRd
fit_numerical_error <- function(msg) {
  # nlminb ("NA/NaN function evaluation"), RTMB ("NA/NaN gradient
  # evaluation") and optim's L-BFGS-B finiteness checks are the
  # vocabulary of an objective that came back undefined
  grepl("NA/NaN|NaN|non-?finite|not finite|infinite|Inf\\b", msg,
        ignore.case = TRUE)
}

#' Runs the optimizer call and rewrites any failure into a message that
#' names the model and what to try next.
#'
#' Whatever the optimizer throws reaches the user raw otherwise. RTMB's
#' "NA/NaN gradient evaluation" and nlminb's PORT strings name neither
#' the model nor anything to try, and they arrive with the call stack of
#' stats::nlminb, so the user cannot even tell which of their models
#' failed. Wrap the optimizer call once and say both.
#'
#' A nonlinear model gets the sharper advice, because there the cause is
#' nearly always the start: make_start() can only seed intercepts through
#' a family's init_dpars, and a nonlinear mu has no design of its own, so
#' an nl fit begins at zero, where most nonlinear forms are flat,
#' singular, or undefined. `[brms#734 doctrine]`
#'
#' The condition carries a class of its own, `frmtmb_fit_error`, so the
#' autoscale pre-fit (an inner `fit_assembled()` that has already been
#' through here) is rethrown rather than wrapped a second time.
#'
#' @noRd
fit_error_context <- function(spec, start, REML, control, quadrature,
                              prior, expr) {
  nl <- any(vapply(spec$responses,
                   function(r) length(r$nlpars) > 0L, TRUE))
  nl_start <- nl && is.null(start)
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      if (inherits(e, "frmtmb_fit_error")) stop(e)
      emsg <- conditionMessage(e)
      numerical <- fit_numerical_error(emsg)
      msg <- if (nl_start && numerical) {
        paste0("The nonlinear fit failed from its default starting ",
               "values (", emsg, "). Nonlinear models ",
               "need starting values in the right region: ",
               "par_template(formula, data) names the parameters, and ",
               "the names go straight back, e.g. ",
               "start = list(beta = c(ult = 5000)). A normal() or ",
               "student_t() prior on a nonlinear parameter places its ",
               "start from the prior location instead")
      } else if (numerical) {
        paste0("The optimizer failed on this model (",
               vb_fit_detail(spec, REML, control, quadrature, prior),
               "): ", emsg,
               ". The likelihood was undefined or unbounded somewhere ",
               "the optimizer stepped. Refit with verbose = TRUE to see ",
               "which stage broke, try another optimizer ",
               "(frmtmb_control(optimizer = \"optim\")), or start the ",
               "fit nearer the optimum with the `start` argument of ",
               "frm()")
      } else {
        # the real message first, then only what is certainly true:
        # which model it came from
        paste0(emsg, " (raised while fitting: ",
               vb_fit_detail(spec, REML, control, quadrature, prior),
               ")")
      }
      frm_stop(msg, call. = FALSE, class = "frmtmb_fit_error")
    }),
    warning = function(w) {
      # nlminb's own "NA/NaN function evaluation" is the optimizer
      # noticing the same undefined objective the error above names;
      # letting both through would report the failure twice
      if (nl_start && grepl("NA/NaN", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' What a `start` component name means when a NONLINEAR PARAMETER of
#' the same model carries the same name.
#'
#' `start` and `newparams` are keyed by par-template COMPONENT
#' (`beta`, `betad`, `b`, `theta`, `thetaac`, `thetar`, `miss`), and a
#' nonlinear parameter may be named after any of them and fit
#' perfectly: its coefficients live inside `beta` as `b_(Intercept)`,
#' not in the component called `b`. So `start = list(b = 1)` on such a
#' model sets the RANDOM-EFFECT VECTOR, and reported the length of that
#' vector ("start$b must have length 120") - a message about an object
#' the caller never meant. Neither reading can be dropped, so name both.
#'
#' `arg` is the caller's own spelling, so a `frm_simulate()` message
#' does not answer a `newparams$b` error by talking about `start$b`.
#'
#' @noRd
nl_start_collision_msg <- function(nm, tpl, arg = "start") {
  what <- switch(nm,
    beta = "the fixed-effect coefficients of the location predictors",
    betad = "the coefficients of the distributional parameters",
    b = "the random-effect vector (the conditional modes)",
    theta = "the covariance parameters",
    thetaac = "the autocorrelation parameters",
    thetar = "the residual-correlation parameters",
    miss = "the imputed missing values",
    "a parameter-template component")
  pfx <- paste0(nm, "_")
  own <- character(0)
  for (cp in c("beta", "betad")) {
    cn <- names(tpl[[cp]]) %||% character(0)
    hit <- cn[startsWith(cn, pfx)]
    if (length(hit)) {
      own <- c(own, paste0(cp, " = c(`", hit[1L], "` = ...)"))
      break
    }
  }
  paste0("`", arg, "$", nm, "` sets ", what, ", not the nonlinear ",
         "parameter '", nm, "' of the same name: `", arg, "` is keyed ",
         "by parameter-template component, and a nonlinear parameter's ",
         "coefficients sit INSIDE one of those components",
         if (length(own)) {
           paste0(" - here as ", arg, " = list(", own[1L], ")")
         },
         ". par_template() lists both")
}

#' The response with its primary predictor taken out, by least squares
#' on that predictor's own design matrix.
#'
#' A family whose starting value describes the SHAPE of the residual
#' cannot read that shape off the response. skew_normal()'s alpha is
#' the case: a covariate large enough to dominate the marginal
#' distribution gives the raw response the opposite skew, alpha then
#' starts on the wrong side of zero, and zero is a stationary point
#' with singular information, so the optimizer stops there with a
#' clean convergence code (dev/skewinit-findings.md).
#'
#' Least squares, not the model: this runs before any parameter has a
#' value, the answer only has to get a sign right, and a fit is not
#' affordable here. Grouping structure is not removed either, so a
#' random effect's variance stays in the residual; it widens the
#' residual without biasing its third moment.
#'
#' An intercept-only predictor, and anything that does not resolve to a
#' design at all, gives back the RESPONSE rather than the centred
#' response. There is nothing to take out in either case, and an
#' initializer then computes the same floating-point number it computed
#' from `y` before, so a model without a mu covariate does not move at
#' all. Centring moves `sd()` in its last bit on about one sample in a
#' thousand, and a start that moves in its last bit moves every iterate
#' after it.
#'
#' @noRd
mu_residuals <- function(frame, resp, cache = NULL) {
  # linpreds key a response by name on a multivariate model and by
  # position on a univariate one; an environment takes only a string
  key <- as.character(resp)
  if (!is.null(cache) && !is.null(cache[[key]])) return(cache[[key]])
  y <- frame[["y"]][[resp]]
  out <- NULL
  if (is.numeric(y) && is.null(dim(y))) {
    fam <- frame[["spec"]]$responses[[resp]]$family
    prim <- (fam[["primary_dpars"]] %||% "mu")[1L]
    for (l in frame[["linpreds"]]) {
      if (!identical(l[["resp"]], resp) ||
            !identical(l[["dpar"]], prim)) next
      X <- if (is.null(l[["constant"]])) l[["X"]] else NULL
      # offset() is part of the predictor and is not in X, so a model
      # that enters its covariate as an offset kept the covariate's
      # skew in the "residual" and started alpha from the wrong side
      # exactly as the raw response did
      off <- l[["offset"]]
      has_off <- is.numeric(off) && length(off) == length(y) &&
        all(is.finite(off))
      yo <- if (has_off) as.numeric(y) - as.numeric(off) else as.numeric(y)
      if (design_ok(X) && nrow(X) == length(y)) {
        out <- if (!has_off && ncol(X) == 1L &&
                     all(as.numeric(X[, 1L]) == 1)) {
          # an intercept has nothing to take out; see this function's
          # header for why the response and not the centred response
          as.numeric(y)
        } else {
          # qr() pivots, so a rank-deficient design still gives
          # residuals; as.matrix() because a sparse() design is a
          # Matrix, which base::qr() does not take
          tryCatch(as.numeric(base::qr.resid(base::qr(as.matrix(X)), yo)),
                   error = function(e) NULL)
        }
      } else if (has_off) {
        out <- yo                    # an offset with no design of its own
      }
      break
    }
  }
  # nothing resolved to a design, so nothing is taken out and the
  # answer is the response, which is what an initializer saw before
  if (is.null(out) || length(out) != length(y) || anyNA(out)) {
    out <- as.numeric(y)
  }
  if (!is.null(cache)) cache[[key]] <- out
  out
}

#' The cold starting values: the parameter template with each linear
#' predictor's intercept seeded from the family's own initializer, then
#' any nonlinear parameter a prior's location places, then whatever the
#' user gave in `start` written over the result. Component names are
#' checked here, where the error can still name the template.
#'
#' `prior_entries` are the resolved prior entries of the fit; see
#' prior_nl_starts() for why only nonlinear coefficients read them.
#' `announce` is set at the one call site per user-level fit, so the
#' autoscale pre-fit and the recovery restarts do not repeat the
#' message.
#'
#' @noRd
make_start <- function(frame, start, prior_entries = NULL,
                       announce = FALSE) {
  tpl <- frame[["par_template"]]
  # computed at most once per response, and only for a family that asks
  resid_cache <- new.env(parent = emptyenv())
  for (lp in frame[["linpreds"]]) {
    if (!is.null(lp[["constant"]])) next   # mapped; keep link(constant)
    resp <- frame[["spec"]]$responses[[lp[["resp"]]]]
    init_fn <- resp$family[["init_dpars"]][[lp[["dpar"]]]]
    if (is.null(init_fn)) next
    icpt <- match("(Intercept)", colnames(lp[["X"]]))
    # ONLY a dpar that declares a stationary point gets its start
    # spread over a design with no intercept. Every other dpar keeps
    # the old rule, an intercept or nothing, and that restriction is
    # the point: placing a least-squares start for every family read
    # as an improvement and is NOT safe. The value is right on the
    # PREDICTOR scale and can be ruinous on the PARAMETER scale,
    # because nlminb judges its step relative to the parameter. On
    # `ypois ~ 0 + xt` with xt at 1e-6 it starts the coefficient at
    # 1.6e6, nlminb reports X-convergence before moving, and the fit
    # lands 1967 log-likelihood units below glm() without a warning.
    # dev/test-backlog.md has the measurement and what a safe general
    # version would need.
    st_all <- resp$family[["post"]][["stationary"]]
    declares <- is.list(st_all) && is.list(st_all[[lp[["dpar"]]]])
    if (is.na(icpt) && !declares) next
    raw <- if (length(formals(init_fn)) >= 3L) {
      init_fn(frame[["y"]][[lp[["resp"]]]],
              frame[["aterm_values"]][[lp[["resp"]]]],
              mu_residuals(frame, lp[["resp"]], resid_cache))
    } else {
      init_fn(frame[["y"]][[lp[["resp"]]]],
              frame[["aterm_values"]][[lp[["resp"]]]])
    }
    val <- lp[["link"]]$linkfun(raw)
    b <- if (is.finite(val)) {
      if (is.na(icpt)) predictor_at(lp[["X"]], val) else NULL
    }
    if (!is.null(b)) {
      tpl[[lp[["par"]]]][lp[["idx"]]] <- b
    } else if (is.finite(val)) {
      # the intercept path, byte for byte what it was before this lane
      if (!is.na(icpt)) tpl[[lp[["par"]]]][lp[["idx"]][icpt]] <- val
    } else if (announce) {
      # A bounded link sends an init at or past its bound to Inf, and
      # the value used to be dropped without a word: the fit then
      # started from the template default and looked like a slow or
      # failed optimization with nothing pointing at the cause.
      # Announced on the same flag as the prior-start message, so the
      # autoscale pre-fit and the recovery restarts do not repeat it.
      frm_warning("Starting value ", format(raw[1L]), " for ", lp[["dpar"]],
                  " is ", format(val[1L]), " through its ", lp[["link"]]$name,
                  " link, so it was ignored and that intercept starts from ",
                  "zero on the link scale. Give init_dpars a value inside ",
                  "the link's range, or pass start =", call. = FALSE)
    }
  }
  placed <- prior_nl_starts(frame, prior_entries)
  # [[ ]] throughout: $beta would partial-match nothing here today, but
  # the template's component names are a moving set
  for (p in placed) tpl[["beta"]][p$idx] <- p$value
  claimed <- integer(0)
  if (!is.null(start)) {
    nl_named <- unique(unlist(lapply(frame[["spec"]]$responses,
                                     function(r) r$nlpars %||% character(0))))
    for (nm in names(start)) {
      if (!nm %in% names(tpl)) {
        frm_stop("Unknown start component: '", nm, "' (template has: ",
                 paste(names(tpl), collapse = ", "), ")", call. = FALSE)
      }
      if (nm %in% nl_named) {
        # the component reading is the one that applies, but it is not
        # the one the caller is likely to have meant
        msg <- nl_start_collision_msg(nm, tpl)
        tpl[[nm]] <- tryCatch(
          resolve_start_component(tpl[[nm]], start[[nm]], nm),
          # the length and name errors describe the COMPONENT, so on
          # this model they describe the wrong object; carry the
          # collision into the error rather than leaving it to a
          # warning the error would outrun
          error = function(e) {
            frm_stop(conditionMessage(e), ". ", msg, call. = FALSE)
          })
        frm_warning(msg, call. = FALSE)
        next
      }
      tpl[[nm]] <- resolve_start_component(tpl[[nm]], start[[nm]], nm)
    }
    claimed <- start_claimed_idx(frame[["par_template"]][["beta"]],
                                 start[["beta"]], "beta")
  }
  if (announce && length(placed)) {
    # what `start` set is the user's doing, not the prior's
    kept <- Filter(function(p) !p$idx %in% claimed, placed)
    if (length(kept)) {
      frm_message("Nonlinear starting values placed at the prior locations: ",
                  paste(paste0(names(kept), " = ",
                               vapply(kept, function(p) format(p$value), "")),
                        collapse = ", "),
                  ". Give `start` to choose your own.")
    }
  }
  tpl
}

#' The post-fit verdict: the optimizer status, the maximum absolute
#' gradient, and, when a report is already there, the Hessian and the
#' standard errors. Each failure is a separate warning, so
#' `suppressWarnings()` silences them and the fit is still returned.
#'
#' @noRd
check_convergence <- function(fit, control) {
  # collected first, warned second, so the verbose summary line can
  # report how many diagnostics the fit raised
  msgs <- character(0)
  if (fit$opt$convergence != 0) {
    # a nonlinear model that started at zero is the likeliest cause, and
    # `start` is the only lever, so name it here [brms#734]
    nl_hint <- if (any(vapply(fit$spec$responses,
                              function(r) length(r$nlpars) > 0L, TRUE))) {
      paste0(". This is a nonlinear model; unless `start` or a ",
             "located prior placed them, the fit began at zero, which ",
             "is rarely in the right region (see par_template())")
    } else ""
    msgs <- c(msgs, paste0("Optimizer did not report convergence: ",
                           fit$opt$message, nl_hint))
  }
  # under autoscale the gradient is judged in the same natural units
  # the optimizer used (a 1e6-scale column bounds its coefficient's
  # absolute gradient near machine noise times 1e6)
  g <- if (!length(fit$opt$par)) NA_real_ else {
    try(max(abs(fit$obj$gr(fit$opt$par) * (fit$par_units %||% 1))),
        silent = TRUE)
  }
  if (inherits(g, "try-error")) g <- NA_real_
  if (!is.null(fit$importance)) {
    # An importance-corrected objective is a Monte Carlo estimate, and
    # its gradient carries an O(N^-1/2) error that is exactly zero only
    # where the estimator itself is exact (a gaussian response at its
    # own anchor). No optimizer drives that below `grad_tol`, so
    # judging this fit by it would warn on every correct fit. The fit
    # has its OWN convergence criterion, the round-to-round parameter
    # move, and its own accuracy report, the effective sample sizes and
    # the Monte Carlo standard error; the gradient stays visible in
    # `fit$importance$grad`.
    if (isTRUE(fit$importance$capped)) {
      # Two different fits reach `capped`, and only one of them wants
      # more rounds. A stalled iteration takes the SAME step every
      # round, so the shift it reports is that step times the round
      # count and more rounds buy a proportionally larger wrong
      # answer; telling that user to raise the cap is advice about a
      # number that is not an estimate. imp_stalled() separates them.
      msgs <- c(msgs, if (imp_stalled(fit$importance$moves)) {
        paste0("The importance correction moved by the same amount, ",
               format(fit$importance$moved, digits = 3), ", in every ",
               "one of its ", fit$importance$rounds, " rounds, so it ",
               "has not converged and its total shift is that step ",
               "times the round count rather than an estimate. The ",
               "step is a property of the draws and not of the data, ",
               "so raising frmtmb_control(importance_rounds =) would ",
               "only move the estimates proportionally further. A ",
               "variance component the Laplace fit has already ",
               "collapsed does this, because nothing is left to ",
               "reweight: check VarCorr() before reading the ",
               "corrected estimates")
      } else {
        paste0("The importance correction used all ",
               fit$importance$rounds, " of its rounds and ",
               "the estimates were still moving by ",
               format(fit$importance$moved, digits = 3),
               " at the last one. Raise ",
               "frmtmb_control(importance_rounds =), or ",
               "raise the draw count so each round lands ",
               "in the same place")
      })
    }
  } else if (is.finite(g) && g > control$grad_tol) {
    msgs <- c(msgs, paste0("Large maximum absolute gradient at the ",
                           "optimum (", format(g, digits = 3),
                           "); the fit may not have converged. ",
                           "diagnose() names the offending parameter; ",
                           "see the 'Convergence problems' section of ",
                           "vignette('diagnostics') for the remedies"))
  }
  # Covariance verdicts are only known once sdreport has run (se =
  # TRUE); the lazy path surfaces them through vcov()/summary()/
  # diagnose() instead. pdHess does not imply usable standard errors:
  # the Cholesky it comes from succeeds on a Hessian LAPACK's solver
  # then refuses as computationally singular, and cov.fixed is NaN.
  sdr <- fit$cache$sdr
  if (!is.null(sdr) && !is.null(sdr$pdHess) && !isTRUE(sdr$pdHess)) {
    # "overparameterized" is one of two causes and the wrong one when a
    # direction is flat; flat_par_note() names the parameters when it is
    msgs <- c(msgs, paste0("Hessian is not positive definite; standard ",
                           "errors are unreliable. The model may be ",
                           "overparameterized", flat_par_note(fit)))
  } else if (!is.null(sdr) && length(sdr$cov.fixed) &&
             any(!is.finite(sdr$cov.fixed))) {
    msgs <- c(msgs, paste0("Some standard errors are not finite: the ",
                           "covariance could not be recovered from the ",
                           "Hessian", flat_par_note(fit),
                           ". diagnose() names the offending ",
                           "parameters; see the 'Convergence problems' ",
                           "section of vignette('diagnostics')"))
  }
  for (m in msgs) frm_warning(m, call. = FALSE)
  invisible(list(grad = g, warnings = msgs))
}
