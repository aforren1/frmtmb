# brms-style prior specification: the vocabulary (set_prior, prior,
# the prior_* constructors, brmsprior translation) and the fit-route
# machinery that resolves it for frm()'s MAP penalty, frm_simulate()
# and par_template(). The sampling-only side (default priors,
# non-centering, the tmbstan bridge) is in the frmtmb.sample package;
# the resolution machinery below is exported for it to reach.

#' Set up priors brms-style
#'
#' Builds prior specifications with brms spelling:
#' `set_prior("normal(0, 5)", class = "b")`. Combine several with `+` or
#' `c()`. Distributions: `normal(mu, sd)`, `student_t(df, mu, sd)`,
#' `cauchy(mu, sd)`, `exponential(rate)`, `logistic(mu, s)`,
#' `gamma(shape, rate)`, `lkj(eta)`; an empty string sets bounds only.
#'
#' Classes and their scales:
#' - `"b"`: population-level coefficients of `dpar` (default: the
#'   location parameters), excluding the intercept; narrow to one
#'   coefficient with `coef`. Link scale.
#' - `"Intercept"`: the intercept of `dpar`, at the MEAN of that
#'   sub-formula's predictors, which is the intercept brms's own
#'   `Intercept` prior constrains. Link scale. See Where an intercept
#'   prior lands. On an ordinal family the thresholds are the intercept,
#'   here as in brms, so this class addresses the whole threshold
#'   vector; see Ordinal thresholds.
#' - `"sd"`: random-effect standard deviations (and smoothing SDs), on
#'   the NATURAL sd scale with the log-Jacobian applied, so
#'   `set_prior("exponential(1)", class = "sd")` means what it says;
#'   narrow with `group`. As in brms, `resp`, `dpar` and `nlpar` select
#'   the standard deviations of that response, distributional parameter
#'   and nonlinear parameter only, and leaving one empty selects the
#'   empty one: with no `dpar` the specification reaches the location
#'   parameter's blocks and not a `phi ~ (1 | g)` block, and in a
#'   multivariate model a specification with no `resp` reaches nothing
#'   and is refused. A block that spans several predictors,
#'   `(1 | q | g)` in both `mu` and `sigma`, is the one exception, as it
#'   is in brms: there a specification that leaves a field empty reaches
#'   every standard deviation of the block, and a more specific one
#'   takes over its own.
#' - `"cor"`: the CORRELATION of a random-effect block, as a whole.
#'   `lkj(eta)` only, and it addresses a BLOCK the way class `"sd"`
#'   does, by `group`; `set_prior("lkj(2)", class = "cor")` covers every
#'   correlated block of the model, which is brms's spelling. See
#'   The LKJ prior below.
#' - `"ar"`, `"ma"`, `"cosy"`, `"cortime"`: the R-side residual
#'   correlation of an `ar()`, `ma()`, `arma()`, `cosy()` or `unstr()`
#'   term, under brms's own class names. `"ar"`, `"ma"` and `"cosy"`
#'   take an ordinary density on the NATURAL coefficient with that
#'   map's Jacobian applied, as class `"sd"` does; `"cortime"` takes
#'   `lkj(eta)` on an `unstr()` time correlation as a whole. Narrow to
#'   one response with `resp`. See Residual correlation below.
#' - `"rescor"`: the residual correlation BETWEEN responses of a
#'   multivariate model (`set_rescor(TRUE)`), as a whole. `lkj(eta)`
#'   only, as brms spells it, and no `resp`.
#' - `"theta"`: raw internal covariance parameters (escape hatch).
#'   `coef` names one by its internal name and spans all three
#'   covariance components: `"theta_2"` for a random-effect block,
#'   `"thetaac_1"` for a residual autocorrelation, `"thetar_1"` for a
#'   residual correlation. This is the one spelling that reaches a
#'   single parameter of a structure whose natural coefficients are not
#'   free of one another.
#' - a DISTRIBUTIONAL parameter's own name (`"sigma"`, `"shape"`,
#'   `"phi"`, `"zi"`, `"nu"`, ...): a density on that parameter ITSELF,
#'   on its NATURAL scale, with the inverse link's log-Jacobian
#'   applied, so `set_prior("student_t(3, 0, 2.5)", class = "sigma")`
#'   means what it says. Available where the parameter has no predictor
#'   of its own; see A distributional parameter's own class.
#'
#' In a MULTIVARIATE model, every specification of class `"b"`,
#' `"Intercept"`, `"sd"`, `"ar"`, `"ma"`, `"cosy"` or `"cortime"`, and
#' every distributional parameter's own class, names its response with
#' `resp`. A specification without `resp` is refused, as brms refuses
#' it: write one specification per response. Classes `"rescor"` and
#' `"theta"` take no `resp`, and `"cor"` needs none. On a model with ONE
#' response, `resp` is refused, as in brms.
#'
#' Two spellings are frmtmb's own, and brms refuses both: `"cor"` with
#' `resp`, which narrows the prior to that response's blocks, and
#' `"Intercept"` with `nlpar`, which addresses the intercept of one
#' nonlinear parameter (see Nonlinear parameters). A class `"b"`
#' specification that reaches no coefficient, because the predictor it
#' names has no population-level slope (`y ~ 1`), is refused, as in
#' brms. In the same way, a class
#' `"b"` or `"Intercept"` specification on a NONLINEAR location names
#' its nonlinear parameter with `nlpar`; see Nonlinear parameters.
#' And where a family's location is SEVERAL distributional parameters,
#' as in [categorical()], [multinomial()] and [mixture()], a class
#' `"b"`, `"Intercept"` or `"sd"` specification names one of them with
#' `dpar` (`"mub"`, `"mu1"`), and [default_prior()] lists the rows that
#' way, as brms does.
#'
#' Two specifications for the same slot (the same class, coef, group,
#' resp, dpar and nlpar) are refused wherever a prior is passed in, as
#' brms refuses them; that includes a density followed by a bounds-only
#' specification, so write the density and its bounds in one call.
#' When specifications for DIFFERENT slots reach the same parameter,
#' the more specific one applies whatever order they are written in, as
#' in brms: a `coef` specification over its class, and a `group`
#' specification over a class-wide `"sd"` or `"cor"` one. Between
#' classes the later one applies: a class `"theta"` prior on a position
#' a `"cor"` prior covers replaces that whole LKJ term if it comes later,
#' and the other way round.
#' `lb`/`ub` become hard bounds. See Hard bounds.
#'
#' @section Where an intercept prior lands:
#' brms centers its design matrix and constrains the intercept at the
#' MEAN of the predictors, recovering the reported one as
#' `b_Intercept = Intercept - dot_product(means_X, b)`. frmtmb
#' parameterizes by the intercept at zero and evaluates a class
#' `"Intercept"` density at `b0 + means_X'b`, so the two packages
#' constrain the same quantity. The map between the two
#' parameterizations is unit triangular and carries no Jacobian.
#'
#' This matters whenever a predictor is not centered. On
#' `Reaction ~ Days`, where `mean(Days)` is 4.5, the intercept at zero
#' is strongly correlated with the slope and the intercept at the mean
#' is orthogonal to it, so a prior on the first biases the slope and a
#' prior on the second does not. Earlier releases used the intercept
#' at zero, and `brms::get_prior()`'s own default there moved the slope
#' by 0.068 standard errors; it now moves it by 0.00003.
#'
#' Every sub-formula with an intercept is centered separately, as brms
#' does (`means_X`, `means_X_sigma`, ...). A NONLINEAR parameter's
#' sub-formula is not centered on either side, and neither are a
#' smooth's unpenalized columns or a `mo()` term, which sit outside
#' brms's `Xc` as well. To put a density on the intercept at zero, name
#' it as a coefficient instead: `class = "b", coef = "Intercept"`,
#' which is also how a `brms::bf(center = FALSE)` model's prior arrives.
#'
#' @section Ordinal thresholds:
#' `cumulative()`, `sratio()`, `cratio()` and `acat()` have no intercept
#' column: the thresholds replace it. brms priors them as its
#' `Intercept` class and so does frmtmb, so
#' `set_prior("student_t(3, 0, 2.5)", class = "Intercept")` on an
#' ordinal model addresses the whole threshold vector. It addresses the
#' THRESHOLDS, at the mean of the predictors, with the log-Jacobian of
#' the map from frmtmb's internal storage; `cumulative()` and
#' `sratio()` hold `(tau_1, log increments)`, which is the same map
#' Stan's `ordered` type applies, and `cratio()` and `acat()` hold the
#' thresholds themselves. `lb`/`ub` are refused there, because one
#' number cannot box a whole vector of ordered thresholds.
#'
#' `prior = list(tau_raw = prior_normal(0, 5))` reaches the same
#' parameters on the INTERNAL scale, one entry per threshold, which is
#' the escape hatch to use when the increments rather than the
#' thresholds are what a prior is about.
#'
#' @section Hard bounds:
#' `lb`/`ub` are how a box constraint is written. A specification may
#' carry bounds alone (`prior = ""`), a distribution alone, or both. A
#' bounds-only specification for a slot another specification already
#' gives a density is refused as a duplicate, as in brms; a bound on a
#' coefficient beside a class-wide density is a different slot and
#' boxes that coefficient under the class density.
#'
#' A bound is addressed exactly like the distribution beside it, so
#' `set_prior("", nlpar = "guess", lb = 0, ub = 1)` bounds the nonlinear
#' parameter `guess`, and `dpar`, `resp`, `group` and `coef` narrow a
#' bound the same way they narrow a density. As with a distribution,
#' class `"b"` with `nlpar` covers EVERY coefficient of that parameter;
#' `coef` picks out one.
#'
#' The scale is the parameter's own: class `"sd"` bounds a standard
#' deviation on the sd scale (frmtmb stores its log), classes `"ar"`,
#' `"ma"` and `"cosy"` bound the natural coefficient of a residual
#' structure, and class `"theta"` bounds an internal covariance
#' parameter itself, one position at a time with `coef = "theta_2"`,
#' `"thetaac_1"` or `"thetar_1"`. Everything else is bounded on the
#' internal (link) scale, so a bound on a log-linked dispersion is a
#' bound on its logarithm.
#'
#' In [frm()] a bound is a box constraint handed to the optimizer; in
#' `frmtmb.sample::frm_sample()` it becomes one of Stan's constrained
#' transforms. Both take this spelling and no other: the `lower`/`upper`
#' arguments of releases before 0.49 are gone rather than aliased, and a
#' call still using them fails as an unused argument. Every outer
#' parameter they could reach has a class here, down to a single
#' internal covariance parameter.
#'
#' Where the two spellings differed, this one broadcasts: a bound
#' carried by `nlpar =` covers every coefficient of that parameter, the
#' way a prior does, and `coef` narrows it to one. When a class-wide
#' and a coefficient-specific specification both bound a parameter, the
#' coefficient-specific one applies.
#'
#' @section Residual correlation:
#' frmtmb holds an `ar()`, `ma()`, `arma()`, `cosy()` or `unstr()`
#' residual block in one unconstrained vector, chosen so the optimizer
#' cannot step outside the stationary and invertible region. A prior is
#' still written about the parameter brms names, and carried onto that
#' vector with the log Jacobian of the map, exactly as class `"sd"`
#' carries a density on a standard deviation onto its logarithm. So
#' `set_prior("normal(0, 0.5)", class = "ar")` is a density on the AR
#' coefficient itself, and `summary()` reports the parameter the prior
#' was written about.
#'
#' Bounds behave the same way where the map allows it. A first-order
#' `ar`, `ma` or `cosy` coefficient is a monotone function of one
#' internal parameter, so `lb`/`ub` map exactly onto a box. At order two
#' and above they do not: `ar[1]` is a function of every internal
#' parameter of the block at once, so no box in internal space is the
#' box asked for, and `lb`/`ub` are refused rather than approximated.
#' Little is lost, because the parameterization already guarantees
#' stationarity and invertibility, which is what such a bound is usually
#' for; where a hard box really is wanted, `class = "theta"` with
#' `coef = "thetaac_1"` bounds one internal parameter.
#'
#' `cosy` is bounded below at `-1/(d - 1)` for `d` time points, where a
#' compound-symmetric matrix stops being positive definite, and a bound
#' outside that window is refused rather than clamped. brms bounds
#' `cosy` on `[0, 1]` instead, so a negative estimate here has no brms
#' counterpart.
#'
#' @section The LKJ prior:
#' `lkj(eta)` is the density `det(C)^(eta - 1)` over a block's
#' correlation matrix `C`, normalized: `eta = 1` is uniform over
#' correlation matrices, larger `eta` concentrates toward the identity.
#' frmtmb holds a correlation as an unconstrained row-normalized
#' Cholesky parameter rather than as `C`, so the density is carried onto
#' those parameters with the exact Jacobian of that map (the derivation
#' is in the source of `R/priors.R`; `tests/testthat/test-lkj.R` checks
#' the sampled correlations against the closed-form LKJ marginals). The
#' prior a FLAT correlation parameter carries instead is
#' `(1 - rho^2)^(-3/2)`, which is improper.
#'
#' It fits `us()` and `gr(cov = )` blocks of two or more terms, which
#' hold a whole correlation matrix, and the one-parameter structures
#' `cs()`, `ar1()` and `hetar1()`, whose single bounded correlation
#' takes the LKJ marginal `(1 - rho^2)^(eta - 1)` with that structure's
#' own Jacobian. A `cs()` correlation is bounded below at `-1/(d - 1)`,
#' where a compound-symmetric matrix stops being positive definite, and
#' the density is renormalized over that window. `toep()` is refused:
#' its parameterization is not positive definite everywhere, so it has
#' no correlation matrix to put a density on.
#'
#' @section Nonlinear parameters:
#' `nlpar` addresses one parameter of an `nl = TRUE` formula, brms's
#' spelling: `set_prior("normal(5000, 1000)", nlpar = "ult")`, or
#' `prior(normal(5000, 1000), nlpar = "ult")`. Class `"b"` there covers
#' EVERY coefficient of that parameter, its intercept included, because
#' a nonlinear parameter's sub-formula is not centered and brms holds
#' its intercept in the same coefficient vector as its slopes. That is
#' why the vignette spelling above lands on `ult_(Intercept)` rather
#' than on nothing. Narrow to one column with `coef` (`"Intercept"` and
#' `"(Intercept)"` both name the intercept), or write
#' `class = "Intercept", nlpar = "ult"`, which is frmtmb's spelling of
#' the same slot; brms does not take that spelling. `nlpar` narrows
#' classes `"sd"` and `"cor"` to the random-effect blocks of that
#' parameter as well.
#'
#' A class `"b"` or `"Intercept"` specification without `nlpar` does not
#' reach a nonlinear parameter. On a model whose location is nonlinear,
#' such a specification addresses no parameter, and it is refused, as
#' in brms. A distributional parameter keeps its own spelling there:
#' `set_prior("student_t(3, 0, 2.5)", class = "sigma")` needs no
#' `nlpar`.
#'
#' A prior with a location places [frm()]'s `start` for a nonlinear
#' parameter. `normal()`, `student_t()` and `cauchy()` all carry one,
#' and where `start` does not set a nonlinear coefficient, that
#' coefficient begins at the prior's location, reported in a message.
#' Other parameters keep their usual starts: a prior is a penalty, not
#' a claim about where to begin. Without a located prior a nonlinear
#' model still needs `start`, because `frm()` evaluates the objective
#' AT the starting values; [par_template()] names them.
#'
#' `resp` picks one response of a multivariate model, and a nonlinear
#' parameter of one response is addressed with both `resp` and
#' `nlpar`.
#'
#' @section Translating a brms prior:
#' `frm(prior = )` takes a `brmsprior` object directly, whether it came
#' from `brms::set_prior()`, `brms::prior()` or `brms::get_prior()`. A
#' row applies whatever its `prior` string says, which is brms's own
#' rule; an empty string is brms's flat default and applies nothing.
#' The `source` column is not read. Earlier releases dropped rows marked
#' `source == "default"` were dropped, which lost a prior the user had
#' edited into a `get_prior()` table in place, because brms does not
#' update `source` after that edit.
#'
#' **A specification means here what the same words mean in brms.**
#' Every class brms writes that frmtmb can honor is a class
#' `set_prior()` takes under the same name and with the same meaning,
#' so a table and a hand-written specification reach one code path:
#' - `b`, `Intercept`, `sd`, `cor`, `ar`, `ma`, `cosy`, `cortime` and
#'   `rescor` are frmtmb's own class names and keep their meaning.
#' - a DISTRIBUTIONAL parameter's own class (`sigma`, `shape`, `phi`,
#'   `nu`, `kappa`, `sigma1`, ...) is a density on that parameter
#'   ITSELF, through the dpar's inverse link with that map's
#'   log-Jacobian, which is the change of variables class `"sd"`
#'   performs. A bound travels with it, so `lb = 0` on a log-linked
#'   dispersion becomes no constraint rather than a floor of 1.
#' - `theta`/`theta1`/`theta2`, `simo`, `sds`, `sdgp`, `lscale`,
#'   `sdcar` and `car` are refused by name, each saying where frmtmb
#'   keeps that quantity instead. A refusal is deliberate: translating
#'   one of them would produce a different model rather than no model.
#' - a `coef` on a `sd` or `cor` row is refused for the same reason.
#'   brms narrows such a row to one coefficient of a block; frmtmb
#'   resolves those classes per BLOCK, so applying the row without its
#'   `coef` would put the density on every standard deviation of the
#'   block, which is a wider prior than the one written. `coef` on
#'   class `"b"` narrows as it does in brms and is unaffected.
#'
#' Every refused row of a table is named in ONE message, because a table
#' is edited as a whole and stopping at the first bad row costs a round
#' trip per bad row.
#'
#' brms's `tag` and `check` have no counterpart: `tag` names a prior
#' for reuse inside a Stan program, and `check` passes an unchecked
#' string through to one. frmtmb compiles no Stan program, so both are
#' omitted rather than accepted and ignored.
#'
#' @section A distributional parameter's own class:
#' A distributional parameter has TWO spellings, and which one applies
#' is decided by the model rather than by taste. They are brms's own
#' two, and brms accepts each only on the model the other does not:
#'
#' \describe{
#'   \item{`class = "sigma"`}{a density on sigma ITSELF, on its natural
#'     scale. Available when sigma has no predictor of its own, which
#'     is the model where sigma is a single number.}
#'   \item{`class = "Intercept", dpar = "sigma"`}{a density on the
#'     LINK-scale intercept of sigma's linear predictor. Available when
#'     the model gives sigma a formula, `sigma ~ 1` included.}
#' }
#'
#' Writing a formula for a parameter is what replaces the parameter
#' with a linear predictor, so `bf(y ~ x, sigma ~ 1)` has no `sigma` to
#' put a density on, and `bf(y ~ x)` has no `Intercept_sigma`. Each
#' spelling is refused by name on the other's model, naming the one
#' that applies, because a silently retargeted prior is a different
#' model rather than no model.
#'
#' `class = "b", dpar = "sigma"` addresses that predictor's SLOPES and
#' is on the link scale, in both packages. It needs slopes to exist: on
#' `sigma ~ 1` the predictor is an intercept only, so the row addresses
#' nothing, and it is refused rather than accepted as a silent no-op.
#' brms refuses it there too.
#'
#' A distributional class names ONE parameter, so it takes `resp` and
#' neither `coef` nor `group`; both are refused rather than dropped.
#' On a multivariate model it needs `resp`, as it does in brms: without
#' it, the specification is refused.
#'
#' [get_prior()] lists whichever of the two spellings a model offers.
#' Where a parameter has NEITHER, because frmtmb refuses its class by
#' name, the table leaves it out rather than advertising a slot nothing
#' can fill: a mixture proportion with no predictor of its own is one,
#' and `?set_prior`'s refusal says where that quantity lives.
#'
#' @param prior Distribution string, e.g. `"normal(0, 5)"`, or a
#'   [prior_normal()]/[prior_t()]/[prior_lkj()] object, or `""` for
#'   bounds only.
#' @param class `"b"`, `"Intercept"`, `"sd"`, `"cor"`, `"theta"`, one
#'   of the residual-structure classes `"ar"`, `"ma"`, `"cosy"`,
#'   `"cortime"` and `"rescor"`, or a distributional parameter's own
#'   name (`"sigma"`, `"shape"`, ...). See A distributional parameter's
#'   own class.
#' @param coef Restrict to one coefficient (classes `"b"`/`"Intercept"`).
#' @param group Restrict class `"sd"` or `"cor"` to one grouping factor.
#' @param resp Response of a multivariate model. Required there for
#'   every class except `"cor"`, `"rescor"` and `"theta"`; refused on a
#'   model with one response and on class `"rescor"`.
#' @param dpar Distributional parameter (default: the location
#'   parameters).
#' @param nlpar Nonlinear parameter of an `nl = TRUE` formula. See
#'   Nonlinear parameters.
#' @param lb,ub Optional hard bounds, on the scale described in Hard
#'   bounds.
#' @return A `frmtmb_priorlist`.
#'
#' @srrstats {G2.0,G2.1} `prior` is asserted to be a length-one character
#'   vector before it is parsed, and the parsed distribution's arguments
#'   are asserted to have the arity that distribution requires (two for
#'   `normal`, three for `student_t`, one for `exponential`). A call that
#'   supplies neither a distribution nor bounds errors instead of
#'   producing an empty prior.
#' @srrstats {G2.3a} `class` is checked in two stages, because half of
#'   the vocabulary is a property of the model rather than of the
#'   package: a name that is neither one of frmtmb's own classes nor a
#'   usable parameter name errors here and lists the permitted values,
#'   a name frmtmb refuses on principle (`sds`, `simo`, `car`, ...)
#'   errors here saying where that quantity lives, and any other name
#'   is read as a distributional parameter's and checked against the
#'   model's own parameters when the prior is resolved, which is where
#'   they can be listed. The one distribution that belongs to a single
#'   class, `lkj()`, is checked against it in both directions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), z = rnorm(100),
#'                  g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # `+` combines specifications; the class-wide one goes first so the
#' # coefficient-specific one can override it
#' pr <- set_prior("normal(0, 1)", class = "b") +
#'   set_prior("normal(0, 0.2)", class = "b", coef = "z") +
#'   set_prior("exponential(1)", class = "sd", group = "g")
#' pr
#'
#' # the priors penalize the likelihood: the fit is a MAP estimate
#' fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, prior = pr)
#' fixef_by_dpar(fit)$mu
#' # the tight prior on z shrinks it toward zero
#' fixef_by_dpar(frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd))$mu
#'
#' # an empty distribution string sets a hard bound only
#' set_prior("", class = "b", coef = "x", lb = 0)
#'
#' # a distributional parameter's own class is a density on the
#' # parameter itself, on its own scale, where the model gives that
#' # parameter no predictor. This model does not, so this is the
#' # spelling it offers, and get_prior() lists it
#' fit_s <- frm(bf(y ~ x + z) + gaussian(), data = dd,
#'              prior = set_prior("student_t(3, 0, 2.5)",
#'                                class = "sigma"))
#' sigma(fit_s)
#' # with sigma ~ 1 the model has a log-scale intercept instead, and
#' # that is the slot the prior addresses. brms draws the same line
#' set_prior("student_t(3, 0, 2.5)", class = "Intercept",
#'           dpar = "sigma")
#'
#' # bounds address a nonlinear parameter the way a distribution does,
#' # so a guessing rate is held in [0, 1]
#' set_prior("", nlpar = "guess", lb = 0, ub = 1)
#'
#' # the residual-correlation classes are brms's own names
#' set_prior("normal(0, 0.5)", class = "ar")
#' set_prior("lkj(2)", class = "rescor")
#' # and one internal covariance parameter, by its name
#' set_prior("", class = "theta", coef = "thetaac_1", lb = -2, ub = 2)
#'
#' # class "cor" addresses a correlated block as a whole, brms's
#' # spelling; eta > 1 pulls the correlation toward zero
#' dd$z <- rnorm(100)
#' dd$y2 <- dd$y + rnorm(10, 0, 0.6)[dd$g] * dd$z
#' fitc <- frm(bf(y2 ~ x + z + (z | g)) + gaussian(), data = dd,
#'             prior = set_prior("lkj(4)", class = "cor"))
#' VarCorr(fitc)
#'
#' # get_prior() shows which rows a design offers
#' get_prior(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#'
#' # prior() quotes its first argument, brms's spelling, and reaches
#' # the same machinery
#' prior(normal(0, 1), class = "b")
#' @export
set_prior <- function(prior = "", class = "b", coef = "", group = "",
                      resp = "", dpar = "", nlpar = "", lb = NA,
                      ub = NA) {
  # brms reads every argument as a column and builds one row per
  # element, recycled the way data.frame() recycles, so
  # class = c("b", "sd") is two specifications. A prior object is one
  # value even though it is a list
  args <- list(prior = if (inherits(prior, "frmtmb_prior")) list(prior)
                       else prior,
               class = class, coef = coef, group = group, resp = resp,
               dpar = dpar, nlpar = nlpar, lb = lb, ub = ub)
  n <- lengths(args)
  if (any(n == 0L) || any(max(n) %% n != 0L)) {
    frm_stop("set_prior() recycles its arguments into rows, so every ",
             "argument needs a length that divides the longest; got ",
             paste0(names(n), " = ", n, collapse = ", "), call. = FALSE)
  }
  rows <- lapply(seq_len(max(n)), function(i) {
    a <- lapply(args, function(v) v[[(i - 1L) %% length(v) + 1L]])
    unclass(set_prior_one(a$prior, a$class, a$coef, a$group, a$resp,
                          a$dpar, a$nlpar, a$lb, a$ub))[[1L]]
  })
  structure(rows, class = "frmtmb_priorlist")
}

#' One row of [set_prior()]: every argument is a scalar here.
#'
#' @noRd
set_prior_one <- function(prior, class, coef, group, resp, dpar, nlpar,
                          lb, ub) {
  for (nm in c("coef", "group", "resp", "dpar", "nlpar")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || is.na(v)) {
      frm_stop("`", nm, "` must be a string, not ", arg_desc(v),
               call. = FALSE)
    }
  }
  dist <- parse_prior_dist(prior)
  # brms carries lb/ub as STRINGS (its prior frame is all character),
  # and prior() deparses everything it is given, so a bound arrives
  # here as "0" as often as 0. Normalizing once means the comparisons
  # downstream are numeric, where `"0" > 0` would have been a string
  # comparison that quietly answered FALSE
  lb <- parse_prior_bound(lb, "lb")
  ub <- parse_prior_bound(ub, "ub")
  if (is.null(dist) && is.na(lb) && is.na(ub)) {
    frm_stop("set_prior() needs a distribution, bounds, or both",
             call. = FALSE)
  }
  if (!is.character(class) || length(class) != 1L || is.na(class) ||
      !nzchar(class)) {
    frm_stop("`class` must be non-empty strings", call. = FALSE)
  }
  natural <- !class %in% frmtmb_prior_classes
  if (natural) check_dpar_prior_class(class)
  # lkj is a density over a whole correlation matrix, so it has no
  # meaning on a single coefficient or standard deviation, and the
  # matrix-valued classes have no meaning without it: neither mistake
  # can produce a silently different model
  lkj_classes <- c("cor", "cortime", "rescor")
  is_lkj <- identical(dist$kind, "lkj")
  if (is_lkj && !class %in% lkj_classes) {
    frm_stop("lkj() is a density over a whole correlation matrix; it ",
             "belongs to class = ",
             paste(paste0("\"", lkj_classes, "\""), collapse = ", "),
             " (got class = \"", class, "\")", call. = FALSE)
  }
  if (class %in% lkj_classes && !is_lkj) {
    frm_stop("class = \"", class, "\" takes an lkj() prior, e.g. ",
             "set_prior(\"lkj(2)\", class = \"", class, "\"): it addresses ",
             "a whole correlation matrix, which no per-parameter ",
             "distribution describes", call. = FALSE)
  }
  if (class %in% lkj_classes && (!is.na(lb) || !is.na(ub))) {
    # a bound belongs to ONE parameter, and these classes each name a
    # whole correlation matrix whose entries are not free of one
    # another; accepting it here would silently drop it
    frm_stop("class = \"", class, "\" takes no lb/ub: the bound would apply ",
             "to a whole correlation matrix. Bound one parameter at a time ",
             "with class = \"theta\", whose coef names an internal ",
             "parameter (\"theta_1\", \"thetaac_1\", \"thetar_1\")",
             call. = FALSE)
  }
  # a nonlinear parameter is addressed by nlpar and a distributional
  # one by dpar; frmtmb's frame gives each its own linear predictor, so
  # naming both at once says one thing twice and could say two
  # different things, which is a question about intent rather than a
  # setting to resolve
  if (nzchar(nlpar) && nzchar(dpar)) {
    frm_stop("set_prior() takes `dpar` or `nlpar`, not both: each ",
             "nonlinear parameter has its own linear predictor here, so ",
             "nlpar = \"", nlpar, "\" already names one slot", call. = FALSE)
  }
  ch <- unhonored_coef_refusal(class, coef, group)
  if (!is.null(ch)) frm_stop(ch, call. = FALSE)
  if (natural) {
    # the class IS the parameter, so `dpar` would name it twice and
    # `nlpar` would name something else entirely
    if (nzchar(dpar) && !identical(dpar, class)) {
      frm_stop("class = \"", class, "\" already names the distributional ",
               "parameter this prior is about, so dpar = \"", dpar,
               "\" names a second one. Write one or the other",
               call. = FALSE)
    }
    if (nzchar(nlpar)) {
      frm_stop("class = \"", class, "\" names a distributional parameter ",
               "and nlpar = \"", nlpar, "\" names a nonlinear one. A ",
               "nonlinear parameter's coefficients are class = \"b\" with ",
               "nlpar =, as they are in brms", call. = FALSE)
    }
    dpar <- class
    class <- "Intercept"
  }
  # the string as written is kept beside the parsed density, because
  # brms prints and returns it verbatim: "normal(0,1)" stays
  # "normal(0,1)", and "cauchy(0, 1)" does not come back as the
  # student_t it parses into
  spec <- list(dist = dist, class = class, coef = coef, dpar = dpar,
               group = group, resp = resp, nlpar = nlpar, lb = lb,
               ub = ub,
               prior = if (is.character(prior)) prior else
                 prior_dist_string(dist))
  # written only when TRUE. The link-scale spelling carries no such
  # field at all, and frmtmb.sample's test-sample-direct.R reads its
  # absence, so a default of FALSE would itself be a visible change
  if (natural) spec$natural <- TRUE
  structure(list(spec), class = "frmtmb_priorlist")
}

#' The class names that are frmtmb's own, in the order `?set_prior`
#' documents them. Everything else `set_prior()` accepts is the name of
#' a distributional parameter, which is a class in brms too.
#'
#' @noRd
frmtmb_prior_classes <- c("b", "Intercept", "sd", "cor", "theta",
                          "ar", "ma", "cosy", "cortime", "rescor")

#' Refuse a class name that is neither frmtmb's own nor usable as a
#' distributional parameter's name.
#'
#' There is no model here, so the name cannot be checked against the
#' dpars the model actually has: `resolve_priorlist()` does that, where
#' it can list them. What CAN be checked is the vocabulary brms and
#' frmtmb disagree about, and it is checked with the same messages the
#' translated route uses, so `set_prior(class = "sds")` and a brms
#' `sds` row now say the same thing.
#'
#' @noRd
check_dpar_prior_class <- function(cls) {
  hint <- brms_prior_class_refusal(cls)
  if (!is.null(hint)) {
    frm_stop("class = \"", cls, "\" has no faithful frmtmb spelling. ",
             hint, call. = FALSE)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9_.]*$", cls)) {
    frm_stop("class = \"", cls, "\" is neither one of frmtmb's classes (",
             paste(frmtmb_prior_classes, collapse = ", "),
             ") nor the name of a distributional parameter", call. = FALSE)
  }
  invisible(cls)
}

#' A hard bound as one number or `NA`, from the number, the string
#' brms's prior frame stores, or the deparsed constant [prior()]
#' produces.
#'
#' @noRd
parse_prior_bound <- function(x, arg) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  ok <- length(x) == 1L && (is.character(x) || is.numeric(x) ||
                              is.logical(x))
  if (!ok) {
    frm_stop("`", arg, "` must be a single number or NA, not ", arg_desc(x),
             ": a bound here is a hard box constraint on one parameter, ",
             "not a Stan expression", call. = FALSE)
  }
  if (is.na(x)) return(NA_real_)
  if (is.character(x)) {
    if (!nzchar(x)) return(NA_real_)
    v <- suppressWarnings(as.numeric(x))
    if (is.na(v)) {
      frm_stop("`", arg, "` = ", encodeString(x, quote = "\""),
               " is not a number: a bound here is a hard box constraint ",
               "on one parameter, and only a constant can be one",
               call. = FALSE)
    }
    return(v)
  }
  as.numeric(x)
}

#' Set up priors with brms's quoting spelling
#'
#' `prior()` is [set_prior()] with the distribution given UNQUOTED, as
#' brms's `prior()` takes it: `prior(normal(5000, 1000), nlpar = "ult")`
#' is `set_prior("normal(5000, 1000)", nlpar = "ult")`. Every argument
#' is deparsed rather than evaluated, so `class = b` and `class = "b"`
#' mean the same thing, and a variable holding a distribution is
#' deparsed to its NAME rather than its value: use [prior_string()] to
#' build a prior from strings computed at run time.
#'
#' `prior_()` takes one-sided formulas, calls, names or constants
#' (`prior_(~normal(0, 10), class = ~b)`) and `prior_string()` takes
#' plain strings; both exist so that priors can be built
#' programmatically, and both are brms's.
#'
#' A frmtmb prior and a brms prior are different objects, and with brms
#' attached after frmtmb its `prior()` masks this one. Nothing breaks:
#' [frm()] and `frmtmb.sample::frm_sample()` accept a `brmsprior` object
#' and translate its rows, so `c(prior(...), prior(...))` copied out of a
#' brms script
#' works whichever `prior()` was in scope.
#'
#' @inheritParams set_prior
#' @param ... Any of [set_prior()]'s remaining arguments: `class`,
#'   `coef`, `group`, `resp`, `dpar`, `nlpar`, `lb`, `ub`.
#' @return A `frmtmb_priorlist`.
#' @examples
#' # the brms nonlinear vignette's spelling
#' prior(normal(5000, 1000), nlpar = "ult")
#'
#' # combine with c() or `+`, as with set_prior()
#' c(prior(normal(1, 2), nlpar = "omega"),
#'   prior(normal(45, 10), nlpar = "theta"))
#'
#' # the programmatic spellings
#' prior_(~normal(0, 10), class = ~b)
#' prior_string(paste0("normal(0, ", 2 * 5, ")"), class = "b")
#' @export
prior <- function(prior, ...) {
  cl <- as.list(match.call()[-1L])
  do.call(set_prior, lapply(cl, deparse_prior_arg),
          envir = parent.frame())
}

#' @rdname prior
#' @export
prior_ <- function(prior, ...) {
  cl <- c(list(prior = prior), list(...))
  do.call(set_prior, lapply(cl, deparse_prior_value))
}

#' @rdname prior
#' @export
prior_string <- function(prior, ...) set_prior(prior, ...)

#' brms's `deparse_no_string()`: a character argument is already the
#' string `set_prior()` wants, and anything else is the user's
#' unevaluated code, which is deparsed rather than evaluated.
#'
#' @noRd
deparse_prior_arg <- function(x) {
  if (is.character(x)) x else paste(deparse(x), collapse = "")
}

#' The `prior_()` variant, which reads VALUES rather than unevaluated
#' arguments: a one-sided formula gives up its right-hand side, and a
#' call, name or constant is deparsed as it stands.
#'
#' @noRd
deparse_prior_value <- function(x) {
  if (inherits(x, "formula") && length(x) == 2L) {
    return(paste(deparse(x[[2L]]), collapse = ""))
  }
  if (is.character(x)) return(x)
  if (is.call(x) || is.name(x) || is.atomic(x)) {
    return(paste(deparse(x), collapse = ""))
  }
  frm_stop("prior_() takes one-sided formulas, calls, names or constants; ",
           "got ", arg_desc(x), ". prior_string() takes plain strings",
           call. = FALSE)
}

#' Whatever a `prior =` argument turned out to be, as
#' something the resolver understands.
#'
#' A brms `brmsprior` is a data frame of prior/class/coef/group/resp/
#' dpar/nlpar/lb/ub strings, which is exactly `set_prior()`'s vocabulary
#' written down, so it is TRANSLATED rather than refused: with brms
#' attached its `prior()` masks frmtmb's, and a ported script's
#' `c(prior(...), prior(...))` then arrives here as a brms object
#' through no fault of the caller. The legacy named list and the
#' `"flat"` string pass through untouched.
#'
#' A row applies whatever its `prior` string says, which is brms's own
#' rule, and an empty string is brms's flat default and applies nothing.
#' The `source` column is NOT read: it records who BUILT the row, not
#' who wrote the density in it, and brms does not update it when a user
#' edits the `prior` cell of a `get_prior()` table in place. Keying the
#' drop on it lost the user's own prior and reported it as one brms had
#' filled in.
#'
#' @noRd
as_priorlist <- function(x) {
  # a default_prior() or validate_prior() table is already in this
  # package's vocabulary, so it is read as written rather than
  # translated; brms passes its own table back to brm() the same way
  if (inherits(x, "frmtmb_prior_rows")) {
    return(priorlist_from_rows(as.data.frame(x), what = "the prior table"))
  }
  if (!inherits(x, "brmsprior")) return(x)
  rows <- as.data.frame(x, stringsAsFactors = FALSE)
  chr <- function(nm, i) {
    v <- if (nm %in% names(rows)) rows[[nm]][i] else ""
    if (is.na(v)) "" else as.character(v)
  }
  bnd <- function(nm, i) {
    if (!nm %in% names(rows)) return(NA)
    v <- rows[[nm]][i]
    # brms's frame is all character, and an unset bound is spelled as
    # the empty string in some rows and as NA in others
    if (is.na(v) || (is.character(v) && !nzchar(v))) NA else v
  }
  out <- list()
  # Every refusal in the table is collected and reported together. A
  # table is edited as a whole, so stopping at the first bad row costs
  # one round trip per bad row: y ~ gp(x) carries both an `lscale` row
  # and an `sdgp` row, and naming only the first made it two edits.
  bad <- character(0)
  refuse <- function(i, dist, cls, why) {
    bad[[length(bad) + 1L]] <<-
      paste0("row ", i, " (",
             if (nzchar(dist)) dist else "bounds only",
             ", class = \"", cls, "\"): ", why)
  }
  for (i in seq_len(nrow(rows))) {
    dist <- chr("prior", i)
    lb <- bnd("lb", i)
    ub <- bnd("ub", i)
    cls <- chr("class", i)
    # get_prior()/default_prior() rows with an empty `prior` are
    # "this slot exists and is flat", not a prior to apply. brms echoes
    # the parameter's DECLARED bounds onto those rows as well, so a
    # bound on a flat row is brms restating the model rather than the
    # caller asking for a box; it is carried only where frmtmb can name
    # the class at all, which keeps a flat sds/sdgp/simo row out of a
    # refusal it did not ask for
    if (!nzchar(dist)) {
      if (is.na(lb) && is.na(ub)) next
      if (!is.null(brms_prior_class_refusal(cls))) next
    }
    if (nzchar(chr("tag", i))) {
      refuse(i, dist, cls,
             paste0("tag = \"", chr("tag", i), "\" names a prior for ",
                    "reuse inside a Stan program, which frmtmb does ",
                    "not build. Drop the tag. "))
      next
    }
    hint <- brms_prior_class_refusal(cls)
    if (!is.null(hint)) {
      refuse(i, dist, cls, hint)
      next
    }
    rt <- brms_prior_route(cls, dist)
    ch <- unhonored_coef_refusal(rt$class, chr("coef", i),
                                 chr("group", i))
    if (!is.null(ch)) {
      refuse(i, dist, cls, ch)
      next
    }
    one <- tryCatch(
      set_prior(dist, class = rt$class, coef = chr("coef", i),
                group = chr("group", i), resp = chr("resp", i),
                dpar = rt$dpar %||% chr("dpar", i),
                nlpar = chr("nlpar", i),
                lb = lb, ub = ub),
      error = function(e) {
        refuse(i, dist, cls, paste0(conditionMessage(e), ". "))
        NULL
      })
    if (is.null(one)) next
    # `natural` is set by set_prior(), which read the same class this
    # row carries: the translation is now a pass-through rather than a
    # second rule that could drift from the first
    out[[length(out) + 1L]] <- unclass(one)[[1L]]
  }
  if (length(bad)) {
    frm_stop("A brms prior table has ", length(bad),
             if (length(bad) == 1L) " row" else " rows",
             " with no faithful frmtmb spelling:\n",
             paste0("  ", bad, collapse = "\n"),
             "\nWrite the prior you mean with set_prior() directly",
             call. = FALSE)
  }
  if (!length(out)) return(NULL)
  structure(out, class = "frmtmb_priorlist")
}

#' Why a `coef` or `group` on a row cannot be honored, or `NULL`.
#'
#' brms narrows a `sd` row to one coefficient of a block and writes
#' `exponential_lpdf(sd_1[2] | 1)`, keeping its default on the rest.
#' frmtmb's class `"sd"` addresses a BLOCK: `resolve_priorlist()` reads
#' `group`, `nlpar` and `resp` and never reads `coef`, so a narrowed row
#' would silently apply to every standard deviation of the block. That
#' is a wider prior than the one asked for, which is exactly the failure
#' D1 exists to remove, so it is refused instead. Class `"cor"` reads no
#' `coef` either, and a correlation is not per-coefficient at all.
#'
#' A distributional parameter's own class is one parameter, so it reads
#' neither: it has no coefficients to narrow to and belongs to no
#' random-effect block. brms refuses the same rows, naming a parameter
#' that does not exist (`sigma_x`). This is called with the class the
#' user WROTE, before a natural class is rewritten to its storage pair,
#' which is the only point where that word is still visible.
#'
#' @noRd
unhonored_coef_refusal <- function(cls, coef, group = "") {
  nat <- !cls %in% frmtmb_prior_classes
  if (nat && (nzchar(coef) || nzchar(group))) {
    arg <- if (nzchar(coef)) {
      paste0("coef = \"", coef, "\"")
    } else {
      paste0("group = \"", group, "\"")
    }
    return(paste0("class = \"", cls, "\" is a density on ", cls,
                  " itself, which is one parameter, so ", arg,
                  " names nothing it can narrow to and would be ",
                  "dropped. Drop it; `resp` is the only narrowing a ",
                  "distributional class takes, and `dpar` predictors ",
                  "are addressed with class = \"b\" or ",
                  "class = \"Intercept\". "))
  }
  if (!nzchar(coef) || !cls %in% c("sd", "cor")) return(NULL)
  if (identical(cls, "cor")) {
    return(paste0("class = \"cor\" addresses a whole correlation ",
                  "matrix, so coef = \"", coef, "\" names nothing it ",
                  "can narrow to. Drop the coef; `group` selects the ",
                  "block. "))
  }
  paste0("frmtmb's class = \"sd\" addresses a whole random-effect ",
         "BLOCK, not one coefficient of it, so coef = \"", coef,
         "\" cannot be honored and applying the row without it would ",
         "put the density on every standard deviation of the block. ",
         "Drop the coef to prior the whole block (`group` selects ",
         "which one), or address the one parameter with ",
         "class = \"theta\" and the coef get_prior() lists for it. ")
}

#' brms class names frmtmb spells with the same word and the same
#' meaning, so a translated row keeps its class untouched.
#'
#' @noRd
brms_direct_prior_classes <- c("b", "Intercept", "sd", "cor", "ar",
                               "ma", "cosy", "cortime", "rescor")

#' Why a brms class cannot be carried over, or `NULL` when it can.
#'
#' Each of these names a structure frmtmb holds differently, or does not
#' hold at all. Translating one would produce a DIFFERENT model rather
#' than no model, which is why the refusal is by name and says where the
#' quantity actually lives.
#'
#' @noRd
brms_prior_class_refusal <- function(cls) {
  # brms spells a mixture proportion theta1, theta2, ...; frmtmb's
  # "theta" is the raw internal covariance vector. The whole prefix is
  # refused, because the bare word is the only spelling that used to be
  # caught and "theta2" then fell through to advice that does not work
  if (startsWith(cls, "theta")) {
    return(paste0("brms's \"", cls, "\" is a mixture proportion, and ",
                  "frmtmb's \"theta\" is the raw internal covariance ",
                  "vector: the word names two unrelated sets of ",
                  "parameters, so the row cannot be carried over. A ",
                  "mixture proportion that has its OWN predictor is ",
                  "addressed as class = \"Intercept\" with dpar = that ",
                  "proportion's name, on the link scale where brms puts ",
                  "it too, and get_prior() lists the ones this model ",
                  "has. The proportion brms holds as the reference ",
                  "component is a parameter on neither side, and its ",
                  "row applies nothing in brms either. "))
  }
  switch(cls,
    simo = paste0(
      "brms's \"simo\" is the Dirichlet on a mo() simplex. frmtmb holds ",
      "that simplex as its free softmax coordinates and puts no density ",
      "on it at all, so there is no slot to carry the row into and no ",
      "class that would rename it: drop the row. A flat simo row from ",
      "get_prior() is dropped for you, and only an explicit one reaches ",
      "here. "),
    sds = paste0(
      "brms's \"sds\" is the wiggliness standard deviation of a smooth. ",
      "frmtmb holds a smooth as a random-effect block, so its frmtmb ",
      "spelling is class = \"sd\" with group = the smooth's label, e.g. ",
      "group = \"s(x)\"; get_prior() lists the label this model has. "),
    sdgp = paste0(
      "brms's \"sdgp\" is the marginal standard deviation of a gp(). ",
      "frmtmb holds a gp() as a random-effect block, so its frmtmb ",
      "spelling is class = \"sd\" with group = the term's label, e.g. ",
      "group = \"gp(x)\"; get_prior() lists the label this model has. "),
    lscale = paste0(
      "brms's \"lscale\" is a gp() length-scale, which frmtmb keeps in ",
      "the raw internal covariance vector rather than in a class of its ",
      "own. Address it as class = \"theta\" with the coef get_prior() ",
      "lists for the block (\"theta_1\", \"theta_2\", ...), remembering ",
      "that a theta prior is on the INTERNAL scale. "),
    sdcar = paste0(
      "brms's \"sdcar\" is the standard deviation of a car() term. ",
      "frmtmb holds a car() as a random-effect block, so its frmtmb ",
      "spelling is class = \"sd\" with group = the term's grouping ",
      "variable; get_prior() lists the label this model has. "),
    car = paste0(
      "brms's \"car\" is the spatial dependence parameter of a car() ",
      "term, which frmtmb keeps in the raw internal covariance vector ",
      "rather than in a class of its own. Address it as ",
      "class = \"theta\" with the coef get_prior() lists for the block ",
      "(\"theta_1\", \"theta_2\", ...), remembering that a theta prior ",
      "is on the INTERNAL scale. "),
    NULL)
}

#' Where one brms prior row lands: the frmtmb class and the dpar it
#' needs.
#'
#' Every class brms writes that frmtmb can honor is now a class
#' `set_prior()` takes under the same name, so a translated row keeps
#' its class and `set_prior()` decides the placement. That is the whole
#' point of the flip: there is one rule, written once, and the two
#' routes cannot drift apart. Everything `brms_prior_class_refusal()`
#' names is refused instead.
#'
#' @noRd
brms_prior_route <- function(cls, dist) {
  if (cls %in% brms_direct_prior_classes) {
    return(list(class = cls, dpar = NULL, natural = FALSE))
  }
  hint <- brms_prior_class_refusal(cls)
  if (is.null(hint)) {
    return(list(class = cls, dpar = NULL, natural = TRUE))
  }
  frm_stop("A brms prior with class = \"", cls, "\" (", dist, ") has no ",
           "faithful frmtmb spelling. ", hint,
           "Write the prior you mean with set_prior() directly",
           call. = FALSE)
}

#' The gate on its own, for callers that only want the verdict.
#'
#' @noRd
check_brms_prior_class <- function(cls, dist) {
  invisible(brms_prior_route(cls, dist)$class)
}

#' Turn a brms-style prior string such as `"normal(0, 5)"` into a
#' `frmtmb_prior` object, or `NULL` for the empty string. It errors on an
#' unparsable string, on non-numeric arguments, and on the wrong number
#' of arguments for the named distribution.
#'
#' @noRd
parse_prior_dist <- function(prior) {
  if (inherits(prior, "frmtmb_prior")) return(prior)
  if (!is.character(prior) || length(prior) != 1L || is.na(prior)) {
    frm_stop("A prior must be one string such as \"normal(0, 5)\" or a ",
             "prior object, not ", arg_desc(prior), call. = FALSE)
  }
  if (prior == "") return(NULL)
  # the NAME is matched case-insensitively and may carry digits, so that
  # brms's shrinkage priors reach the unsupported-density message below
  # rather than the generic parse failure. R2D2() is the one that needs
  # both: an upper-case name and an empty argument list
  m <- regmatches(
    prior,
    regexec("^\\s*([A-Za-z_][A-Za-z_0-9]*)\\s*\\(([^)]*)\\)\\s*$",
            prior))[[1]]
  if (length(m) != 3) {
    frm_stop("Cannot parse prior '", prior,
             "'; expected e.g. \"normal(0, 5)\"", call. = FALSE)
  }
  kind <- m[2]
  pars <- as.numeric(strsplit(m[3], ",", fixed = TRUE)[[1]])
  if (anyNA(pars)) {
    frm_stop("Non-numeric arguments in prior '", prior, "'", call. = FALSE)
  }
  # one arity per density, checked before the switch so that every
  # density refuses the same way and names its own parameters
  takes <- prior_dist_params[[kind]]
  if (!is.null(takes) && length(pars) != length(takes)) {
    frm_stop("Prior '", prior, "' gives ", length(pars), " argument",
             if (length(pars) != 1L) "s", "; ", kind, "() takes ",
             length(takes), ": ", kind, "(", paste(takes, collapse = ", "),
             ")", call. = FALSE)
  }
  switch(kind,
    normal = prior_normal(pars[1], pars[2]),
    student_t = prior_t(pars[1], pars[2], pars[3]),
    cauchy = prior_t(1, pars[1], pars[2]),
    exponential = {
      if (!(pars[1] > 0)) {
        frm_stop("Prior '", prior, "': the rate of exponential() must ",
                 "be positive, not ", pars[1], call. = FALSE)
      }
      structure(list(kind = "exponential", rate = pars[1]),
                class = "frmtmb_prior")
    },
    lkj = prior_lkj(pars[1]),
    logistic = prior_logistic(pars[1], pars[2]),
    gamma = prior_gamma(pars[1], pars[2]),
    inv_gamma = prior_inv_gamma(pars[1], pars[2]),
    beta = prior_beta(pars[1], pars[2]),
    frm_stop("Unsupported prior distribution '", kind,
             "' (supported: normal, student_t, cauchy, exponential, ",
             "logistic, gamma, inv_gamma, beta, lkj)", call. = FALSE)
  )
}

#' The parameters each prior density string takes, in order, for the
#' arity refusal of parse_prior_dist(). The names are brms's.
#'
#' @noRd
prior_dist_params <- list(
  normal = c("mu", "sigma"), student_t = c("nu", "mu", "sigma"),
  cauchy = c("mu", "sigma"), exponential = "beta", lkj = "eta",
  logistic = c("mu", "sigma"), gamma = c("alpha", "beta"),
  inv_gamma = c("alpha", "beta"), beta = c("alpha", "beta"))

#' The location of a prior distribution, or `NA` where it has none.
#' `normal` and `student_t` (which `cauchy` parses into) carry one;
#' `exponential` and `lkj` do not. Read by the nonlinear start
#' placement in R/par-template.R.
#'
#' @noRd
prior_dist_location <- function(dist) {
  loc <- switch(dist$kind %||% "", normal = , t = , logistic = dist$location,
                NULL)
  if (is.null(loc) || length(loc) != 1L || !is.finite(loc)) {
    return(NA_real_)
  }
  as.numeric(loc)
}

#' @export
"+.frmtmb_priorlist" <- function(e1, e2) {
  for (e in list(e1, e2)) {
    if (!inherits(e, "frmtmb_priorlist")) {
      frm_stop("`+` combines prior specifications made by set_prior(), ",
               "as in set_prior(\"normal(0, 1)\") + set_prior(\"normal(0, ",
               "2)\", class = \"sd\"); one side is ", arg_desc(e),
               call. = FALSE)
    }
  }
  structure(c(unclass(e1), unclass(e2)), class = "frmtmb_priorlist")
}

#' @export
c.frmtmb_priorlist <- function(...) {
  structure(do.call(c, lapply(list(...), unclass)),
            class = "frmtmb_priorlist")
}

#' Print a prior specification
#'
#' brms's layout: one specification prints as the parameter it
#' addresses followed by its density, `b_x ~ normal(0, 1)`, with any
#' bounds in front as `<lower=0>`; several print as a table with one row
#' each, the table `as.data.frame()` returns, where a row with no density
#' of its own shows its class row's density as `"(vectorized)"`. An
#' empty prior prints nothing, as in brms.
#'
#' @param x A `frmtmb_priorlist`.
#' @param show_df Print as a table (`TRUE`) or as one line per
#'   specification (`FALSE`). `NULL`, the default, prints a table when
#'   there is more than one specification, as brms does.
#' @param ... Must be empty.
#' @return `x`, invisibly.
#' @examples
#' set_prior("normal(0,1)", coef = "x")
#' set_prior("cauchy(0,1)", class = "sd", group = "g")
#' set_prior("normal(0, 2)", class = c("b", "sd"))
#' @export
print.frmtmb_priorlist <- function(x, show_df = NULL, ...) {
  frm_check_dots(...)
  n <- length(unclass(x))
  if (is.null(show_df)) show_df <- n != 1L
  check_flag(show_df, "show_df")
  df <- as.data.frame(x)
  if (!n) {
    # brms prints nothing for an empty prior, not an empty table
  } else if (show_df) {
    shown <- df
    # a row with no density of its own shows the density of the class
    # row it sits under, marked "(vectorized)", as brms prints it
    for (i in which(!nzchar(shown$prior))) {
      up <- setdiff(which(df$class == df$class[i] &
                            df$resp == df$resp[i] & df$dpar == df$dpar[i] &
                            df$nlpar == df$nlpar[i] & !nzchar(df$coef) &
                            !nzchar(df$group)), i)
      if (length(up)) shown$source[i] <- "(vectorized)"
      shown$prior[i] <- if (length(up) && nzchar(df$prior[up[1L]])) {
        df$prior[up[1L]]
      } else {
        "(flat)"
      }
    }
    print.data.frame(shown, row.names = FALSE)
  } else {
    cat(prior_row_lines(df), sep = "\n")
  }
  ov <- attr(x, "overrides")
  if (length(ov)) {
    cat("plus internal-scale overrides on: ",
        paste(names(ov), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' One line per row in brms's `.print_prior()` layout: bounds, then
#' the parameter name built from class, group, resp, dpar, nlpar and
#' coef, then the density.
#'
#' @noRd
prior_row_lines <- function(df) {
  usc <- function(v) ifelse(nzchar(v), paste0("_", v), "")
  group <- usc(df$group)
  deeper <- nzchar(df$resp) | nzchar(df$dpar) | nzchar(df$nlpar) |
    nzchar(df$coef)
  group <- ifelse(deeper & nzchar(group), paste0(group, "_"), group)
  lb <- ifelse(is.na(df$lb), "", paste0("lower=", df$lb))
  ub <- ifelse(is.na(df$ub), "", paste0("upper=", df$ub))
  bound <- ifelse(nzchar(lb) & nzchar(ub), paste0("<", lb, ",", ub, "> "),
                  ifelse(nzchar(lb), paste0("<", lb, "> "),
                         ifelse(nzchar(ub), paste0("<", ub, "> "), "")))
  prior <- ifelse(nzchar(df$prior), df$prior, "(flat)")
  paste0(bound, df$class, group, usc(df$resp), usc(df$dpar),
         usc(df$nlpar), usc(df$coef), " ~ ", prior)
}

#' The string a parsed density is written as, for a specification built
#' from a prior object rather than a string. It parses back to the same
#' object.
#'
#' @noRd
prior_dist_string <- function(dist) {
  if (is.null(dist)) return("")
  kind <- if (identical(dist$kind, "t")) "student_t" else dist$kind
  paste0(kind, "(", paste(unlist(dist[-1L]), collapse = ", "), ")")
}

#' @export
as.data.frame.frmtmb_priorlist <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  frm_check_dots(...)
  specs <- unclass(x)
  chr <- function(f) {
    vapply(specs, function(s) {
      v <- s[[f]]
      if (is.null(v) || is.na(v)) "" else as.character(v)
    }, "")
  }
  bnd <- function(f) {
    vapply(specs, function(s) {
      v <- s[[f]]
      if (is.null(v) || is.na(v)) NA_character_ else as.character(v)
    }, "")
  }
  # the spelling a specification is WRITTEN with, so a density on a
  # distributional parameter itself reads as that parameter's class
  sp <- lapply(specs, spec_spelling)
  data.frame(
    prior = vapply(specs, function(s) {
      s[["prior"]] %||% prior_dist_string(s[["dist"]])
    }, ""),
    class = vapply(sp, `[[`, "", "class"),
    coef = chr("coef"), group = chr("group"), resp = chr("resp"),
    dpar = vapply(sp, `[[`, "", "dpar"),
    nlpar = chr("nlpar"), lb = bnd("lb"), ub = bnd("ub"),
    source = rep("user", length(specs)),
    row.names = row.names, stringsAsFactors = FALSE
  )
}

#' Column access on a prior specification
#'
#' A `frmtmb_priorlist` holds parsed specifications, and brms's prior
#' object is a data frame of strings. `$` reads the columns brms's
#' object has, from the table `as.data.frame()` builds, so
#' `set_prior("normal(0, 2)", class = c("b", "sd"))$class` is
#' `c("b", "sd")` in both packages. The name must match a column
#' exactly; anything else is `NULL`, as it is on a data frame.
#' Assigning a column rebuilds every specification from the edited
#' table with [set_prior()], so an edit is checked as a new call would
#' be.
#'
#' @param x A `frmtmb_priorlist`.
#' @param name A column: `prior`, `class`, `coef`, `group`, `resp`,
#'   `dpar`, `nlpar`, `lb`, `ub` or `source`.
#' @param value The new column, recycled as a data frame column is.
#' @return For `$`, a character vector with one element per
#'   specification. For `$<-`, the rebuilt `frmtmb_priorlist`.
#' @examples
#' pr <- set_prior("normal(0, 2)", class = c("b", "sd"))
#' pr$class
#' pr$prior[2] <- "exponential(1)"
#' pr
#' @name frmtmb_priorlist-columns
#' @export
`$.frmtmb_priorlist` <- function(x, name) {
  as.data.frame(x)[[name, exact = TRUE]]
}

#' @rdname frmtmb_priorlist-columns
#' @export
`$<-.frmtmb_priorlist` <- function(x, name, value) {
  df <- as.data.frame(x)
  cols <- setdiff(names(df), "source")
  if (!name %in% cols) {
    frm_stop("A prior specification has no column '", name, "' to assign. ",
             "The columns are ", paste(cols, collapse = ", "), call. = FALSE)
  }
  df[[name]] <- value
  out <- priorlist_from_rows(df, what = "the edited prior specification") %||%
    empty_prior()
  attr(out, "overrides") <- attr(x, "overrides")
  out
}

#' An empty prior specification
#'
#' brms's `empty_prior()`: a prior object with no specifications, to add
#' specifications to with `+` or `c()`. [frm()] given it applies no
#' prior.
#'
#' @return A `frmtmb_priorlist` of length zero.
#' @examples
#' pr <- empty_prior()
#' pr <- pr + set_prior("normal(0, 1)", class = "b")
#' pr
#' @export
empty_prior <- function() {
  structure(list(), class = "frmtmb_priorlist")
}

#' Transform an object into a prior specification
#'
#' brms's `as.brmsprior()`. A data frame (or anything
#' `as.data.frame()` accepts) with a `prior` column becomes the object
#' [set_prior()] returns, one specification per row. Missing columns
#' take `set_prior()`'s defaults (`class = "b"`, empty `coef`, `group`,
#' `resp`, `dpar` and `nlpar`, no bounds) and columns `set_prior()` does
#' not take are dropped.
#'
#' The name is brms's, and so is the job: turn a table into the prior
#' object this package fits with. That object is a `frmtmb_priorlist`
#' rather than a `brmsprior`, because this package parses a density
#' when it is written rather than when the model is compiled. A
#' `brmsprior` built by brms is translated row by row, as [frm()]
#' translates one, and a table from [default_prior()] or
#' [validate_prior()] is read as written.
#'
#' A row whose `prior` is empty or `"(flat)"` and that carries no bound
#' applies nothing, and is left out: it is the flat default a
#' [default_prior()] table lists for a slot nobody has set. So is a row
#' whose `source` is `"(vectorized)"`, which repeats the density of the
#' class row above it.
#'
#' @param x A data frame with a `prior` column, a `brmsprior`, a
#'   [default_prior()] table or a `frmtmb_priorlist`.
#' @return A `frmtmb_priorlist`.
#' @examples
#' as.brmsprior(data.frame(prior = "normal(0,1)", coef = c("a", "b")))
#' @export
as.brmsprior <- function(x) {
  if (inherits(x, "frmtmb_priorlist")) return(x)
  if (inherits(x, "brmsprior")) {
    return(as_priorlist(x) %||% empty_prior())
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"prior" %in% names(x)) {
    frm_stop("as.brmsprior() needs a `prior` column", call. = FALSE)
  }
  priorlist_from_rows(x, what = "the table") %||% empty_prior()
}

#' Build a `frmtmb_priorlist` from a table in `set_prior()`'s own
#' vocabulary, one call per row. Every refused row is named in one
#' message, as `as_priorlist()` does for a brms table, because a table is
#' edited as a whole.
#'
#' @noRd
priorlist_from_rows <- function(x, what) {
  defaults <- list(class = "b", coef = "", group = "", resp = "",
                   dpar = "", nlpar = "", lb = NA, ub = NA)
  cell <- function(nm, i) {
    if (!nm %in% names(x)) return(defaults[[nm]])
    v <- x[[nm]][[i]]
    if (nm %in% c("lb", "ub")) return(v)
    if (is.na(v)) "" else as.character(v)
  }
  unset <- function(v) {
    is.null(v) || is.na(v) || (is.character(v) && !nzchar(v))
  }
  if ("tag" %in% names(x) && any(nzchar(x[["tag"]][!is.na(x[["tag"]])]))) {
    frm_stop("A tag names a prior for reuse inside a Stan program, which ",
             "frmtmb does not build. Drop the `tag` column", call. = FALSE)
  }
  out <- list()
  bad <- character(0)
  for (i in seq_len(nrow(x))) {
    if ("source" %in% names(x) &&
          identical(as.character(x[["source"]][[i]]), "(vectorized)")) {
      next
    }
    pr <- cell("prior", i)
    if (identical(pr, "(flat)")) pr <- ""
    lb <- cell("lb", i)
    ub <- cell("ub", i)
    if (!nzchar(pr) && unset(lb) && unset(ub)) next
    one <- tryCatch(
      set_prior(pr, class = cell("class", i), coef = cell("coef", i),
                group = cell("group", i), resp = cell("resp", i),
                dpar = cell("dpar", i), nlpar = cell("nlpar", i),
                lb = if (unset(lb)) NA else lb,
                ub = if (unset(ub)) NA else ub),
      error = function(e) {
        bad[[length(bad) + 1L]] <<- paste0("row ", i, ": ",
                                            conditionMessage(e))
        NULL
      })
    if (!is.null(one)) out <- c(out, unclass(one))
  }
  if (length(bad)) {
    frm_stop("Cannot read ", length(bad),
             if (length(bad) == 1L) " row" else " rows", " of ", what, ":\n",
             paste0("  ", bad, collapse = "\n"), call. = FALSE)
  }
  if (!length(out)) return(NULL)
  structure(out, class = "frmtmb_priorlist")
}

#' Default priors: the slots a prior can target
#'
#' brms's `default_prior()`, and `get_prior()`, which brms has kept as
#' its alias since 2.20.14: one row per slot a prior can target, with
#' the class/coef/dpar/group values to pass to [set_prior()]. Classes
#' `"sd"` and `"cor"` are targeted by `group` and `nlpar`; the
#' residual-correlation classes (`"ar"`, `"ma"`, `"cosy"`, `"cortime"`)
#' by `resp`; class `"rescor"` has one row and takes no `resp`; and
#' class `"theta"` rows name the raw internal covariance parameters
#' (escape hatch, including correlations one at a time, across all
#' three covariance components). Where the location is several
#' distributional parameters (a categorical or mixture model), each
#' one's `b`, `Intercept` and `sd` rows carry its `dpar`, as in brms.
#'
#' A nonlinear parameter's coefficients are listed under class `"b"`
#' with its name in the `nlpar` column, the intercept among them, which
#' is how brms lists them and what [set_prior()] addresses (see its
#' Nonlinear parameters section).
#'
#' A distributional parameter is listed under whichever of its two
#' spellings the model offers: its OWN class where the model gives it
#' no predictor, and class `"Intercept"` with its name in the `dpar`
#' column where the model gives it a formula. brms lists the same two,
#' the same way round. See A distributional parameter's own class in
#' [set_prior()].
#'
#' An ordinal family has no intercept column, so its class `"Intercept"`
#' row names the THRESHOLD vector, which is what the same row means in
#' brms. See the Ordinal thresholds section of [set_prior()].
#'
#' @section Where the table differs from brms's:
#' The class `"theta"` rows are frmtmb's own. They name the internal
#' covariance parameters, which are real parameters of the fit, and for
#' a `gp()` length-scale, a `car()` dependence parameter or one entry of
#' a structured covariance they are the only spelling that reaches one
#' parameter. brms has no such class, so a brms table for the same model
#' has fewer rows.
#'
#' brms lists a class `"sd"` row per coefficient of each block
#' (`sd_patient__Intercept`). frmtmb's class `"sd"` addresses a whole
#' block and refuses a `coef`, so those rows are not listed; class
#' `"theta"` reaches one standard deviation on its own.
#'
#' A flat slot reads `"(flat)"` in the `prior` column, where brms stores
#' an empty string and prints `(flat)`.
#'
#' @section Which route the defaults describe:
#' Every column but `prior` is a property of the design, and the design
#' does not change with what is attached. The `prior` column is a
#' property of a ROUTE, and frmtmb has two of them with different
#' defaults: [frm()] is maximum likelihood and is flat in every slot
#' until a prior is set, while `frmtmb.sample::frm_sample()` applies
#' brms's weakly-informative defaults on both of its routes. `route`
#' makes the caller say which one is being asked about, so that the
#' answer is a property of the question rather than of the search path.
#'
#' `route = "fit"`, the default, reads no registry at all. Its table is
#' identical whatever extension packages are loaded.
#'
#' `route = "sample"` reads the defaults `frm_sample()` would apply,
#' which only frmtmb.sample can state. Without that package loaded the
#' call is refused rather than answered `(flat)`, because a flat table
#' would be a wrong answer about the sampling route and not a missing
#' one.
#'
#' brms's `default_prior()` describes what `brm()` would use, so the
#' brms reading of this function is `route = "sample"`. `route = "fit"`
#' has no brms counterpart: it describes `frm()`.
#'
#' With frmtmb.sample attached, the example below with
#' `route = "sample"` added returns the same rows as the example below with brms's densities in
#' the `prior` column instead of `(flat)`: a Student-t on the intercept
#' centered on the response, a half-Student-t on `sigma` and on each
#' standard deviation, and `lkj(1)` on each correlation. Population-level
#' slopes stay `(flat)`, as they are in brms.
#'
#' @param object A `bf()` formula (with family), a plain formula, or
#'   an already fitted `frmtmb_fit`.
#' @param formula The same as `object`, under the name brms's
#'   `get_prior()` gives it.
#' @param data A data frame of model data (ignored when `object` is a
#'   fit).
#' @param family Family, when `object` does not carry one.
#' @param data2 Structural objects, as in [frm()] (ignored when
#'   `object` is a fit, which carries its own).
#' @param route Which route's defaults the `prior` column reports.
#'   `"fit"` (the default) reports the defaults [frm()] applies, which
#'   are flat in every slot; it consults no registry, so its answer
#'   does not depend on which packages are attached. `"sample"` reports
#'   the defaults `frmtmb.sample::frm_sample()` applies, and refuses
#'   when no package has registered any. The returned object records
#'   the route and `print()` names it on its first line.
#' @param ... For `get_prior()`, the arguments of `default_prior()`.
#' @return A data frame of class `frmtmb_prior_rows` with brms's columns
#'   `prior`, `class`, `coef`, `group`, `resp`, `dpar`, `nlpar`, `lb`,
#'   `ub` and `source`, and a `route` attribute. `source` is
#'   `"default"` on every row; [validate_prior()] marks the rows a user
#'   set.
#' @seealso [validate_prior()] for the table with a prior filled in.
#' @examples
#' dd <- data.frame(y = rnorm(60), x = rnorm(60),
#'                  g = factor(rep(1:6, 10)))
#' # what frm() applies: flat, whatever else is loaded. For what
#' # frm_sample() applies, see "Which route the defaults describe"
#' default_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' # the same table under brms's older name
#' get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' @export
default_prior <- function(object, data = NULL, family = NULL,
                          data2 = list(), route = c("fit", "sample")) {
  route <- frm_match_arg(route)
  # refused before the frame is assembled: a route nothing can answer is
  # unanswerable for every model, so the work would be thrown away
  if (identical(route, "sample")) require_prior_defaults()
  # check_trials = FALSE because brms's default_prior() answers for a
  # binomial model written without trials(), which frm() refuses
  design <- prior_design(object, data, family, data2, check_trials = FALSE)
  prior_table(design$spec, design$frame, route)
}

#' @rdname default_prior
#' @export
get_prior <- function(formula, ...) {
  default_prior(formula, ...)
}

#' The spec and frame a prior table and a prior check are read from.
#'
#' @noRd
prior_design <- function(object, data, family, data2,
                         check_trials = TRUE) {
  if (inherits(object, "frmtmb_fit")) {
    return(list(spec = object$spec, frame = object$frame))
  }
  bform <- resolve_deferred_families(as_bform(object, family), data)
  spec <- parse_spec(bform)
  frame <- assemble_frame(spec, data, data2 = data2,
                         check_trials = check_trials)
  # the table is built from the spec's primary_dpars and nlpars, so
  # it needs the finalized families: a family that derived its dpars
  # from the response would otherwise be tabled under the vocabulary
  # it was written with, and the fit route would disagree with the
  # formula route. No in-repo family does that today.
  spec <- carry_finalized_responses(spec, frame)
  list(spec = spec, frame = frame)
}

#' Check a prior against a model
#'
#' brms's `validate_prior()`: resolve a prior specification against the
#' model it is meant for, refuse it if any part of it addresses nothing
#' the model has, and return the whole [default_prior()] table with the
#' prior filled in. Fitting runs the same check, so a prior this
#' function accepts is one [frm()] accepts, and a prior it refuses is
#' refused here with the message the fit would give, before any fitting
#' work.
#'
#' A row the prior sets reads `source = "user"`. A row that takes its
#' density from a class row above it, the coefficients under a class
#' `"b"` prior or the blocks under a class `"sd"` prior, reads
#' `source = "(vectorized)"` and repeats that density, as brms prints
#' it. Every other row keeps its default.
#'
#' Two specifications for the same slot are refused, as [frm()] and
#' brms refuse them. Where specifications for different slots reach one
#' parameter, the table shows the more specific one, which is the one
#' [frm()] applies; see [set_prior()].
#'
#' The table can be passed back as `prior =`: [frm()] reads a
#' `frmtmb_prior_rows` table the way [as.brmsprior()] does, applying the
#' rows with a density or a bound and skipping the flat and
#' vectorized ones.
#'
#' @inheritParams default_prior
#' @param prior A `frmtmb_priorlist` from [set_prior()], a `brmsprior`,
#'   or `NULL` for none.
#' @param formula A `bf()` formula (with family) or a plain formula.
#' @return A `frmtmb_prior_rows` table, as [default_prior()] returns.
#' @examples
#' dd <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60),
#'                  g = factor(rep(1:6, 10)))
#' validate_prior(prior(normal(0, 10), class = b) +
#'                  prior(cauchy(0, 2), class = sd),
#'                y ~ x + z + (1 | g), data = dd)
#' # a prior on a coefficient the model does not have is refused
#' try(validate_prior(prior(normal(0, 1), coef = w), y ~ x, data = dd))
#' @export
validate_prior <- function(prior, formula, data, family = NULL,
                           data2 = list(), route = c("fit", "sample")) {
  route <- frm_match_arg(route)
  if (identical(route, "sample")) require_prior_defaults()
  pl <- as_priorlist(prior) %||% empty_prior()
  if (!inherits(pl, "frmtmb_priorlist")) {
    frm_stop("validate_prior() takes a prior built by set_prior() or ",
             "prior(), or a brms prior table; got ", arg_desc(prior),
             call. = FALSE)
  }
  check_prior_slots(pl)
  design <- prior_design(formula, data, family, data2)
  # the fit's own resolution, so every refusal is the one frm() gives
  resolve_priorlist(design, pl)
  tab <- prior_table(design$spec, design$frame, route)
  fill_prior_table(tab, pl)
}

#' Write a prior specification into a [default_prior()] table: the
#' rows it addresses by name first, in the order the resolver applies
#' them so the more specific one wins, and then the rows below a class
#' row it set. Reading the table must never show a density the fit does
#' not apply.
#'
#' @noRd
fill_prior_table <- function(tab, pl) {
  cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
  bare <- function(v) par_name_bare(v)
  # which rows a specification's resp, dpar and nlpar reach, read the
  # way the resolver reads them. No resp means every response; in a
  # multivariate model only a class with no response key gets that far,
  # because resp_missing_refusal() stops the others. A
  # class cor specification with no dpar or nlpar reaches the blocks of
  # every predictor. A class sd one reaches its own prefix only, as in
  # brms (sd_spec_reach())
  covers <- function(resp, dpar, nlpar, class, rows) {
    wide <- identical(class, "cor")
    sd <- identical(class, "sd")
    (rows$resp == resp | (!sd & !nzchar(resp))) &
      (rows$dpar == dpar | (wide & !nzchar(dpar))) &
      (rows$nlpar == nlpar | (wide & !nzchar(nlpar)))
  }
  user <- as.data.frame(structure(prior_specificity_order(pl),
                                  class = "frmtmb_priorlist"))
  for (i in seq_len(nrow(user))) {
    u <- user[i, ]
    hit <- tab$class == u$class & tab$group == u$group &
      covers(u$resp, u$dpar, u$nlpar, u$class, tab) &
      (tab$coef == u$coef | bare(tab$coef) == bare(u$coef))
    if (!any(hit)) {
      # the resolver accepted it, so it addresses something: a slot the
      # table lists under a wider row, such as a response of a
      # multivariate model. It is appended rather than lost
      extra <- tab[rep(1L, 1L), , drop = FALSE]
      extra[1L, cols] <- u[cols]
      extra$prior <- "(flat)"
      extra$lb <- NA_real_
      extra$ub <- NA_real_
      tab <- rbind(tab, extra)
      hit <- c(logical(nrow(tab) - 1L), TRUE)
    }
    if (nzchar(u$prior)) {
      tab$prior[hit] <- u$prior
      tab$source[hit] <- "user"
    }
    if (!is.na(u$lb)) {
      tab$lb[hit] <- as.numeric(u$lb)
      tab$source[hit] <- "user"
    }
    if (!is.na(u$ub)) {
      tab$ub[hit] <- as.numeric(u$ub)
      tab$source[hit] <- "user"
    }
  }
  # a row the prior did not name takes its density from the user row
  # the resolver applies over it: of the class rows that reach it, the
  # most specific, as prior_specificity_order() ranks them. Parents are
  # looked for anywhere in the table, because a dpar- or resp-specific
  # class row the prior names may have been appended at its end
  nz <- function(v) as.integer(nzchar(v))
  rank <- 10L * nz(tab$group) + nz(tab$resp) + nz(tab$dpar) + nz(tab$nlpar)
  is_user <- tab$source == "user" & tab$prior != "(flat)"
  for (i in seq_len(nrow(tab))) {
    if (!identical(tab$source[i], "default") ||
        !identical(tab$prior[i], "(flat)")) {
      next
    }
    if (!nzchar(tab$coef[i]) && !nzchar(tab$group[i]) &&
        !nzchar(tab$dpar[i]) && !nzchar(tab$resp[i]) &&
        !nzchar(tab$nlpar[i])) {
      next
    }
    parent <- vapply(seq_len(nrow(tab)), function(j) {
      covers(tab$resp[j], tab$dpar[j], tab$nlpar[j], tab$class[j],
             tab[i, ])
    }, TRUE)
    up <- which(tab$class == tab$class[i] & parent & is_user &
                  !nzchar(tab$coef) &
                  (if (nzchar(tab$coef[i])) {
                    tab$group == tab$group[i] | !nzchar(tab$group)
                  } else if (nzchar(tab$group[i])) {
                    !nzchar(tab$group)
                  } else {
                    !nzchar(tab$group)
                  }) &
                  seq_len(nrow(tab)) != i)
    if (!length(up)) next
    best <- up[rank[up] == max(rank[up])]
    tab$prior[i] <- tab$prior[best[length(best)]]
    tab$source[i] <- "(vectorized)"
  }
  rownames(tab) <- NULL
  tab
}

#' The [default_prior()] table for an assembled model.
#'
#' @noRd
prior_table <- function(spec, frame, route) {
  multi <- length(spec$responses) > 1L
  rows <- list()
  # The fit route reads no registry, deliberately. It once did, and a
  # loaded sampling package then changed the answer of a call that had
  # asked about frm(): same call, same fit, a different table depending
  # on what was attached. Asking the caller which route they mean is
  # what makes the answer theirs rather than the search path's.
  defs <- if (identical(route, "sample")) {
    registered_prior_defaults(spec, frame)
  } else {
    list()
  }
  add <- function(class, coef = "", group = "", dpar = "", nlpar = "",
                  resp = "") {
    # a default speaks for a CLASS, not for one coefficient of it: brms
    # reports the class row and leaves the per-coefficient rows flat, and
    # a per-coefficient row here would claim a default nothing applies.
    # A default on cor is written class-wide, as brms's cor rows carry no
    # prefix; one on sd is written per prefix, as brms's sd rows are
    d <- if (nzchar(coef)) NULL else
      defs[[prior_slot_key(class,
                           if (identical(class, "cor")) "" else dpar,
                           nlpar, resp)]]
    # brms's column order, so a table read by position lines up with
    # the one brms returns
    rows[[length(rows) + 1L]] <<- data.frame(
      prior = d %||% "(flat)", class = class, coef = coef, group = group,
      resp = resp, dpar = dpar, nlpar = nlpar, lb = NA_real_,
      ub = NA_real_, source = "default"
    )
  }

  for (lp in frame[["linpreds"]]) {
    if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]])) next
    rspec <- spec$responses[[lp[["resp"]]]]
    # location dpars are the default target (dpar = ""), matching
    # set_prior()'s resolution
    dpar_lab <- lp_prior_dpar(rspec, lp[["dpar"]])
    resp_lab <- if (multi) lp[["resp"]] else ""
    nl_lab <- if (lp[["dpar"]] %in% (rspec$nlpars %||% character(0))) {
      lp[["dpar"]]
    } else {
      ""
    }
    cn <- colnames(lp[["X"]])
    # a nonlinear parameter's sub-formula is not centered, so its whole
    # coefficient vector including the intercept is class "b" - brms's
    # own listing, and what set_prior(nlpar =) addresses
    if (nzchar(nl_lab)) {
      if (length(cn)) {
        add("b", dpar = dpar_lab, nlpar = nl_lab, resp = resp_lab)
        for (co in cn) {
          add("b", coef = co, dpar = dpar_lab, nlpar = nl_lab,
              resp = resp_lab)
        }
      }
      next
    }
    # a distributional parameter with no predictor of its own is
    # addressed by its own class, on its own scale: brms lists the same
    # row, and set_prior() refuses the link-scale spelling there
    if (nzchar(dpar_lab) && !lp[["dpar"]] %in% rspec$primary_dpars &&
          !dpar_has_predictor(spec, lp)) {
      # unless the class is one this package refuses BY NAME. A mixture
      # proportion is a predictor-free dpar whose word collides with
      # frmtmb's own "theta", so `set_prior()` refuses it and points at
      # a spelling the shape gate then refuses back. Advertising a slot
      # nothing can fill is worse than leaving it out, and a flat row
      # for an unnameable class is already left out one branch above
      if (!is.null(brms_prior_class_refusal(lp[["dpar"]]))) next
      if (length(cn)) add(lp[["dpar"]], resp = resp_lab)
      next
    }
    if ("(Intercept)" %in% cn) {
      add("Intercept", dpar = dpar_lab, resp = resp_lab)
    } else if (!nzchar(dpar_lab) &&
                 identical(rspec$family[["type"]], "ordinal") &&
                 length(frame[["par_template"]][["tau_raw"]] %||%
                          numeric(0))) {
      # an ordinal family has no intercept column: the thresholds
      # replace it, and class "Intercept" is what addresses them here as
      # it does in brms
      add("Intercept", resp = resp_lab)
    }
    # cs() terms are class "b" rows under their own coef, as brms lists
    # them; the resolver reaches them through the same class
    others <- c(setdiff(cn, "(Intercept)"),
                vapply(lp[["cs"]] %||% list(), cs_term_coef, ""))
    if (length(others)) {
      add("b", dpar = dpar_lab, resp = resp_lab)
      for (co in others) {
        add("b", coef = co, dpar = dpar_lab, resp = resp_lab)
      }
    }
  }

  sd_rows <- list()
  cor_rows <- list()
  for (bk in frame[["re_blocks"]]) {
    key <- list(group = bk[["group_name"]], nlpar = block_nlpar(spec, frame,
                                                                bk),
                dpar = block_dpar(spec, frame, bk),
                resp = if (multi) block_resp(frame, bk) else "")
    # brms lists sd rows per prefix, so a block spanning mu and sigma
    # has a row for each (sd_spec_reach() is the rule the rows describe)
    sp <- unique(block_sd_prefix(spec, frame, bk)[, c("resp", "dpar",
                                                      "nlpar")])
    for (i in seq_len(nrow(sp))) {
      sd_rows[[length(sd_rows) + 1L]] <- list(
        group = bk[["group_name"]], nlpar = sp$nlpar[i], dpar = sp$dpar[i],
        resp = sp$resp[i])
    }
    if (identical(block_cor_prior(bk), "lkj")) {
      cor_rows[[length(cor_rows) + 1L]] <- key
    }
  }
  if (length(sd_rows)) {
    # the class-wide row, once per prefix
    for (k in sd_rows) add("sd", dpar = k$dpar, nlpar = k$nlpar,
                           resp = k$resp)
    for (k in sd_rows) {
      add("sd", group = k$group, dpar = k$dpar, nlpar = k$nlpar,
          resp = k$resp)
    }
  }
  if (length(cor_rows)) {
    add("cor")
    for (k in cor_rows) {
      add("cor", group = k$group, dpar = k$dpar, nlpar = k$nlpar,
          resp = k$resp)
    }
  }
  # the R-side residual structures, under the class names brms shows for
  # them. `resp` narrows a class to one response, as it does above
  for (rs in names(frame[["autocor"]] %||% list())) {
    ac <- frame[["autocor"]][[rs]]
    resp_lab <- if (multi) rs else ""
    for (cl in names(autocor_prior_classes)) {
      if (length(autocor_class_idx(ac, cl))) add(cl, resp = resp_lab)
    }
  }
  if (length(frame[["par_template"]][["thetar"]] %||% numeric(0))) add("rescor")

  # class "theta" is the raw internal escape hatch, and it spans all
  # three covariance components: the residual-correlation ones carry no
  # per-parameter natural-scale class, so these rows are the only way to
  # name one of them on its own.
  #
  # [[ ]] rather than $ here and below: `$theta` partial-matches
  # `thetaac` in a model that has a residual structure and no random
  # effects, which offered a random-effect row for a model with none
  if (length(frame[["par_template"]][["theta"]] %||% numeric(0))) add("theta")
  for (cmp in theta_components) {
    v <- frame[["par_template"]][[cmp]] %||% numeric(0)
    if (!length(v)) next
    for (nm in par_template_names(v, cmp)) add("theta", coef = nm)
  }

  out <- unique(do.call(rbind, rows))
  rownames(out) <- NULL
  attr(out, "route") <- route
  class(out) <- c("frmtmb_prior_rows", "data.frame")
  out
}

#' @export
print.frmtmb_prior_rows <- function(x, ...) {
  # a table of default priors that does not say which route it
  # describes is the ambiguity `route` exists to remove, so the label
  # travels with the object rather than with the call that made it. A
  # ROW subset keeps it, correctly: fewer rows of a fit-route table are
  # still a fit-route table. A COLUMN subset does not, because
  # `[.data.frame` drops the attribute there, and a table missing the
  # `prior` column has no route to claim anyway. The NULL branch prints
  # such a table, and any assembled by other means, without a label.
  route <- attr(x, "route")
  if (identical(route, "fit")) {
    cat("route = \"fit\": the prior defaults frm() applies\n")
  } else if (identical(route, "sample")) {
    cat("route = \"sample\": the prior defaults frm_sample() applies\n")
  }
  print(as.data.frame(x), ...)
  invisible(x)
}

#' @export
as.data.frame.frmtmb_prior_rows <- function(x, row.names = NULL,
                                            optional = FALSE, ...) {
  # the label belongs to the class, so coercion drops both together. A
  # coerced table that kept a stray `route` attribute would not be
  # identical() to the plain data frame a reader would write by hand,
  # which is the one thing coercion is asked for.
  attr(x, "route") <- NULL
  class(x) <- "data.frame"
  as.data.frame(x, row.names = row.names, optional = optional, ...)
}

#' The coef name brms gives a cs() term's row: the term with its `cs`
#' prefix dropped, spelled as brms spells a column name (`cs(x)` is
#' `"x"`).
#'
#' @noRd
cs_term_coef <- function(ct) {
  brms_rename(sub("^cs", "", ct[["label"]]))
}

#' Whether a response's location is several distributional parameters,
#' each with its own prior rows.
#'
#' brms names them (`mub`, `muc` for a categorical response; `mu1`,
#' `mu2` for a mixture) and lists every `b`, `Intercept` and `sd` row
#' under that `dpar`, never under an empty one (dev/mvprior-log/
#' brms-probe.txt, cases cat, mix, catre, mixre). A categorical or
#' multinomial response keeps that naming with ONE non-reference
#' category, as brms does. Every other family with several location
#' dpars (the latent-class, hidden Markov, mixture_mvn() and
#' accumulator families) follows the same rule, which is frmtmb's for
#' families brms does not have: a prior with no dpar does not say which
#' location it means.
#'
#' @noRd
multi_location <- function(rspec) {
  fam <- rspec$family
  loc <- setdiff(rspec$primary_dpars, rspec$nlpars %||% character(0))
  length(loc) > 1L || !is.null(fam[["mix"]]) ||
    isTRUE(fam[["family"]] %in% c("categorical", "multinomial"))
}

#' Whether brms has this several-location family, so a refusal may say
#' "as in brms" rather than name frmtmb's own rule.
#'
#' @noRd
brms_multi_location <- function(rspec) {
  fam <- rspec$family
  fn <- fam[["family"]] %||% ""
  fn %in% c("categorical", "multinomial") ||
    (startsWith(fn, "mixture(") && !is.null(fam[["mix"]]) &&
       !identical(fam[["mix"]][["type"]], "categorical"))
}

#' The `dpar` a linear predictor's prior rows carry, brms's label: `""`
#' for the location parameter of a family that has one, and for a
#' nonlinear parameter (which is addressed by `nlpar`); the parameter's
#' own name otherwise. `primary_dpars` alone is not the location,
#' because it also holds the nonlinear parameters REML integrates and
#' every location of a family that has several.
#'
#' @noRd
lp_prior_dpar <- function(rspec, dp) {
  if (dp %in% (rspec$nlpars %||% character(0))) return("")
  if (dp %in% rspec$primary_dpars && !multi_location(rspec)) return("")
  dp
}

#' The linear predictors a random-effect block draws its columns from.
#' A block merges the components that share an `|ID|` key, and those can
#' come from different predictors, so this is a set rather than one
#' value; `bk$dpar` is the first component's and is the fallback for a
#' block whose components predate the key.
#'
#' @noRd
block_linpreds <- function(frame, bk) {
  keys <- vapply(bk[["components"]] %||% list(),
                 function(cp) cp$lp_key %||% "", "")
  lps <- frame[["linpreds"]][keys[nzchar(keys)]]
  Filter(Negate(is.null), lps)
}

#' The nonlinear parameter a block belongs to, or `""`. Blocks that
#' straddle several are left unlabeled rather than assigned to one.
#'
#' @noRd
block_nlpar <- function(spec, frame, bk) {
  lps <- block_linpreds(frame, bk)
  np <- unique(vapply(lps, function(lp) {
    nl <- spec$responses[[lp[["resp"]]]]$nlpars %||% character(0)
    if (lp[["dpar"]] %in% nl) lp[["dpar"]] else ""
  }, ""))
  np <- np[nzchar(np)]
  if (length(np) == 1L) np else ""
}

#' The distributional parameter a block belongs to, or `""` for the
#' location parameter of a family that has one, a nonlinear parameter,
#' and a block that straddles several. brms lists a block of
#' `phi ~ (1 | g)` as class `"sd"` with `dpar = "phi"`, and without the
#' column two blocks on the same factor in different predictors printed
#' as one row.
#'
#' @noRd
block_dpar <- function(spec, frame, bk) {
  dp <- unique(vapply(block_linpreds(frame, bk), function(lp) {
    lp_prior_dpar(spec$responses[[lp[["resp"]]]], lp[["dpar"]])
  }, ""))
  if (length(dp) == 1L) dp else ""
}

#' The response a block belongs to, or `""` when it straddles several.
#'
#' @noRd
block_resp <- function(frame, bk) {
  rs <- unique(vapply(block_linpreds(frame, bk), `[[`, "", "resp"))
  if (length(rs) == 1L) rs else ""
}

#' Whether a class `"cor"` specification addresses this block.
#' `group` has always narrowed here; `resp`, `dpar` and `nlpar` narrow
#' the same way. Class `"sd"` is resolved per standard deviation by
#' `sd_spec_reach()` instead.
#'
#' @noRd
block_addressed <- function(spec, frame, bk, s) {
  if (nzchar(s$group) && !identical(bk[["group_name"]], s$group)) return(FALSE)
  want_np <- nzchar(s$nlpar %||% "")
  want_rs <- nzchar(s$resp %||% "")
  want_dp <- nzchar(s$dpar)
  if (!want_np && !want_rs && !want_dp) return(TRUE)
  lps <- block_linpreds(frame, bk)
  if (want_rs && !s$resp %in% vapply(lps, `[[`, "", "resp")) return(FALSE)
  dp <- vapply(lps, `[[`, "", "dpar")
  if (want_np && !identical(block_nlpar(spec, frame, bk), s$nlpar)) {
    return(FALSE)
  }
  if (want_dp && !s$dpar %in% dp) return(FALSE)
  TRUE
}

#' The prefix brms gives each COLUMN of a block: the response (in a
#' multivariate model), the distributional parameter (`""` for the
#' location parameters) and the nonlinear parameter (`""` for none).
#' brms keys its class `"sd"` prior rows by this prefix.
#'
#' One row per column of the block, in column order. A block with no
#' component record takes the block-level labels for every column.
#'
#' @noRd
block_col_prefix <- function(spec, frame, bk) {
  multi <- length(spec$responses) > 1L
  d <- bk[["dim"]] %||% 0L
  out <- data.frame(resp = rep(if (multi) block_resp(frame, bk) else "", d),
                    dpar = rep(block_dpar(spec, frame, bk), d),
                    nlpar = rep(block_nlpar(spec, frame, bk), d),
                    stringsAsFactors = FALSE)
  for (cp in bk[["components"]] %||% list()) {
    lp <- frame[["linpreds"]][[cp[["lp_key"]] %||% ""]]
    if (is.null(lp) || is.null(cp[["offset"]]) || is.null(cp[["dim"]])) next
    rspec <- spec$responses[[lp[["resp"]]]]
    nl <- rspec$nlpars %||% character(0)
    dp <- lp[["dpar"]]
    cols <- cp[["offset"]] + seq_len(cp[["dim"]])
    cols <- cols[cols <= d]
    out$resp[cols] <- if (multi) lp[["resp"]] else ""
    out$dpar[cols] <- lp_prior_dpar(rspec, dp)
    out$nlpar[cols] <- if (dp %in% nl) dp else ""
  }
  out
}

#' The prefix of each standard deviation of a block, one row per entry
#' of `block_sd_idx()`, which `k` numbers. A structure with one standard
#' deviation per column takes the column's prefix; one whose standard
#' deviations are shared across columns takes every prefix its columns
#' have, one row each under the same `k`.
#'
#' @noRd
block_sd_prefix <- function(spec, frame, bk) {
  sd_i <- block_sd_idx(bk)
  cp <- block_col_prefix(spec, frame, bk)
  if (!length(sd_i)) return(cbind(cp[0L, , drop = FALSE], k = integer(0)))
  if (length(sd_i) == nrow(cp)) return(cbind(cp, k = seq_along(sd_i)))
  u <- unique(cp)
  do.call(rbind, lapply(seq_along(sd_i), function(j) cbind(u, k = j)))
}

#' Which standard deviations of a block a class `"sd"` specification
#' reaches, by brms 2.23.0's rule (measured in
#' dev/correct-log/brms-priors.txt and brms-priors2.txt).
#'
#' brms lists class `"sd"` rows per prefix (response, distributional
#' parameter, nonlinear parameter), and a specification applies where
#' its prefix EQUALS the standard deviation's: with no `dpar` it is the
#' location parameter's, not every predictor's; in a multivariate model
#' with no `resp` it is no response's, and brms refuses it. The one
#' exception is a block that spans several prefixes, `(1 | q | g)` in
#' `mu` and `sigma`: brms reads its prior rows per prefix with each
#' field also matching an empty one, so a specification with fewer
#' fields reaches every column there, and the more specific one wins,
#' which the caller's specificity order already arranges.
#'
#' Returns a list: `reach`, a logical per entry of `block_sd_idx()`, and
#' `exact`, whether the specification's own prefix is one this block
#' has, which is what brms requires for the row to exist at all.
#'
#' @noRd
sd_spec_reach <- function(spec, frame, bk, s) {
  sp <- block_sd_prefix(spec, frame, bk)
  n <- length(block_sd_idx(bk))
  none <- list(reach = logical(n), exact = FALSE)
  if (!n) return(none)
  if (nzchar(s$group) && !identical(bk[["group_name"]], s$group)) {
    return(none)
  }
  want <- c(resp = s$resp %||% "", dpar = s$dpar %||% "",
            nlpar = s$nlpar %||% "")
  eq <- sp$resp == want[["resp"]] & sp$dpar == want[["dpar"]] &
    sp$nlpar == want[["nlpar"]]
  spans <- nrow(unique(sp[, c("resp", "dpar", "nlpar")])) > 1L
  fits <- if (spans) {
    (sp$resp == want[["resp"]] | !nzchar(want[["resp"]])) &
      (sp$dpar == want[["dpar"]] | !nzchar(want[["dpar"]])) &
      (sp$nlpar == want[["nlpar"]] | !nzchar(want[["nlpar"]]))
  } else {
    eq
  }
  reach <- vapply(seq_len(n), function(j) any(fits[sp$k == j]), TRUE)
  list(reach = reach, exact = any(eq))
}

#' What a refused class `"sd"` specification could have named: the
#' prefixes this model's standard deviations carry, as `set_prior()`
#' arguments.
#'
#' @noRd
sd_prefix_hint <- function(spec, frame) {
  sp <- unique(do.call(rbind, c(
    list(data.frame(resp = character(0), dpar = character(0),
                    nlpar = character(0))),
    lapply(frame[["re_blocks"]], function(bk) {
      block_sd_prefix(spec, frame, bk)[, c("resp", "dpar", "nlpar")]
    }))))
  if (!NROW(sp)) return("This model has no random-effect standard deviations")
  lab <- apply(sp, 1L, function(r) {
    f <- c(if (nzchar(r[["resp"]])) paste0("resp = \"", r[["resp"]], "\""),
           if (nzchar(r[["dpar"]])) paste0("dpar = \"", r[["dpar"]], "\""),
           if (nzchar(r[["nlpar"]])) paste0("nlpar = \"", r[["nlpar"]], "\""))
    if (length(f)) paste(f, collapse = ", ") else "no resp, dpar or nlpar"
  })
  paste0("As in brms, class \"sd\" is addressed per response, ",
         "distributional parameter and nonlinear parameter, and a ",
         "specification without one of these names none of the others. ",
         "This model's standard deviations take: ",
         paste0("(", lab, ")", collapse = "; "))
}

#' Where a specification was pointing, for the refusals that report a
#' target no design offers. Written once so both refusals name the same
#' fields in the same order.
#'
#' @noRd
spec_target <- function(s) {
  sp <- spec_spelling(s)
  paste0("class=", sp$class,
         if (nzchar(s$coef)) paste0(", coef=", s$coef),
         if (nzchar(s$group)) paste0(", group=", s$group),
         if (nzchar(s$resp %||% "")) paste0(", resp=", s$resp),
         if (nzchar(sp$dpar)) paste0(", dpar=", sp$dpar),
         if (nzchar(s$nlpar %||% "")) paste0(", nlpar=", s$nlpar))
}

#' The class and dpar a specification is WRITTEN with, which is not
#' always the pair it is stored as.
#'
#' A density on a distributional parameter itself is stored as that
#' parameter's intercept slot plus the `natural` flag, because that is
#' the slot the resolver assigns to. It is written, in frmtmb and in
#' brms alike, as the parameter's own class. Printing it back as
#' `class = "Intercept", dpar = "sigma"` would print a spelling this
#' model refuses, which is the opposite of what
#' `print.frmtmb_priorlist()` promises: that what it prints can be
#' pasted back into `set_prior()`.
#'
#' @noRd
spec_spelling <- function(s) {
  dp <- s$dpar %||% ""
  if (isTRUE(s$natural) && identical(s$class, "Intercept") &&
        nzchar(dp)) {
    return(list(class = dp, dpar = ""))
  }
  list(class = s$class, dpar = dp)
}

#' The nonlinear parameters a model declares, for the refusal that has
#' to say what `nlpar` could have named.
#'
#' @noRd
model_nlpars <- function(spec) {
  unique(unlist(lapply(spec$responses, function(r) {
    r$nlpars %||% character(0)
  }))) %||% character(0)
}

#' Copy the bounds of a prior specification onto an existing entry. This
#' is how a coefficient's bounds-only specification boxes the entry a
#' class-wide distribution created, instead of being lost.
#'
#' @noRd
entry_bounds <- function(entry, s) {
  if (!is.na(s$lb)) entry$lb <- s$lb
  if (!is.na(s$ub)) entry$ub <- s$ub
  entry
}

#' The specifications of a priorlist in the order they are applied: the
#' positions each class holds are kept, and within a class the less
#' specific specifications move before the more specific ones.
#'
#' brms applies a coefficient's prior over its class's and a group's
#' over the class-wide one whatever the order written, and so does this.
#' Before, the later specification won, so `coef = "x"` written before a
#' class `"b"` prior was silently overridden by it. Keeping each class
#' in its own positions leaves the order between classes alone, which
#' is the only precedence `"cor"` and the `"theta"` hatch have; a stable
#' sort keeps equally specific specifications in written order, which
#' is the order frmtmb.sample stacks its defaults in.
#'
#' @noRd
prior_specificity_order <- function(pl) {
  specs <- unclass(pl)
  if (length(specs) < 2L) return(specs)
  nz <- function(s, f) as.integer(nzchar(s[[f]] %||% ""))
  cls <- vapply(specs, function(s) spec_spelling(s)$class, "")
  rank <- vapply(specs, function(s) {
    100L * nz(s, "coef") + 10L * nz(s, "group") + nz(s, "resp") +
      nz(s, "dpar") + nz(s, "nlpar")
  }, 0L)
  out <- specs
  for (k in unique(cls)) {
    at <- which(cls == k)
    out[at] <- specs[at[order(rank[at])]]
  }
  out
}

#' The slot one specification addresses, brms's key: class as written,
#' coef, group, resp, dpar and nlpar.
#'
#' @noRd
prior_slot_label <- function(s) {
  sp <- spec_spelling(s)
  df <- data.frame(prior = "", class = sp$class,
                   coef = par_name_bare(s$coef %||% ""),
                   group = s$group %||% "", resp = s$resp %||% "",
                   dpar = sp$dpar, nlpar = s$nlpar %||% "",
                   lb = NA_character_, ub = NA_character_,
                   stringsAsFactors = FALSE)
  sub(" ~ .*$", "", prior_row_lines(df))
}

#' The slot one specification addresses as a comparable key: class as
#' written, coef without the parentheses of `(Intercept)`, group, resp,
#' dpar and nlpar, each kept apart.
#'
#' @noRd
prior_spec_slot <- function(s) {
  sp <- spec_spelling(s)
  paste(sp$class, par_name_bare(s$coef %||% ""), s$group %||% "",
        s$resp %||% "", sp$dpar, s$nlpar %||% "", sep = "\r")
}

#' Refuse two specifications for the same slot
#'
#' brms refuses them, whether the two are identical, carry different
#' densities, or are a density and a bounds-only specification. Earlier
#' releases applied the later one, and a later bounds-only specification
#' tightened an earlier density. Called where a user passes a prior:
#' [frm()], [validate_prior()], [frm_simulate()], [par_template()] and
#' `frmtmb.sample::frm_sample()`. Not inside the resolver, because
#' frmtmb.sample stacks its own specifications there.
#'
#' @param prior A `frmtmb_priorlist`; anything else is returned
#'   unchecked.
#' @return `prior`, invisibly.
#' @noRd
check_prior_slots <- function(prior) {
  if (!inherits(prior, "frmtmb_priorlist")) return(invisible(prior))
  specs <- unclass(prior)
  # the fields themselves, not the name brms prints: b with dpar = "sigma"
  # and b with coef = "sigma" both print b_sigma and are two slots
  keys <- vapply(specs, prior_spec_slot, "")
  dup_at <- which(duplicated(keys))
  if (length(dup_at)) {
    first <- dup_at[!duplicated(keys[dup_at])]
    labels <- vapply(specs[first], prior_slot_label, "")
    counts <- vapply(keys[first], function(k) sum(keys == k), 0L)
    frm_stop("Duplicated prior specifications are not allowed: ",
             paste0("'", labels, "' is given ", counts, " times",
                    collapse = "; "),
             ". Write one specification per slot, with its density and ",
             "both bounds in the same call, e.g. ",
             "set_prior(\"normal(0, 1)\", class = \"b\", lb = 0)",
             call. = FALSE)
  }
  invisible(prior)
}

#' The covariance components class `"theta"` addresses, in the order it
#' searches them for a `coef` name.
#'
#' @noRd
theta_components <- c("theta", "thetaac", "thetar")

#' The classes whose prior rows brms keys by response in a multivariate
#' model, so that a specification without `resp` names no row there.
#' Class `"sd"` is keyed the same way and is resolved per standard
#' deviation by `sd_spec_reach()`. Classes `"cor"` and `"rescor"` carry
#' no prefix in brms, and `"theta"` is frmtmb's own.
#'
#' @noRd
resp_keyed_prior_classes <- c("b", "Intercept", "ar", "ma", "cosy",
                              "cortime")

#' Why a specification with no `resp` cannot apply to this model, or
#' `NULL`.
#'
#' brms 2.23.0 refuses a class `"b"`, `"Intercept"`, distributional
#' parameter or residual-correlation prior with no `resp` in a
#' multivariate model: its rows there are keyed by response, and an
#' empty `resp` matches none of them (measured in
#' dev/mvprior-log/brms-probe.txt). frmtmb applied such a specification
#' to every response that had the slot, which gave a different model
#' from the same call.
#'
#' @noRd
resp_missing_refusal <- function(spec, frame, s) {
  if (length(spec$responses) < 2L || nzchar(s$resp %||% "")) return(NULL)
  if (!s$class %in% resp_keyed_prior_classes) return(NULL)
  sp <- spec_spelling(s)
  # the responses the same specification resolves on once it names
  # them, found by resolving it: reading them off the default_prior()
  # table offered calls the resolver then refused, and missed ones it
  # accepts (a dpar the table did not list, class Intercept with nlpar)
  design <- list(spec = spec, frame = frame)
  resolves <- function(s1) {
    one <- structure(list(s1), class = "frmtmb_priorlist")
    !inherits(tryCatch(resolve_priorlist(design, one), error = identity),
              "error")
  }
  # (resp, dpar) pairs. A response whose location is several dpars
  # takes a location prior only with its dpar, so those are offered as
  # well, one per location
  cand <- list()
  for (r in names(spec$responses)) {
    s1 <- s
    s1$resp <- r
    if (resolves(s1)) {
      cand[[length(cand) + 1L]] <- c(r, sp$dpar)
      next
    }
    rspec <- spec$responses[[r]]
    if (nzchar(sp$dpar) || nzchar(s$nlpar %||% "") ||
          !s$class %in% c("b", "Intercept") || !multi_location(rspec)) {
      next
    }
    locs <- setdiff(rspec$primary_dpars, rspec$nlpars %||% character(0))
    for (dp in locs) {
      s2 <- s1
      s2$dpar <- dp
      if (resolves(s2)) cand[[length(cand) + 1L]] <- c(r, dp)
    }
  }
  arg <- function(rd) {
    r <- rd[[1L]]
    dp <- rd[[2L]]
    paste0("set_prior(",
           if (is.null(s$dist)) "\"\"" else
             encodeString(s$prior %||% "", quote = "\""),
           ", class = \"", sp$class, "\"",
           if (nzchar(s$coef %||% "")) paste0(", coef = \"", s$coef, "\""),
           if (nzchar(dp)) paste0(", dpar = \"", dp, "\""),
           if (nzchar(s$nlpar %||% "")) paste0(", nlpar = \"", s$nlpar,
                                               "\""),
           ", resp = \"", r, "\"",
           if (!is.na(s$lb)) paste0(", lb = ", format(s$lb)),
           if (!is.na(s$ub)) paste0(", ub = ", format(s$ub)), ")")
  }
  rs <- cand
  paste0("A prior with no resp names no parameter of a multivariate ",
         "model (", spec_target(s), "). As in brms, class \"", sp$class,
         "\" is addressed per response here, and a specification ",
         "without resp reaches none of them. ",
         if (length(rs)) {
           paste0("Write one specification per slot this model has: ",
                  paste(vapply(rs, arg, ""), collapse = " + "))
         } else {
           paste0("No response of this model has this slot; ",
                  "default_prior() lists the ones it has")
         })
}

#' Which component and positions a class `"theta"` specification names.
#' An empty `coef` means the whole `theta` component, which is what the
#' class has always meant; a `coef` is matched against the internal
#' names of all three covariance components, so `"thetaac_1"` and
#' `"thetar_2"` reach the residual-correlation parameters that carry no
#' natural-scale class of their own.
#'
#' @noRd
theta_coef_target <- function(frame, coef) {
  nms_of <- function(cmp) {
    v <- frame[["par_template"]][[cmp]] %||% numeric(0)
    if (!length(v)) character(0) else par_template_names(v, cmp)
  }
  if (!nzchar(coef %||% "")) {
    n <- length(nms_of("theta"))
    if (!n) {
      frm_stop("class = \"theta\" names no parameter: this model has no ",
               "random-effect covariance parameters. The residual-",
               "correlation ones are addressed by name, e.g. ",
               "coef = \"thetaac_1\"", call. = FALSE)
    }
    return(list(comp = "theta", idx = seq_len(n)))
  }
  for (cmp in theta_components) {
    k <- match(coef, nms_of(cmp))
    if (!is.na(k)) return(list(comp = cmp, idx = k))
  }
  # the bare position keeps the older `coef = "2"` spelling working
  if (grepl("^[0-9]+$", coef)) {
    k <- as.integer(coef)
    if (k >= 1L && k <= length(nms_of("theta"))) {
      return(list(comp = "theta", idx = k))
    }
  }
  have <- unlist(lapply(theta_components, nms_of), use.names = FALSE)
  frm_stop("class = \"theta\" coef = ", encodeString(coef, quote = "\""),
           " names no covariance parameter of this model. It has ",
           if (length(have)) paste(have, collapse = ", ") else "none",
           call. = FALSE)
}

#' Resolve a priorlist against a fit: per-parameter prior entries (the
#' more specific specification of a class applies whatever the order,
#' and between classes the later one) plus named bound vectors. Same-slot
#' duplicates are refused where users pass priors, not here, because
#' frmtmb.sample stacks its own defaults under a fit's prior and a
#' call's prior before resolving them. An
#' entry holds `comp`, a scalar `idx`, `dist`, `scale` ("internal" or
#' "sd"), and `lb`/`ub` on the entry's own scale. `frm_simulate()`
#' rejects draws outside those bounds; `frm_sample()` uses the named
#' bound vectors instead.
#'
#' @noRd
resolve_priorlist <- function(fit, pl) {
  frame <- fit$frame
  assigned <- list()   # key "comp.idx" -> entry
  lower <- c()
  upper <- c()
  # a class "cor" entry covers SEVERAL theta positions at once, so the
  # key names them all; a one-position entry keys exactly as before
  nm_of <- function(comp, idx) paste0(comp, ".", paste(idx, collapse = ","))

  # Assigning over a position that a joint (class "cor") entry already
  # covers has to RETIRE that entry, and so does a joint entry covering
  # positions that per-parameter entries claimed. Otherwise both
  # densities would be added. With that, "later wins" reads the same
  # across the two spellings, which are different classes and so keep
  # their written order.
  claim <- function(comp, idx) {
    drop <- vapply(assigned, function(e) {
      identical(e$comp, comp) && length(intersect(e$idx, idx)) > 0L
    }, TRUE)
    assigned <<- assigned[!drop]
  }

  nlpars <- model_nlpars(fit$spec)

  target_coefs <- function(s) {
    # (comp, idx, name) triplets for classes b / Intercept
    out <- list()
    want_np <- nzchar(s$nlpar %||% "")
    if (want_np && !s$nlpar %in% nlpars) {
      frm_stop("nlpar = \"", s$nlpar, "\" names no nonlinear parameter of ",
               "this model. It has ",
               if (length(nlpars)) paste(nlpars, collapse = ", ") else
                 "none (write nl = TRUE in bf() to declare them)",
               ". A distributional parameter is addressed with dpar =",
               call. = FALSE)
    }
    nl_loc <- character(0)
    multi_loc <- character(0)
    loc_brms <- TRUE
    for (lp in frame[["linpreds"]]) {
      if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]])) next
      rspec <- fit$spec$responses[[lp[["resp"]]]]
      if (nzchar(s$resp %||% "") && !identical(lp[["resp"]], s$resp)) next
      # primary_dpars carries the nonlinear parameters a nonlinear
      # location is built from, because REML integrates them. For a
      # prior they are not the location: brms keys their rows by nlpar
      # and refuses a specification without it (brms-probe.txt, case nl)
      is_nlpar <- lp[["dpar"]] %in% (rspec$nlpars %||% character(0))
      is_prim <- lp[["dpar"]] %in% rspec$primary_dpars
      # a location that is one of several is addressed by its dpar, as
      # brms lists it; without one the specification used to reach all
      is_loc <- is_prim && !is_nlpar &&
        !nzchar(lp_prior_dpar(rspec, lp[["dpar"]]))
      if (want_np) {
        if (!identical(lp[["dpar"]], s$nlpar)) next
      } else if (nzchar(s$dpar)) {
        if (!identical(lp[["dpar"]], s$dpar)) next
      } else if (!is_loc) {
        if (is_nlpar && is_prim) {
          nl_loc <- c(nl_loc, lp[["dpar"]])
        } else if (is_prim) {
          multi_loc <- c(multi_loc, lp[["dpar"]])
          loc_brms <- loc_brms && brms_multi_location(rspec)
        }
        next
      }
      cn <- colnames(lp[["X"]])
      pick <- if (s$class == "Intercept") {
        which(cn == "(Intercept)")
      } else if (nzchar(s$coef)) {
        # brms writes an intercept as "Intercept"; the design matrix
        # spells it "(Intercept)", and both name the same column
        which(cn == s$coef | par_name_bare(cn) == par_name_bare(s$coef))
      } else if (want_np) {
        # a nonlinear parameter's sub-formula is NOT centered, so its
        # intercept sits in the same coefficient vector as its slopes
        # and class "b" covers it. This is the whole reason brms's
        # prior(normal(5000, 1000), nlpar = "ult") lands on an
        # intercept-only nonlinear parameter
        seq_along(cn)
      } else {
        which(cn != "(Intercept)")
      }
      # `name` is what a BOUND is keyed by, and resolve_bounds() matches
      # against outer_par_names(): the template spelling, which carries
      # the dpar/nlpar/resp prefix ("guess_(Intercept)"). The design
      # matrix column alone ("(Intercept)") names no outer parameter and
      # collides across sub-formulas, so read the name off the template
      # position rather than off the column
      pnm <- par_template_names(frame[["par_template"]][[lp[["par"]]]],
        lp[["par"]])
      # `link` answers where a `natural` density belongs, and `center`
      # is brms's centering offset for THIS sub-formula's intercept; both
      # are properties of the linear predictor rather than of the
      # specification, so they are read here and carried on the target
      ctr <- if (length(pick) && s$class == "Intercept" && !want_np) {
        lp_center_offset(frame, lp)
      }
      for (k in pick) {
        out[[length(out) + 1L]] <- list(comp = lp[["par"]],
                                        idx = lp[["idx"]][k],
                                        name = pnm[lp[["idx"]][k]],
                                        link = lp[["link"]],
                                        center = ctr)
      }
      # a cs() term's threshold-specific coefficients are class "b" in
      # brms, under the term's own coef name: class "b" puts the density
      # on every one of them (to_vector(bcs)), coef = "x" on the row of x
      # (bcs[1]) (dev/mvprior-log/cs-brms.txt). They live in their own
      # component, outside the design matrix, which is why class "b"
      # used to miss them: nothing on ord ~ cs(x), the ordinary slopes
      # only on ord ~ z + cs(x)
      if (identical(s$class, "b") && !want_np) {
        for (ct in lp[["cs"]] %||% list()) {
          lab <- cs_term_coef(ct)
          if (nzchar(s$coef) &&
                !identical(par_name_bare(lab), par_name_bare(s$coef))) {
            next
          }
          v <- frame[["par_template"]][[ct[["par"]]]]
          cnm <- par_template_names(v, ct[["par"]])
          for (k in seq_along(v)) {
            out[[length(out) + 1L]] <- list(comp = ct[["par"]], idx = k,
                                            name = cnm[k], link = NULL,
                                            center = NULL)
          }
        }
      }
    }
    if (!length(out) && length(nl_loc)) {
      # before, class "b" here matched nothing and was silently dropped,
      # and class "Intercept" or a coef reached every nonlinear
      # parameter at once
      frm_stop("Prior target not found (", spec_target(s), "): this ",
               "model's location is nonlinear, and as in brms its ",
               "coefficients are class = \"b\" with nlpar = one of ",
               paste(unique(nl_loc), collapse = ", "), ", the intercept ",
               "among them. A class \"Intercept\" or coef = \"Intercept\" ",
               "prior names one nonlinear parameter's intercept with nlpar = ",
               "as well; without it, it reaches none",
               call. = FALSE)
    }
    if (!length(out) && length(multi_loc)) {
      locs <- unique(multi_loc)
      one <- length(locs) == 1L
      frm_stop("Prior target not found (", spec_target(s), "): this ",
               "model's location is ",
               if (one) {
                 paste0("the distributional parameter ", locs,
                        ", and a class \"", s$class, "\" prior names it ",
                        "with dpar = \"", locs, "\"")
               } else {
                 paste0("several distributional parameters, ",
                        paste(locs, collapse = ", "), ", and a class \"",
                        s$class, "\" prior names one of them with dpar =, ",
                        "one specification per parameter")
               },
               if (loc_brms) ", as in brms" else
                 paste0(". That is frmtmb's rule for a family whose ",
                        "location is several distributional parameters"),
               call. = FALSE)
    }
    if (!length(out) && identical(s$class, "b") && !nzchar(s$coef) &&
          !want_np) {
      # before, this matched nothing and was dropped without a word,
      # while prior_summary() still listed it (user decision, 2026-09-24)
      frm_stop("Prior target not found (", spec_target(s), "): the ",
               if (nzchar(s$dpar)) paste0(s$dpar, " predictor") else
                 "location predictor",
               if (nzchar(s$resp %||% "")) paste0(" of response ", s$resp),
               " has no population-level slope, so a class \"b\" prior ",
               "reaches no parameter. brms refuses it too. A prior on the ",
               "intercept is class = \"Intercept\"", call. = FALSE)
    }
    if (!length(out) &&
          (nzchar(s$coef) || s$class == "Intercept" || want_np ||
             nzchar(s$resp %||% ""))) {
      # an ordinal model has no intercept column at all: the thresholds
      # replace it, and a bare class = "Intercept" reaches them. Saying
      # so here is the difference between "this model has no such slot"
      # and "you narrowed the row past the slot it has"
      hint <- if (s$class == "Intercept" && has_ordinal_thresholds(fit)) {
        paste0(". On an ordinal family the thresholds ARE the ",
               "intercept, and a bare class = \"Intercept\" with no ",
               "coef, dpar or nlpar addresses the whole threshold ",
               "vector; prior = list(tau_raw = ) reaches the same ",
               "parameters on the internal scale")
      } else {
        ""
      }
      frm_stop("Prior target not found (", spec_target(s), ")", hint,
               call. = FALSE)
    }
    out
  }

  # One residual-autocorrelation class over every block that carries it.
  # `resp` narrows to one response, as it does everywhere else; the
  # refusal names what the model actually has, so a class aimed at the
  # wrong structure says so rather than silently matching nothing.
  resolve_ac_class <- function(s) {
    acs <- frame[["autocor"]] %||% list()
    hit <- FALSE
    for (rs in names(acs)) {
      ac <- acs[[rs]]
      if (nzchar(s$resp %||% "") && !identical(rs, s$resp)) next
      within <- autocor_class_idx(ac, s$class)
      if (is.null(within) || !length(within)) next
      hit <- TRUE
      idx <- ac[["theta_idx"]][within]
      tr <- autocor_trans(ac, s$class)
      nms <- par_template_names(frame[["par_template"]][["thetaac"]], "thetaac")
      if (!is.null(s$dist)) {
        dst <- if (identical(s$class, "cortime")) {
          lkj_dist(s$dist$eta, list(kind = "chol", d = ac[["d"]]))
        } else {
          trans_dist(s$dist, tr)
        }
        claim("thetaac", idx)
        assigned[[nm_of("thetaac", idx)]] <<-
          list(comp = "thetaac", idx = idx, dist = dst,
               scale = "internal", lb = NA, ub = NA)
      }
      # set_prior() has already refused lb/ub on the matrix-valued
      # classes, so only the transformed scalar maps reach this
      if (!is.na(s$lb) || !is.na(s$ub)) {
        # a box on unconstrained coefficients IS a box on the internal
        # ones, so only the transformed maps are limited to order one
        if (length(idx) > 1L && !identical(tr$map, "identity")) {
          frm_stop("class = \"", s$class, "\" takes no lb/ub at order ",
                   length(idx), ": coefficient ", s$class,
                   "[1] is a function of every one of this block's ",
                   length(idx), " internal parameters, so a bound on it is ",
                   "not a bound on any of them. The parameterization already ",
                   "keeps the process ",
                   if (identical(s$class, "ar")) "stationary" else "invertible",
                   "; bound an internal parameter with class = \"theta\", ",
                   "coef = \"", nms[idx[1L]], "\" if a box is really wanted",
                   call. = FALSE)
        }
        if (!is.na(s$lb)) {
          lower[nms[idx]] <<- ac_bound_theta(s$lb, tr, ac, s$class, "lb")
        }
        if (!is.na(s$ub)) {
          upper[nms[idx]] <<- ac_bound_theta(s$ub, tr, ac, s$class, "ub")
        }
      }
    }
    if (!hit) {
      have <- vapply(names(acs), function(rs) {
        paste0(acs[[rs]]$label, " [", rs, "]")
      }, "")
      frm_stop("No residual autocorrelation matches ", spec_target(s), ". ",
               if (length(have)) {
                 paste0("This model's residual structure is ",
                        paste(have, collapse = ", "),
                        ", which carries no \"", s$class, "\" parameter")
               } else {
                 paste0("This model has no residual autocorrelation term ",
                        "(write one with ar(), ma(), arma(), cosy() or ",
                        "unstr() in the formula)")
               }, call. = FALSE)
    }
  }

  # The residual correlation of a multivariate model: one unstructured
  # matrix held as the same row-normalized Cholesky a `us` block uses,
  # so the LKJ density and its Jacobian carry over unchanged.
  resolve_rescor <- function(s) {
    # one correlation matrix across the responses, so a resp names a part
    # of it that has no prior of its own; brms refuses it ("Lrescor_y1"),
    # and 0.62.0 applied it to the whole matrix with the resp ignored
    if (nzchar(s$resp %||% "")) {
      frm_stop("class = \"rescor\" takes no resp: the residual ",
               "correlation is one set of parameters across all responses, ",
               "and lkj() is a density on the whole matrix. Drop resp = \"",
               s$resp, "\", as brms requires", call. = FALSE)
    }
    n_r <- length(frame[["par_template"]][["thetar"]] %||% numeric(0))
    if (!n_r) {
      frm_stop("No residual correlation matches ", spec_target(s),
               ". This model has none: it needs two or more responses and ",
               "set_rescor(TRUE)", call. = FALSE)
    }
    idx <- seq_len(n_r)
    claim("thetar", idx)
    assigned[[nm_of("thetar", idx)]] <<-
      list(comp = "thetar", idx = idx,
           dist = lkj_dist(s$dist$eta, list(kind = "chol",
                                            d = length(fit$spec$responses))),
           scale = "internal", lb = NA, ub = NA)
  }

  # brms's class "Intercept" on an ordinal family names the THRESHOLDS:
  # its ordinal program declares them as the `Intercept` vector and puts
  # the default student_t there. frmtmb holds them in `tau_raw`, which
  # had no class spelling at all, so the row used to reach the resolver
  # and stop with a bare "Prior target not found".
  ordinal_threshold_entry <- function(s) {
    raw <- frame[["par_template"]][["tau_raw"]] %||% numeric(0)
    if (!length(raw) || length(fit$spec$responses) != 1L) return(NULL)
    rspec <- fit$spec$responses[[1L]]
    if (!identical(rspec$family[["type"]], "ordinal")) return(NULL)
    if (nzchar(s$coef) || nzchar(s$dpar) || nzchar(s$nlpar %||% "")) {
      return(NULL)
    }
    if (nzchar(s$resp %||% "") &&
          !identical(s$resp, rspec$resp_name)) {
      return(NULL)
    }
    # cumulative() and sratio() hold (tau_1, log increments), which is
    # the same map Stan's `ordered` type applies, so the density on the
    # thresholds carries that map's log-Jacobian. cratio() and acat()
    # hold the thresholds themselves and brms declares them unordered,
    # so neither side has a Jacobian there
    ordered <- rspec$family[["family"]] %in% c("cumulative", "sratio")
    list(comp = "tau_raw", idx = seq_along(raw), dist = s$dist,
         scale = if (ordered) "ordthres" else "internal",
         link = NULL, offset = ordinal_center_offset(frame, rspec),
         lb = s$lb, ub = s$ub)
  }

  for (s in prior_specificity_order(pl)) {
    # before the shape gate: in a multivariate model that gate would
    # judge the spelling against whichever response it met first
    no_resp <- resp_missing_refusal(fit$spec, frame, s)
    if (!is.null(no_resp)) frm_stop(no_resp, call. = FALSE)
    if (length(fit$spec$responses) > 1L && nzchar(s$resp %||% "") &&
          !s$resp %in% names(fit$spec$responses)) {
      # said first, so that no later refusal describes a response that
      # does not exist
      frm_stop("Prior target not found (", spec_target(s), "): resp = \"",
               s$resp, "\" names no response of this model. It has ",
               paste(names(fit$spec$responses), collapse = ", "),
               call. = FALSE)
    }
    if (length(fit$spec$responses) == 1L && nzchar(s$resp %||% "")) {
      # brms keys a univariate model's rows by no response, so it refuses
      # resp = "y" there ("b_y"); frmtmb applied it (user decision,
      # 2026-09-24)
      frm_stop("resp = \"", s$resp, "\" (", spec_target(s), "): resp ",
               "applies only to a multivariate model, and this model has ",
               "one response. Drop resp, as brms requires", call. = FALSE)
    }
    bad_shape <- dpar_shape_refusal(fit, s)
    if (!is.null(bad_shape)) frm_stop(bad_shape, call. = FALSE)
    ord_th <- if (s$class == "Intercept") ordinal_threshold_entry(s)
    if (!is.null(ord_th)) {
      if (!is.null(s$dist)) {
        claim("tau_raw", ord_th$idx)
        assigned[[nm_of("tau_raw", ord_th$idx)]] <- ord_th
      }
      if (!is.na(s$lb) || !is.na(s$ub)) {
        frm_stop("class = \"Intercept\" on an ordinal family addresses the ",
                 "whole threshold vector, so lb/ub would box every ",
                 "threshold with one number. Bound one at a time with ",
                 "class = \"theta\", or write the prior through ",
                 "prior = list(tau_raw = ) on the internal scale",
                 call. = FALSE)
      }
    } else if (s$class %in% c("b", "Intercept")) {
      for (tg in target_coefs(s)) {
        key <- nm_of(tg$comp, tg$idx)
        pm <- coef_placement(s, tg)
        if (!is.null(s$dist)) {
          assigned[[key]] <- list(comp = tg$comp, idx = tg$idx,
                                  dist = s$dist, scale = pm$scale,
                                  link = pm$link, offset = pm$offset,
                                  lb = s$lb, ub = s$ub)
        } else if (!is.null(assigned[[key]])) {
          assigned[[key]] <- entry_bounds(assigned[[key]], s)
        }
        # a bound belongs to the quantity the DENSITY is about, so a
        # `natural` placement has to carry it back through the link
        # before it can box an internal parameter. brms writes lb = 0 on
        # every dispersion default, and on a log-linked sigma that is
        # log(0) = -Inf, not a floor of 1 (see R3 of the punch re-check)
        if (!is.na(s$lb)) lower[tg$name] <- internal_bound(s$lb, pm, "lb")
        if (!is.na(s$ub)) upper[tg$name] <- internal_bound(s$ub, pm, "ub")
      }
    } else if (s$class == "sd") {
      # brms's rule: a row exists only for a prefix some block has, and
      # it reaches the standard deviations of that prefix (sd_spec_reach)
      exact <- any(vapply(frame[["re_blocks"]], function(bk) {
        sd_spec_reach(fit$spec, frame, bk, s)$exact
      }, TRUE))
      hit <- FALSE
      for (bk in if (exact) frame[["re_blocks"]]) {
        sd_i <- covstruct_registry[[bk[["covstruct"]]]]$sd_idx(bk[["dim"]])
        sd_i <- sd_i[sd_spec_reach(fit$spec, frame, bk, s)$reach]
        for (k in sd_i) {
          hit <- TRUE
          i <- bk[["theta_idx"]][k]
          key <- nm_of("theta", i)
          if (!is.null(s$dist)) {
            assigned[[key]] <- list(comp = "theta", idx = i,
                                    dist = s$dist, scale = "sd",
                                    lb = s$lb, ub = s$ub)
          } else if (!is.null(assigned[[key]])) {
            assigned[[key]] <- entry_bounds(assigned[[key]], s)
          }
          nm_theta <- paste0("theta_", i)
          if (!is.na(s$lb)) {
            lower[nm_theta] <- if (s$lb > 0) log(s$lb) else -Inf
          }
          if (!is.na(s$ub)) upper[nm_theta] <- log(s$ub)
        }
      }
      if (!hit) {
        frm_stop("No random-effect SDs match ", spec_target(s), ". ",
                 sd_prefix_hint(fit$spec, frame), call. = FALSE)
      }
    } else if (s$class == "cor") {
      hit <- FALSE
      refused <- character(0)
      for (bk in frame[["re_blocks"]]) {
        if (!block_addressed(fit$spec, frame, bk, s)) next
        cs <- bk[["covstruct"]]
        if (cs %in% names(lkj_refusals)) {
          refused <- c(refused, paste0(bk[["term_label"]], " [", cs, "]: ",
                                       unname(lkj_refusals[[cs]])))
          next
        }
        spec <- block_cor_spec(bk)
        if (is.null(spec)) next
        hit <- TRUE
        idx <- bk[["theta_idx"]][spec$idx]
        claim("theta", idx)
        assigned[[nm_of("theta", idx)]] <-
          list(comp = "theta", idx = idx,
               dist = lkj_dist(s$dist$eta, spec), scale = "internal",
               lb = NA, ub = NA)
      }
      if (!hit) {
        # a refused structure is named with its reason; otherwise the
        # model simply has no correlation to prior, and the message says
        # what it does have
        have <- unique(vapply(frame[["re_blocks"]], function(bk) {
          paste0(bk[["term_label"]], " [", bk[["covstruct"]], "]")
        }, ""))
        frm_stop("No random-effect correlations match ", spec_target(s),
                 ". ",
                 if (length(refused)) {
                   paste0("No LKJ density fits ",
                          paste(refused, collapse = "; "))
                 } else if (length(have)) {
                   paste0("These blocks have no correlation parameter: ",
                          paste(have, collapse = ", "))
                 } else {
                   "This model has no random-effect blocks"
                 }, call. = FALSE)
      }
    } else if (s$class == "theta") {
      # The raw internal escape hatch. All three covariance components
      # share one internal naming scheme, and their names are distinct,
      # so `coef` names the COMPONENT as well as the position:
      # "thetaac_1" reaches a residual-correlation parameter that no
      # natural-scale class can express a box constraint on.
      tgt <- theta_coef_target(frame, s$coef)
      cmp <- tgt$comp
      nms <- par_template_names(frame[["par_template"]][[cmp]] %||% numeric(0),
                                cmp)
      for (i in tgt$idx) {
        key <- nm_of(cmp, i)
        if (!is.null(s$dist)) {
          claim(cmp, i)
          assigned[[key]] <-
            list(comp = cmp, idx = i, dist = s$dist,
                 scale = "internal", lb = s$lb, ub = s$ub)
        } else if (!is.null(assigned[[key]])) {
          assigned[[key]] <- entry_bounds(assigned[[key]], s)
        }
        if (!is.na(s$lb)) lower[nms[i]] <- s$lb
        if (!is.na(s$ub)) upper[nms[i]] <- s$ub
      }
    } else if (s$class %in% names(autocor_prior_classes)) {
      resolve_ac_class(s)
    } else if (s$class == "rescor") {
      resolve_rescor(s)
    }
  }
  list(entries = unname(assigned), lower = lower, upper = upper)
}

#' brms's centering offset for one linear predictor's intercept, or
#' `NULL` where brms would not center.
#'
#' brms constrains the intercept at the MEAN of the predictors: its Stan
#' program declares `Intercept` against a centered design and recovers
#' the reported one as `b_Intercept = Intercept - dot_product(means_X,
#' b)`. frmtmb parameterizes by the raw intercept, so the quantity
#' brms's prior is about is `b0 + means_X'b`, and that is what this
#' offset supplies. The map between the two parameterizations is unit
#' triangular, so it carries no Jacobian, which is why the measured S7
#' residual is a pure density difference.
#'
#' Only the PARAMETRIC columns are centered. That is brms's own rule,
#' read off its generated code rather than assumed: `Xs` (a smooth's
#' unpenalized part) and `Xmo` (a mo() term) sit outside `Xc` and
#' outside `means_X`, and frmtmb appends those columns after
#' `n_param_cols`. A sub-formula with no intercept column (an ordinal
#' one, where thresholds replace it) or with no other parametric column
#' has nothing to center.
#'
#' @noRd
lp_center_offset <- function(frame, lp) {
  X <- lp[["X"]]
  np <- lp[["n_param_cols"]] %||% 0L
  if (is.null(X) || !nrow(X) || np < 2L) return(NULL)
  cn <- colnames(X)[seq_len(np)]
  if (!"(Intercept)" %in% cn) return(NULL)
  k <- which(cn != "(Intercept)")
  if (!length(k)) return(NULL)
  m <- as.numeric(Matrix::colMeans(X[, k, drop = FALSE]))
  keep <- which(m != 0)
  # a mean-zero design is brms's own arithmetic with nothing in it, and
  # dropping those columns keeps them off the tape
  if (!length(keep)) return(NULL)
  list(comp = lp[["par"]], idx = lp[["idx"]][k[keep]], w = m[keep])
}

#' Does this model hold ordinal thresholds a prior can address?
#'
#' @noRd
has_ordinal_thresholds <- function(fit) {
  raw <- fit$frame[["par_template"]][["tau_raw"]] %||% numeric(0)
  length(raw) > 0L && length(fit$spec$responses) == 1L &&
    identical(fit$spec$responses[[1L]]$family[["type"]], "ordinal")
}

#' brms's centering offset for an ordinal threshold vector, or `NULL`.
#'
#' The ordinal program centers every column of its design (there is no
#' intercept column to leave out) and recovers the reported thresholds
#' as `b_Intercept = Intercept + dot_product(means_X, b)`, so the
#' quantity its prior is about is `tau - means_X'b`. The sign is the
#' opposite of the ordinary intercept's, because an ordinal linear
#' predictor enters the density as `tau - eta`.
#'
#' @noRd
ordinal_center_offset <- function(frame, rspec) {
  for (lp in frame[["linpreds"]]) {
    if (!identical(lp[["resp"]], rspec$resp_name)) next
    if (!lp[["dpar"]] %in% rspec$primary_dpars) next
    X <- lp[["X"]]
    np <- lp[["n_param_cols"]] %||% 0L
    if (is.null(X) || !nrow(X) || np < 1L) return(NULL)
    k <- seq_len(np)
    m <- as.numeric(Matrix::colMeans(X[, k, drop = FALSE]))
    keep <- which(m != 0)
    if (!length(keep)) return(NULL)
    return(list(comp = lp[["par"]], idx = lp[["idx"]][k[keep]],
                w = -m[keep]))
  }
  NULL
}

#' Where the density of one class b / Intercept specification sits:
#' `scale`, the `link` a natural placement needs, and the centering
#' `offset`.
#'
#' The log link is the case class `"sd"` already implements exactly
#' (`exp()` with the coefficient as its log-Jacobian), so it reuses that
#' path rather than a parallel one; the identity link makes natural and
#' link scale the same quantity; every other link goes through the link
#' object's own `linkinv` and `mu_eta`.
#'
#' @noRd
coef_placement <- function(s, tg) {
  ctr <- tg$center
  if (!isTRUE(s$natural)) {
    return(list(scale = "internal", link = NULL, offset = ctr))
  }
  nm <- tg$link$name %||% "identity"
  if (identical(nm, "log")) {
    return(list(scale = "sd", link = NULL, offset = ctr))
  }
  if (identical(nm, "identity")) {
    return(list(scale = "internal", link = NULL, offset = ctr))
  }
  list(scale = "natural", link = tg$link, offset = ctr)
}

#' Does this distributional parameter have a PREDICTOR, in the sense
#' that decides which of the two brms spellings applies to it?
#'
#' brms's rule is "did the model write a formula for it", not "does the
#' design have more than an intercept": `bf(y ~ x, sigma ~ 1)` declares
#' `Intercept_sigma` and retires `sigma`, exactly as `sigma ~ x` does.
#' Its design is a single intercept column, identical to the one an
#' unwritten `sigma` gets, so the frame cannot tell them apart. The
#' spec can: `plain_dpar()` (R/parse.R) builds an unwritten dpar out of
#' four slots, where a written formula goes through `parse_linpred()`
#' and carries its special-term slots as well.
#'
#' The structural tests come first so that a dpar which is predicted on
#' any reading stays predicted even if that record ever changes shape.
#'
#' @noRd
dpar_has_predictor <- function(spec, lp) {
  X <- lp[["X"]]
  if (!is.null(X) && ncol(X) > 1L) return(TRUE)
  if (!is.null(lp[["Z"]])) return(TRUE)
  if (length(lp[["smooths"]] %||% list())) return(TRUE)
  if (length(lp[["gps"]] %||% list())) return(TRUE)
  if (length(lp[["mo"]] %||% list())) return(TRUE)
  dp <- spec$responses[[lp[["resp"]]]]$dpars[[lp[["dpar"]]]]
  "acterms" %in% names(dp %||% list())
}

#' Why a dpar-addressed specification does not fit this model's shape,
#' or `NULL`.
#'
#' The two brms spellings for a distributional parameter are mutually
#' exclusive BY MODEL SHAPE, and frmtmb now says the same:
#' `class = "sigma"` is a density on sigma itself and needs a sigma
#' with no predictor; `class = "Intercept", dpar = "sigma"` is a
#' density on the link-scale intercept of sigma's linear predictor and
#' needs a sigma that has one. brms refuses each on the other's model
#' ("The following priors do not correspond to any model parameter"),
#' and refusing by name is what keeps a ported script from silently
#' meaning something else.
#'
#' @noRd
dpar_shape_refusal <- function(fit, s) {
  dp <- s$dpar %||% ""
  if (!nzchar(dp) || !s$class %in% c("b", "Intercept")) return(NULL)
  if (nzchar(s$nlpar %||% "")) return(NULL)
  frame <- fit$frame
  lps <- Filter(function(lp) {
    identical(lp[["dpar"]], dp) &&
      (!nzchar(s$resp %||% "") || identical(lp[["resp"]], s$resp))
  }, frame[["linpreds"]])
  nat <- isTRUE(s$natural)
  # a nonlinear parameter has a linear predictor like a dpar, and would
  # otherwise be told it "has a predictor of its own", which is true and
  # not the sentence a user who wrote class = "ult" needs
  if (nat && dp %in% model_nlpars(fit$spec)) {
    return(paste0("class = \"", dp, "\" names a NONLINEAR parameter of ",
                  "this model, whose coefficients are class = \"b\" ",
                  "with nlpar = \"", dp, "\", the intercept among them, ",
                  "as brms lists them. A class of its own belongs to a ",
                  "distributional parameter"))
  }
  spelling <- if (nat) {
    paste0("class = \"", dp, "\"")
  } else {
    paste0("class = \"", s$class, "\", dpar = \"", dp, "\"")
  }
  if (!length(lps)) {
    have <- unique(vapply(frame[["linpreds"]], function(lp) {
      lp[["dpar"]]
    }, ""))
    return(paste0(spelling, " names no distributional parameter of ",
                  "this model. It has ", paste(have, collapse = ", "),
                  ". frmtmb's own classes are ",
                  paste(frmtmb_prior_classes, collapse = ", ")))
  }
  for (lp in lps) {
    if (!is.null(lp[["constant"]])) {
      return(paste0(spelling, ": ", dp, " is fixed at ",
                    format(lp[["constant"]]), " in this model, so it ",
                    "holds no parameter a prior can reach. Drop the ",
                    "constant from bf() to estimate it"))
    }
    if (!is.null(lp[["nl_body"]])) {
      return(paste0(spelling, ": ", dp, " is computed by an nlf() body ",
                    "here, so it has no coefficient of its own. Prior ",
                    "the parameters the body is written from instead"))
    }
    # a location parameter always has the main formula, even when that
    # formula is an intercept only, so it takes the Intercept spelling
    pred <- lp[["dpar"]] %in%
      fit$spec$responses[[lp[["resp"]]]]$primary_dpars ||
      dpar_has_predictor(fit$spec, lp)
    link <- lp[["link"]]$name %||% "identity"
    # the third spelling brms decides by model shape. `dpar ~ 1` gives
    # the parameter a predictor, so it takes the Intercept spelling,
    # but it gives it no population-level SLOPES, so class "b" there
    # addresses an empty set. An empty match is silent: the penalty was
    # bit-identical to no prior at all, which is the one failure this
    # gate exists to prevent
    if (identical(s$class, "b") && pred &&
          !length(setdiff(colnames(lp[["X"]]) %||% character(0),
                          "(Intercept)"))) {
      return(paste0("class = \"b\", dpar = \"", dp, "\" addresses the ",
                    "population-level slopes of ", dp,
                    "'s linear predictor, and this model's ", dp,
                    " predictor is an intercept only, so it has none. ",
                    "The prior brms takes there is ",
                    "class = \"Intercept\", dpar = \"", dp,
                    "\", a density on the ", link, "-scale intercept; ",
                    "give ", dp, " a predictor to have slopes to prior"))
    }
    if (nat && pred) {
      return(paste0("class = \"", dp, "\" is a density on ", dp,
                    " itself, which brms accepts only where ", dp,
                    " has no predictor of its own. This model gives ",
                    dp, " one, so the prior brms takes there is ",
                    "class = \"Intercept\", dpar = \"", dp,
                    "\", a density on the ", link, "-scale intercept ",
                    "of that predictor"))
    }
    if (!nat && !pred) {
      if (identical(s$class, "Intercept")) {
        return(paste0("class = \"Intercept\", dpar = \"", dp, "\" is a ",
                      "density on the ", link, "-scale intercept of ",
                      dp, "'s linear predictor, which brms accepts ",
                      "only where ", dp, " has one. This model writes ",
                      "no ", dp, " formula, so the prior brms takes ",
                      "there is class = \"", dp, "\", a density on ",
                      dp, " itself. Write ", dp, " ~ 1 in bf() to have ",
                      "a link-scale intercept to prior"))
      }
      return(paste0("class = \"b\", dpar = \"", dp, "\" addresses ", dp,
                    "'s slopes, and this model gives ", dp,
                    " no predictor, so it has none. Write class = \"",
                    dp, "\" for a density on ", dp,
                    " itself, or give ", dp, " a formula in bf()"))
    }
  }
  NULL
}

#' A user-facing bound carried onto the internal parameter the box
#' actually constrains.
#'
#' @noRd
internal_bound <- function(v, pm, which) {
  if (identical(pm$scale, "internal")) return(v)
  lf <- if (identical(pm$scale, "sd")) log else pm$link$linkfun
  b <- suppressWarnings(as.numeric(lf(v)))
  if (is.finite(b)) return(b)
  # the bound lies outside the parameter's own support. Below it a lower
  # bound constrains nothing, which is exactly brms's lb = 0 on a
  # log-linked dispersion parameter; an upper bound there would be an
  # empty box and is a question about intent rather than a number to
  # invent
  if (identical(which, "lb") && (is.nan(b) || b == -Inf)) return(-Inf)
  # isTRUE(), because a bound outside the support maps to NaN and
  # `NaN == Inf` is NA, which would make the `if` itself the error
  # instead of the sentence below
  if (identical(which, "ub") && isTRUE(b == Inf)) return(Inf)
  frm_stop(which, " = ", v, " is outside the support of the parameter this ",
           "prior is about, so it describes an empty box", call. = FALSE)
}

#' Log density of one prior entry value (AD-safe), with the change of
#' variables and the centering offset the entry asks for.
#'
#' @noRd
prior_logdens <- function(x, dist, scale, link = NULL, offset = 0) {
  jac <- 0
  if (identical(scale, "ordthres")) {
    # (tau_1, log increments) -> the thresholds, which is the map Stan's
    # `ordered` type applies, so its log-Jacobian is the sum of the
    # increments. The centering shift comes after the map, because the
    # quantity brms priors is the CENTERED threshold vector
    "[<-" <- RTMB::ADoverload("[<-")
    K1 <- length(x)
    tau <- rep(x[1], K1)
    if (K1 > 1L) {
      for (k in 2:K1) tau[k] <- tau[k - 1] + exp(x[k])
      jac <- sum(x[-1])
    }
    return(sum(prior_base_logdens(tau + offset, dist)) + jac)
  }
  # brms's Intercept prior is about the intercept at the predictor
  # MEANS; the tape carries the intercept at zero. The map between them
  # is unit triangular, so the shift enters the density and nothing
  # enters the Jacobian
  if (length(offset) && !identical(offset, 0)) x <- x + offset
  if (identical(scale, "sd")) {
    jac <- x          # theta = log sd; add the Jacobian
    x <- exp(x)
  } else if (identical(scale, "natural")) {
    # the same change of variables for a link that is not the log: the
    # density belongs to the dpar, and the tape carries its link-scale
    # coefficient
    jac <- log(abs(link$mu_eta(x)))
    x <- link$linkinv(x)
  }
  # a density written about a TRANSFORMED parameter (an AR coefficient,
  # a cosy correlation): evaluate it at the natural value and add the
  # map's log Jacobian, the same change of variables "sd" performs
  # inline above. One number for the whole segment, because the map is
  # not elementwise.
  if (identical(dist$kind, "trans")) {
    return(sum(prior_base_logdens(ac_trans_value(x, dist$trans),
                                  dist$inner)) +
             ac_trans_logjac(x, dist$trans))
  }
  prior_base_logdens(x, dist) + jac
}

#' The density itself, with no change of variables.
#'
#' @noRd
prior_base_logdens <- function(x, dist) {
  switch(dist$kind,
    normal = RTMB::dnorm(x, dist$location, dist$scale, log = TRUE),
    t = RTMB::dt((x - dist$location) / dist$scale, df = dist$df,
                 log = TRUE) - log(dist$scale),
    exponential = log(dist$rate) - dist$rate * x,
    # written through logspace_add rather than log1p(exp(-z)), which
    # overflows in the lower tail exactly where brms's logistic(0, 1) on
    # a mixture proportion is doing its work
    logistic = {
      z <- (x - dist$location) / dist$scale
      -z - log(dist$scale) - 2 * RTMB::logspace_add(0 * z, -z)
    },
    gamma = dist$shape * log(dist$rate) - lgamma(dist$shape) +
      (dist$shape - 1) * log(x) - dist$rate * x,
    # brms's inv_gamma(shape, scale) is Stan's, so the second argument
    # is the SCALE of the inverse gamma, which is the rate of the gamma
    # on 1/x. Written out rather than through dgamma(1/x) so the tape
    # carries one expression instead of a reciprocal and a correction
    inv_gamma = dist$shape * log(dist$scale) - lgamma(dist$shape) -
      (dist$shape + 1) * log(x) - dist$scale / x,
    # log1p rather than log(1 - x): a zi or hu near one is exactly where
    # brms's beta(1, 1) is doing its work, and that is where log(1 - x)
    # loses its last digits
    beta = (dist$shape1 - 1) * log(x) +
      (dist$shape2 - 1) * log1p(-x) -
      (lgamma(dist$shape1) + lgamma(dist$shape2) -
         lgamma(dist$shape1 + dist$shape2)),
    # a JOINT density over a whole correlation, so `x` is the block's
    # correlation segment and the value is one number, not one per
    # element (see lkj_logdens)
    lkj = lkj_logdens(x, dist)
  )
}

#' Prior objects, addressed by internal parameter name
#'
#' The named-list prior spelling, as opposed to [set_prior()]'s classes.
#' Priors written this way apply on the INTERNAL parameter scale:
#' coefficients are on their link scale, and covariance parameters
#' (`theta_*`) are the unconstrained parameterization (log-SDs,
#' scaled-Cholesky terms), so `prior_normal(0, 1)` on `theta_1` is a
#' lognormal prior on that standard deviation.
#'
#' [frm()] takes them as a MAP penalty, and so does
#' `frmtmb.sample::frm_sample()`, where they take over exactly the
#' parameters they name and leave the rest of the prior stack in place.
#' [par_template()] and [get_prior()] name the addressable slots.
#'
#' @param location,scale,df Prior parameters. `scale` is also the
#'   second argument of `prior_inv_gamma()`, brms's
#'   `inv_gamma(shape, scale)`.
#' @return A `frmtmb_prior` object.
#' @seealso [set_prior()] for the class-based spelling, which is the one
#'   most models want.
#' @examples
#' # the objects themselves are cheap descriptions
#' prior_normal(0, 2)
#' prior_t(df = 3, location = 0, scale = 1)
#'
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#'
#' # names are internal parameter names, or whole components. theta_1
#' # is a log-SD, so a normal there is a lognormal on the SD.
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'            prior = list(beta = prior_normal(0, 5),
#'                         theta_1 = prior_t(3, 0, 1)))
#' prior_summary(fit)
#' @name frmtmb-priors
NULL

#' @rdname frmtmb-priors
#' @export
prior_normal <- function(location = 0, scale = 1) {
  # a length-2 location built a prior that recycled against the whole
  # parameter block it was attached to, and `stopifnot(scale > 0)`
  # reported a negative scale as "scale > 0 is not TRUE", which names
  # the test rather than the argument
  check_number(location, "location")
  check_positive(scale, "scale")
  structure(list(kind = "normal", location = location, scale = scale),
            class = "frmtmb_prior")
}

#' @rdname frmtmb-priors
#' @export
prior_t <- function(df = 3, location = 0, scale = 1) {
  check_positive(df, "df")
  check_number(location, "location")
  check_positive(scale, "scale")
  structure(list(kind = "t", df = df, location = location, scale = scale),
            class = "frmtmb_prior")
}

#' @rdname frmtmb-priors
#' @export
prior_logistic <- function(location = 0, scale = 1) {
  check_number(location, "location")
  check_positive(scale, "scale")
  structure(list(kind = "logistic", location = location, scale = scale),
            class = "frmtmb_prior")
}

#' @param shape,rate Gamma prior parameters, brms's spelling
#'   (`gamma(shape, rate)`), both positive. The density has support on
#'   the positive line, so it belongs to a quantity that lives there: a
#'   brms `shape`, `phi`, `nu` or `kappa` row, which arrives on the
#'   parameter's own scale.
#' @rdname frmtmb-priors
#' @export
prior_gamma <- function(shape = 1, rate = 1) {
  check_positive(shape, "shape")
  check_positive(rate, "rate")
  structure(list(kind = "gamma", shape = shape, rate = rate),
            class = "frmtmb_prior")
}

#' The inverse gamma has support on the positive line, so it belongs to
#' a quantity that lives there: brms's default on a `shape` parameter
#' and on a `gp()` length-scale. Its `scale` is brms's second
#' `inv_gamma(shape, scale)` argument, documented with the other scale
#' parameters above rather than a second time here, which R CMD check
#' reads as a duplicated argument.
#'
#' @rdname frmtmb-priors
#' @export
prior_inv_gamma <- function(shape = 1, scale = 1) {
  check_positive(shape, "shape")
  check_positive(scale, "scale")
  structure(list(kind = "inv_gamma", shape = shape, scale = scale),
            class = "frmtmb_prior")
}

#' @param shape1,shape2 Beta prior parameters, both positive. The
#'   density has support on `(0, 1)`, so it belongs to a quantity that
#'   lives there: brms's default on a zero-inflation `zi` or a hurdle
#'   `hu`, which reach it through their logit link.
#' @rdname frmtmb-priors
#' @export
prior_beta <- function(shape1 = 1, shape2 = 1) {
  check_positive(shape1, "shape1")
  check_positive(shape2, "shape2")
  structure(list(kind = "beta", shape1 = shape1, shape2 = shape2),
            class = "frmtmb_prior")
}

#' @param eta LKJ shape. `1` is uniform over correlation matrices,
#'   larger values concentrate toward the identity, and `0 < eta < 1`
#'   pushes toward the boundary.
#' @rdname frmtmb-priors
#' @export
prior_lkj <- function(eta = 1) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) ||
      eta <= 0) {
    frm_stop("prior_lkj(eta =) takes one finite positive number; eta = 1 ",
             "is uniform over correlation matrices", call. = FALSE)
  }
  structure(list(kind = "lkj", eta = eta), class = "frmtmb_prior")
}

# ------------------------------------------------- the LKJ prior ------
#
# THE DENSITY, ON FRMTMB'S OWN PARAMETERS. LKJ(eta) is the density
# `p(C) = c_d(eta)^-1 det(C)^(eta - 1)` over d x d correlation matrices;
# `eta = 1` is uniform over them. frmtmb never holds `C`. It holds `t`,
# the strictly-lower entries of a unit-diagonal lower-triangular matrix
# whose rows are then normalized (us_chol_L()), so the prior has to be
# carried onto `t` with the Jacobian of that map, and the result is what
# is implemented here. The derivation, in three steps:
#
# 1. ON THE CHOLESKY FACTOR. With `C = L L'`, `L` lower-triangular with
#    positive diagonal and unit-norm rows, the Jacobian of the map from
#    the strict lower triangle of `C` to that of `L` is
#    `prod_{i>=2} L_ii^(d - i)`, so
#      p(L) = c_d(eta)^-1 prod_{i>=2} L_ii^(d - i + 2 eta - 2),
#    using `det(C) = prod L_ii^2`. That is Stan's
#    `lkj_corr_cholesky_lpdf`, and it is the standard route.
#
# 2. FROM L TO t, ROW BY ROW. Row `i` of frmtmb's unnormalized matrix is
#    `(t_i, 1, 0, ...)` with `t_i` of length `m = i - 1`, so row `i` of
#    `L` is that vector over `sqrt(1 + ||t_i||^2)` and
#      L_ii = (1 + ||t_i||^2)^(-1/2).
#    The free entries of row `i` of `L` are `u_i = t_i / sqrt(1 + s)`,
#    `s = ||t_i||^2`, whose Jacobian matrix is
#    `(1 + s)^(-1/2) (I - t_i t_i' / (1 + s))`, with determinant
#      (1 + s)^(-m/2) * (1 - s/(1 + s)) = (1 + s)^(-(m + 2)/2)
#                                       = L_ii^(m + 2) = L_ii^(i + 1).
#
# 3. THE PRODUCT. Row `i` therefore contributes
#    `L_ii^(d - i + 2 eta - 2) * L_ii^(i + 1) = L_ii^(2 eta + d - 1)`:
#    the exponent is the SAME for every row, and with
#    `log L_ii = -log(1 + ||t_i||^2)/2` the whole log density is
#
#      log p(t) = -(eta + (d - 1)/2) * sum_i log(1 + ||t_i||^2)
#                 - sum_i log Z_i.
#
#    `Z_i` normalizes row `i` on its own, because the rows are
#    independent under LKJ: `p(u_i) ∝ (1 - ||u_i||^2)^(a_i - 1)` on the
#    unit ball of R^m with `a_i = eta + (d - i)/2`, and
#    `int (1 - ||u||^2)^(a-1) du = pi^(m/2) Gamma(a) / Gamma(a + m/2)`,
#    which gives `Z_i = pi^((i-1)/2) Gamma(eta + (d-i)/2) /
#    Gamma(eta + (d-1)/2)`. Their product IS the published `c_d(eta)`
#    (checked to 1e-15 against the LKJ 2009 closed form for d = 2..5,
#    tests/testthat/test-lkj.R).
#
# THE d = 2 CHECK, closed form: `rho = t / sqrt(1 + t^2)`, so
# `p(t) = (1 + t^2)^-(eta + 1/2) / (2^(2 eta - 1) B(eta, eta))`, which at
# `eta = 1` is `(1 + t^2)^(-3/2)` and is uniform on rho. Flat on `t`, by
# the same change of variables, is `(1 - rho^2)^(-3/2)`: improper, all
# its mass at |rho| = 1, which is what the LKJ prior replaces.
#
# THE ONE-PARAMETER STRUCTURES (cs, homcs, ar1, hetar1) hold a single
# bounded correlation instead of a whole matrix, so there is no matrix
# for the density above to be about. They take the d = 2 form,
# `p(rho) ∝ (1 - rho^2)^(eta - 1)`, which is the LKJ marginal, with each
# structure's own Jacobian: `rho = t/sqrt(1 + t^2)` for ar1 (identical to
# the d = 2 case above), and the scaled logistic onto `(-1/(d-1), 1)` for
# cs, whose normalizer picks up the mass LKJ puts below `-1/(d-1)`.

#' `sum_i log Z_i`: the log of the LKJ normalizing constant for
#' dimension `d`, assembled from the per-row constants of the derivation
#' above.
#'
#' @noRd
lkj_lognorm <- function(eta, d) {
  if (d < 2L) return(0)
  i <- seq_len(d - 1L) + 1L
  sum((i - 1) / 2 * log(pi) + lgamma(eta + (d - i) / 2) -
        lgamma(eta + (d - 1) / 2))
}

#' The positions of the correlation segment belonging to each ROW of the
#' Cholesky factor. `L[lower.tri(L)] <- t` fills column-major, so the
#' entries of one row are not contiguous.
#'
#' @noRd
lkj_rows <- function(d) {
  ii <- row(diag(d))[lower.tri(diag(d))]
  lapply(seq_len(d - 1L) + 1L, function(i) which(ii == i))
}

#' The internal prior object the objective evaluates: the user's `eta`
#' plus everything about the block's map that does not depend on the
#' parameters, computed once here so the taped density is arithmetic
#' only.
#'
#' @noRd
lkj_dist <- function(eta, spec) {
  d <- spec$d
  out <- list(kind = "lkj", eta = eta, map = spec$kind, d = d)
  if (identical(spec$kind, "chol")) {
    out$rows <- lkj_rows(d)
    out$pow <- eta + (d - 1) / 2
    out$lognorm <- lkj_lognorm(eta, d)
  } else if (identical(spec$kind, "ar1")) {
    out$pow <- eta + 0.5
    out$lognorm <- lkj_lognorm(eta, 2L)
  } else {
    a <- 1 / (d - 1)
    out$a <- a
    # the marginal restricted to (-a, 1), renormalized over that window
    out$lognorm <- lkj_lognorm(eta, 2L) +
      log(1 - stats::pbeta((1 - a) / 2, eta, eta))
  }
  structure(out, class = "frmtmb_prior")
}

#' The LKJ log density at one block's correlation parameters (AD-safe).
#' One number for the whole block, whatever its map.
#'
#' @noRd
lkj_logdens <- function(t, dist) {
  if (identical(dist$map, "chol")) {
    q <- 0
    for (r in dist$rows) q <- q + log(1 + sum(t[r] * t[r]))
    return(-dist$pow * q - dist$lognorm)
  }
  if (identical(dist$map, "ar1")) {
    return(-dist$pow * log(1 + t[1] * t[1]) - dist$lognorm)
  }
  a <- dist$a
  rho <- -a + (1 + a) / (1 + exp(-t[1]))
  # d rho / d t = (rho + a)(1 - rho)/(1 + a), written in rho so that the
  # logistic is computed once
  ld <- log(rho + a) + log(1 - rho) - log(1 + a) - dist$lognorm
  # eta = 1 has no density factor at all, and writing the term anyway
  # would evaluate 0 * log(0) = NaN where the logistic saturates. `eta`
  # is a constant, so this branch is resolved when the tape is built.
  if (dist$eta != 1) ld <- ld + (dist$eta - 1) * log(1 - rho * rho)
  ld
}

# --------------------------------------------------------------------
# Residual-correlation priors: the R-side autocorrelation classes and
# `rescor`. brms names these surfaces `ar`, `ma`, `cosy`, `cortime` and
# `rescor`, and frmtmb holds all of them in two internal vectors
# (`thetaac`, `thetar`) on unconstrained scales. A prior is therefore
# written on the NATURAL parameter and carried inward with the log
# Jacobian of the map, which is what class "sd" already does for a
# log standard deviation.
# --------------------------------------------------------------------

#' Which brms class names an R-side autocorrelation surface, and which
#' `struct` values can carry each one.
#'
#' @noRd
autocor_prior_classes <- list(
  ar      = c("ar", "arma"),
  ma      = c("ma", "arma"),
  cosy    = "cosy",
  cortime = "unstr"
)

#' Positions WITHIN one residual block's `thetaac` segment that a class
#' addresses, or `NULL` when the block's structure has no such
#' parameter. `arma` lays its AR coefficients out first.
#'
#' @noRd
autocor_class_idx <- function(ac, cls) {
  st <- ac[["struct"]]
  if (!st %in% autocor_prior_classes[[cls]]) return(NULL)
  switch(cls,
    ar = if (ac[["p"]]) seq_len(ac[["p"]]),
    ma = if (ac[["q"]]) {
      (if (identical(st, "arma")) ac[["p"]] else 0L) + seq_len(ac[["q"]])
    },
    cosy = 1L,
    cortime = seq_len(autocor_n_cor(ac[["d"]]))
  )
}

#' The map one autocorrelation class carries its prior through. `cosy`
#' is a single logistic onto the positive-definite window; the ARMA
#' coefficients are a partial-autocorrelation transform composed with
#' the Levinson-Durbin recursion; `cortime` is the row-normalized
#' Cholesky the LKJ density already knows, so it carries no `trans`.
#'
#' @noRd
autocor_trans <- function(ac, cls) {
  # brms's cov = FALSE coefficients are unconstrained and ARE the
  # internal parameters, so a prior on them needs no change of variables
  if (autocor_is_cond(ac)) return(list(map = "identity"))
  if (identical(cls, "cosy")) return(list(map = "cosy",
                                          a = 1 / (ac[["d"]] - 1)))
  if (identical(cls, "cortime")) return(NULL)
  list(map = "levinson")
}

#' A prior on a transformed parameter: the user's density, plus the map
#' that takes the internal vector to the parameter the density is
#' written about.
#'
#' @noRd
trans_dist <- function(inner, trans) {
  structure(list(kind = "trans", inner = inner, trans = trans),
            class = "frmtmb_prior")
}

#' Internal vector -> the natural parameters a transformed prior is
#' written about. AD-safe.
#'
#' @noRd
ac_trans_value <- function(th, tr) {
  if (identical(tr$map, "identity")) return(th)
  if (identical(tr$map, "cosy")) {
    return(-tr$a + (1 + tr$a) / (1 + exp(-th[1])))
  }
  autocor_levinson(autocor_pacf(th))
}

#' Log absolute determinant of that map's Jacobian, so the density the
#' user wrote on the natural scale stays a density once carried onto the
#' internal one.
#'
#' `th -> pacf` is elementwise, contributing `-3/2 log(1 + th^2)` each.
#' The Levinson step that extends order `k - 1` to `k` rewrites the
#' earlier coefficients as `(I - pac_k J) phi`, with `J` the exchange
#' matrix; `J` has eigenvalue `+1` with multiplicity `ceiling(m/2)` and
#' `-1` with multiplicity `floor(m/2)` at size `m`, which gives that
#' step's determinant in closed form. The new coordinate is `pac_k`
#' itself and contributes 1. `test-priors-autocor-classes.R` checks the
#' whole expression against a numeric Jacobian.
#'
#' @noRd
ac_trans_logjac <- function(th, tr) {
  if (identical(tr$map, "identity")) return(0)
  if (identical(tr$map, "cosy")) {
    s <- 1 / (1 + exp(-th[1]))
    return(log(1 + tr$a) + log(s) + log(1 - s))
  }
  pac <- autocor_pacf(th)
  lj <- -1.5 * sum(log(1 + th * th))
  n <- length(th)
  if (n > 1L) {
    for (k in 2:n) {
      m <- k - 1L
      lj <- lj + ceiling(m / 2) * log(1 - pac[k]) +
        floor(m / 2) * log(1 + pac[k])
    }
  }
  lj
}

#' A natural-scale bound as a bound on the internal parameter. Only a
#' MONOTONE scalar map can carry one: a higher-order ARMA coefficient is
#' a function of several internal parameters at once, so no box in
#' internal space is the box the user asked for.
#'
#' @noRd
ac_bound_theta <- function(v, tr, ac, cls, what) {
  if (identical(tr$map, "identity")) return(v)
  if (identical(tr$map, "cosy")) {
    a <- tr$a
    if (v <= -a || v >= 1) {
      frm_stop("class = \"cosy\" ", what, " = ", v, " is outside the window ",
               "a compound-symmetric correlation of ", ac[["d"]],
               " time points can occupy, (", format(-a), ", 1)",
               call. = FALSE)
    }
    return(stats::qlogis((v + a) / (1 + a)))
  }
  if (abs(v) >= 1) {
    frm_stop("class = \"", cls, "\" ", what, " = ", v,
             " is outside (-1, 1), which is the whole range a first-order ",
             cls, " coefficient can take", call. = FALSE)
  }
  v / sqrt(1 - v * v)
}

#' Accepts the legacy named list of prior objects OR a priorlist; returns
#' entries plus bounds.
#'
#' @noRd
resolve_prior_input <- function(fit, prior) {
  # the argument boundaries translate a brms prior object already; this
  # covers the internal callers that reach the resolver directly
  prior <- as_priorlist(prior)
  if (inherits(prior, "frmtmb_priorlist")) {
    return(resolve_priorlist(fit, prior))
  }
  legacy <- resolve_priors(fit, prior)
  entries <- list()
  for (e in legacy) {
    for (i in e$idx) {
      entries[[length(entries) + 1L]] <-
        list(comp = e$comp, idx = i, dist = e$prior, scale = "internal")
    }
  }
  list(entries = entries, lower = c(), upper = c())
}

# --- the fit route: prior entries and bounds onto the parameters ------
#
# Reached from frm() (MAP / regularized ML), par_template() and
# frm_simulate() without going through frm_sample(), so this stays in
# core when the sampling surface leaves.

#' Resolve a named prior list to per-component index/parameter vectors.
#' Names may be individual parameters (as in outer_par_names()) or whole
#' components ("beta", "betad", "theta", "thetar", "thetaac").
#'
#' @noRd
resolve_priors <- function(fit, prior) {
  if (!is.list(prior) || is.null(names(prior)) ||
        !all(nzchar(names(prior)))) {
    frm_stop("`prior` must be made by set_prior() or be a list with a ",
             "name for each element, the parameter it applies to, as in ",
             "list(x = prior_normal(0, 1)); got ", arg_desc(prior),
             if (is.list(prior)) " without names for all elements",
             call. = FALSE)
  }
  tpl <- fit$frame[["par_template"]]
  comp_names <- list()
  for (cp in setdiff(names(tpl), c("b", "miss"))) {
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    if (cp == "betad" && length(fit$frame[["betad_fixed_idx"]])) {
      v[fit$frame[["betad_fixed_idx"]]] <- NA   # mapped: no prior
    }
    comp_names[[cp]] <- v
  }
  entries <- list()
  add <- function(comp, idx, pr) {
    entries[[length(entries) + 1L]] <<- list(comp = comp, idx = idx,
                                             prior = pr)
  }
  for (nm in names(prior)) {
    pr <- prior[[nm]]
    if (!inherits(pr, "frmtmb_prior")) {
      frm_stop("prior[['", nm, "']] must be a prior object ",
               "(prior_normal(), prior_t())", call. = FALSE)
    }
    if (identical(pr$kind, "lkj")) {
      # this spelling addresses parameters one at a time; the LKJ
      # density is over a block's whole correlation and needs the
      # structure's map, which only the class spelling carries
      frm_stop("prior_lkj() addresses a block's whole correlation, so it ",
               "cannot be given by parameter name; write ",
               "set_prior(\"lkj(", format(pr$eta), ")\", class = \"cor\")",
               call. = FALSE)
    }
    if (nm %in% names(comp_names)) {
      idx <- which(!is.na(comp_names[[nm]]))
      add(nm, idx, pr)
      next
    }
    hit <- FALSE
    for (cp in names(comp_names)) {
      # both spellings: the template's own `(Intercept)` and the
      # parenthesis-free `Intercept` the draws, variables() and
      # hypothesis() all use. Priors are written against names the user
      # read off one of those surfaces
      i <- which(comp_names[[cp]] == nm |
                   par_name_bare(comp_names[[cp]]) == par_name_bare(nm))
      if (length(i)) {
        add(cp, i, pr)
        hit <- TRUE
        break
      }
    }
    if (!hit) {
      frm_stop("Unknown parameter in prior: '", nm, "'. Available: ",
               paste(par_name_bare(unlist(comp_names))[
                 !is.na(unlist(comp_names))], collapse = ", "),
               " or component names ",
               paste(names(comp_names), collapse = ", "), call. = FALSE)
    }
  }
  entries
}

#' The centering offset of one entry, evaluated on the tape: a single
#' number, `means_X'b` over the coefficients the offset names.
#'
#' @noRd
entry_offset <- function(e, pars) {
  o <- e$offset
  if (is.null(o)) return(0)
  sum(pars[[o$comp]][o$idx] * o$w)
}

#' AD-safe negative log prior over resolved per-parameter entries
#' (each: comp, idx, dist, scale, and optionally link/offset; see
#' prior_logdens).
#'
#' @noRd
neg_log_prior_fn <- function(entries) {
  function(pars) {
    nlp <- 0
    for (e in entries) {
      nlp <- nlp - sum(prior_logdens(pars[[e$comp]][e$idx], e$dist,
                                     e$scale, e$link,
                                     entry_offset(e, pars)))
    }
    nlp
  }
}

#' Named bound specs -> full-length vectors over the outer parameters.
#'
#' @noRd
resolve_bounds <- function(fit, lower, upper) {
  nm <- outer_par_names(fit)
  mk <- function(x, fill) {
    out <- rep(fill, length(nm))
    if (is.null(x)) return(out)
    if (is.null(names(x)) || any(names(x) == "")) {
      frm_stop("Bounds must be named numeric vectors over parameter ",
               "names, e.g. c(x = 0)", call. = FALSE)
    }
    # the paren-tolerant addressing of confint(parm =), so a name copied
    # out of a hypothesis() expression works here too, plus the bare
    # name of an intercept-only nonlinear parameter: bounds on an ODE
    # model are written against the parameters of the dynamics (la),
    # not against their design-matrix spelling (la_(Intercept))
    pos <- apply_nlpar_alias(fit, names(x), match_par_name(names(x), nm))
    if (anyNA(pos)) {
      frm_stop("Unknown parameter(s) in bounds: ",
               paste(names(x)[is.na(pos)], collapse = ", "), ". Available: ",
               paste(nm, collapse = ", "),
               " (parentheses may be dropped, and intercept-only nonlinear ",
               "parameters may be named bare)", call. = FALSE)
    }
    out[pos] <- as.numeric(x)
    out
  }
  list(lower = mk(lower, -Inf), upper = mk(upper, Inf))
}
