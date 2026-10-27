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
#' `gamma(shape, rate)`, `inv_gamma(shape, scale)`, `beta(a, b)`,
#' `uniform(lb, ub)`, `lognormal(mu, sd)`, `lkj(eta)`, `dirichlet(alpha)`;
#' an empty string sets bounds only.
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
#'   narrow with `group`, and to one coefficient's standard deviation
#'   with `group` and `coef` (`coef = "x"` for the `x` column of
#'   `(1 + x | g)`, as [default_prior()] lists it). See Group-level
#'   standard deviations. As in brms, `resp`, `dpar` and `nlpar` select
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
#' - `"lncor"`: the correlation matrix of the log ratios of a
#'   [logistic_normal()] response with three or more categories, as a
#'   whole, brms's class of the same name. `lkj(eta)` only.
#' - `"delta"`: the distance between neighboring thresholds of an
#'   ordinal family with `threshold = "equidistant"`, brms's class of
#'   the same name, on its own scale; `group` selects the vector of one
#'   level of `thres(gr = )`. See Ordinal thresholds.
#' - `"simo"`: the simplex of one `mo()` term, `dirichlet(alpha)` only,
#'   addressed by brms's coef: the term's label with `1` appended
#'   (`"mox1"`, `"mox:fb1"`, as [get_prior()] lists them), and by `dpar`
#'   for a term outside the location parameter. `alpha` is one number
#'   or one per weight. On a MAP fit the density carries the Jacobian of
#'   the coordinates the simplex is held in, so the mode under
#'   `dirichlet(a)` is the mode of `dirichlet(a + 1/2)` on the weights,
#'   pulled further toward the uniform simplex; see "Monotonic effects"
#'   in [frm()].
#' - `"sbhaz"`: the baseline-hazard simplex of a `cox()` response
#'   (frmtmb.survival), or the simplex a family declares under that name
#'   (`frmtmb_family(simplex_pars =)`); `dirichlet(alpha)` only, brms's
#'   class of the same name; no `coef`.
#' - `"theta"`: raw internal covariance parameters (escape hatch).
#'   `coef` names one by its internal name and spans all three
#'   covariance components: `"theta_2"` for a random-effect block,
#'   `"thetar_1"` for a residual correlation, and `"thetaac_1"` for a
#'   parameter of a residual term that has no class of its own. This is
#'   the one spelling that reaches a single parameter of a structure
#'   whose natural coefficients are not free of one another. A
#'   residual-correlation parameter that a class reaches (`"ar"`,
#'   `"ma"`, `"cosy"`, `"cortime"`, frmtmb.spatial's `"lagsar"`) has
#'   that class as its one slot, as in brms: its `thetaac_` row is not
#'   listed, and a prior on it is refused with the class to use.
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
#' @section Group-level standard deviations:
#' `set_prior(class = "sd", group = "g", coef = "x")` puts the density
#' on one standard deviation: that of the `x` column of the `g` block.
#' brms puts the same prior on `sd_g__x`. The other standard deviations
#' of the block keep the prior of the group or of the class. A `coef`
#' without a `group` is refused, because brms has no such row.
#'
#' The density is on the standard deviation itself, with the
#' log-Jacobian of the log map. frmtmb holds a correlated block as log
#' standard deviations and separate correlation parameters. It does not
#' hold a Cholesky factor. So each standard deviation is one parameter,
#' and the prior on it is exact. The prior does not change the
#' correlations. A class `"cor"` prior addresses those. In a Student-t
#' block (`gr(g, dist = "student")`), the parameter is the scale. That
#' is also what brms's `sd` is in that term.
#'
#' Some structures give their columns one shared standard deviation,
#' for example `homdiag()`, `ar1()` and `homcs()`. They have no standard
#' deviation per coefficient. A `coef` that names one of their columns
#' is refused by name. Drop the `coef` to put the prior on the shared
#' standard deviation. A block of another response, distributional
#' parameter or nonlinear parameter does not cause this refusal.
#'
#' The group labels are brms's. A nested `(1 | g1/g2)` has the groups
#' `g1` and `g1:g2`. A multi-membership `(1 | mm(m1, m2))` has the group
#' `mmm1m2`. [default_prior()] lists them so.
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
#' it as a coefficient instead: `class = "b", coef = "Intercept"`.
#' A formula written `0 + Intercept + x`, or `bf(center = FALSE)`, makes
#' that the intercept's only prior slot, as brms does: the intercept is
#' then class `"b"` and a class `"b"` prior without a coef reaches it
#' too, while class `"Intercept"` has no slot and is refused.
#'
#' @section Ordinal thresholds:
#' `cumulative()`, `sratio()`, `cratio()` and `acat()` have no intercept
#' column: the thresholds replace it. brms priors them as its
#' `Intercept` class and so does frmtmb, so
#' `set_prior("student_t(3, 0, 2.5)", class = "Intercept")` on an
#' ordinal model addresses the whole threshold vector. It addresses the
#' THRESHOLDS, at the mean of the predictors (at zero under
#' `bf(center = FALSE)`, as in brms), with the log-Jacobian of
#' the map from frmtmb's internal storage. `cumulative()` holds
#' `(tau_1, log increments)`, which is the same map Stan's `ordered`
#' type applies to brms's ordered thresholds. `sratio()`, `cratio()` and
#' `acat()` hold the thresholds themselves, which brms declares as an
#' unconstrained vector, so there is no Jacobian. `lb`/`ub` are refused
#' there, because one number cannot box a whole vector of thresholds.
#'
#' One threshold is addressed by its number, as brms lists it:
#' `set_prior("normal(0, 1)", class = "Intercept", coef = "2")` is a
#' density on the second threshold, and it replaces a class-wide
#' density there and only there. The vector stays one entry, so the
#' log-Jacobian of `cumulative()`'s map enters once. Under
#' `thres(gr = )` a threshold is named by its level and its number,
#' `group = "a", coef = "2"`. brms has no per-threshold parameter
#' under `threshold = "equidistant"`, and the row is refused there, as
#' in brms.
#'
#' `prior = list(tau_raw = prior_normal(0, 5))` reaches the same
#' parameters on the INTERNAL scale, one entry per threshold, which is
#' the escape hatch to use when the increments rather than the
#' thresholds are what a prior is about.
#'
#' The threshold structure changes what class `"Intercept"` addresses,
#' as it does in brms. Under `threshold = "equidistant"` it is the
#' first threshold alone, at the mean of the predictors, and class
#' `"delta"` is the distance between thresholds: on `delta` itself,
#' with the log-Jacobian of the log that `cumulative()` and
#' `hurdle_cumulative()` hold it on, since brms bounds it below by 0
#' there. Under `threshold = "sum_to_zero"` class `"Intercept"` is
#' refused: brms's prior there is on thresholds it declares before
#' centering them, whose common location only that prior places, and
#' frmtmb estimates the centered thresholds directly, so the parameter
#' does not exist.
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
#' parameter itself, one position at a time with `coef = "theta_2"` or
#' `"thetar_1"` (or `"thetaac_1"` where no class reaches that
#' parameter, see class `"theta"` above). Everything else is bounded on the
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
#' for, and brms's `cov = TRUE` form is of order one only. Under
#' `cov = FALSE` the coefficients are unconstrained and are the internal
#' parameters, so `lb`/`ub` apply at any order, as in brms. A bound at
#' an end of the parameter's own interval, such as brms's `lb = -1` on a
#' first-order `ar`, bounds nothing and is accepted, so a
#' [default_prior()] table can be passed back as `prior =`.
#'
#' `cosy` is bounded below at `-1/(d - 1)` for `d` time points, where a
#' compound-symmetric matrix stops being positive definite, and a bound
#' outside that window is refused rather than clamped. brms bounds
#' `cosy` on `[0, 1]` instead, so a negative estimate here has no brms
#' counterpart. `frmtmb.sample::frm_sample()` samples `cosy` on brms's
#' `[0, 1]`: its default prior there is brms's flat one, and a `cosy`
#' prior that sets no `lb` keeps `lb = 0`, as brms keeps the bound it
#' declares. `default_prior(route = "sample")` lists that interval, and
#' the default route lists frm()'s window.
#'
#' Under a first-order `cov = TRUE` `ar()`, `ma()` or `arma()`, frmtmb's
#' `sigma` (in `summary()`, `sigma()` and `VarCorr()`) is the MARGINAL
#' sd of a residual, and brms's `sigma` is the innovation sd, smaller by
#' `c = sqrt(1 - ar^2)`, `1 / sqrt(1 + ma^2)` or
#' `sqrt((1 - ar^2) / (1 + ma^2 + 2 ar ma))`: the same model, another
#' name. A density on class `"sigma"` is read as brms reads it, at
#' `sigma * c`, with the Jacobian of `c`, so `set_prior("student_t(3, 0,
#' 2.5)", class = "sigma")` is brms's prior on its sigma. A `lb` above 0
#' or a `ub` on `sigma` is refused there, because a box on the
#' innovation sd moves with `ar` and is no box on the marginal sd frmtmb
#' holds. Under `sigma ~ x` with such a term, a prior on the intercept of
#' `sigma` (`class = "Intercept", dpar = "sigma"`) is still read at
#' frmtmb's log marginal sd, where brms reads it at the log innovation
#' sd; the slopes (`class = "b", dpar = "sigma"`) are the same in both.
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
#' - `simo` and `sbhaz` keep their meaning too: a `dirichlet()` on a
#'   `mo()` simplex and on a `cox()` baseline.
#' - `theta`/`theta1`/`theta2`, `sds`, `sdgp`, `lscale`,
#'   `sdcar`, `car`, `meanme`, `sdme` and `corme` are refused by name,
#'   each saying where frmtmb keeps that quantity instead, or that it
#'   has no prior slot. A refusal is deliberate: translating
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
#'   `"cortime"` and `"rescor"`, the simplex classes `"simo"` and
#'   `"sbhaz"`, or a distributional parameter's own
#'   name (`"sigma"`, `"shape"`, ...). See A distributional parameter's
#'   own class.
#' @param coef Restrict to one coefficient (classes `"b"`/`"Intercept"`),
#'   to the standard deviation of one coefficient (class `"sd"`, with
#'   `group`), or name the simplex of class `"simo"`. A coefficient is
#'   named as brms names it and as [default_prior()] lists it: `logxp`
#'   for the column of `log(xp)`, `IxE2` for `I(x^2)`. The name of the
#'   design column, `log(xp)`, also works.
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
#'   a name frmtmb refuses on principle (`sds`, `car`, ...)
#'   errors here saying where that quantity lives, and any other name
#'   is read as a distributional parameter's and checked against the
#'   model's own parameters when the prior is resolved, which is where
#'   they can be listed. The two distributions that belong to some classes
#'   only, `lkj()` and `dirichlet()`, are checked against them in both
#'   directions.
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
#' set_prior("", class = "theta", coef = "theta_2", lb = -2, ub = 2)
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

#' The bounds a class has of its own, brms 2.23.0's: `lb = 0` on the
#' standard deviations and length scales, the bounds of `dpar_bounds()`
#' on a distributional parameter's own class, none on a coefficient, an
#' intercept or frmtmb's internal `theta`.
#'
#' @noRd
prior_class_bounds <- function(s) {
  if (isTRUE(s$natural)) {
    dc <- sub("[0-9]+$", "", s$dpar %||% "")
    # a logistic_normal() category's sigma<k> is brms's class sigma
    if (startsWith(dc, "sigma")) dc <- "sigma"
    return(switch(dc,
      sigma = , shape = , phi = , kappa = , beta = , bs = , disc = ,
      nu = c(lb = TRUE, ub = FALSE),
      zi = , hu = , zoi = , coi = , bias = , quantile = ,
      ndt = c(lb = TRUE, ub = TRUE),
      c(lb = FALSE, ub = FALSE)))
  }
  switch(s$class,
    sd = , sdgp = , lscale = , delta = c(lb = TRUE, ub = FALSE),
    cosy = c(lb = TRUE, ub = TRUE),
    c(lb = FALSE, ub = FALSE))
}

#' brms's `check_prior_content()` warnings: a density with a lower bound
#' (or both bounds) on a parameter that has no such bound of its own,
#' with no `lb` (`ub`) given. The parameter is then not constrained, and
#' the density is a wall a sampler meets rather than a support it is
#' transformed onto: brms warns, and so does this, in brms's words
#' (dev/sampleprior-brmsunif.R).
#'
#' It reads the WHOLE prior a call is given, as brms reads its combined
#' prior frame: a coefficient's or a group's row is bounded by the
#' bounds of the class row above it (same class, dpar, nlpar and resp,
#' no coef or group), which frmtmb applies to it too. So it runs where
#' the prior is used (check_prior_slots()), not in each set_prior()
#' call, and the defaults, which never pass through there, are not
#' checked: brms checks a user's rows only.
#'
#' @noRd
prior_bound_warning <- function(rows) {
  lb_kind <- c("lognormal", "exponential", "gamma", "inv_gamma")
  ulb_kind <- c("beta", "uniform")
  lb_msg <- ub_msg <- character(0)
  cls_key <- function(s) {
    paste(s$class, s$dpar %||% "", s$nlpar %||% "", s$resp %||% "",
          sep = "|")
  }
  top <- Filter(function(s) {
    !nzchar(s$coef %||% "") && !nzchar(s$group %||% "")
  }, rows)
  top_lb <- vapply(top, function(s) !is.na(s$lb), NA)
  top_ub <- vapply(top, function(s) !is.na(s$ub), NA)
  top_key <- vapply(top, cls_key, "")
  for (s in rows) {
    k <- s$dist$kind %||% ""
    if (!k %in% c(lb_kind, ulb_kind)) next
    own <- prior_class_bounds(s)
    up <- top_key == cls_key(s)
    if (any(top_lb[up])) own[["lb"]] <- TRUE
    if (any(top_ub[up])) own[["ub"]] <- TRUE
    sp <- spec_spelling(s)
    lab <- paste0(sp$class, if (nzchar(sp$dpar)) paste0("_", sp$dpar),
                  if (nzchar(s$coef %||% "")) paste0("_", s$coef), " ~ ",
                  s$prior %||% "")
    if (!own[["lb"]] && is.na(s$lb)) lb_msg <- c(lb_msg, lab)
    if (k %in% ulb_kind && !own[["ub"]] && is.na(s$ub)) {
      ub_msg <- c(ub_msg, lab)
    }
  }
  if (length(lb_msg)) {
    frm_warning("It appears as if you have specified a lower bounded ",
                "prior on a parameter that has no natural lower bound.\n",
                "If this is really what you want, please specify ",
                "argument 'lb' of 'set_prior' appropriately.\n",
                "Warning occurred for prior \n",
                paste(lb_msg, collapse = "\n"), call. = FALSE,
                class = "frmtmb_prior_unbounded")
  }
  if (length(ub_msg)) {
    frm_warning("It appears as if you have specified an upper bounded ",
                "prior on a parameter that has no natural upper bound.\n",
                "If this is really what you want, please specify ",
                "argument 'ub' of 'set_prior' appropriately.\n",
                "Warning occurred for prior \n",
                paste(ub_msg, collapse = "\n"), call. = FALSE,
                class = "frmtmb_prior_unbounded")
  }
  invisible(NULL)
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
  lkj_classes <- c("cor", "cortime", "rescor", "lncor")
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
    # an unstr() correlation has one prior slot, its class, as in brms
    # (autocor_class_coefs()), so class "theta" is no route for it
    frm_stop("class = \"", class, "\" takes no lb/ub: the bound would apply ",
             "to a whole correlation matrix. ",
             if (identical(class, "cortime")) {
               paste0("brms takes none either; an unstr() correlation has ",
                      "one prior slot, this class, and lkj() is its density")
             } else {
               paste0("Bound one parameter at a time with class = ",
                      "\"theta\", whose coef names an internal parameter ",
                      "(\"theta_1\", \"thetar_1\")")
             }, call. = FALSE)
  }
  # brms's rules for a simplex, with its words: the Dirichlet and the
  # simplex classes go together, a simo row names one simplex, and a
  # simplex takes no bound (dev/moparity-brms-simo.R). sbhaz is a cox()
  # baseline, one per response, so it takes no coef. A class that is
  # not frmtmb's own may name a simplex a family declares
  # (frmtmb_family(simplex_pars =)), which only the fit can tell, so
  # the same rules are applied to it there (prior_simplex_check()).
  # A dirichlet() on a class that is no simplex is refused there too,
  # for every class: `class = "sigma"` can be told from a declared
  # simplex only by the model, and refusing its other spelling,
  # `class = "Intercept", dpar = "sigma"`, here gave one prior two
  # refusal points
  is_dir <- identical(dist$kind, "dirichlet")
  simplex_classes <- c("simo", "sbhaz")
  if (class %in% simplex_classes) {
    if (!is_dir) {
      frm_stop("Currently 'dirichlet' is the only valid prior for simplex ",
               "parameters: class = \"", class, "\" takes dirichlet(alpha), ",
               "e.g. set_prior(\"dirichlet(c(2, 1, 1))\", class = \"simo\", ",
               "coef = \"mox1\")", call. = FALSE)
    }
    if (!is.na(lb) || !is.na(ub)) {
      frm_stop("class = \"", class, "\" takes no lb/ub: a simplex is ",
               "bounded by its own constraint", call. = FALSE)
    }
    if (identical(class, "sbhaz") && (nzchar(coef) || nzchar(group))) {
      frm_stop("class = \"sbhaz\" is the whole baseline simplex of a cox() ",
               "response (frmtmb.survival), so coef and group name nothing ",
               "it can narrow to", call. = FALSE)
    }
  }
  if (identical(class, "simo")) {
    if (!nzchar(coef)) {
      frm_stop("class = \"simo\" addresses one simplex, named by its coef ",
               "as get_prior() lists it (\"mox1\" for mo(x), \"mox:z1\" ",
               "for mo(x):z); brms refuses a row without one too",
               call. = FALSE)
    }
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
frmtmb_prior_classes <- c("b", "Intercept", "sd", "sdgp", "lscale", "cor",
                          "theta", "ar", "ma", "cosy", "cortime", "rescor",
                          "delta", "simo", "sbhaz", "lncor")

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
    # brms fills its corme row with lkj(1), which is flat over the
    # correlation matrices and so applies nothing; only a row somebody
    # edited is a prior frmtmb cannot honor
    if (identical(cls, "corme") &&
        identical(gsub("[[:space:]]", "", dist), "lkj(1)")) {
      next
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
#' Class `"cor"` reads no `coef`: a correlation is not per-coefficient
#' at all. A class `"sd"` row narrowed to one coefficient is honored
#' (`sd_spec_reach()`), as brms honors it.
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
  if (cls %in% c("sdgp", "lscale") && nzchar(group)) {
    return(paste0("class = \"", cls, "\" has no group: a gp() term is not ",
                  "a grouping term, and brms's rows for it carry none. ",
                  "`coef` names the term, as get_prior() lists it. "))
  }
  if (!nzchar(coef) || !identical(cls, "cor")) return(NULL)
  paste0("class = \"cor\" addresses a whole correlation ",
         "matrix, so coef = \"", coef, "\" names nothing it ",
         "can narrow to. Drop the coef; `group` selects the ",
         "block. ")
}

#' brms class names frmtmb spells with the same word and the same
#' meaning, so a translated row keeps its class untouched.
#'
#' @noRd
brms_direct_prior_classes <- c("b", "Intercept", "sd", "sdgp", "lscale",
                               "cor", "ar", "ma", "cosy", "cortime",
                               "rescor", "delta", "simo", "sbhaz", "lncor")

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
    meanme = paste0(
      "brms's \"meanme\" is the mean of a me() term's latent values. ",
      "frmtmb estimates it by maximum likelihood and has no prior slot ",
      "for it, so there is no class to carry the row into: drop the ",
      "row. A flat meanme row from get_prior() is dropped for you. "),
    sdme = paste0(
      "brms's \"sdme\" is the SD of a me() term's latent values. frmtmb ",
      "estimates it by maximum likelihood and has no prior slot for it, ",
      "so there is no class to carry the row into: drop the row. A flat ",
      "sdme row from get_prior() is dropped for you. "),
    corme = paste0(
      "brms's \"corme\" is the correlation of the latent values of ",
      "several me() terms. frmtmb estimates it by maximum likelihood and ",
      "has no prior slot for it: drop the row. brms's default lkj(1) is ",
      "flat and is dropped for you, and set_mecor(FALSE) fixes the ",
      "correlation at zero. "),
    sds = paste0(
      "brms's \"sds\" is the wiggliness standard deviation of a smooth. ",
      "frmtmb holds a smooth as a random-effect block, so its frmtmb ",
      "spelling is class = \"sd\" with group = the smooth's label, e.g. ",
      "group = \"s(x)\"; get_prior() lists the label this model has. "),
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
  # brms's Dirichlet takes a vector, dirichlet(c(2, 1, 1)), whose
  # parentheses the one-call pattern below cannot read
  if (grepl("^\\s*dirichlet\\s*[(]", prior)) {
    return(prior_dirichlet(parse_dirichlet_arg(prior)))
  }
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
    # brms's default on a non-decision time without a formula is
    # uniform(0, min_Y), and on an ordinal discrimination under the
    # identity link lognormal(0, 1)
    uniform = {
      if (!all(is.finite(pars)) || !(pars[1] < pars[2])) {
        frm_stop("Prior '", prior, "': uniform(lb, ub) needs finite ",
                 "bounds with lb < ub", call. = FALSE)
      }
      structure(list(kind = "uniform", lower = pars[1], upper = pars[2]),
                class = "frmtmb_prior")
    },
    lognormal = {
      if (!(pars[2] > 0)) {
        frm_stop("Prior '", prior, "': the scale of lognormal() must be ",
                 "positive, not ", pars[2], call. = FALSE)
      }
      structure(list(kind = "lognormal", location = pars[1],
                     scale = pars[2]), class = "frmtmb_prior")
    },
    frm_stop("Unsupported prior distribution '", kind,
             "' (supported: normal, student_t, cauchy, exponential, ",
             "logistic, gamma, inv_gamma, beta, uniform, lognormal, lkj, ",
             "dirichlet)", call. = FALSE)
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
  inv_gamma = c("alpha", "beta"), beta = c("alpha", "beta"),
  uniform = c("alpha", "beta"), lognormal = c("mu", "sigma"))

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
  if (identical(dist$kind, "dirichlet")) {
    a <- dist$alpha
    return(paste0("dirichlet(", if (length(a) == 1L) format(a) else {
      paste0("c(", paste(format(a), collapse = ", "), ")")
    }, ")"))
  }
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
#' `"sd"` and `"cor"` are targeted by `group` and `nlpar`, and class
#' `"sd"` lists a row per coefficient of each group, as brms does; the
#' residual-correlation classes (`"ar"`, `"ma"`, `"cosy"`, `"cortime"`)
#' by `resp`, with `lb` and `ub` where the parameter lives on an
#' interval (brms's -1 and 1 on a first-order `cov = TRUE` `ar` or `ma`;
#' for `cosy`, frm()'s window `-1/(d - 1)` to 1 on the default route and
#' brms's 0 to 1 on `route = "sample"`); class `"rescor"` has one row
#' and takes no `resp`; and class `"theta"` rows name the raw internal
#' covariance parameters (escape hatch, including correlations one at a
#' time, across all three covariance components), except a
#' residual-correlation parameter that one of those classes reaches,
#' which has that class as its one row, as in brms. Where the location
#' is several distributional parameters (a categorical or mixture
#' model), each one's `b`, `Intercept` and `sd` rows carry its `dpar`,
#' as in brms.
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
#' A covariance structure brms does not have, and whose columns share
#' one standard deviation (`homdiag()`, `ar1()`, ...), has no class
#' `"sd"` row per coefficient. Every other block has one per column,
#' as in brms (`sd_patient__Intercept`).
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
  # check_response = FALSE because brms's default_prior() never reads
  # the response: it answers for a binomial model written without
  # trials() and for a response the family refuses; frm() refuses both
  design <- prior_design(object, data, family, data2,
                         check_response = FALSE)
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
                         check_response = TRUE) {
  if (inherits(object, "frmtmb_fit")) {
    return(list(spec = object$spec, frame = object$frame))
  }
  bform <- expand_dot_bform(as_bform(object, family), data)
  bform <- resolve_deferred_families(bform, data)
  spec <- parse_spec(bform)
  frame <- assemble_frame(spec, data, data2 = data2,
                         check_response = check_response)
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
  # brms's validate_prior() never reads the response either, so it
  # answers the models default_prior() answers
  design <- prior_design(formula, data, family, data2,
                         check_response = FALSE)
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
  # both sides in brms's spelling: the table writes sx_1 where a prior
  # written against it before may say s(x).fx1
  bare <- function(v) par_name_bare(prior_coef_brms(v))
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
  tab <- prior_rows_sd_order(tab)
  rownames(tab) <- NULL
  tab
}

#' The rows of a filled table in brms's order.
#' brms's validate_prior() and prior_summary() order the rows by
#' `order(class, group, resp, dpar, nlpar, coef)`, not by default_prior()'s
#' `order(resp, dpar, ...)`: class by class, and within class `"sd"`
#' every class-wide row first, then group by group, the responses inside
#' each (dev/prilabel-rev-log/brmsfit.txt, dev/brmsnames-1024-order.R).
#' frmtmb's own class `"theta"`, which brms does not have, stays last.
#'
#' @noRd
prior_rows_sd_order <- function(tab) {
  if (nrow(tab) < 2L) return(tab)
  keep <- which(tab$class != "theta")
  s <- as.data.frame(tab)[keep, , drop = FALSE]
  idx <- c(keep[with(s, order(class, group, resp, dpar, nlpar, coef))],
           which(tab$class == "theta"))
  # assigned in place, so the table keeps its class and attributes
  tab[seq_len(nrow(tab)), ] <- tab[idx, , drop = FALSE]
  tab
}

#' brms's spelling of a coefficient in the prior table.
#'
#' A design column is named by `model.matrix()`, and brms renames it
#' (`rename()`, `brms_rename()` here): `(Intercept)` is `Intercept`,
#' `log(xp)` is `logxp`, `I(x^2)` is `IxE2` and `poly(x, 2)1` is
#' `polyx21`, the names `variables()` gives after `b_`. The unpenalized
#' column of a smooth, `s(x).fx1`, is brms's `sx_1`. The resolver takes
#' either spelling (`prior_coef_match()`), so a prior written against
#' the design column still lands.
#'
#' @noRd
prior_coef_brms <- function(cn) {
  out <- brms_rename(cn)
  fx <- grepl("[.]fx[0-9]+$", cn)
  if (any(fx)) {
    lab <- brms_rename(sub("[.]fx[0-9]+$", "", cn[fx]))
    out[fx] <- paste0(lab, "_", sub("^.*[.]fx([0-9]+)$", "\\1", cn[fx]))
  }
  out
}

#' A block's coefficient rows in brms's spelling and brms's order.
#'
#' brms sorts its prior table by coefficient within a class
#' (`.default_prior()`, `order(..., class, group, coef)`), so a slope
#' named `z` follows a smooth's `sx_1` whatever order the formula put
#' them in. `order()` is brms's own call here, locale and all.
#'
#' @noRd
prior_coef_rows <- function(cn) {
  out <- prior_coef_brms(cn)
  out[order(out)]
}

#' Which design columns a written `coef` names: the column name, the same
#' with parentheses dropped (`Intercept`), or brms's spelling of it.
#'
#' @noRd
prior_coef_match <- function(cn, coef) {
  cn == coef | par_name_bare(cn) == par_name_bare(coef) |
    prior_coef_brms(cn) == coef
}

#' Which of brms's coefficient names `bn` (a block's standard deviations,
#' a `cs()` term's coef) a written `coef` names: brms's name itself, or
#' the design column it renames, so `coef = "I(x^2)"` reaches `IxE2`
#' here as it does for class `"b"` (`prior_coef_match()`).
#'
#' @noRd
prior_coef_match_brms <- function(bn, coef) {
  !is.na(bn) & (bn == coef | bn == par_name_bare(coef) |
                  bn == brms_rename(coef))
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
                  resp = "", lb = NA_real_, dflt = NULL, ub = NA_real_) {
    # a default speaks for a CLASS, not for one coefficient of it: brms
    # reports the class row and leaves the per-coefficient rows flat, and
    # a per-coefficient row here would claim a default nothing applies.
    # A default on cor is written class-wide, as brms's cor rows carry no
    # prefix; one on sd is written per prefix, as brms's sd rows are
    # The exception is a default WRITTEN per coefficient, which a gp()
    # length scale's is: brms tunes its inv_gamma to each sub-GP's own
    # distances and lists it on that sub-GP's row
    slot <- prior_slot_key(class, if (identical(class, "cor")) "" else dpar,
                           nlpar, resp)
    d <- if (nzchar(coef)) defs[[paste0(slot, "|", coef)]] else defs[[slot]]
    # brms's column order, so a table read by position lines up with
    # the one brms returns
    rows[[length(rows) + 1L]] <<- data.frame(
      prior = d %||% dflt %||% "(flat)", class = class, coef = coef,
      group = group,
      resp = resp, dpar = dpar, nlpar = nlpar, lb = lb,
      ub = ub, source = "default"
    )
  }

  simo_rows <- list()
  for (lp in frame[["linpreds"]]) {
    # an equated dpar has no parameter of its own; its target's row is
    # the one brms lists
    if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]]) ||
          !is.null(lp[["equate"]])) next
    rspec <- spec$responses[[lp[["resp"]]]]
    # location dpars are the default target (dpar = ""), matching
    # set_prior()'s resolution
    dpar_lab <- lp_prior_dpar(rspec, lp[["dpar"]])
    # the shared nu of a Student-t rescor model has no response
    resp_lab <- if (multi && !isTRUE(lp[["shared"]])) lp[["resp"]] else ""
    nl_lab <- if (lp[["dpar"]] %in% (rspec$nlpars %||% character(0))) {
      lp[["dpar"]]
    } else {
      ""
    }
    # one simo row per mo() simplex, as brms lists them; added after the
    # other classes, where brms's table puts them
    for (mi in lp[["mo"]] %||% list()) {
      simo_rows[[length(simo_rows) + 1L]] <- list(
        coef = paste0(mi$label, "1"),
        dpar = if (nzchar(nl_lab)) "" else dpar_lab, nlpar = nl_lab,
        resp = resp_lab)
    }
    cn <- colnames(lp[["X"]])
    # a nonlinear parameter's sub-formula is not centered, so its whole
    # coefficient vector including the intercept is class "b" - brms's
    # own listing, and what set_prior(nlpar =) addresses
    if (nzchar(nl_lab)) {
      if (length(cn)) {
        add("b", dpar = dpar_lab, nlpar = nl_lab, resp = resp_lab)
        for (co in prior_coef_rows(cn)) {
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
    icpt_b <- isFALSE(lp[["center"]])
    if ("(Intercept)" %in% cn && !icpt_b) {
      add("Intercept", dpar = dpar_lab, resp = resp_lab)
    } else if (!is.null(ob <- ord_lp_block(frame, spec, lp)) &&
                 (!ob$shared || identical(ob$k, 1L)) &&
                 length(frame[["par_template"]][[ob$comp]] %||%
                          numeric(0))) {
      # an ordinal family has no intercept column: the thresholds
      # replace it, and class "Intercept" is what addresses them here as
      # it does in brms. Under sum_to_zero brms's Intercept is a vector
      # frmtmb does not have (ordinal_threshold_entry()), and under
      # equidistant the distance between thresholds is class "delta",
      # listed once for the model and once per level of thres(gr = ).
      # An ordinal mixture's component k lists its own under dpar
      # "mu<k>"; the vector an order = "mu" mixture shares is listed
      # once, with no dpar, as brms lists its fixed_Intercept
      th <- ob$fam[["thres"]]
      type <- th[["type"]] %||% "flexible"
      grouped <- isTRUE(th[["grouped"]])
      tdp <- ob$prior_dpar
      if (identical(type, "equidistant")) {
        # brms bounds delta below by 0 where the thresholds are ordered,
        # and says so in its table
        lb_d <- if (ob$ordered) 0 else NA_real_
        add("delta", dpar = tdp, resp = resp_lab, lb = lb_d)
        if (grouped) {
          for (g in th[["groups"]]) {
            add("delta", group = g, dpar = tdp, resp = resp_lab, lb = lb_d)
          }
        }
      }
      if (!identical(type, "sum_to_zero")) {
        add("Intercept", dpar = tdp, resp = resp_lab)
        # under flexible thresholds brms lists one row per threshold,
        # coef = "1", "2", ..., under each level of thres(gr = ) when
        # they are grouped; equidistant thresholds have none. The count
        # is the block's own component: in a mixture the plain
        # "tau_raw" component does not exist, and each row carries the
        # block's dpar, "mu<k>" under order = "none" and none for the
        # vector an order = "mu" mixture shares, as brms lists them
        flex <- identical(type, "flexible")
        if (!grouped && flex) {
          K1 <- length(frame[["par_template"]][[ob$comp]])
          for (k in seq_len(K1)) {
            add("Intercept", coef = as.character(k), dpar = tdp,
                resp = resp_lab)
          }
        }
        # one row per level of thres(gr = ), each a threshold vector of
        # its own that group = addresses, as brms lists them
        if (grouped) {
          for (gi in seq_along(th[["groups"]])) {
            g <- th[["groups"]][gi]
            add("Intercept", group = g, dpar = tdp, resp = resp_lab)
            if (flex) {
              for (k in seq_len(th[["nthres"]][gi])) {
                add("Intercept", coef = as.character(k), group = g,
                    dpar = tdp, resp = resp_lab)
              }
            }
          }
        }
      }
    }
    # cs() terms are class "b" rows under their own coef, as brms lists
    # them; the resolver reaches them through the same class
    # brms spells the uncentered intercept's coef "Intercept"
    others <- c(if (icpt_b) sub("^[(]Intercept[)]$", "Intercept", cn) else
                  setdiff(cn, "(Intercept)"),
                vapply(lp[["cs"]] %||% list(), cs_term_coef, ""))
    if (length(others)) {
      add("b", dpar = dpar_lab, resp = resp_lab)
      for (co in prior_coef_rows(others)) {
        add("b", coef = co, dpar = dpar_lab, resp = resp_lab)
      }
    }
  }

  sd_rows <- list()
  cor_rows <- list()
  gp_rows <- list()
  for (bk in frame[["re_blocks"]]) {
    key <- list(group = bk[["group_name"]], nlpar = block_nlpar(spec, frame,
                                                                bk),
                dpar = block_dpar(spec, frame, bk),
                resp = if (multi) block_resp(frame, bk) else "")
    if (bk[["covstruct"]] %in% c("gp", "hsgp")) {
      # brms's classes sdgp and lscale, keyed by the prefix class "sd"
      # uses, coef the sub-GP (frame_gp()'s sfx1 and sfx2)
      sp <- unique(block_sd_prefix(spec, frame, bk)[, c("resp", "dpar",
                                                        "nlpar")])
      for (i in seq_len(nrow(sp))) {
        gp_rows[[length(gp_rows) + 1L]] <- list(
          nlpar = sp$nlpar[i], dpar = sp$dpar[i], resp = sp$resp[i],
          sfx1 = bk[["gp_brms"]][["sfx1"]],
          sfx2 = bk[["gp_brms"]][["sfx2"]])
      }
      next
    }
    # brms lists sd rows per prefix, so a block spanning mu and sigma
    # has a row for each, and below each group row one per coefficient
    # (sd_spec_reach() is the rule the rows describe)
    sp <- block_sd_prefix(spec, frame, bk)
    cf <- block_sd_coefs(frame, bk)
    sd_rows[[length(sd_rows) + 1L]] <- data.frame(
      group = bk[["group_name"]], nlpar = sp$nlpar, dpar = sp$dpar,
      resp = sp$resp, coef = if (nrow(sp)) cf[sp$k] else character(0),
      stringsAsFactors = FALSE)
    if (identical(block_cor_prior(bk), "lkj")) {
      cor_rows[[length(cor_rows) + 1L]] <- key
    }
  }
  sd_rows <- do.call(rbind, sd_rows)
  if (NROW(sd_rows)) {
    # the class-wide row of each prefix, its group rows and their
    # coefficient rows, in brms's order(resp, dpar, nlpar, class, group,
    # coef) (.default_prior()); `""` sorts first, so each wider row
    # comes before the rows it covers
    wide <- unique(rbind(
      transform(sd_rows, group = "", coef = ""),
      transform(sd_rows, coef = ""),
      sd_rows[!is.na(sd_rows$coef), , drop = FALSE]))
    wide <- wide[with(wide, order(resp, dpar, nlpar, group, coef)), ,
                 drop = FALSE]
    for (i in seq_len(nrow(wide))) {
      add("sd", coef = wide$coef[i], group = wide$group[i],
          dpar = wide$dpar[i], nlpar = wide$nlpar[i], resp = wide$resp[i])
    }
  }
  if (length(cor_rows)) {
    add("cor")
    for (k in cor_rows) {
      add("cor", group = k$group, dpar = k$dpar, nlpar = k$nlpar,
          resp = k$resp)
    }
  }
  # brms lists lscale before sdgp, each class-wide and then per sub-GP,
  # bounded below by 0
  for (cl in c("lscale", "sdgp")) {
    pre <- unique(lapply(gp_rows, `[`, c("dpar", "nlpar", "resp")))
    for (k in pre) add(cl, dpar = k$dpar, nlpar = k$nlpar, resp = k$resp,
                       lb = 0)
    for (k in gp_rows) {
      for (co in if (identical(cl, "sdgp")) k$sfx1 else k$sfx2) {
        add(cl, coef = co, dpar = k$dpar, nlpar = k$nlpar, resp = k$resp,
            lb = 0)
      }
    }
  }
  # the R-side residual structures, under the class names brms shows for
  # them. `resp` narrows a class to one response, as it does above
  for (rs in names(frame[["autocor"]] %||% list())) {
    ac <- frame[["autocor"]][[rs]]
    resp_lab <- if (multi) rs else ""
    for (cl in names(autocor_prior_classes)) {
      if (!length(autocor_class_idx(ac, cl))) next
      bd <- autocor_class_bounds(ac, cl, route)
      add(cl, resp = resp_lab, lb = bd[[1L]], ub = bd[[2L]])
    }
    # a registry term's bounded parameters (sar()'s lagsar), under
    # their natural names and on their interval, as brms lists them
    for (pr in autocor_seam_bounded(ac)) {
      add(pr$name, resp = resp_lab, lb = pr$lower, ub = pr$upper)
    }
  }
  if (length(frame[["par_template"]][["thetar"]] %||% numeric(0))) add("rescor")
  # a logistic_normal() response's correlations, brms's class lncor;
  # frmtmb.sample supplies brms's lkj(1) default, as for rescor. One row
  # per such response, which a multivariate model names, as brms does
  for (b in lncor_blocks(list(spec = spec, frame = frame))) {
    add("lncor", resp = if (multi) b$resp else "")
  }
  # brms's default on a simplex is dirichlet(1), the flat density on it;
  # frm() applies no density there, and frm_sample() applies brms's
  for (sr in simo_rows) {
    add("simo", coef = sr$coef, dpar = sr$dpar, nlpar = sr$nlpar,
        resp = sr$resp,
        dflt = if (identical(route, "sample")) "dirichlet(1)")
  }
  # and on a simplex a family declares, one row per response under its
  # brms name: class sbhaz for a cox() baseline
  for (tm in family_simplex_terms(frame, spec)) {
    add(tm$name, resp = if (multi) tm$resp else "",
        dflt = if (identical(route, "sample")) "dirichlet(1)")
  }

  # class "theta" is the raw internal escape hatch, and it spans all
  # three covariance components, so these rows name one internal
  # parameter on its own.
  #
  # [[ ]] rather than $ here and below: `$theta` partial-matches
  # `thetaac` in a model that has a residual structure and no random
  # effects, which offered a random-effect row for a model with none
  if (length(frame[["par_template"]][["theta"]] %||% numeric(0))) add("theta")
  # a residual-correlation parameter that a class reaches (ar, cosy,
  # cortime, a registry term's lagsar) is listed once, under that class
  # above, as brms lists it: a second row would be a second prior slot
  # on the same parameter
  own_class <- autocor_class_coefs(frame)$coef
  for (cmp in theta_components) {
    v <- frame[["par_template"]][[cmp]] %||% numeric(0)
    if (!length(v)) next
    for (nm in par_template_names(v, cmp)) {
      if (identical(cmp, "thetaac") && nm %in% own_class) next
      add("theta", coef = nm)
    }
  }

  out <- unique(do.call(rbind, rows))
  # brms's row order (.default_prior(): resp, dpar, nlpar, class, group,
  # coef, with R's order(), locale and all), so b comes before cor,
  # Intercept and sd; frmtmb's own class "theta", which brms does not
  # have, stays last in its own order
  if (!is.null(out)) {
    own <- out$class == "theta"
    brm <- out[!own, , drop = FALSE]
    brm <- brm[with(brm, order(resp, dpar, nlpar, class, group, coef)), ,
               drop = FALSE]
    out <- rbind(brm, out[own, , drop = FALSE])
  }
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
#' dpars (the latent-class, hidden Markov, multivariate normal mixture
#' and accumulator families) follows the same rule, which is frmtmb's for
#' families brms does not have: a prior with no dpar does not say which
#' location it means.
#'
#' @noRd
multi_location <- function(rspec) {
  fam <- rspec$family
  loc <- setdiff(rspec$primary_dpars, rspec$nlpars %||% character(0))
  length(loc) > 1L || !is.null(fam[["mix"]]) ||
    isTRUE(fam[["family"]] %in% c("categorical", multinomial_families))
}

#' Whether brms has this several-location family, so a refusal may say
#' "as in brms" rather than name frmtmb's own rule.
#'
#' @noRd
brms_multi_location <- function(rspec) {
  fam <- rspec$family
  fn <- fam[["family"]] %||% ""
  fn %in% c("categorical", multinomial_families) ||
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
#' A `coef` narrows the specification to the standard deviation of that
#' one coefficient, as brms's `sd_1[2]` row does, and then the prefix
#' must EQUAL the standard deviation's even in a block that spans
#' several: brms puts `set_prior(class = "sd", group = "g", coef =
#' "Intercept")` on mu's intercept alone in `(1 + x | p | g)` with
#' `sigma ~ (1 | p | g)` (dev/prilabel-log/brms-ref3.txt). A structure
#' whose columns share one standard deviation has none per coefficient
#' (`block_sd_coefs()`); the caller refuses a `coef` that names one of
#' its columns.
#'
#' Returns a list: `reach`, a logical per entry of `block_sd_idx()`, and
#' `exact`, whether the specification's own prefix (and coefficient) is
#' one this block has, which is what brms requires for the row to exist
#' at all.
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
  if (nzchar(s$coef %||% "")) {
    cf <- block_sd_coefs(frame, bk)
    hit <- prior_coef_match_brms(cf, s$coef)
    reach <- vapply(seq_len(n), function(j) hit[j] && any(eq[sp$k == j]),
                    TRUE)
    return(list(reach = reach, exact = any(reach)))
  }
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

#' The structures whose standard deviations no brms coefficient names,
#' beside the fields (`cs_is_field()`: a smooth's, brms's class `"sds"`,
#' a gp() term's, `"sdgp"`, a CAR or SPDE field's): a reduced-rank or
#' fixed block, which have none of their own per column.
#'
#' @noRd
sd_coef_unaddressable <- c("rr", "equalto")

#' brms's coefficient name for each standard deviation of a block, one
#' per entry of `block_sd_idx()`: `Intercept`, `x`, `fb`, `moxo`, as
#' brms's class `"sd"` rows name them (`brms_rename()` of the column
#' name, which is what brms's `get_rnames()` gives before the prefix).
#' `NA` for every entry of a block whose columns share their standard
#' deviations (`homdiag()`, `ar1()`, `homcs()`, ...), where no entry is
#' one coefficient's, and of a block in `sd_coef_unaddressable`.
#'
#' @noRd
block_sd_coefs <- function(frame, bk) {
  n <- length(block_sd_idx(bk))
  d <- bk[["dim"]] %||% 0L
  if (!n || n != d || bk[["covstruct"]] %in% sd_coef_unaddressable ||
        cs_is_field(bk) ||
        length(bk[["cnms"]]) != d) {
    return(rep(NA_character_, n))
  }
  block_col_coefs(bk)
}

#' brms's coefficient name for each COLUMN of a block, in column order:
#' each component's column names, `brms_rename()`d.
#'
#' @noRd
block_col_coefs <- function(bk) {
  d <- bk[["dim"]] %||% 0L
  cf <- bk[["cnms"]] %||% rep(NA_character_, d)
  for (cp in bk[["components"]] %||% list()) {
    if (is.null(cp[["offset"]]) || is.null(cp[["cnms"]])) next
    pos <- cp[["offset"]] + seq_along(cp[["cnms"]])
    ok <- pos <= d
    cf[pos[ok]] <- cp[["cnms"]][ok]
  }
  brms_rename(cf)
}

#' A block as a message names it: a `||` or `gr(cor = FALSE)` term as
#' it was written, `(1 + f || g)`, and any other term by its label and
#' structure, `1 + f | g [us]`. The label alone would read a `||` term
#' as the correlated `1 + f | g`, which is the term it is not.
#'
#' @noRd
block_written_label <- function(bk) {
  wt <- unique(unlist(lapply(bk[["components"]] %||% list(), function(cp) {
    cp[["written_text"]]
  })))
  if (length(wt)) return(paste(wt, collapse = " + "))
  paste0(bk[["term_label"]], " [", bk[["covstruct"]], "]")
}

#' What a refused per-coefficient class `"sd"` specification could have
#' named: the coefficients each group's standard deviations carry, as
#' `default_prior()` lists them.
#'
#' @noRd
sd_coef_hint <- function(frame, s) {
  bks <- Filter(function(bk) {
    !nzchar(s$group) || identical(bk[["group_name"]], s$group)
  }, frame[["re_blocks"]])
  grp <- vapply(bks, function(bk) bk[["group_name"]] %||% "", "")
  cf <- lapply(bks, function(bk) {
    v <- block_sd_coefs(frame, bk)
    v[!is.na(v)]
  })
  lab <- vapply(unique(grp), function(g) {
    v <- unique(unlist(cf[grp == g]))
    if (!length(v)) return(NA_character_)
    paste0("group \"", g, "\": ", paste0("\"", v, "\"", collapse = ", "))
  }, "")
  lab <- lab[!is.na(lab)]
  if (!length(lab)) {
    return(paste0("No standard deviation of ",
                  if (nzchar(s$group)) paste0("group \"", s$group, "\"")
                  else "this model",
                  " is one coefficient's. "))
  }
  paste0("Per-coefficient standard deviations: ",
         paste(lab, collapse = "; "), ". ")
}

#' The columns of a block a `coef` names whose standard deviation is
#' shared with other columns, so that a prior on it cannot be one
#' coefficient's. `character(0)` when there is none. A column counts
#' only where the specification's prefix (resp, dpar, nlpar) EQUALS the
#' column's, the rule of `sd_spec_reach()`'s coef branch: a
#' `sigma ~ homdiag(1 + x | g)` block does not stop a coef prior with no
#' dpar from reaching mu's `(1 + x | g)`, as brms puts it there.
#'
#' @noRd
sd_coef_shared <- function(spec, frame, bk, s) {
  if (!nzchar(s$coef %||% "") || !length(block_sd_idx(bk)) ||
        bk[["covstruct"]] %in% sd_coef_unaddressable) {
    return(character(0))
  }
  if (nzchar(s$group) && !identical(bk[["group_name"]], s$group)) {
    return(character(0))
  }
  if (!all(is.na(block_sd_coefs(frame, bk)))) return(character(0))
  cf <- block_col_coefs(bk)
  cp <- block_col_prefix(spec, frame, bk)
  if (length(cf) != nrow(cp)) return(character(0))
  eq <- cp$resp == (s$resp %||% "") & cp$dpar == (s$dpar %||% "") &
    cp$nlpar == (s$nlpar %||% "")
  cf[eq & prior_coef_match_brms(cf, s$coef)]
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
    # a gp() term's sd is class "sdgp", which the caller names
    lapply(Filter(function(bk) !bk[["covstruct"]] %in% c("gp", "hsgp"),
                  frame[["re_blocks"]]), function(bk) {
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
#' written, coef in brms's spelling (`prior_coef_brms()`: `(Intercept)`
#' is `Intercept`, the design column `I(x^2)` is `IxE2`, a smooth's
#' `s(x).fx1` is `sx_1`, each pair naming one slot), group, resp, dpar
#' and nlpar, each kept apart.
#'
#' @noRd
prior_spec_slot <- function(s) {
  sp <- spec_spelling(s)
  paste(sp$class, prior_coef_brms(s$coef %||% ""), s$group %||% "",
        s$resp %||% "", sp$dpar, s$nlpar %||% "", sep = "\r")
}

#' Drop the specifications of a stored prior that the updated model
#' cannot take, as brms's `update()` drops them.
#'
#' brms's `update.brmsfit()` passes the stored prior back to `brm()`
#' marked `allow_invalid_prior`, and a specification that matches no
#' parameter of the new model is then dropped with no message: on
#' `brms_overview`'s fit2, which removes the last correlation, brms
#' recompiles without its class-wide `lkj(2)` and samples. frmtmb's
#' `update()` refused the same call, because frm() refuses such a
#' prior when it is written directly, as brms does. `update()` sets the
#' same mark on the prior it carries over, and this drops what does not
#' resolve on the new model. It says what it dropped, where brms is
#' silent: a refusal is replaced by a drop, and the drop must not be
#' less diagnostic than the refusal was, so the user still reads which
#' of their specifications is gone.
#'
#' A specification is kept exactly when it resolves alone on the new
#' model, which is the test `resp_missing_refusal()` uses. The
#' specifications the user passes to `update()` itself are marked
#' `update_strict` (`update_merge_prior()`), are never dropped, and are
#' checked as frm() checks any prior.
#'
#' @noRd
prior_drop_unmatched <- function(prior, spec, frame) {
  if (!isTRUE(attr(prior, "allow_invalid_prior"))) return(prior)
  attr(prior, "allow_invalid_prior") <- NULL
  strict <- attr(prior, "update_strict")
  attr(prior, "update_strict") <- NULL
  # drop1() refits through the same door, and the message names it
  who <- attr(prior, "drop_who") %||% "update()"
  attr(prior, "drop_who") <- NULL
  if (!inherits(prior, "frmtmb_priorlist")) return(prior)
  design <- list(spec = spec, frame = frame)
  specs <- unclass(prior)
  ok <- vapply(specs, function(s) {
    one <- structure(list(s), class = "frmtmb_priorlist")
    !inherits(tryCatch(resolve_priorlist(design, one),
                       error = keep_missing(identity)),
              "error")
  }, NA)
  if (length(strict) == length(ok)) ok[strict] <- TRUE
  if (all(ok)) return(prior)
  lab <- vapply(specs[!ok], function(s) {
    what <- s[["prior"]] %||% prior_dist_string(s[["dist"]])
    # a specification with bounds and no density has an empty prior
    if (!length(what) || !nzchar(what)) what <- "bounds"
    paste0(what, " on ", prior_slot_label(s))
  }, "")
  one <- sum(!ok) == 1L
  frm_message(who, " drops ",
              if (one) "the prior specification" else
                paste(sum(!ok), "prior specifications"),
              " of the original fit that ", if (one) "matches" else "match",
              " no parameter of the ",
              if (identical(who, "update()")) {
                "updated model, as brms drops "
              } else {
                "refitted model, as update() drops "
              },
              if (one) "it" else "them", ": ", paste(lab, collapse = "; "))
  if (!any(ok)) return(NULL)
  structure(specs[ok], class = class(prior))
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
  prior_bound_warning(specs)
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
                              "cortime", "delta")

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
  # the nu a Student-t rescor model shares across its responses is one
  # parameter with no response in it, and brms's row for it has none;
  # it is the one dpar such a model addresses without resp
  if (!is.null(spec[["rescor_nu"]]) && identical(s$dpar %||% "", "nu")) {
    return(NULL)
  }
  sp <- spec_spelling(s)
  # the responses the same specification resolves on once it names
  # them, found by resolving it: reading them off the default_prior()
  # table offered calls the resolver then refused, and missed ones it
  # accepts (a dpar the table did not list, class Intercept with nlpar)
  design <- list(spec = spec, frame = frame)
  resolves <- function(s1) {
    one <- structure(list(s1), class = "frmtmb_priorlist")
    !inherits(tryCatch(resolve_priorlist(design, one),
                       error = keep_missing(identity)),
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
               "random-effect covariance parameters. A residual-",
               "correlation parameter is addressed by its class (\"ar\", ",
               "\"ma\", \"cosy\", \"cortime\"), as in brms, and one with ",
               "no class by its internal name in coef, as ",
               "get_prior() lists it", call. = FALSE)
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
      if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]]) ||
            !is.null(lp[["equate"]])) next
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
      # `0 + Intercept` or center = FALSE: the intercept is a class "b"
      # coefficient like any other, and class "Intercept" has no slot here
      icpt_b <- isFALSE(lp[["center"]])
      pick <- if (s$class == "Intercept") {
        if (!icpt_b) which(cn == "(Intercept)") else integer(0)
      } else if (nzchar(s$coef)) {
        # brms writes an intercept as "Intercept" and a smooth's
        # unpenalized column as sx_1; the design matrix spells them
        # "(Intercept)" and s(x).fx1, and both name the same column
        which(prior_coef_match(cn, s$coef))
      } else if (want_np) {
        # a nonlinear parameter's sub-formula is NOT centered, so its
        # intercept sits in the same coefficient vector as its slopes
        # and class "b" covers it. This is the whole reason brms's
        # prior(normal(5000, 1000), nlpar = "ult") lands on an
        # intercept-only nonlinear parameter
        seq_along(cn)
      } else {
        which(cn != "(Intercept)" | icpt_b)
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
                                        center = ctr,
                                        resp = lp[["resp"]],
                                        dpar = lp[["dpar"]])
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
                !prior_coef_match_brms(lab, s$coef) &&
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
      hint <- if (s$class == "Intercept" &&
                    has_ordinal_thresholds(fit, s$resp)) {
        paste0(". On an ordinal family the thresholds ARE the ",
               "intercept, and a bare class = \"Intercept\" with no ",
               "coef, dpar or nlpar addresses the whole threshold ",
               "vector; prior = list(tau_raw = ) reaches the same ",
               "parameters on the internal scale")
      } else if (s$class == "Intercept" && any(vapply(
        frame[["linpreds"]], function(lp) isFALSE(lp[["center"]]), NA))) {
        paste0(". An intercept written 0 + Intercept, or under ",
               "center = FALSE, is an ordinary coefficient, as in brms: ",
               "address it with class = \"b\", coef = \"Intercept\"")
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
  # brms's classes "sdgp" and "lscale" on a gp() term's theta: the sd
  # and the length scales on brms's natural scales, a sub-GP picked by
  # brms's coef ("gpxfa"). The length scale brms priors is on its scaled
  # inputs, exp(theta) / div, so the density sits at theta - log(div)
  # (`shift`) with the same log-Jacobian class "sd" takes. resp, dpar
  # and nlpar select by class "sd"'s rule (sd_spec_reach()).
  resolve_gp_class <- function(s) {
    hit <- FALSE
    # the prefix only: a gp() coef names a sub-GP, read below, not a
    # column's standard deviation
    s_pre <- s
    s_pre$coef <- ""
    for (bk in frame[["re_blocks"]]) {
      if (!bk[["covstruct"]] %in% c("gp", "hsgp")) next
      meta <- bk[["gp_brms"]]
      if (is.null(meta)) next
      if (!any(sd_spec_reach(fit$spec, frame, bk, s_pre)$reach)) next
      if (identical(s$class, "sdgp")) {
        pick <- 1L
        cf <- meta$sfx1
      } else {
        pick <- 1L + seq_along(meta$sfx2)
        cf <- meta$sfx2
      }
      keep <- !nzchar(s$coef) | cf == s$coef
      shift <- if (identical(s$class, "lscale")) {
        -log(bk[["gp_lscale_div"]] %||% 1)
      } else 0
      for (k in pick[keep]) {
        hit <- TRUE
        i <- bk[["theta_idx"]][k]
        key <- nm_of("theta", i)
        if (!is.null(s$dist)) {
          claim("theta", i)
          assigned[[key]] <<- list(comp = "theta", idx = i, dist = s$dist,
                                   scale = "sd", shift = shift,
                                   lb = s$lb, ub = s$ub)
        } else if (!is.null(assigned[[key]])) {
          assigned[[key]] <<- entry_bounds(assigned[[key]], s)
        }
        nm_theta <- paste0("theta_", i)
        if (!is.na(s$lb)) {
          lower[nm_theta] <<- if (s$lb > 0) log(s$lb) - shift else -Inf
        }
        if (!is.na(s$ub)) upper[nm_theta] <<- log(s$ub) - shift
      }
    }
    if (!hit) {
      have <- unlist(lapply(frame[["re_blocks"]], function(bk) {
        m <- bk[["gp_brms"]]
        if (is.null(m)) NULL else if (identical(s$class, "sdgp")) m$sfx1 else
          m$sfx2
      }))
      frm_stop("No gp() ", if (identical(s$class, "sdgp")) {
        "standard deviations"
      } else "length scales", " match ", spec_target(s), ". ",
      if (length(have)) {
        paste0("This model's are coef = ",
               paste0("\"", unique(have), "\"", collapse = ", "))
      } else "This model has no gp() term", call. = FALSE)
    }
    invisible(NULL)
  }

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
                   ". brms has no such bound either: its cov = TRUE form ",
                   "is of order one only. A bound needs order one, or the ",
                   "cov = FALSE form, whose coefficients are unconstrained",
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

  # A logistic_normal() response's correlations, held in the same chart
  # as the residual correlations above, so the same LKJ density applies.
  resolve_lncor <- function(s) {
    bl <- lncor_blocks(fit)
    # resp narrows the class to one response of a multivariate model,
    # named as frmtmb names it or as brms does (its Stan name, Ya for
    # Y_a), since brms's prior table reads the second
    if (nzchar(s$resp %||% "")) {
      bl <- Filter(function(b) {
        s$resp %in% c(b$resp, brms_stan_name(b$resp))
      }, bl)
    }
    if (!length(bl)) {
      frm_stop("No logistic-normal correlation matches ", spec_target(s),
               ". Class \"lncor\" is the correlation matrix of a ",
               "logistic_normal() response with three or more categories, ",
               "and this model has none", call. = FALSE)
    }
    for (b in bl) {
      idx <- seq_along(fit$frame[["par_template"]][[b$comp]])
      claim(b$comp, idx)
      assigned[[nm_of(b$comp, idx)]] <<-
        list(comp = b$comp, idx = idx,
             dist = lkj_dist(s$dist$eta,
                             list(kind = "chol", d = length(b$sx$pred))),
             scale = "internal", lb = NA, ub = NA)
    }
  }

  # brms's class "Intercept" on an ordinal family names the THRESHOLDS:
  # its ordinal program declares them as the `Intercept` vector and puts
  # the default student_t there. frmtmb holds them in `tau_raw`, which
  # had no class spelling at all, so the row used to reach the resolver
  # and stop with a bare "Prior target not found".
  # A multivariate model holds each ordinal response's thresholds under
  # a name of its own, and there the specification always carries resp
  # (resp_missing_refusal() has refused it otherwise).
  ordinal_threshold_entry <- function(s) {
    rs <- fit$spec$responses
    rspec <- if (length(rs) == 1L) rs[[1L]] else rs[[s$resp %||% ""]]
    if (is.null(rspec)) return(NULL)
    if (!identical(rspec$family[["type"]], "ordinal")) return(NULL)
    # coef = "1", "2", ... names one threshold (below), so only an
    # nlpar leaves the thresholds
    if (nzchar(s$nlpar %||% "")) return(NULL)
    # the block the dpar names: "" for a plain family's thresholds and
    # for the vector an order = "mu" mixture shares, mu<k> for one
    # component's own under order = "none", as brms keys the rows
    ob <- ord_prior_block(frame, fit$spec, rspec, s$dpar)
    if (is.null(ob)) return(NULL)
    comp <- ob$comp
    raw <- frame[["par_template"]][[comp]] %||% numeric(0)
    if (!length(raw)) return(NULL)
    if (nzchar(s$resp %||% "") &&
          !identical(s$resp, rspec$resp_name)) {
      return(NULL)
    }
    # cumulative() holds (tau_1, log increments), which is the same map
    # Stan's `ordered` type applies, so the density on the thresholds
    # carries that map's log-Jacobian. sratio(), cratio() and acat()
    # hold the thresholds themselves and brms declares them unordered
    # (brms:::has_ordered_thres() is FALSE for all three), so neither
    # side has a Jacobian there
    ordered <- ob$ordered
    th <- ob$fam[["thres"]]
    grouped <- isTRUE(th[["grouped"]])
    type <- th[["type"]] %||% "flexible"
    if (identical(type, "sum_to_zero")) {
      frm_stop("Prior target not found (", spec_target(s), "): with ",
               "threshold = 'sum_to_zero' there is no threshold vector ",
               "for class \"Intercept\" to address. brms puts that prior ",
               "on thresholds it declares before centering them, whose ",
               "common location the likelihood cannot see and only the ",
               "prior places; frmtmb estimates the centered thresholds ",
               "directly, so that parameter does not exist here. ",
               "prior = list(", comp, " = ) reaches the free directions ",
               "on the internal scale", call. = FALSE)
    }
    if (nzchar(s$group) && !grouped) {
      frm_stop("Prior target not found (", spec_target(s), "): group = ",
               "on class \"Intercept\" names one threshold vector of a ",
               "model with grouped thresholds, thres(gr = ), and this ",
               "model has one threshold vector. Drop group", call. = FALSE)
    }
    # equidistant: brms's class "Intercept" is about the first threshold
    # alone, `first_Intercept`, centered as the whole vector would be;
    # the distance is class "delta"
    equi <- identical(type, "equidistant")
    # brms's per-threshold rows, coef = "1", "2", ...: under flexible
    # thresholds each is one element of the vector, under each level of
    # thres(gr = ) when the thresholds are grouped. brms has no such
    # parameter under equidistant thresholds and refuses the row
    k <- NULL
    if (nzchar(s$coef)) {
      if (equi) {
        frm_stop("Prior target not found (", spec_target(s), "): with ",
                 "threshold = 'equidistant' class \"Intercept\" is the ",
                 "first threshold alone, so it takes no coef; brms ",
                 "refuses the row too (\"The following priors do not ",
                 "correspond to any model parameter\"). The distance ",
                 "between thresholds is class \"delta\"", call. = FALSE)
      }
      if (grouped && !nzchar(s$group)) {
        frm_stop("Prior target not found (", spec_target(s), "): the ",
                 "thresholds are grouped by thres(gr = ), so one ",
                 "threshold is named by its level and its number, ",
                 "group = and coef =, as brms lists them", call. = FALSE)
      }
      k <- suppressWarnings(as.integer(s$coef))
      if (is.na(k) || !identical(as.character(k), s$coef) || k < 1L) {
        frm_stop("Prior target not found (", spec_target(s), "): coef = ",
                 "on class \"Intercept\" of an ordinal family numbers a ",
                 "threshold, \"1\", \"2\", ..., as brms lists them",
                 call. = FALSE)
      }
    }
    check_k <- function(n_k) {
      if (!is.null(k) && k > n_k) {
        frm_stop("Prior target not found (", spec_target(s), "): the ",
                 "threshold vector has ", n_k, " threshold",
                 if (n_k != 1L) "s", ", so coef = \"", k, "\" names none ",
                 "of them", call. = FALSE)
      }
    }
    if (!grouped) {
      check_k(length(raw))
      return(list(list(comp = comp, idx = if (equi) 1L else seq_along(raw),
                       dist = s$dist,
                       scale = if (ordered && !equi) "ordthres" else
                         "internal",
                       link = NULL,
                       # brms centers no design under order = "mu"
                       # (stan_center_X() is FALSE with fixed thresholds)
                       offset = if (!ob$shared) {
                         ordinal_center_offset(frame, rspec,
                                               ob$lp[["dpar"]])
                       },
                       lb = s$lb, ub = s$ub, coef_k = k)))
    }
    # one entry per group: each slice is a vector of its own, so an
    # ordered map and its Jacobian are per slice. brms does not center
    # the design of a model with grouped thresholds, so no offset
    lay <- thres_layout(th[["nthres"]], type)
    gs <- seq_len(lay$G)
    if (nzchar(s$group)) {
      gs <- which(th[["groups"]] == s$group)
      if (!length(gs)) {
        frm_stop("Prior target not found (", spec_target(s), "): the ",
                 "thresholds are grouped by thres(gr = ), whose levels ",
                 "are ", paste0("\"", th[["groups"]], "\"", collapse = ", "),
                 call. = FALSE)
      }
    }
    lapply(gs, function(g) {
      check_k(lay$nthres[g])
      list(comp = comp,
           idx = if (equi) lay$rstart[g] else lay$rstart[g]:lay$rend[g],
           dist = s$dist,
           scale = if (ordered && !equi) "ordthres" else "internal",
           link = NULL, offset = NULL, lb = s$lb, ub = s$ub, coef_k = k)
    })
  }

  # brms's class "delta": the distance between neighboring thresholds
  # of an equidistant threshold vector, one per vector. brms bounds it
  # below by 0 for the ordered families, where frmtmb holds log(delta),
  # so the density there is on exp() of the internal value with the log
  # Jacobian, the "sd" placement; the other families hold delta itself
  ordinal_delta_entry <- function(s) {
    rs <- fit$spec$responses
    rspec <- if (length(rs) == 1L) rs[[1L]] else rs[[s$resp %||% ""]]
    # component k of an ordinal mixture: dpar = "mu<k>", brms's row
    ob <- if (!is.null(rspec)) {
      ord_prior_block(frame, fit$spec, rspec, s$dpar)
    }
    th <- ob$fam[["thres"]]
    if (is.null(ob) || !identical(th[["type"]], "equidistant")) {
      frm_stop("Prior target not found (", spec_target(s), "): class ",
               "\"delta\" is the distance between the thresholds of an ",
               "ordinal family with threshold = 'equidistant', and this ",
               "model has none",
               if (nzchar(s$dpar)) paste0(" under dpar = \"", s$dpar, "\""),
               call. = FALSE)
    }
    if (nzchar(s$coef) || nzchar(s$nlpar %||% "")) {
      frm_stop("Prior target not found (", spec_target(s), "): class ",
               "\"delta\" is one parameter per threshold vector, so it ",
               "takes no coef or nlpar. group = selects the vector ",
               "of one level of thres(gr = ), and dpar = \"mu<k>\" the ",
               "component of an ordinal mixture", call. = FALSE)
    }
    comp <- ob$comp
    lay <- thres_layout(th[["nthres"]], "equidistant")
    gs <- seq_len(lay$G)
    if (nzchar(s$group)) {
      gs <- if (isTRUE(th[["grouped"]])) {
        which(th[["groups"]] == s$group)
      } else integer(0)
      if (!length(gs)) {
        frm_stop("Prior target not found (", spec_target(s), "): ",
                 if (isTRUE(th[["grouped"]])) {
                   paste0("the thresholds are grouped by thres(gr = ), ",
                          "whose levels are ",
                          paste0("\"", th[["groups"]], "\"",
                                 collapse = ", "))
                 } else {
                   paste0("group = names one threshold vector of a model ",
                          "with grouped thresholds, thres(gr = ), and this ",
                          "model has one threshold vector. Drop group")
                 }, call. = FALSE)
      }
    }
    ordered <- ob$ordered
    lapply(gs, function(g) {
      list(comp = comp, idx = lay$rstart[g] + 1L, dist = s$dist,
           scale = if (ordered) "sd" else "internal", link = NULL,
           offset = NULL, lb = s$lb, ub = s$ub)
    })
  }

  fam_simplex <- unique(vapply(family_simplex_terms(frame, fit$spec),
                               function(tm) tm$name, ""))
  for (s in prior_specificity_order(pl)) {
    # brms names a response by its Stan name (Ya for Y_a), which reaches
    # the response frmtmb names Y_a
    if (nzchar(s$resp %||% "")) {
      s$resp <- resp_from_stan_name(fit$spec, s$resp)
    }
    s <- prior_simplex_check(s, fam_simplex)
    # before the shape gate: in a multivariate model that gate would
    # judge the spelling against whichever response it met first
    no_resp <- resp_missing_refusal(fit$spec, frame, s)
    if (!is.null(no_resp)) frm_stop(no_resp, call. = FALSE)
    if (!is.null(fit$spec[["rescor_nu"]]) && nzchar(s$resp %||% "") &&
          identical(s$dpar %||% "", "nu")) {
      # the same rule as class = "rescor": one parameter across the
      # responses, which brms names nu with no response in it
      frm_stop("A prior on nu (", spec_target(s), ") takes no resp in ",
               "this model: with rescor = TRUE the Student-t responses ",
               "share one nu. Drop resp = \"", s$resp, "\", as brms ",
               "requires", call. = FALSE)
    }
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
    # brms's class of a registry term's parameter (sar()'s lagsar): a
    # density on the natural parameter, carried inward by the map the
    # term's bounds declare; a parameter with no declared map has only
    # its internal spelling
    sc <- spec_spelling(s)$class
    if (isTRUE(s$natural) && sc %in% autocor_seam_classes(frame)) {
      s$class <- sc
      s$dpar <- NULL
      resolve_ac_class(s)
      next
    }
    if (isTRUE(s$natural) &&
          sc %in% autocor_seam_classes(frame, bounded = FALSE)) {
      ap <- frmtmb_autocor_pars(fit)
      frm_stop("class = \"", sc, "\" names a parameter of ",
               ap$term[ap$name == sc][1L], ", whose term declares no ",
               "interval, so a density on it cannot be carried to the ",
               "internal parameter. Write it there: class = \"theta\", ",
               "coef = \"", ap$coef[ap$name == sc][1L], "\"",
               call. = FALSE)
    }
    bad_shape <- dpar_shape_refusal(fit, s)
    if (!is.null(bad_shape)) frm_stop(bad_shape, call. = FALSE)
    ord_th <- if (s$class == "Intercept") ordinal_threshold_entry(s)
    if (!is.null(ord_th)) {
      if (!is.null(s$dist)) {
        for (e in ord_th) {
          key <- nm_of(e$comp, e$idx)
          if (!is.null(e$coef_k)) {
            # one threshold of the vector: its density replaces the
            # class-wide one there and only there. The vector stays one
            # entry, so an ordered map's log-Jacobian enters once
            # (prior_logdens()), as Stan's `ordered` type adds it once
            old <- assigned[[key]]
            dists <- if (is.null(old)) {
              rep(list(NULL), length(e$idx))
            } else if (identical(old$dist$kind, "vec")) {
              old$dist$dists
            } else {
              rep(list(old$dist), length(e$idx))
            }
            dists[e$coef_k] <- list(e$dist)
            e$dist <- list(kind = "vec", dists = dists)
          }
          e$coef_k <- NULL
          claim(e$comp, e$idx)
          assigned[[key]] <- e
        }
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
        # a density on sigma itself under a first-order cov = TRUE ar(),
        # ma() or arma() is brms's density on ITS sigma, the innovation
        # sd, whoever writes it (a default, a user, a sampling table
        # passed back), on frm()'s route as on frm_sample()'s
        innov <- if (isTRUE(s$natural) && identical(tg$dpar, "sigma") &&
                       !is.null(tg$resp)) {
          autocor_innov(frame[["autocor"]][[tg$resp]])
        }
        # and a box on the innovation sd moves with ar, so it is no box
        # on frmtmb's marginal sd; lb = 0 is the same box on either
        if (!is.null(innov) && (isTRUE(s$lb > 0) || !is.na(s$ub))) {
          frm_stop("class = \"sigma\" takes no ",
                   if (isTRUE(s$lb > 0)) "lb" else "ub", " here: under ",
                   frame[["autocor"]][[tg$resp]][["label"]], " brms's ",
                   "sigma, which a prior on class \"sigma\" is about, is ",
                   "the innovation sd, frmtmb's marginal sd times a factor ",
                   "of the correlation parameters, so a bound on it is not ",
                   "a bound on the parameter frmtmb holds. Write the ",
                   "density alone", call. = FALSE)
        }
        if (!is.null(s$dist)) {
          assigned[[key]] <- list(comp = tg$comp, idx = tg$idx,
                                  dist = s$dist, scale = pm$scale,
                                  link = pm$link, offset = pm$offset,
                                  lb = s$lb, ub = s$ub)
          if (!is.null(innov)) assigned[[key]]$innov <- innov
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
      if (nzchar(s$coef %||% "") && !nzchar(s$group)) {
        # brms keys a coefficient's sd row by its group and has no
        # row sd_x: "The following priors do not correspond to any
        # model parameter" (dev/prilabel-log/brms-ref.txt)
        frm_stop("class = \"sd\" with coef = \"", s$coef, "\" needs a ",
                 "`group`: brms lists a coefficient's standard deviation ",
                 "under its grouping factor (sd_<group>__", s$coef,
                 ") and has no row for coef = \"", s$coef, "\" alone. ",
                 "Add group = , as default_prior() lists it",
                 call. = FALSE)
      }
      for (bk in frame[["re_blocks"]]) {
        shared <- sd_coef_shared(fit$spec, frame, bk, s)
        if (length(shared)) {
          frm_stop("class = \"sd\", coef = \"", s$coef, "\" names a ",
                   "column of ", bk[["term_label"]], ", whose ",
                   bk[["covstruct"]], "() structure gives its columns ",
                   "one shared standard deviation, so no standard ",
                   "deviation is that coefficient's alone. Drop the coef ",
                   "to put the prior on the shared one", call. = FALSE)
        }
      }
      # brms's rule: a row exists only for a prefix some block has, and
      # it reaches the standard deviations of that prefix (sd_spec_reach)
      exact <- any(vapply(frame[["re_blocks"]], function(bk) {
        sd_spec_reach(fit$spec, frame, bk, s)$exact
      }, TRUE))
      hit <- FALSE
      for (bk in if (exact) frame[["re_blocks"]]) {
        # a gp() term's sd is brms's class "sdgp", not "sd"
        if (bk[["covstruct"]] %in% c("gp", "hsgp")) next
        sd_i <- cs_entry(bk[["covstruct"]])$sd_idx(bk)
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
                 if (nzchar(s$coef %||% "")) sd_coef_hint(frame, s),
                 sd_prefix_hint(fit$spec, frame),
                 if (any(vapply(frame[["re_blocks"]], function(bk) {
                   bk[["covstruct"]] %in% c("gp", "hsgp")
                 }, NA))) {
                   paste0(". A gp() term's standard deviation is class ",
                          "\"sdgp\", as in brms")
                 }, call. = FALSE)
      }
    } else if (s$class %in% c("sdgp", "lscale")) {
      resolve_gp_class(s)
    } else if (s$class == "cor") {
      hit <- FALSE
      refused <- character(0)
      for (bk in frame[["re_blocks"]]) {
        if (!block_addressed(fit$spec, frame, bk, s)) next
        cs <- bk[["covstruct"]]
        if (cs %in% names(lkj_refusals)) {
          refused <- c(refused, paste0(block_written_label(bk), ": ",
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
        have <- unique(vapply(frame[["re_blocks"]], block_written_label,
                              ""))
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
      # one slot per parameter: the internal name of a parameter that a
      # class reaches would be a second prior on it, and which of the
      # two applied followed the order they were written in
      own <- autocor_class_coefs(frame)
      hit <- if (identical(cmp, "thetaac")) {
        intersect(nms[tgt$idx], own$coef)
      } else character(0)
      if (length(hit)) {
        k <- match(hit[1L], own$coef)
        frm_stop("class = \"theta\", coef = \"", hit[1L], "\" is the ",
                 "internal parameter of ", own$class[k], ", which has a ",
                 "class of its own, as in brms. Write the prior with ",
                 "class = \"", own$class[k], "\"",
                 if (length(fit$spec$responses) > 1L) {
                   paste0(", resp = \"", own$resp[k], "\"")
                 },
                 "; one parameter takes one prior", call. = FALSE)
      }
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
    } else if (s$class == "lncor") {
      resolve_lncor(s)
    } else if (s$class %in% c("simo", "sbhaz", fam_simplex)) {
      # one joint density over the simplex's chart coordinates, in the
      # chart the frame reads them in (dirichlet_logdens())
      tgts <- if (identical(s$class, "simo")) {
        simo_targets(fit, s)
      } else {
        family_simplex_targets(fit, s)
      }
      for (tm in tgts) {
        idx <- seq_len(tm$D - 1L)
        claim(tm$zeta, idx)
        assigned[[nm_of(tm$zeta, idx)]] <-
          list(comp = tm$zeta, idx = idx,
               dist = dirichlet_dist(s$dist$alpha, tm$D, tm$chart),
               scale = "internal", lb = NA, ub = NA)
      }
    } else if (s$class == "delta") {
      for (e in ordinal_delta_entry(s)) {
        key <- nm_of(e$comp, e$idx)
        if (!is.null(s$dist)) {
          claim(e$comp, e$idx)
          assigned[[key]] <- e
        } else if (!is.null(assigned[[key]])) {
          assigned[[key]] <- entry_bounds(assigned[[key]], s)
        }
        # a bound is on delta itself, carried through the log that an
        # ordered family holds it on
        nm_b <- par_template_names(frame[["par_template"]][[e$comp]],
                                   e$comp)[e$idx]
        pm <- list(scale = e$scale)
        if (!is.na(s$lb)) lower[nm_b] <- internal_bound(s$lb, pm, "lb")
        if (!is.na(s$ub)) upper[nm_b] <- internal_bound(s$ub, pm, "ub")
      }
    }
  }
  list(entries = unname(assigned), lower = lower, upper = upper)
}

#' The `mo()` simplexes a class `"simo"` specification names: brms's
#' coef, the term's label with `1` appended (`mox1`, `mox:fb1`), in the
#' linear predictor its `dpar` (or `nlpar`) and `resp` select, `mu`'s
#' when it names none. A coef that names no simplex is refused with the
#' ones the model has, and a concentration vector of the wrong length
#' with brms's words.
#'
#' @noRd
simo_targets <- function(fit, s) {
  out <- list()
  have <- character(0)
  for (lp in fit$frame[["linpreds"]]) {
    if (!length(lp[["mo"]])) next
    rspec <- fit$spec$responses[[lp[["resp"]]]]
    dp <- if (lp[["dpar"]] %in% (rspec$nlpars %||% character(0))) {
      lp[["dpar"]]
    } else {
      lp_prior_dpar(rspec, lp[["dpar"]])
    }
    for (mi in lp[["mo"]]) {
      cf <- paste0(mi$label, "1")
      have <- c(have, paste0(cf, if (nzchar(dp)) paste0(" (dpar ", dp, ")")))
      if (nzchar(s$resp %||% "") && !identical(lp[["resp"]], s$resp)) next
      want <- if (nzchar(s$nlpar %||% "")) s$nlpar else s$dpar %||% ""
      if (!identical(dp, want) || !identical(cf, s$coef)) next
      a <- s$dist$alpha
      if (length(a) != 1L && length(a) != mi$D) {
        frm_stop("Invalid Dirichlet prior. Expected input of length ",
                 mi$D, " for simo coef \"", cf, "\", the ", mi$D,
                 " weights of ", mi$label, "; got ", length(a),
                 call. = FALSE)
      }
      out[[length(out) + 1L]] <- list(zeta = mi$zeta, D = mi$D,
                                      chart = mi[["chart"]])
    }
  }
  if (!length(out)) {
    frm_stop("No mo() simplex matches ", spec_target(s), ". ",
             if (length(have)) {
               paste0("This model's simplexes are coef ",
                      paste0("\"", have, "\"", collapse = ", "))
             } else {
               "This model has no mo() term"
             }, call. = FALSE)
  }
  out
}

#' The simplexes a family declares under the name a specification's
#' class gives (`frmtmb_family(simplex_pars =)`): one per response
#' whose family declares it, narrowed by `resp`, in the chart the
#' family reads it in (the `<name>_chart` addition term). For class
#' `"sbhaz"`, a `cox()` baseline, this is brms's `dirichlet(con_sbhaz)`,
#' whose default is `dirichlet(1)`.
#'
#' @noRd
family_simplex_targets <- function(fit, s) {
  out <- list()
  frame <- fit$frame
  for (tm in family_simplex_terms(frame, fit$spec)) {
    if (!identical(tm$name, s$class)) next
    if (nzchar(s$resp %||% "") && !identical(tm$resp, s$resp)) next
    a <- s$dist$alpha
    if (length(a) != 1L && length(a) != tm$D) {
      frm_stop("Invalid Dirichlet prior. Expected input of length ", tm$D,
               " for class \"", s$class, "\", the ", tm$D, " weights of ",
               "the ", tm$family, "() simplex; got ", length(a),
               call. = FALSE)
    }
    out[[length(out) + 1L]] <- list(
      zeta = tm$zeta, D = tm$D,
      chart = frame[["aterm_values"]][[tm$resp]][[tm$chart_key]])
  }
  if (!length(out)) {
    frm_stop("No simplex matches ", spec_target(s), ": class \"", s$class,
             "\" names a simplex a family declares under that name, and ",
             "no family of this model declares one",
             if (identical(s$class, "sbhaz")) {
               paste0(". \"sbhaz\" is the baseline hazard of cox(), in ",
                      "frmtmb.survival")
             }, call. = FALSE)
  }
  out
}

#' brms's rules for a simplex, applied at the fit to a class that only
#' the fit can tell is a simplex: one a family declares
#' (`frmtmb_family(simplex_pars =)`). `set_prior()` applies the same
#' rules to `"simo"` and `"sbhaz"` when the row is written, so they are
#' not repeated here for those two. A `dirichlet()` on any class that
#' is no simplex of the model is refused here, however it was written.
#'
#' `set_prior()` reads a class that is not frmtmb's own as a
#' distributional parameter's name and stores it as that parameter's
#' `"Intercept"`. Only the fit knows when the name is a declared simplex
#' instead, so the row is given its class back here, and the spec is
#' returned.
#'
#' @noRd
prior_simplex_check <- function(s, fam_simplex) {
  cls <- s$class %||% ""
  dp <- s$dpar %||% ""
  # the name the user wrote, which a refusal must repeat
  said <- if (identical(cls, "Intercept") && nzchar(dp)) dp else cls
  # `class = "Intercept", dpar = "sigma"` as written, not the natural
  # `class = "sigma"` that set_prior() stores under the same two fields
  link_scale <- identical(cls, "Intercept") && nzchar(dp) &&
    !isTRUE(s$natural)
  written <- if (nzchar(dp) && (link_scale || !identical(cls, "Intercept"))) {
    paste0("class = \"", cls, "\", dpar = \"", dp, "\"")
  } else {
    paste0("class = \"", said, "\"")
  }
  is_dir <- identical(s$dist$kind, "dirichlet")
  if (cls %in% c("simo", "sbhaz")) return(s)
  if (said %in% fam_simplex) {
    if (!is_dir) {
      frm_stop("Currently 'dirichlet' is the only valid prior for simplex ",
               "parameters: class = \"", said, "\" is a simplex of this ",
               "model and takes dirichlet(alpha)", call. = FALSE)
    }
    if (link_scale) {
      frm_stop(written, " names the intercept of a linear predictor, and \"",
               dp, "\" is a simplex this model's family declares, which has ",
               "none. Write class = \"", dp, "\"", call. = FALSE)
    }
    if (!is.na(s$lb %||% NA) || !is.na(s$ub %||% NA)) {
      frm_stop("class = \"", said, "\" is a simplex a family declares, ",
               "which takes no lb/ub: its own constraint bounds it",
               call. = FALSE)
    }
    if (nzchar(s$coef %||% "") || nzchar(s$group %||% "")) {
      frm_stop("class = \"", said, "\" is a whole simplex a family ",
               "declares, one per response, so coef and group name ",
               "nothing it can narrow to", call. = FALSE)
    }
    s$class <- said
    s$dpar <- ""
  } else if (is_dir) {
    # brms's class "theta" of a mixture is the Dirichlet on its mixing
    # proportions, a different set of parameters from frmtmb's theta
    mix <- if (startsWith(said, "theta")) {
      paste0(". In brms, class = \"", said, "\" of a mixture() takes a ",
             "dirichlet() on the mixing proportions; ",
             trimws(brms_prior_class_refusal(said)))
    } else ""
    frm_stop("dirichlet() is the density of a simplex; it belongs to ",
             "class = \"simo\" (a mo() simplex), class = \"sbhaz\" (the ",
             "baseline of a cox() model, frmtmb.survival) or a simplex a ",
             "family declares, not ", written, mix, call. = FALSE)
  }
  s
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
has_ordinal_thresholds <- function(fit, resp = NULL) {
  rs <- fit$spec$responses
  rspec <- if (length(rs) == 1L) rs[[1L]] else rs[[resp %||% ""]]
  if (is.null(rspec)) return(FALSE)
  nm <- extra_tpl_name(fit$frame, rspec$resp_name, "tau_raw")
  raw <- fit$frame[["par_template"]][[nm]] %||% numeric(0)
  length(raw) > 0L && identical(rspec$family[["type"]], "ordinal")
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
ordinal_center_offset <- function(frame, rspec, dpar = NULL) {
  for (lp in frame[["linpreds"]]) {
    if (!identical(lp[["resp"]], rspec$resp_name)) next
    if (!lp[["dpar"]] %in% rspec$primary_dpars) next
    # an ordinal mixture's component k is centered on its own mu<k>
    if (!is.null(dpar) && !identical(lp[["dpar"]], dpar)) next
    # bf(center = FALSE): brms's thresholds are then the uncentered ones
    if (isFALSE(lp[["center"]])) return(NULL)
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
  # `sigma ~ 0 + Intercept` is a predictor with one class "b"
  # coefficient in brms, not the scalar a class "sigma" prior names
  if (isFALSE(lp[["center"]])) return(TRUE)
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
    ac_cls <- autocor_seam_classes(frame)
    return(paste0(spelling, " names no distributional parameter of ",
                  "this model. It has ", paste(have, collapse = ", "),
                  if (length(ac_cls)) {
                    paste0(". Its residual correlation terms take ",
                           paste0("class = \"", ac_cls, "\"",
                                  collapse = ", "))
                  },
                  ". frmtmb's own classes are ",
                  paste(frmtmb_prior_classes, collapse = ", ")))
  }
  for (lp in lps) {
    if (!is.null(lp[["equate"]])) {
      return(paste0(spelling, ": ", dp, " is equated to ", lp[["equate"]],
                    " in this model, so it holds no parameter of its own. ",
                    "Address ", lp[["equate"]], " instead"))
    }
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
prior_logdens <- function(x, dist, scale, link = NULL, offset = 0,
                          lc = 0) {
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
    jac <- natural_log_jacobian(x, link)
    x <- link$linkinv(x)
  }
  # a density about c * the parameter (brms's innovation sd, `lc` its
  # log factor): read there, with the factor's Jacobian. `lc` is a
  # function of other parameters, so it is not a constant on the tape
  if (!identical(lc, 0)) {
    x <- x * exp(lc)
    jac <- jac + lc
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

#' The log Jacobian of a natural-scale density's link, `log |mu_eta|`.
#' A logit link's `mu_eta` is `p (1 - p)`, which is exactly 0 once
#' `plogis(eta)` rounds to 1, past `eta` near 37, so its log was -Inf
#' there: a mode on that side (an `rdm()` ndt fraction at 1 under
#' `ndt_group()`, dev/sampleprior-rdm2.R) had an infinite objective and
#' no chain could start. The logit and frmtmb.eam's scaled logit take it
#' as `log(top) - log1p(exp(-eta)) - log1p(exp(eta))`, which is finite
#' for every `eta`, as Stan's `log_inv_logit() + log1m_inv_logit()` is.
#'
#' @noRd
natural_log_jacobian <- function(x, link) {
  nm <- link$name %||% ""
  if (identical(nm, "logit") || startsWith(nm, "scaled_logit")) {
    top <- if (identical(nm, "logit")) 1 else link$linkinv(Inf)
    return(log(top) - RTMB::logspace_add(0 * x, -x) -
             RTMB::logspace_add(0 * x, x))
  }
  log(abs(link$mu_eta(x)))
}

#' The density itself, with no change of variables.
#'
#' @noRd
prior_base_logdens <- function(x, dist) {
  # NaN where `v` is below 0 and exactly 0 at and above it: the support
  # is CLOSED, as Stan's is, so a value a link rounds onto the edge
  # (plogis(40) is 1) or a mode on it keeps the density's value there.
  # The 1e-300 is absorbed by any |v| above 1e-284, which leaves the
  # guard nan_unless_pos() itself except at the edge, where its gradient
  # is still 0 (rtmb-pitfalls item 22: NaN only for subnormal arguments)
  support_guard <- function(v) nan_unless_pos(v + 1e-300)
  if (identical(dist$kind, "vec")) {
    # one density per element of a threshold vector, NULL for an
    # element with none (a coef row without a class-wide row)
    out <- 0
    for (k in seq_along(dist$dists)) {
      if (!is.null(dist$dists[[k]])) {
        out <- out + prior_base_logdens(x[k], dist$dists[[k]])
      }
    }
    return(out)
  }
  switch(dist$kind,
    normal = RTMB::dnorm(x, dist$location, dist$scale, log = TRUE),
    t = RTMB::dt((x - dist$location) / dist$scale, df = dist$df,
                 log = TRUE) - log(dist$scale),
    # each density with a bounded support is NaN outside it on the tape
    # too (support_guard()): (shape - 1) * log(x) at shape 1 is a
    # constant-zero operand the tape folds, and exponential() had no
    # guard at all (dev/rtmb-pitfalls.md, item 22)
    exponential = log(dist$rate) - dist$rate * x + support_guard(x),
    # written through logspace_add rather than log1p(exp(-z)), which
    # overflows in the lower tail exactly where brms's logistic(0, 1) on
    # a mixture proportion is doing its work
    logistic = {
      z <- (x - dist$location) / dist$scale
      -z - log(dist$scale) - 2 * RTMB::logspace_add(0 * z, -z)
    },
    gamma = dist$shape * log(dist$rate) - lgamma(dist$shape) +
      (dist$shape - 1) * log(x) - dist$rate * x + support_guard(x),
    # brms's inv_gamma(shape, scale) is Stan's, so the second argument
    # is the SCALE of the inverse gamma, which is the rate of the gamma
    # on 1/x. Written out rather than through dgamma(1/x) so the tape
    # carries one expression instead of a reciprocal and a correction
    inv_gamma = dist$shape * log(dist$scale) - lgamma(dist$shape) -
      (dist$shape + 1) * log(x) - dist$scale / x + support_guard(x),
    # log1p rather than log(1 - x): a zi or hu near one is exactly where
    # brms's beta(1, 1) is doing its work, and that is where log(1 - x)
    # loses its last digits
    beta = (dist$shape1 - 1) * log(x) +
      (dist$shape2 - 1) * log1p(-x) -
      (lgamma(dist$shape1) + lgamma(dist$shape2) -
         lgamma(dist$shape1 + dist$shape2)) +
      support_guard(x) + support_guard(1 - x),
    # NaN outside the interval, where brms's lb and ub reject the value.
    # The guard has no constant operand, which the tape would fold away
    # (dev/rtmb-pitfalls.md, item 22)
    uniform = -log(dist$upper - dist$lower) +
      support_guard(x - dist$lower) + support_guard(dist$upper - x),
    lognormal = RTMB::dnorm(log(x), dist$location, dist$scale,
                            log = TRUE) - log(x) + support_guard(x),
    # a JOINT density over a whole correlation, so `x` is the block's
    # correlation segment and the value is one number, not one per
    # element (see lkj_logdens)
    lkj = lkj_logdens(x, dist),
    # a JOINT density over a mo() simplex's chart coordinates, one number
    dirichlet = dirichlet_logdens(x, dist)
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

#' @param alpha Dirichlet concentrations, positive: one number for every
#'   weight of the simplex, or one per weight. brms's
#'   `dirichlet(alpha)`, the density of a `mo()` simplex (class
#'   `"simo"`); `alpha = 1` is uniform on the simplex.
#' @rdname frmtmb-priors
#' @export
prior_dirichlet <- function(alpha = 1) {
  if (!is.numeric(alpha) || !length(alpha) || anyNA(alpha) ||
        any(!is.finite(alpha)) || any(alpha <= 0)) {
    frm_stop("prior_dirichlet(alpha =) takes finite positive ",
             "concentrations, one for every weight or one per weight; ",
             "alpha = 1 is uniform on the simplex", call. = FALSE)
  }
  structure(list(kind = "dirichlet", alpha = as.numeric(alpha)),
            class = "frmtmb_prior")
}

#' The concentrations of a `"dirichlet(...)"` string: a number or a
#' vector written with `c()`, as brms takes them. The argument is read
#' as a number, never evaluated as R code.
#'
#' @noRd
parse_dirichlet_arg <- function(prior) {
  inner <- sub("^\\s*dirichlet\\s*[(](.*)[)]\\s*$", "\\1", prior)
  inner <- sub("^\\s*c\\s*[(](.*)[)]\\s*$", "\\1", inner)
  vals <- suppressWarnings(as.numeric(strsplit(inner, ",",
                                               fixed = TRUE)[[1L]]))
  if (identical(inner, prior) || !length(vals) || anyNA(vals)) {
    frm_stop("Cannot parse prior '", prior, "'; a Dirichlet is written ",
             "\"dirichlet(1)\" or \"dirichlet(c(2, 1, 1))\"", call. = FALSE)
  }
  vals
}

# ------------------------------------------------- the Dirichlet ------
#
# brms puts `dirichlet(alpha)` on the weights `w` of a mo() simplex,
#   log p(w) = lgamma(sum a) - sum lgamma(a) + sum (a_k - 1) log w_k,
# a density on the D - 1 free weights. frmtmb holds the simplex in a
# chart of D - 1 coordinates `z` (mo_simplex()), so the prior on `z` is
# that density times the Jacobian of the chart, and each chart has its
# own:
#
# SOFTMAX (frmtmb.sample's sampler): w = exp(c(0, z)) / sum(.). The
#   Jacobian of z -> (w_2, ..., w_D) is prod(w), so
#     log p(z) = lgamma(sum a) - sum lgamma(a) + sum a_k log w_k,
#   which at a = 1 is mo_simplex_nlp()'s lgamma(D) + sum(log w).
#
# STEREOGRAPHIC (the fit's chart): z is the stereographic projection of
#   a point u of the unit sphere in R^D from the pole -(1, ..., 1) /
#   sqrt(D), and w = u^2. The projection is conformal with factor
#   2 / (1 + |z|^2), so the sphere's surface element is
#   (2 / (1 + |z|^2))^(D - 1) dz. On one orthant of the sphere, the
#   map u -> w has dw = 2^(D - 1) prod |u_k| dS: the uniform density on
#   the sphere carried to the simplex is Dirichlet(1/2, ..., 1/2), whose
#   normalizer pi^(D/2) / Gamma(D/2) against the sphere's area
#   2 pi^(D/2) / Gamma(D/2) over 2^D orthants gives the constant. The
#   chart covers the sphere once and the simplex 2^D times, one sheet
#   per orthant, so the density on z that integrates to one is
#     log p(z) = log p(w) + sum log|u_k| + (D - 1) log 4
#                - (D - 1) log(1 + |z|^2) - D log 2
#              = C + sum (2 a_k - 1) log|u_k| - (D - 1) log(1 + |z|^2)
#                + (D - 2) log 2.
#   A weight at a face is u_k = 0, so with a_k < 1/2 the density is
#   unbounded there and a MAP fit has no maximum; frm() refuses it.
#   dev/moparity-dirichlet-jac.R checks both Jacobians against finite
#   differences and the normalization by integration.

#' Refuse a MAP fit whose Dirichlet has no maximum in the fit's chart.
#'
#' In the stereographic chart a weight's factor is `|u_k|^(2 a_k - 1)`
#' (see above), unbounded at the face `u_k = 0` when `a_k < 1/2`: the
#' posterior density of the coordinates has no mode, and the optimizer
#' would run to the face and report the infinite density's location as
#' an estimate. frmtmb.sample's `frm_sample()` samples the same prior
#' in softmax coordinates, where it is proper and bounded.
#'
#' @noRd
check_map_dirichlet <- function(entries) {
  for (e in entries) {
    d <- e$dist
    if (!identical(d$kind, "dirichlet") || !is.null(d$chart)) next
    if (any(d$alpha < 0.5)) {
      frm_stop("A dirichlet() prior with a concentration below 1/2 (",
               prior_dist_string(d), ") has no maximum: its density is ",
               "unbounded where a weight is 0, so a MAP fit has no ",
               "estimate to report. Use concentrations of 1/2 or more, ",
               "or sample the posterior with frmtmb.sample's frm_sample()",
               call. = FALSE)
    }
  }
  invisible(entries)
}

#' The resolved Dirichlet the objective evaluates on one simplex's
#' coordinates: the concentrations, recycled to the simplex's `D`
#' weights, and the chart the coordinates are in.
#'
#' @noRd
dirichlet_dist <- function(alpha, D, chart = NULL) {
  a <- if (length(alpha) == 1L) rep(alpha, D) else alpha
  structure(list(kind = "dirichlet", alpha = a, D = D, chart = chart,
                 lognorm = lgamma(sum(a)) - sum(lgamma(a))),
            class = "frmtmb_prior")
}

#' The Dirichlet log density of one simplex's chart coordinates `z`,
#' Jacobian included (AD-safe). One number for the whole simplex.
#'
#' @noRd
dirichlet_logdens <- function(z, dist) {
  "c" <- RTMB::ADoverload("c")
  a <- dist$alpha
  D <- dist$D
  if (identical(dist$chart, "softmax")) {
    v <- c(0, z)
    lse <- v[1L]
    for (k in seq_along(v)[-1L]) lse <- RTMB::logspace_add(lse, v[k])
    # log w = v - lse, computed without forming w
    return(dist$lognorm + sum(a * (v - lse)))
  }
  P <- -1 / sqrt(D)
  r2 <- sum(z * z)
  u <- (2 * as.vector(mo_sphere_basis(D) %*% z) + (r2 - 1) * P) / (r2 + 1)
  out <- dist$lognorm - (D - 1) * log(1 + r2) + (D - 2) * log(2)
  # a_k = 1/2 has no factor in u_k, and writing 0 * log(0) at a face
  # would be NaN; `a` is a constant, so this is settled at tape build
  for (k in which(2 * a - 1 != 0)) {
    out <- out + (2 * a[k] - 1) * log(abs(u[k]))
  }
  out
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
  if (autocor_is_seam(ac)) {
    hit <- Filter(function(pr) identical(pr$name, cls),
                  autocor_seam_bounded(ac))
    return(if (length(hit)) hit[[1L]]$idx)
  }
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

#' Every residual-correlation parameter that a class reaches (`ar`,
#' `ma`, `cosy`, `cortime`, and a registry term's bounded parameter such
#' as `lagsar`): its internal name (`thetaac_<k>`), the class and the
#' response. Such a parameter has one prior slot, its class, as in
#' brms, so its internal name is not listed in the prior table and a
#' prior on it is refused with the class to use instead. A registry
#' parameter with no declared interval has no class, and keeps its
#' internal spelling.
#'
#' @noRd
autocor_class_coefs <- function(frame) {
  out <- data.frame(coef = character(0), class = character(0),
                    resp = character(0), stringsAsFactors = FALSE)
  for (rs in names(frame[["autocor"]] %||% list())) {
    ac <- frame[["autocor"]][[rs]]
    cls <- c(names(autocor_prior_classes),
             vapply(autocor_seam_bounded(ac), `[[`, "", "name"))
    for (cl in cls) {
      within <- autocor_class_idx(ac, cl)
      if (!length(within)) next
      out <- rbind(out, data.frame(
        coef = paste0("thetaac_", ac[["theta_idx"]][within]), class = cl,
        resp = rs, stringsAsFactors = FALSE))
    }
  }
  out
}

#' The interval a class's natural parameter lives on, `c(lb, ub)`, as
#' the prior table shows it; `NA` where it is not an interval. These are
#' brms's bounds where brms has the term: -1 and 1 on a `cov = TRUE` AR
#' or MA coefficient of order one, which is also the image of frmtmb's
#' map, and none on brms's unconstrained `cov = FALSE` coefficients.
#' `cosy` depends on the route. brms bounds it on `[0, 1]`, and
#' `frm_sample()` samples it there; `frm()` fits it on its whole
#' positive-definite window, `(-1/(d - 1), 1)`, which has no brms
#' counterpart. A higher order is not a box, and brms does not have it;
#' `cortime` is a matrix.
#'
#' @noRd
autocor_class_bounds <- function(ac, cls, route = "fit") {
  none <- c(NA_real_, NA_real_)
  if (autocor_is_seam(ac)) {
    hit <- Filter(function(pr) identical(pr$name, cls),
                  autocor_seam_bounded(ac))
    return(if (length(hit)) c(hit[[1L]]$lower, hit[[1L]]$upper) else none)
  }
  if (autocor_is_cond(ac) || identical(cls, "cortime")) return(none)
  if (identical(cls, "cosy")) {
    if (identical(route, "sample")) return(c(0, 1))
    # the same expression as autocor_trans()'s window, so a table passed
    # back as prior = names its end exactly (ac_bound_theta())
    return(if (isTRUE(ac[["d"]] >= 2)) c(-(1 / (ac[["d"]] - 1)), 1) else none)
  }
  if (length(autocor_class_idx(ac, cls)) == 1L) c(-1, 1) else none
}

#' brms's default prior on each residual-correlation class that frmtmb
#' holds through a transform: flat on the NATURAL parameter, as brms's
#' `(flat)` is. Written for the sampling extension as a uniform density
#' on the natural parameter, which core carries inward with the map's
#' log Jacobian, so the density on the internal parameter is the
#' Jacobian alone. Without it the sampler is flat on the internal
#' parameter instead, which for a first-order AR coefficient is
#' `(1 - rho^2)^(-3/2)` on `rho`, with its mass at the ends.
#'
#' One row per class and response: `lower` and `upper` are the ends of
#' the uniform, which spans the whole image of the map so that it is
#' constant wherever the parameter can go; `lb` and `ub` are a bound
#' brms declares inside that image (`cosy`'s 0), `NA` where there is
#' none; `interval` is `FALSE` for a covariance-form AR or MA of order
#' `k > 1`, which brms does not have. Its image is the stationary (or
#' invertible) region, which is not a box but lies inside the one
#' `|phi_j| <= choose(k, j)` gives, so a uniform on
#' `(-choose(k, k %/% 2), choose(k, k %/% 2))` is flat on it.
#'
#' A `cov = FALSE` term has no row: its coefficients are unconstrained
#' and are the internal parameters, so flat there is brms's flat
#' already. Neither has `cortime`, whose default is `lkj(1)`.
#'
#' @noRd
autocor_flat_priors <- function(fit) {
  out <- data.frame(resp = character(0), class = character(0),
                    lower = numeric(0), upper = numeric(0),
                    lb = numeric(0), ub = numeric(0),
                    interval = logical(0), stringsAsFactors = FALSE)
  # a registry term's bounds live in its package, which a saved fit
  # may come back without; refused there with the package to load
  if (length(fit$frame[["autocor"]])) require_fit_packages(fit)
  for (rs in names(fit$frame[["autocor"]] %||% list())) {
    ac <- fit$frame[["autocor"]][[rs]]
    if (autocor_is_cond(ac)) next
    cls <- c(setdiff(names(autocor_prior_classes), "cortime"),
             vapply(autocor_seam_bounded(ac), `[[`, "", "name"))
    for (cl in cls) {
      within <- autocor_class_idx(ac, cl)
      if (!length(within)) next
      img <- autocor_class_bounds(ac, cl, "fit")
      bd <- autocor_class_bounds(ac, cl, "sample")
      interval <- !anyNA(img)
      if (!interval) {
        if (!cl %in% c("ar", "ma")) next
        k <- length(within)
        img <- c(-1, 1) * choose(k, k %/% 2)
        bd <- c(NA_real_, NA_real_)
      }
      out <- rbind(out, data.frame(
        resp = rs, class = cl, lower = img[1L], upper = img[2L],
        lb = if (isTRUE(bd[1L] > img[1L])) bd[1L] else NA_real_,
        ub = if (isTRUE(bd[2L] < img[2L])) bd[2L] else NA_real_,
        interval = interval, stringsAsFactors = FALSE))
    }
  }
  out
}

#' A chart for each bounded scalar residual-correlation parameter, so
#' that a sampler can draw it through ONE logistic onto its natural
#' interval `(L, U)`, as Stan's own bounded transform does in brms. The
#' alternative, the sampler's log bound on the internal parameter
#' stacked on the term's own map, gave 5 or 6 divergent transitions per
#' 4 chains of 2000 on a cosy near its bound at 0, where brms has none
#' and this chart has none (dev/acprior-cosydiv.R).
#'
#' `lower` and `upper` are the named internal bounds a prior resolved to
#' (`resolve_prior_input()$lower`). One term per first-order parameter
#' with a finite bound on at least one side, whose map is not the
#' identity (a `cov = FALSE` coefficient is the sampled parameter, and
#' the sampler's own bound already is brms's transform there). Each term
#' has `comp`, `idx` (within `comp`), `name`, `class`, `resp`, the ends
#' `L` and `U` on the natural scale, and three maps of the chart
#' coordinate `v`: `to_internal(v)` (tape-safe), `logjac(v)`, the log of
#' `|d internal / d v|` (tape-safe), and `from_internal(th)` (numeric,
#' for a start). A start outside `(L, U)` is moved just inside it.
#' `to_natural(th)` gives the natural value of an internal one and
#' `image` the ends of the map's own interval, for a message about a
#' start outside the bound.
#'
#' @noRd
autocor_bound_chart <- function(fit, lower, upper) {
  lower <- unlist(lower %||% list())
  upper <- unlist(upper %||% list())
  get_b <- function(b, nm, fill) if (nm %in% names(b)) b[[nm]] else fill
  out <- list()
  if (length(fit$frame[["autocor"]])) require_fit_packages(fit)
  for (rs in names(fit$frame[["autocor"]] %||% list())) {
    ac <- fit$frame[["autocor"]][[rs]]
    cls <- c(setdiff(names(autocor_prior_classes), "cortime"),
             vapply(autocor_seam_bounded(ac), `[[`, "", "name"))
    for (cl in cls) {
      within <- autocor_class_idx(ac, cl)
      if (length(within) != 1L) next
      tr <- autocor_trans(ac, cl)
      if (is.null(tr) || identical(tr$map, "identity")) next
      nm <- paste0("thetaac_", ac[["theta_idx"]][within])
      lo <- get_b(lower, nm, -Inf)
      hi <- get_b(upper, nm, Inf)
      if (!is.finite(lo) && !is.finite(hi)) next
      img <- unname(switch(tr$map, bounded = c(tr$lo, tr$hi),
                           cosy = c(-tr$a, 1), c(-1, 1)))
      L <- if (is.finite(lo)) as.numeric(ac_trans_value(lo, tr)) else img[1L]
      U <- if (is.finite(hi)) as.numeric(ac_trans_value(hi, tr)) else img[2L]
      # the internal bound is the image of the natural one the prior
      # wrote, and its round trip carries a few ulps (5.6e-17 for an
      # lb = 0), which would let a draw sit a hair past it
      snap <- function(b) {
        r <- round(b, 10)
        if (abs(b - r) <= 64 * .Machine$double.eps * max(1, abs(b))) r else b
      }
      L <- snap(unname(L))
      U <- snap(unname(U))
      out[[length(out) + 1L]] <- ac_chart_term(
        nm, cl, rs, ac[["theta_idx"]][within], L, U, tr)
    }
  }
  out
}

#' brms's `sigma` under a first-order `cov = TRUE` `ar()`, `ma()` or
#' `arma()` is the INNOVATION sd: its `cholesky_cor_ar1()`,
#' `cholesky_cor_ma1()` and `cholesky_cor_arma1()` are covariances with
#' diagonal `1 / (1 - ar^2)`, `1 + ma^2` and
#' `(1 + ma^2 + 2 ar ma) / (1 - ar^2)`, where frmtmb's are correlations
#' (dev/acprior-armacheck.R: the same correlation to 2.2e-16). So
#' brms's sigma is frmtmb's `sigma * c`, `c` one over the square root of
#' that diagonal. `autocor_innov(ac)` says where `c` comes from, `NULL`
#' for every other block (cosy and unstr are correlations in brms too).
#'
#' @noRd
autocor_innov <- function(ac) {
  if (is.null(ac) || autocor_is_seam(ac) || autocor_is_cond(ac)) {
    return(NULL)
  }
  st <- ac[["struct"]]
  p <- if (st %in% c("ar", "arma")) ac[["p"]] else 0L
  q <- if (st %in% c("ma", "arma")) ac[["q"]] else 0L
  if (!st %in% c("ar", "ma", "arma") || p > 1L || q > 1L || p + q < 1L) {
    return(NULL)
  }
  list(comp = "thetaac", idx = ac[["theta_idx"]], p = p, q = q)
}

#' `log c` of `autocor_innov()` at the block's internal parameters,
#' tape-safe: `-0.5 log` of brms's diagonal.
#'
#' @noRd
innov_logc <- function(th, innov) {
  # 1 - ar^2 is 1 / (1 + th^2) exactly, which keeps its digits where ar
  # rounds to 1
  if (!innov$q) return(-0.5 * log(1 + th[1] * th[1]))
  ma <- autocor_pacf(th[innov$p + 1L])
  if (!innov$p) return(-0.5 * log(1 + ma * ma))
  ar <- autocor_pacf(th[1])
  -0.5 * log(1 + th[1] * th[1]) - 0.5 * log(1 + ma * ma + 2 * ar * ma)
}

#' The responses whose sigma brms reads as an innovation sd
#' (`autocor_innov()`), for the sampling extension's default on sigma.
#'
#' @noRd
autocor_innov_sd <- function(fit) {
  acs <- fit$frame[["autocor"]] %||% list()
  names(acs)[vapply(acs, function(ac) !is.null(autocor_innov(ac)), NA)]
}

#' One term of `autocor_bound_chart()`.
#'
#' @noRd
ac_chart_term <- function(nm, cl, rs, idx, L, U, tr) {
  # natural -> internal, the inverse of ac_trans_value() for a scalar map
  inv <- function(r) {
    if (identical(tr$map, "bounded")) {
      q <- (r - tr$lo) / (tr$hi - tr$lo)
      return(log(q) - log(1 - q))
    }
    if (identical(tr$map, "cosy")) {
      q <- (r + tr$a) / (1 + tr$a)
      return(log(q) - log(1 - q))
    }
    r / sqrt(1 - r * r)
  }
  list(
    comp = "thetaac", idx = idx, name = nm, class = cl, resp = rs,
    L = L, U = U,
    to_internal = function(v) inv(L + (U - L) / (1 + exp(-v))),
    # d internal / d v = (d natural / d v) / (d natural / d internal)
    logjac = function(v) {
      s <- 1 / (1 + exp(-v))
      th <- inv(L + (U - L) * s)
      # log s and log(1 - s) through logspace_add, finite at any v
      log(U - L) - RTMB::logspace_add(0 * v, -v) -
        RTMB::logspace_add(0 * v, v) - ac_trans_logjac(th, tr)
    },
    from_internal = function(th) {
      q <- (as.numeric(ac_trans_value(th, tr)) - L) / (U - L)
      stats::qlogis(min(max(q, 1e-8), 1 - 1e-8))
    },
    to_natural = function(th) as.numeric(ac_trans_value(th, tr)),
    # the ends of the map's own image, which a bound narrows
    image = unname(switch(tr$map, bounded = c(tr$lo, tr$hi),
                          cosy = c(-tr$a, 1), c(-1, 1))))
}

#' The map one autocorrelation class carries its prior through. `cosy`
#' is a single logistic onto the positive-definite window; the ARMA
#' coefficients are a partial-autocorrelation transform composed with
#' the Levinson-Durbin recursion; `cortime` is the row-normalized
#' Cholesky the LKJ density already knows, so it carries no `trans`.
#'
#' @noRd
autocor_trans <- function(ac, cls) {
  if (autocor_is_seam(ac)) {
    hit <- Filter(function(pr) identical(pr$name, cls),
                  autocor_seam_bounded(ac))[[1L]]
    return(list(map = "bounded", lo = hit$lower, hi = hit$upper))
  }
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
  if (identical(tr$map, "bounded")) {
    return(tr$lo + (tr$hi - tr$lo) / (1 + exp(-th[1])))
  }
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
  if (identical(tr$map, "bounded")) {
    # d rho / d theta = (hi - lo) p (1 - p), the same form as cosy's
    s <- 1 / (1 + exp(-th[1]))
    return(log(tr$hi - tr$lo) + log(s) + log(1 - s))
  }
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
#' A bound AT the end of the parameter's own interval holds nothing back,
#' since the map never reaches that end, so it is no bound at all
#' (`-Inf` or `Inf`). That is what the prior table lists there, as brms's
#' does, and a table passed back as `prior =` must fit.
#'
#' @noRd
ac_bound_theta <- function(v, tr, ac, cls, what) {
  if (identical(tr$map, "identity")) return(v)
  ends <- switch(tr$map, bounded = c(tr$lo, tr$hi), cosy = c(-tr$a, 1),
                 c(-1, 1))
  # an end is computed (1 / an eigenvalue for sar()) and carries a few
  # ulps, so a bound typed at it, ub = 1 for a row-standardized W, is
  # that end within the snap autocor_bound_chart() applies
  at_end <- function(e) {
    abs(v - e) <= 64 * .Machine$double.eps * max(1, abs(e))
  }
  if (identical(what, "lb") && at_end(ends[1L])) return(-Inf)
  if (identical(what, "ub") && at_end(ends[2L])) return(Inf)
  if (identical(tr$map, "bounded")) {
    if (v <= tr$lo || v >= tr$hi) {
      frm_stop("class = \"", cls, "\" ", what, " = ", v, " is outside ",
               "the interval of ", ac[["label"]], ", (",
               format(tr$lo, digits = 15), ", ",
               format(tr$hi, digits = 15), ")", call. = FALSE)
    }
    return(stats::qlogis((v - tr$lo) / (tr$hi - tr$lo)))
  }
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
    if (is.null(v)) {
      v <- paste0(cp, "_", seq_along(tpl[[cp]]), recycle0 = TRUE)
    }
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
    if (identical(pr$kind, "dirichlet")) {
      # the same reason: a density over a whole simplex, in its chart
      frm_stop("prior_dirichlet() addresses a whole mo() simplex, so it ",
               "cannot be given by parameter name; write ",
               "set_prior(\"", prior_dist_string(pr), "\", class = ",
               "\"simo\", coef = ...) with the coef get_prior() lists",
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
  # `shift` is a constant, the change of scale between the parameter
  # held and the one the density is about (a gp() length scale on
  # brms's scaled inputs)
  sh <- e[["shift"]] %||% 0
  o <- e$offset
  if (is.null(o)) return(sh)
  sum(pars[[o$comp]][o$idx] * o$w) + sh
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
      lc <- if (!is.null(e$innov)) {
        innov_logc(pars[[e$innov$comp]][e$innov$idx], e$innov)
      } else 0
      nlp <- nlp - sum(prior_logdens(pars[[e$comp]][e$idx], e$dist,
                                     e$scale, e$link,
                                     entry_offset(e, pars), lc = lc))
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

#' A prior's `resp` as frmtmb names the response: brms's Stan name of a
#' response (`brms_stan_name()`, `Ya` for `Y_a`) maps to it, and a name
#' that is already frmtmb's, or that matches nothing, is kept.
#'
#' @noRd
resp_from_stan_name <- function(spec, resp) {
  rn <- names(spec$responses)
  if (resp %in% rn) return(resp)
  hit <- match(resp, brms_stan_name(rn))
  if (is.na(hit)) resp else rn[hit]
}
