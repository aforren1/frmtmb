#' Define a model family
#'
#' Constructs a family object for use with [frm()]. The log-density
#' function must be vectorized and AD-compatible: it is evaluated on RTMB
#' 'advector' objects during taping, so it must use RTMB-overloaded
#' operations (RTMB and RTMBdist `d*` functions, plain arithmetic) and must
#' not branch on parameter values.
#'
#' @param family Character name of the family.
#' @param dpars Character vector of distributional parameter names. The
#'   first entry must be `"mu"`.
#' @param links Named list mapping each dpar to a link name or a link
#'   object. [frmtmb-links] lists the names and says what a link object
#'   has to carry.
#' @param lpdf Function `(y, dpars, aterms)` returning the vectorized
#'   log-density. `dpars` is a named list of advector vectors; `aterms` is a
#'   named list of numeric addition-term values (for example `trials`).
#'   Read `dpars` by name, never by position: it also carries reserved
#'   entries that are not distributional parameters. Take `log(mu)`,
#'   `log(1 - mu)` and `1 - mu` from [frmtmb-robust-dpars] rather than
#'   writing them out, because an inverse link saturates and the plain
#'   arithmetic returns `NaN` for the value and the gradient alike in
#'   the tail.
#' @param valid_y Optional function `(y, aterms)` that signals an error for
#'   invalid responses. Called once at assembly time.
#' @param init_dpars Optional named list of functions `(y, aterms)` giving a
#'   response-scale starting value per dpar (applied to the intercept
#'   through the link).
#'
#'   A function that declares a THIRD argument, `(y, aterms, resid)`,
#'   is given the response with its `mu` predictor taken out: any
#'   `offset()` subtracted, then least squares on that predictor's own
#'   design. A start that depends on the SHAPE of the residual rather
#'   than on the response cannot be read off `y`: a covariate can give
#'   the raw response the opposite skew, whether it enters as a column
#'   or as an offset (see `skew_normal()`). When `mu` is an intercept
#'   alone and carries no offset there is nothing to take out and
#'   `resid` is the response itself, so such a model starts exactly
#'   where it did. A two-argument function is called with two
#'   arguments, so nothing already written has to change.
#'
#'   The value an initializer returns names the PREDICTOR, and it is
#'   placed in whichever coefficients that design uses to carry one: an
#'   intercept when there is one, and otherwise the least-squares
#'   solution, which puts the value in every cell of `~ 0 + g`. A
#'   design that was skipped for having no intercept used to start
#'   every coefficient at zero.
#' @param type One of `"continuous"`, `"discrete"`, `"ordinal"`,
#'   `"categorical"`. The last two say the modelled response is a
#'   distribution over `1..K` categories rather than a number, so
#'   [fitted()] returns an `n x K` probability matrix; `"ordinal"`
#'   shares one latent predictor across the categories and
#'   `"categorical"` gives each non-reference category its own.
#' @param post Named list of numeric helper functions used by
#'   [fitted()], [predict()] and [residuals()]: `mean_fn(dpars,
#'   aterms)` (the response mean), `var_fn(dpars, aterms)` (for pearson
#'   residuals) and `dev_fn(y, dpars, aterms)` (the unit deviance
#'   `2 * (loglik of the saturated fit - loglik at the fitted value)`,
#'   for `residuals(type = "deviance")`). A family that omits one is
#'   refused by the method that needs it.
#'
#'   `post$fit_check(fit, resp)` is different in kind: it is run once,
#'   when a fit FINISHES, and its return value is discarded. It is where
#'   a family says something about where the optimizer landed, which
#'   nothing else can: `logLik()` reads the optimizer's own value, so a
#'   family whose likelihood is floored or degenerate in some region had
#'   no way to report it. Warn from it rather than stopping; a hook that
#'   throws is caught, reported as a warning naming the family, and the
#'   fit is returned regardless.
#'
#'   `post$stationary` is a declaration rather than a function: a named
#'   list, one entry per dpar, each `list(at =, tol =, from =)`. It says
#'   that the likelihood has a stationary point where that dpar's
#'   PREDICTOR equals `at` on the LINK scale at every row. A fit whose
#'   fitted predictor is within `tol` of that everywhere is refitted
#'   from each value in `from`, with the coefficients solved so the
#'   predictor takes that value, and whichever optimum is best is kept.
#'
#'   Reading the predictor rather than the coefficients is what makes
#'   this work for a design without an intercept: `alpha ~ 0 + g` sits
#'   on the same point as `alpha ~ g` with no coefficient named
#'   `"(Intercept)"` to read or to write a restart into. Declaring a
#'   stationary point is ALSO what earns a dpar a starting value on an
#'   intercept-less design at all: a dpar that declares none keeps the
#'   plain rule, an intercept or nothing, because a least-squares start
#'   is right on the predictor scale and can be ruinous on the
#'   parameter scale. Non-finite entries in `from` are dropped, so
#'   `from = c(Inf, -2)` restarts from one side only. Declare one
#'   only for a point the likelihood has at EVERY sample, such as
#'   `alpha = 0` for the skew normal, where the information is singular
#'   and no convergence test can separate a stall from a maximum. The
#'   extra fits only ever run when the declared point is reached.
#'
#'   `mixture()` carries a component's declaration through under the
#'   dpar name the mixture gives it (`alpha1`, `alpha2`, ...).
#' @param sim Optional numeric simulator `(dpars, aterms, n)` returning `n`
#'   response draws; used by [simulate()], `posterior_predict()` and
#'   [frm_simulate()]. It is stateless and rowwise: it sees the
#'   distributional parameters and nothing else. A family whose extra
#'   parameters (`extra_pars`) enter the draw declares a fourth
#'   argument `extra` instead.
#' @param sim_ctx Optional structured simulator `(ctx)` for a family
#'   whose draws are not rowwise (see Structured simulators). It takes
#'   precedence over `sim`.
#' @param sim_refusal Optional one-sentence reason why the family has no
#'   simulator, appended to the refusal each entry point raises. Use it
#'   when the omission is a decision rather than a gap.
#' @param primary_dpars Which dpars receive the main model formula
#'   (default `"mu"`). Families with several location predictors (for
#'   example multinomial's per-category `mu2`, `mu3`, ...) list them all;
#'   these live in the `beta` parameter vector and are integrated out
#'   under REML.
#' @param lcdf Optional vectorized AD log-safe CDF `(q, dpars, aterms)`
#'   returning probabilities; enables `cens()` and `trunc()` addition
#'   terms. It is the only thing either term asks of a family, whatever
#'   the family's `type`: a discrete family that supplies one is
#'   censored under the inclusive convention (see Censoring a discrete
#'   response), which reads a lower bound as `F(q - 1)` and so calls
#'   this at one below a recorded value.
#' @param lccdf Optional vectorized AD LOG SURVIVOR function
#'   `(q, dpars, aterms)` returning `log(1 - F(q))` directly. A family
#'   that declares it scores a RIGHT-censored row from it instead of
#'   from `log(1 - F)`, which cannot be accurate once `F` rounds to one
#'   (see Right censoring and the representable tail). It is optional
#'   and independent of `lcdf`: a family that supplies only `lccdf`
#'   accepts right censoring and refuses left censoring, interval
#'   censoring and `trunc()`, each by name.
#' @param required_aterms The addition-term values the density cannot do
#'   without, named as they reach `aterms`: `"vint1"`, `"vreal2"`,
#'   `"trials"`. A character vector names the terms it needs ALL of. A
#'   LIST adds the alternative: an element of length one is a term that
#'   must be there, and an element of length more than one is a set of
#'   spellings, ANY one of which will do, so
#'   `required_aterms = list(c("dec", "vint1"), "vreal1")` reads as one
#'   of `dec` or `vint1`, and `vreal1`. A family that reads the same
#'   datum from either of two terms declares it that way instead of
#'   hand-rolling the refusal.
#'
#'   Frame assembly refuses a model that leaves a requirement unmet,
#'   naming the family, the terms and the spelling that supplies them;
#'   for a set of alternatives it names all of them and writes the first
#'   into the example formula. Without the declaration an absent term
#'   reaches the density as `NULL`, arithmetic on it gives `numeric(0)`,
#'   and the log-likelihood becomes a sum over nothing: a fit that
#'   returns, with a log-likelihood of zero. Declare every per-row datum
#'   the density indexes.
#' @param accepts_aterms The addition terms this family reads or lets the
#'   core act on, named as a formula writes them and without
#'   parentheses: `c("weights", "trials", "cens")`. Frame assembly
#'   refuses any other term on the response, by name, and lists the ones
#'   the family takes. `character(0)` declares a family that takes none.
#'   `NULL`, the default, accepts every registered term, which is what a
#'   family written before this argument existed keeps.
#'
#'   `required_aterms` is a conjunction of what the density cannot do
#'   without; this is the complementary allow-list, and the two are read
#'   together, so a required term need not be repeated here. Without a
#'   declaration an unread term is parsed, stored on the fit and
#'   silently ignored: `wiener()` accepted a `vint()` it cannot use, and
#'   `lba()` accepted a `dec()`, both giving a fit bit-identical to the
#'   one without the term.
#'
#'   Naming `"se"` here is more than an allow-list entry: it is how a
#'   family OPTS IN to `se()`. That term is the one whose entire effect
#'   is inside the density - the core hands over `aterms[["se"]]` and
#'   maps out the residual scale it replaces, and does nothing else with
#'   it - so a family that does not declare it is refused the term
#'   rather than given one it would ignore. `NULL`, which accepts every
#'   other term, declares nothing and so does not open this one. Read
#'   `aterms[["se"]]` as the known standard deviation, and
#'   `aterms[["se_sigma"]]` to honor `se(x, sigma = TRUE)`, which asks
#'   for the known and estimated scales in quadrature.
#' @param se_dpar The dpar that a known standard error replaces, named
#'   so that the core can map it out: `se_dpar = "tau"` maps out `tau`
#'   exactly as the convention maps out `sigma`. `NA` declares that
#'   `se()` replaces NO dpar, which is the shape of a family whose whole
#'   scale IS the known standard error. `NULL`, the default, reads the
#'   convention: the dpar named `sigma`, if the family has one.
#'
#'   Only a family that declares `se()` reaches this. `se()` without
#'   `sigma = TRUE` says the residual scale is known, so the dpar it
#'   replaces has to stop being estimated. A dpar the density never
#'   reads is a flat direction and a NaN standard error, so a declaring
#'   family with no `sigma` and another free dpar is refused: the core
#'   cannot tell a second SCALE from a genuine SHAPE. This argument is
#'   how the family says which one it has. `se_dpar = NA` is a promise
#'   that every remaining dpar is read alongside the known standard
#'   error, the way a skew or a tail index is.
#' @param exclusive_aterms Sets of addition-term values that say the SAME
#'   thing to the density, so that at most one of each set may be
#'   supplied. Named as `required_aterms` names them, in values rather
#'   than in terms: `exclusive_aterms = list(c("dec", "vint1"))`. A
#'   character vector is one such set; a list is several. The FIRST name
#'   of a set is the spelling the density reads, and the refusal tells
#'   the user to keep it.
#'
#'   This is what an allow-list cannot express, because both spellings
#'   are legitimately on it. `wiener()` reads its boundary indicator from
#'   `dec()` and falls back to `vint1`, so `rt | dec(u) + vint(1 - u)`
#'   passed every guard and fitted with a log-likelihood bit-identical to
#'   the `dec()`-only model: the user said one thing twice, and
#'   inconsistently, and nothing complained. Declare the alternatives of
#'   an any-of `required_aterms` group here when the density reads only
#'   one of them; leave them undeclared when it reads both, as `gddm()`
#'   does, where `dec()` and `vint1` carry different data.
#'
#'   Read WITH `required_aterms`, not instead of it: an any-of group in
#'   one and the same set in the other together mean "exactly one".
#'   Declaring a set that holds two values `required_aterms` demands
#'   TOGETHER is refused here, because no model could then be fitted.
#' @param family_finalize Optional function `(fam, y, aterms)` returning
#'   a family. It runs once at frame assembly, after the response is
#'   coerced and validated and before any link is used, and whatever it
#'   returns is the family the rest of the fit sees. It is how a family
#'   derives a link bound, a default, or an extra slot FROM the data
#'   instead of asking the user for a quantity the framework already
#'   holds (see Deriving a family from the data).
#' @param extra_pars Optional function `(y, aterms)` returning a named
#'   list of numeric starting vectors for family-level parameters outside
#'   the dpar system (for example ordinal thresholds). They join the
#'   parameter template under their own names and reach `lpdf` as its
#'   fourth argument.
#' @param drop_intercept If `TRUE`, the intercept column is removed from
#'   the main formula's design matrix (ordinal families: thresholds take
#'   its place).
#' @param structure Optional [frmtmb_structure()] for a family whose
#'   likelihood does not factorize over rows (a group-level [mixture()],
#'   a hidden Markov chain). It carries the non-rowwise log-likelihood,
#'   the frame block that likelihood reads, and the capability flags
#'   that say which post-fit methods the family can answer. `lpdf` stays
#'   required even then, for the rowwise contract, and may be a stub
#'   that refuses.
#' @return An object of class `frmtmb_family`. `$` on it also answers
#'   brms's link fields, computed from `links` when read: `link`, the
#'   name of the link for the mean (for an ordinal family, the
#'   distribution function its thresholds are read through;
#'   `"logit"` for a categorical or multinomial family, whose categories
#'   share one joint link; `"identity"` for a family with no `mu`, such
#'   as a [mixture()]), `linkfun` and `linkinv` for that link, and
#'   `link_<dpar>`, the link name of every other distributional
#'   parameter, as in `beta_binomial()$link_phi`. Because they are
#'   read from `links` each time, a `family_finalize()` that replaces a
#'   link, or any other edit, cannot leave them stale.
#'
#'   `$` does not partial-match on a family object. `fam$lpd` is an
#'   error naming the field it would have matched, rather than
#'   silently returning it. Use the full name.
#' @seealso [frmtmb-robust-dpars] for the accessors a density uses to
#'   stay exact where an inverse link saturates,
#'   [frmtmb_structure()] for a likelihood that does not
#'   factorize over rows, [frmtmb_register_aterm()] for giving the
#'   family's per-row data a name of its own instead of `vint()`,
#'   [frmtmb_register_compat()] for telling [frm_compat()] what the
#'   family does and does not combine with, and
#'   [frmtmb-extension-api] for the accessors a family outside frmtmb
#'   may use after a fit
#' @section Structured simulators:
#' Some families cannot draw a response one row at a time: a group-level
#' [mixture()] draws one class per group, a [mixture_mvn()] draw needs
#' the class covariances, which are family-level extras rather than
#' dpars, and a hidden Markov family walks a chain per sequence. Those
#' families supply `sim_ctx(ctx)` instead of `sim(dpars, aterms, n)`.
#'
#' `ctx` is a list with `fit` (any object carrying `spec`, `frame` and
#' `estimates` - a fitted model, one posterior draw, or the de novo
#' shim), `family`, `rspec`, `resp`, `dpars` (the evaluated numeric
#' distributional parameters), `aterms`, `n`, `extra` (the family-level
#' extra parameters) and the frame structures `autocor` and `block`
#' (the structured family's own data; see [frmtmb_structure()]).
#' Read its fields with `[[ ]]`.
#'
#' The same `sim_ctx()` serves [simulate()], `posterior_predict()` and
#' [frm_simulate()]. Because a structured draw covers whole sequences or
#' groups, `trunc()` rejection and `newdata` cannot apply to it and are
#' refused.
#' @examples
#' # a custom family is a plain R log-density over taped parameters
#' dd <- data.frame(y = rbinom(100, 5, 0.4),
#'                  size = 5, x = rnorm(100))
#' fam <- custom_family(
#'   "vbinom", dpars = "mu", links = list(mu = "logit"),
#'   lpdf = function(y, dpars, aterms) {
#'     RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
#'   },
#'   # the density indexes vint1, so a model without it is refused
#'   # rather than fitted against a zero-length log-likelihood
#'   required_aterms = "vint1",
#'   type = "discrete"
#' )
#' fit <- frm(bf(y | vint(size) ~ x) + fam, data = dd)
#' fixef(fit)
#' @srrstats {G2.0,G2.1} The family contract is asserted on both length
#'   and type before anything is built: `family` must be a length-one
#'   character vector, `dpars` a character vector of length at least one,
#'   and `lpdf` a function. A family supplied as a `stats::family()`
#'   object, a family constructor, or a name is dispatched to a single
#'   internal representation, and an unrecognized value errors naming the
#'   supported families.
#' @srrstats {RE4.12} The transform used on each linear predictor and its
#'   inverse are both available. Every distributional parameter carries a
#'   link with `linkfun`, `linkinv`, and `mu_eta` (the derivative), and
#'   `links` selects them per parameter. They are reachable through the
#'   fit with `family()`, through `insight::link_function()` and
#'   `insight::link_inverse()`, and are applied by `frm_linpred(type =)` to
#'   move between the link and response scales. Link functions are
#'   written out over plain arithmetic rather than taken from
#'   `stats::make.link()`, because the latter clamps at C level in ways
#'   the AD tape cannot see.
#'
#' @section Slot call order:
#' The order the slots run in is part of the contract, because a family
#' that derives anything from the data depends on it. Measured on an
#' instrumented family, one fit of one response:
#'
#' 1. `valid_y(y, aterms)`, once, at frame assembly, with the response
#'    coerced and the addition terms evaluated.
#' 2. `family_finalize(fam, y, aterms)`, once, immediately after, and
#'    still before any link function is called.
#' 3. `aterm_data()` and `extra_pars(y, aterms)`, once each, at
#'    assembly.
#' 4. `init_dpars[[dpar]](y, aterms)`, once per dpar, when the starting
#'    values are built; each value goes straight through that dpar's
#'    `linkfun`.
#' 5. `linkinv` and `lpdf`, on the tape, from then on.
#'
#' `post$mean_fn`, `post$var_fn`, `post$dev_fn` and `sim` are never
#' called by [frm()]. They run only when a post-fit method asks for
#' them.
#'
#' So a link may read anything `valid_y` or `family_finalize` computed,
#' and neither of those may read anything a link produced.
#'
#' @section Deriving a family from the data:
#' Some families are not fully determined until the response is in hand.
#' A shifted family whose density is zero below a non-decision time
#' wants a link bounded above by `min(y)`, so that the constraint is
#' structural rather than left to the optimizer. The bound is a property
#' of the data, and the family object is built before [frm()] sees any.
#'
#' `family_finalize` closes that gap. It receives the family, the
#' coerced response and the evaluated addition terms, and returns the
#' family the fit will use:
#'
#' ```r
#' shifted <- frmtmb_family(
#'   "shifted", dpars = c("mu", "ndt"),
#'   links = list(mu = "log", ndt = "log"),
#'   lpdf = function(y, dpars, aterms) { ... },
#'   family_finalize = function(fam, y, aterms) {
#'     # a logit onto (0, min(y)): ndt cannot reach the smallest
#'     # observation, whatever the optimizer tries
#'     ub <- min(y)
#'     fam$links$ndt <- list(
#'       name    = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
#'       linkfun = function(mu) log(mu / (ub - mu)),
#'       linkinv = function(eta) ub / (1 + exp(-eta)),
#'       mu_eta  = function(eta) {
#'         p <- 1 / (1 + exp(-eta))
#'         ub * p * (1 - p)
#'       }
#'     )
#'     fam
#'   }
#' )
#' ```
#'
#' Replacing a link this way replaces it everywhere: the starting
#' values, the tape, and every post-fit method read the finalized
#' family. The alternative an extension reached for before this slot
#' existed was to have `valid_y` write the bound into an environment the
#' link closures read at call time, which works only for as long as the
#' call order happens to hold and leaves the family object lying about
#' what it is.
#'
#' @section Right censoring and the representable tail:
#' Core forms a right-censored row's contribution as `log(Fub - F(y))`,
#' which without truncation is `log(1 - F)`. A double cannot represent
#' the complement of a probability that has rounded to one, so past that
#' point the contribution is not merely inaccurate: it is CONSTANT, and
#' its gradient is exactly zero. An optimizer then prices such a row the
#' same however far it moves, and fits every other row as if the
#' survivor were free. That failure is silent - the fit converges, and
#' `logLik()` and `AIC()` report the floored number.
#'
#' `lccdf` removes the class for right censoring by giving core the
#' quantity it actually needs. Measured on a standard normal tail, one
#' process:
#'
#' \tabular{lrr}{
#'   \strong{z} \tab \strong{log(1 - pnorm(z))} \tab
#'     \strong{pnorm(z, lower.tail = FALSE, log.p = TRUE)} \cr
#'   8.0 \tab -35.013 \tab -35.013 \cr
#'   8.3 \tab -Inf \tab -37.494 \cr
#'   37 \tab -Inf \tab -689.031 \cr
#'   500 \tab -Inf \tab -125007.13
#' }
#'
#' Both the value and the derivative are exact on the right, over the
#' whole range.
#'
#' The built-in families that declare `lccdf` are [gaussian()],
#' [lognormal()], [exponential()], [weibull()] and [cox()]. The last
#' three have `log S` in closed form (`-q/mu`, `-(q/scale)^shape`,
#' `-H0(t) * mu`); the first two use `pnorm()`'s own log upper tail.
#'
#' Two families with a CDF do NOT declare one, for measured reasons.
#' [inverse.gaussian()] gains nothing:
#' `RTMBdist::pinvgauss(lower.tail = FALSE, log.p = TRUE)` is computed
#' on the probability scale and reaches `-Inf` at the same
#' `log S = -34` that `log(1 - F)` does. [poisson()] is censored (see
#' Censoring a discrete response) and would use the slot, but cannot
#' declare one: `RTMB::ppois(lower.tail = FALSE, log.p = TRUE)` is
#' exact in R (`-2773.28` at `q = 700`, `lambda = 5`, where
#' `log(1 - F)` is `-Inf`)
#' and does not TAPE - inside `MakeADFun` it reaches
#' `stats::ppois` and errors with "Non-numeric argument to mathematical
#' function". An exact discrete log survivor has to be written out
#' before poisson can have one.
#'
#' `lccdf` fixes RIGHT censoring and nothing else. Left censoring is
#' still `log(F(y) - Flb)`, interval censoring is still a difference of
#' CDFs, and the truncation normalizer is still `log(Fub - Flb)`, so a
#' LEFT-TRUNCATED survival model - delayed entry, which is routine -
#' meets the identical representability problem from the other side.
#' Closing that needs a windowed log-difference slot, and this is the
#' first step rather than the last one.
#'
#' @section Censoring a discrete response:
#' A censoring bound on a discrete response NAMES a value the response
#' can take, and the value is INCLUDED in the event:
#'
#' \tabular{lll}{
#'   \strong{code} \tab \strong{means} \tab \strong{scored as} \cr
#'   `0` "none" \tab `Y == y` \tab `f(y)` \cr
#'   `-1` "left" \tab `Y <= y` \tab `F(y)` \cr
#'   `1` "right" \tab `Y >= y` \tab `1 - F(y - 1)` \cr
#'   `2` "interval" \tab `y <= Y <= y2` \tab `F(y2) - F(y - 1)`
#' }
#'
#' Equivalently: every LOWER edge enters the CDF as `F(edge - 1)`, and
#' upper edges are unchanged because `F` already includes its argument.
#' It is the rule `trunc(lb = )` has always followed on a discrete
#' response, so one number means one thing however a response is
#' bounded, and it is what a count recorded as "5 or more" means.
#'
#' It DIFFERS from brms for RIGHT and INTERVAL censoring only, where
#' brms emits `poisson_lccdf(y | mu)`, that is `P(Y > y)`, and reads an
#' interval as `(y, y2]`. LEFT censoring is `P(Y <= y)` in both packages
#' and agrees exactly. Migrating a right- or interval-censored count
#' model changes its log-likelihood; on 200 poisson draws at
#' `lambda = 4` right censored at 6, the two readings differ by 20.8 log
#' units and 3.7 percent of the estimate. Subtract one from every
#' right-censored and interval lower bound to reproduce a brms fit.
#'
#' The divergence is deliberate, because **brms is internally
#' inconsistent here and frmtmb cannot be both.** brms's own discrete
#' truncation emits `poisson_lccdf(lb - 1 | mu)`, an INCLUSIVE lower
#' bound, `P(Y >= lb)`. So in brms `trunc(lb = 6)` means `Y >= 6` while
#' a right-censored row recorded at 6 means `Y > 6`: one number, two
#' meanings, on one response. frmtmb's `trunc()` reproduces brms bit for
#' bit (poisson on `y = 3,4,5,6` at `b0 = log 4`, `trunc(lb = 2)` gives
#' 6.9990718955 under both, against 6.2954805379 for the exclusive
#' reading), so its `cens()` had to choose between matching brms's
#' censoring and matching its own truncation. It matches its own.
#'
#' It is also the only reading consistent with the censored SIMULATOR,
#' which predates all of this: `simulate(censored = TRUE)` caps a draw
#' with `pmin(pmax(y, lo), hi)`, recording the value `k` exactly when
#' the latent draw is `>= k`. Measured on 4000 draws from a fit censored
#' at 7, the simulated mass at that point is 0.13250, against
#' `P(Y >= 7) = 0.12190` inclusive and `P(Y > 7) = 0.05763` exclusive.
#' The other reading would silently decouple the likelihood from the
#' simulator that `dharma_residuals()` rests on.
#'
#' Two consequences worth knowing. A one-point interval (`y2 == y`) is
#' legal and is the exact observation `P(Y = y)`, where on a continuous
#' response it is refused as an event of probability zero. And the
#' shift assumes the support is the unit integer lattice, so a
#' non-integer censoring bound is refused rather than moved onto a
#' point the family has no mass at.
#'
#' `residuals(type = "osa")` is refused on a censored discrete fit:
#' inclusive bounds make an uncensored row's support `[lo + 1, hi - 1]`
#' rather than the `[lo, hi]` the one-step window is built on.
#'
#' @section Tape-safe scope:
#' `lpdf` and `lcdf` run with RTMB's tape-safe `c()`, `[<-` and
#' `diag<-` in scope automatically (the
#' `"c" <- RTMB::ADoverload("c")` boilerplate is spliced in unless the
#' function already binds it), so base spellings keep the
#' automatic-differentiation class. A helper the density CALLS still
#' needs its own bindings: lexical scope does not travel into other
#' functions.
#'
#' @export
frmtmb_family <- function(family, dpars, links, lpdf, valid_y = NULL,
                          init_dpars = list(), type = "continuous",
                          post = list(), sim = NULL, sim_ctx = NULL,
                          sim_refusal = NULL,
                          primary_dpars = "mu", lcdf = NULL,
                          lccdf = NULL,
                          required_aterms = character(0),
                          accepts_aterms = NULL, se_dpar = NULL,
                          exclusive_aterms = list(),
                          family_finalize = NULL,
                          extra_pars = NULL, drop_intercept = FALSE,
                          structure = NULL) {
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    frm_stop("frmtmb_family(family =) must be one string, the name of ",
             "the family, not ", arg_desc(family), call. = FALSE)
  }
  if (!is.character(dpars) || !length(dpars) || anyNA(dpars)) {
    frm_stop("frmtmb_family(dpars =) must name at least one distributional ",
             "parameter as a character vector, not ", arg_desc(dpars),
             call. = FALSE)
  }
  if (!is.function(lpdf)) {
    frm_stop("frmtmb_family(lpdf =) must be a function (y, dpars, aterms) ",
             "returning the log density of each row, not ", arg_desc(lpdf),
             call. = FALSE)
  }
  check_required_aterms(required_aterms)
  accepts_aterms <- check_accepts_aterms(accepts_aterms)
  se_dpar <- check_se_dpar(se_dpar, dpars, primary_dpars, family)
  exclusive_aterms <- check_exclusive_aterms(exclusive_aterms,
                                             required_aterms)
  if (!is.null(family_finalize) && !is.function(family_finalize)) {
    frm_stop("frmtmb_family(family_finalize =) must be a function ",
             "(fam, y, aterms) returning the family", call. = FALSE)
  }
  # `type` selects the response check, the residual scale and what
  # "response" means to predict(), so an unrecognized string used to
  # produce a family that silently behaved as none of them
  check_string_choice(type, "type",
                      c("continuous", "discrete", "ordinal", "categorical"))
  check_flag(drop_intercept, "drop_intercept")
  if (!all(primary_dpars %in% dpars)) {
    frm_stop("`primary_dpars` must be a subset of `dpars`", call. = FALSE)
  }
  if (!setequal(names(links), dpars)) {
    frm_stop("`links` must name every dpar exactly once", call. = FALSE)
  }
  links <- Map(function(lk, dp) get_link(lk, dpar = dp), links, names(links))
  # user-written densities run with the AD overloads in scope; a
  # function that already binds them itself is left untouched
  lpdf <- frmtmb_ad_overload(lpdf)
  if (!is.null(lcdf)) lcdf <- frmtmb_ad_overload(lcdf)
  if (!is.null(lccdf)) {
    if (!is.function(lccdf)) {
      frm_stop("frmtmb_family(lccdf =) must be a function (q, dpars, ",
               "aterms) returning log S(q), or NULL", call. = FALSE)
    }
    lccdf <- frmtmb_ad_overload(lccdf)
  }
  if (!is.null(structure) && !inherits(structure, "frmtmb_structure")) {
    frm_stop("frmtmb_family(structure =) must come from frmtmb_structure(), ",
             "which is what declares a likelihood that does not factorize ",
             "over rows", call. = FALSE)
  }
  # base::structure(), spelled out: the `structure` ARGUMENT above masks
  # the base function inside this body
  base::structure(
    list(family = family, dpars = dpars, links = links, lpdf = lpdf,
         valid_y = valid_y, init_dpars = init_dpars, type = type,
         post = post, sim = sim, sim_ctx = sim_ctx,
         sim_refusal = sim_refusal, primary_dpars = primary_dpars,
         lcdf = lcdf, lccdf = lccdf, required_aterms = required_aterms,
         accepts_aterms = accepts_aterms, se_dpar = se_dpar,
         exclusive_aterms = exclusive_aterms,
         family_finalize = family_finalize, extra_pars = extra_pars,
         drop_intercept = isTRUE(drop_intercept),
         structure = structure),
    class = "frmtmb_family",
    # the package whose code built the family: a refusal core raises
    # about this family's declarations is that package's
    # (?frmtmb-conditions). An attribute, so the list's names and the
    # partial matches of `$` on them do not move.
    frmtmb_package = frm_env_package(parent.frame())
  )
}

#' @rdname frmtmb_family
#' @export
custom_family <- frmtmb_family

#' The names brms reads off a family object for its links. They are not
#' stored on a `frmtmb_family`: `$` computes them from `links` at the
#' moment they are read.
#'
#' Storing them was tried first and went stale through every write path
#' that does not dispatch to a replacement method with a character name
#' (`[<-`, a numeric `[[<-`, an unclass-and-rebuild), and a fit saved
#' before they existed errored on `$link`
#' (dev/reviews/20260916-famlink.md, findings 1 and 2). A value computed
#' at read time has no copy to go stale and needs nothing saved.
#'
#' @noRd
family_link_view_names <- function(x) {
  c("link", "linkfun", "linkinv",
    paste0("link_", family_link_secondary(x)))
}

#' The parameters that get a `link_<dpar>`: every one but the mean, and
#' for a family whose categories share one joint link (categorical,
#' multinomial) not the per-category predictors either, whose `identity`
#' in `links` is where the softmax starts rather than their scale.
#'
#' @noRd
family_link_secondary <- function(x) {
  dpars <- .subset2(x, "dpars")
  out <- setdiff(dpars, "mu")
  if (family_joint_link(x)) out <- setdiff(out, .subset2(x, "primary_dpars"))
  out
}

#' brms's `is_polytomous()`, for the families this package has: the
#' response is a category (ordinal, categorical) or a vector of counts
#' over categories (multinomial). brms refuses a predictive error and
#' `pp_check(type = "error_binned")` for every one of them.
#'
#' @noRd
fam_is_polytomous <- function(fam) {
  isTRUE(fam[["type"]] %in% c("ordinal", "categorical")) ||
    identical(fam[["family"]], "multinomial")
}

#' @noRd
family_joint_link <- function(x) {
  identical(.subset2(x, "type"), "categorical") ||
    identical(.subset2(x, "family"), "multinomial")
}

#' The value of one of brms's link fields, read from `links`.
#'
#' `link` is the mean's link name; an ordinal family reports the
#' distribution function its thresholds are read through, a joint-link
#' family "logit", and a family with no `mu` (a mixture) "identity", which
#' are brms's values. A link still spelled as a string, as an extension
#' writes one inside `family_finalize()`, is resolved for its functions;
#' an unknown name gives `NULL` functions and is left to frame assembly,
#' which refuses it with the dpar named.
#'
#' @noRd
family_link_view <- function(x, name) {
  links <- .subset2(x, "links")
  lk_name <- function(lk) {
    nm <- if (is.list(lk)) lk[["name"]] else lk
    if (is.character(nm) && length(nm) == 1L) nm else NA_character_
  }
  if (startsWith(name, "link_")) {
    return(lk_name(links[[substring(name, 6L)]]))
  }
  ol <- .subset2(x, "ord_link")
  main <- if (!is.null(ol)) {
    ol
  } else if (family_joint_link(x)) {
    "logit"
  } else if ("mu" %in% .subset2(x, "dpars")) {
    # NULL when the object has lost its links: reported as NA, not as a
    # link the family does not carry
    links[["mu"]]
  } else {
    "identity"
  }
  if (identical(name, "link")) return(lk_name(main))
  obj <- if (is.list(main)) main else {
    tryCatch(get_link(main), error = function(e) NULL)
  }
  obj[[name]]
}

#' Family object field access
#'
#' A family object is a list, and `$` on a list PARTIAL-matches: before
#' frmtmb answered brms's field names, `student()$link` returned the
#' `links` list rather than the link name brms code expects.
#'
#' `$` on a `frmtmb_family` answers brms's link fields, `link`, `linkfun`,
#' `linkinv` and `link_<dpar>` for every parameter but the mean, by
#' reading `links` at the moment of the call, so they always describe the
#' links the fit uses, however the object was edited and whenever it was
#' saved. Every other name matches exactly. A name that is only the
#' prefix of one field is an error naming that field, and a name that
#' matches nothing is `NULL`, as on any list.
#'
#' The link fields are not elements of the list, so `names()`, `[[` and
#' `str()` do not show them. An element stored under one of those names
#' is refused when it is read, because `$` would otherwise have to choose
#' between it and the links the fit actually applies.
#'
#' @param x A `frmtmb_family`.
#' @param name The field name.
#' @return The field.
#' @examples
#' fam <- beta_binomial()
#' fam$link
#' fam$link_phi
#' try(fam$link_p)   # a prefix of link_phi, refused
#' @name frmtmb_family-access
#' @keywords internal
NULL

#' @rdname frmtmb_family-access
#' @export
`$.frmtmb_family` <- function(x, name) {
  nms <- names(x)
  # Every derived field is spelled `link`, `linkfun`, `linkinv` or
  # `link_*`, so any other exact name is a stored element whatever the
  # dpars are. Returning it before the view names are built from the
  # dpars removes most of the cost of a read.
  if (!(name %in% c("link", "linkfun", "linkinv") ||
        startsWith(name, "link_")) && name %in% nms) {
    return(.subset2(x, name))
  }
  view <- family_link_view_names(x)
  if (name %in% view) {
    if (name %in% nms) {
      frm_stop("This family object stores an element named `", name, "`, ",
               "which is one of the link fields `$` derives from `links`. ",
               "Remove the element and set the link in `links` instead",
               call. = FALSE)
    }
    return(family_link_view(x, name))
  }
  if (name %in% nms) return(.subset2(x, name))
  hit <- c(nms, view)[startsWith(c(nms, view), name)]
  if (length(hit) == 1L) {
    frm_stop("A family object has no field `", name, "`. `$` would have ",
             "partial-matched `", hit, "`, and it does not do that on a ",
             "family: write `$", hit, "` if that is the field you mean",
             call. = FALSE)
  }
  NULL
}

#' Make a response validator that rejects any value less than or equal
#' to zero. The families with support on `(0, Inf)` use it so a bad
#' response fails at assembly time with the family name in the message.
#'
#' @noRd
positive_y <- function(name) {
  force(name)
  function(y, aterms) {
    if (any(y <= 0)) {
      frm_stop(name, ": response must be strictly positive", call. = FALSE)
    }
  }
}

#' Make a response validator that rejects a response which is not a
#' non-negative integer. The count families use it.
#'
#' @noRd
count_y <- function(name) {
  force(name)
  function(y, aterms) {
    if (any(y < 0) || any(y != round(y))) {
      frm_stop(name, ": response must be non-negative integers", call. = FALSE)
    }
  }
}

# --- Robust (linear-predictor scale) density inputs ------------------
#
# The dpar contract hands a log-density the INVERSE-LINKED value, which
# is what every numeric post-fit path wants. On the tape it throws away
# the far tail: plogis(40) is exactly 1, so `1 - mu` is exactly 0 and a
# binomial log-density is -Inf where the truth is -40, with a gradient
# of NaN. `build_objective()` therefore stores the linear predictor
# beside each dpar under `.eta_<dpar>`, and the accessors below recover
# the quantity a robust density needs from it. Off the tape the entries
# are absent and every accessor returns the plain form. That is correct
# there, because nothing off the tape is differentiated.

#' Resolve the `link` argument of a public accessor.
#'
#' A family outside frmtmb usually holds only the NAME it wrote in
#' `links = list(coh = "logit")`, since [get_link()] is not part of the
#' public surface, so a name is resolved here rather than making that
#' family carry a link object it has no way to build.
#'
#' An object goes through [get_link()] too, rather than being taken on
#' trust. A list that is not an frmtmb link carries no `logit_eta` and
#' no `log_eta`, so trusting it means falling back to the plain
#' arithmetic and returning `-Inf` where the accessor exists to return
#' -40, with nothing said. `stats::make.link("logit")` is exactly such
#' a list: it spells the derivative `mu.eta` where the contract here is
#' `mu_eta`. The validation costs about 26 us, and a density runs in R
#' once per tape build rather than once per iteration, so it is paid
#' once per fit.
#'
#' @noRd
robust_link <- function(link, dpar) {
  get_link(link, dpar = dpar)
}

#' The value of one dpar, or a refusal that names it.
#'
#' `dpars[[dpar]]` is `NULL` for a name that is not there, and `NULL`
#' propagates instead of stopping: `1 - NULL` is `numeric(0)`, so
#' [dpar_complement()] answered a misspelled dpar with an empty density
#' term and no message. The name is written in the family author's own
#' source, so the refusal says which one was not found and what was.
#'
#' @noRd
robust_dpar_value <- function(dpars, dpar, what) {
  v <- dpars[[dpar]]
  if (is.null(v)) {
    have <- setdiff(names(dpars), grep("^\\.eta_", names(dpars),
                                       value = TRUE))
    frm_stop(what, "(): `dpars` has no distributional parameter '", dpar,
             "'. It carries: ", paste(have, collapse = ", "), call. = FALSE)
  }
  v
}

#' The log-odds behind a probability dpar, on the linear-predictor
#' scale, or NULL when the objective did not supply the linear predictor
#' or the dpar's link has no exact log-odds form (identity, inverse).
#'
#' PARTIAL on purpose, and internal for the same reason: the NULL is how
#' a core family switches between two whole density EXPRESSIONS, which
#' is more than the public accessors offer. A family that only needs the
#' quantity uses [dpar_log_complement()] or [dpar_complement()], which
#' fold the fallback in.
#'
#' @noRd
robust_logit <- function(dpars, link, name = "mu") {
  e <- dpars[[paste0(".eta_", name)]]
  if (is.null(e) || is.null(link$logit_eta)) return(NULL)
  link$logit_eta(e)
}

#' The log mean behind a positive-mean dpar, on the linear-predictor
#' scale, or NULL as in [robust_logit()]. Only the log link supplies it,
#' where it is the linear predictor itself.
#'
#' @noRd
robust_logmu <- function(dpars, link, name = "mu") {
  e <- dpars[[paste0(".eta_", name)]]
  if (is.null(e) || is.null(link$log_eta)) return(NULL)
  link$log_eta(e)
}

#' Exact distributional parameters for a custom density
#'
#' The `dpars` list a [frmtmb_family()] log-density is handed carries
#' each distributional parameter on its own response scale. That is what
#' a density can usually use directly, and it is not enough in the tail.
#' An inverse link SATURATES: `stats::plogis(eta)` is exactly one in
#' double precision from `eta = 36.7368005696771` upward, so a density
#' that forms `1 - mu` by subtraction gets exactly zero there and
#' returns `NaN` for the value AND for the gradient, where the truth is
#' an ordinary large negative number. Below that point the subtraction
#' is not fatal, only wrong: the complement it forms is off by 1.0e-3
#' relative at `eta = 30`.
#'
#' The linear predictor never saturated. [frm()] stores it beside each
#' distributional parameter while the objective is taped, and these
#' four accessors recover from it the quantity the density needs.
#' Write a density over them instead of over `log(mu)`, `log(1 - mu)`
#' and `1 - mu`.
#'
#' \describe{
#'   \item{`dpar_log()`}{`log(x)`, for a dpar that must stay positive
#'     (a `shape`, a `phi`, a `sigma`) and for one on the unit
#'     interval alike.}
#'   \item{`dpar_log1m()`}{`log(1 - x)` for a dpar on the unit
#'     interval: a `zi` or `hu` gate, a probability, a coherence. This
#'     is the term that dies first, and it stays finite and exactly
#'     differentiable at a value the optimizer has pushed against 1.}
#'   \item{`dpar_log_complement()`}{Both of the above at once, as `l`
#'     and `l1m`, from ONE log odds, for a density that needs the pair.
#'     A mixture gate does. Take a one-sided function when you need
#'     one term. The pair records a second `RTMB::logspace_add()` over
#'     the whole response and the tape REPLAYS it on every gradient
#'     evaluation: at 1000 rows the accessor tapes 14002 AD nodes
#'     against 10002, and one gradient sweep of the `cross_wishart()`
#'     density costs 633 us against 480 us. What you are not paying
#'     for is the R call, which happens once per fit, when the tape is
#'     built.}
#'   \item{`dpar_complement()`}{`x` as `p` and `1 - x` as `q`, for a
#'     density that needs the natural scale. You call it for `q`; `p`
#'     comes back with it because one log-odds gives both, and a pair
#'     taken from one source cannot disagree with itself.}
#' }
#'
#' @section How far this reaches:
#' The relative error of the plain arithmetic on a logit gate,
#' against what these accessors return, at the linear predictors an
#' optimizer can visit. One process, the log-scale form as the
#' reference:
#'
#' \tabular{lrr}{
#'   \strong{eta} \tab \strong{`1 - plogis(eta)`} \tab
#'     \strong{`log(1 - plogis(eta))`} \cr
#'   20 \tab 3.6e-08 \tab 1.8e-09 \cr
#'   30 \tab 1.0e-03 \tab 3.4e-05 \cr
#'   36 \tab 4.3e-02 \tab 1.2e-03 \cr
#'   36.7368005696771 \tab 1 (exactly 0) \tab 1.9e-02 \cr
#'   40 \tab 1 (exactly 0) \tab Inf, the log is -Inf \cr
#'   700 \tab 1 (exactly 0) \tab Inf, the log is -Inf
#' }
#'
#' The plain form does not fail suddenly. It loses digits from
#' `eta = 20`, is down to three at 30 and to one at 36, and only
#' then becomes `-Inf` with a `NaN` gradient. A fit that walks out
#' there converges to a floored likelihood without a word.
#'
#' @section On the tape and off it:
#' The `.eta_` entries exist only while the objective is taped. Your
#' density reads them exactly once, when the tape is built; what runs
#' on every gradient sweep after that is the arithmetic that read
#' RECORDED, not your R code. Anywhere else the dpar values arrive
#' alone, every accessor falls back to the plain arithmetic, and that
#' is correct because nothing outside the tape is differentiated.
#'
#' Anywhere else is a smaller place than it sounds. `fitted()`,
#' `predict()`, `logLik()`, `simulate()` and `vcov()` do not evaluate
#' the density at all: they read `post$mean_fn`, the link, the
#' optimizer's own value and `sim`. What can reach these accessors off
#' the tape is a family's OWN numeric helpers, `post$dev_fn` and
#' `post$var_fn` among them, and [check_custom_family()]. Counted on a
#' fitted `cross_wishart()`, whose deviance helper is the one that
#' calls them, fourteen post-fit entry points reach the arithmetic
#' exactly once between them, through `residuals(type = "deviance")`.
#'
#' **Test your density on both paths.** A family author who exercises
#' only the off-tape path never runs the branch these accessors exist
#' for, and a density that is wrong on the tape still fits: it converges
#' to the wrong place, or dies at `NA/NaN gradient evaluation` with
#' nothing naming the cause. Build the on-tape list by hand to test it,
#' as the example below does, and put the value under
#' `.eta_<dpar>` yourself for that one purpose.
#'
#' @section The `.eta_` entries are not the API:
#' `.eta_<dpar>` is a reserved name, and reading it directly is not
#' supported even though you can see it in `dpars`. It holds the linear
#' predictor on the LINK scale, so what it means depends on the dpar's
#' link, and a density that assumes one link reads a different quantity
#' the moment the family gains a `link_<dpar>` argument. `-log1p(exp(e))`
#' is `log(1 - x)` on a logit and is nothing at all on an identity link;
#' `e` is `log(x)` on a log link and is not on a softplus. The accessors
#' take the link and branch on what it can supply, which is why they
#' take it as an argument and why it is not optional. Whether an entry
#' is present is also a property of the phase, not of the model, and
#' folding that branch in is most of what these do.
#'
#' @section If you need the raw log odds:
#' There is no accessor for it, and there does not need to be:
#' `dpar_log(d, p, lk) - dpar_log1m(d, p, lk)` reconstructs it. Measured
#' on a logit at `eta` in `{-700, -40, 0, 40, 700}` the absolute error
#' is 0, and it peaks at 5.0e-17 near `eta = 0` where the two terms
#' nearly cancel. A density term linear in the natural parameter cares
#' about that absolute error, not the relative one.
#'
#' @param dpars The named list of distributional parameters the density
#'   was handed. Read it by name, never by position.
#' @param dpar The name of the distributional parameter to read, as
#'   `dpars` names it: `"mu"`, `"zi"`, `"shape"`.
#' @param link That dpar's OWN link, as the link name your family gave
#'   in `links = ` (`"logit"`) or as a link object. It is not optional:
#'   the linear predictor is on the link scale, so nothing can be
#'   recovered from it without knowing which link put it there.
#' @return `dpar_log()` and `dpar_log1m()` a numeric or advector
#'   vector. The other two a list of two such vectors: `l` and `l1m`
#'   for `dpar_log_complement()`, which are what the two one-sided
#'   functions return, and `p` and `q` for `dpar_complement()`.
#' @seealso [frmtmb_family()] for the density these serve,
#'   [frmtmb-links] for which links carry an exact form and which fall
#'   back, and [frmtmb-extension-api] for the accessors a family uses
#'   AFTER a fit, at the estimates
#' @examples
#' # a gate at plogis(40), which is exactly 1 in double precision
#' on_tape <- list(zi = stats::plogis(40), .eta_zi = 40)
#' off_tape <- list(zi = stats::plogis(40))
#'
#' # the value a density needs, and what the subtraction gives instead
#' dpar_log1m(on_tape, "zi", "logit")
#' log(1 - off_tape$zi)
#'
#' # off the tape the entry is absent and the plain form comes back
#' dpar_log1m(off_tape, "zi", "logit")
#'
#' # both terms from one log odds, for a density that needs the pair
#' unlist(dpar_log_complement(on_tape, "zi", "logit"))
#'
#' # the log of a positive dpar, through that dpar's own link
#' dpar_log(list(shape = exp(3), .eta_shape = 3), "shape", "log")
#' dpar_log(list(shape = log1p(exp(3)), .eta_shape = 3), "shape",
#'          "softplus")
#'
#' # a density written over the pair rather than over 1 - mu
#' d <- list(mu = stats::plogis(40), .eta_mu = 40)
#' mp <- dpar_complement(d, "mu", "logit")
#' c(p = mp$p, q = mp$q)
#' @name frmtmb-robust-dpars
NULL

#' @rdname frmtmb-robust-dpars
#' @export
dpar_log <- function(dpars, dpar, link) {
  lk <- robust_link(link, dpar)
  lm <- robust_logmu(dpars, lk, dpar)
  if (!is.null(lm)) return(lm)
  # A unit-interval dpar has no log mean; it has a log odds, and log(x)
  # comes out of that exactly. Without this branch the same call that
  # is exact at eta = -800 on a log link returns -Inf on a logit, where
  # plogis(-800) is exactly 0. One function, two answers, would have
  # been a trap of its own.
  lo <- robust_logit(dpars, lk, dpar)
  if (!is.null(lo)) return(log_inv_logit(lo))
  log(robust_dpar_value(dpars, dpar, "dpar_log"))
}

#' @rdname frmtmb-robust-dpars
#' @export
dpar_log1m <- function(dpars, dpar, link) {
  lo <- robust_logit(dpars, robust_link(link, dpar), dpar)
  # log1p(-p), not log(1 - p): see dpar_log_complement() below
  if (is.null(lo)) {
    return(log1p(-robust_dpar_value(dpars, dpar, "dpar_log1m")))
  }
  log1m_inv_logit(lo)
}

#' @rdname frmtmb-robust-dpars
#' @export
dpar_log_complement <- function(dpars, dpar, link) {
  lo <- robust_logit(dpars, robust_link(link, dpar), dpar)
  if (is.null(lo)) {
    # log1p(-p), not log(1 - p). The two are bit-identical from
    # p = 0.25 up, where `1 - p` is exact, so this buys nothing at the
    # saturating end. It buys the OTHER boundary, which a gate reaches
    # just as often: at p = 1e-17 the subtraction rounds `1 - p` to 1
    # and log(1 - p) is exactly 0, where the value is -1e-17. Measured
    # over 800 points spread through (0, 1), the two spellings differ
    # at 397 of them and every one has p < 0.1.
    p <- robust_dpar_value(dpars, dpar, "dpar_log_complement")
    return(list(l = log(p), l1m = log1p(-p)))
  }
  list(l = log_inv_logit(lo), l1m = log1m_inv_logit(lo))
}

#' @rdname frmtmb-robust-dpars
#' @export
dpar_complement <- function(dpars, dpar, link) {
  lo <- robust_logit(dpars, robust_link(link, dpar), dpar)
  if (is.null(lo)) {
    mu <- robust_dpar_value(dpars, dpar, "dpar_complement")
    return(list(p = mu, q = 1 - mu))
  }
  list(p = exp(log_inv_logit(lo)), q = exp(log1m_inv_logit(lo)))
}

#' `P(a < Z < b)` for standard normal bounds. `pnorm(b) - pnorm(a)`
#' loses every significant digit when both bounds sit in the upper tail,
#' which is exactly where truncated means are computed; the mirrored
#' form keeps full precision there.
#'
#' @noRd
pnorm_diff <- function(a, b) {
  ifelse(a > 0,
         stats::pnorm(-a) - stats::pnorm(-b),
         stats::pnorm(b) - stats::pnorm(a))
}

#' Residual SD including a known `se()` component (meta-analysis):
#' `se()` alone replaces `sigma`; `se(x, sigma = TRUE)` adds them in
#' quadrature.
#'
#' @noRd
resid_sd <- function(sigma, aterms) {
  if (is.null(aterms[["se"]])) return(sigma)
  if (isTRUE(aterms[["se_sigma"]])) sqrt(sigma^2 + aterms[["se"]]^
    2) else aterms[["se"]]
}

#' `trunc()` bounds from an aterm-value list as a pair of length-n
#' numeric vectors, or `NULL` when the response is not truncated. An
#' absent bound is the family's unbounded default; the per-family
#' truncated means all handle the infinities.
#'
#' @noRd
trunc_bounds <- function(aterms, n) {
  lb <- aterms[["trunc_lb"]]
  ub <- aterms[["trunc_ub"]]
  if (is.null(lb) && is.null(ub)) return(NULL)
  list(lb = rep(lb %||% -Inf, length.out = n),
       ub = rep(ub %||% Inf, length.out = n))
}

# On a truncated response this is E[Y | lb <= Y <= ub], the quantity
# fitted(), residuals() and frm_linpred(type = "response") report; per-dpar
# predictions stay untruncated, because they describe the latent
# parameter rather than the observable.

#' @rdname frmtmb-extension-api
#' @export
response_mean <- function(fam, dpars, aterms) {
  mu <- if (!is.null(fam[["post"]]$mean_fn)) {
    fam[["post"]]$mean_fn(dpars, aterms)
  } else if ("mu" %in% names(dpars)) {
    # the convention custom_family() writes down: no mean function, the
    # mean is mu
    dpars[["mu"]]
  } else {
    # a family with neither has no mean, and reporting its first
    # parameter as one reported a race model's drift rate
    frm_stop("family '", fam[["family"]], "' declares no mean: it has no dpar ",
             "named mu and no post$mean_fn, so fitted() and ",
             "frm_linpred(type = \"response\") have nothing to return. ",
             "Ask for type = \"link\" or a dpar by name.", call. = FALSE,
             package = frm_family_package(fam))
  }
  tb <- trunc_bounds(aterms, length(mu))
  if (is.null(tb)) return(mu)
  tmf <- fam[["post"]]$trunc_mean_fn
  if (is.null(tmf)) {
    frm_stop("Family '", fam[["family"]], "' has no truncated mean; ",
             "fitted(), residuals() and frm_linpred(type = \"response\") ",
             "would report the untruncated mean", call. = FALSE,
             package = frm_family_package(fam))
  }
  tmf(dpars, aterms, tb$lb, tb$ub)
}

# --- Unit deviances ---
#
# d_i = 2 * (loglik of the saturated fit - loglik at the fitted value),
# with the dispersion parameter held at its estimate, which is the
# quantity glm() calls the unit deviance. Families that are exponential
# dispersion models reproduce stats::glm()'s dev.resids exactly; the
# rest (negbinomial, beta, tweedie) use the same saturated-likelihood
# definition at a fixed shape.

#' `y log(y / mu)` under the `0 log 0 = 0` convention the saturated fit
#' needs at a zero count.
#'
#' @noRd
ylogy_mu <- function(y, mu) ifelse(y > 0, y * log(y / mu), 0)

#' Negative-binomial unit deviance at a given size (shape). nbinom1
#' feeds it the row's own size `mu / phi`.
#'
#' @noRd
nbinom_deviance <- function(y, mu, size) {
  2 * (ylogy_mu(y, mu) - (y + size) * log((y + size) / (mu + size)))
}

#' Binomial unit deviance on COUNTS out of `size` trials; both terms
#' vanish at the boundaries `y = 0` and `y = size`.
#'
#' @noRd
binomial_deviance <- function(y, mu, size) {
  2 * (ylogy_mu(y, size * mu) + ylogy_mu(size - y, size * (1 - mu)))
}

#' Gamma unit deviance; also the exponential one (shape fixed at 1).
#'
#' @noRd
gamma_deviance <- function(y, mu) 2 * ((y - mu) / mu - log(y / mu))

#' Families that define a unit deviance, for the `residuals()` refusal
#' message. Read off the registry so the list cannot drift; the
#' constructors that need arguments (multinomial) simply drop out.
#'
#' @noRd
deviance_family_names <- function() {
  nms <- names(family_registry)
  ok <- vapply(nms, function(nm) {
    fam <- tryCatch(family_registry[[nm]](), error = function(e) NULL)
    !is.null(fam) && !is.null(fam[["post"]]$dev_fn)
  }, TRUE)
  sort(unique(nms[ok]))
}

#' Deviance residuals: `sign(y - E[Y]) * sqrt(w_i d_i)`. Weights
#' multiply the unit deviance, as in `glm()`. A truncated or censored
#' response is refused: the fitted likelihood is not the family's own
#' density there, so the saturated comparison the unit deviance is built
#' on does not describe the model that was estimated.
#'
#' @noRd
deviance_residuals <- function(fam, y, dpars, aterms, n) {
  dev <- fam[["post"]]$dev_fn
  if (is.null(dev)) {
    frm_stop("residuals(type = \"deviance\") is not available for family '",
             fam[["family"]], "': it has no standard unit deviance. Families ",
             "with one: ", paste(deviance_family_names(), collapse = ", "),
             ". Use type = \"osa\" or dharma_residuals() instead.",
             call. = FALSE,
             package = frm_family_package(fam))
  }
  if (!is.null(trunc_bounds(aterms, n))) {
    frm_stop("residuals(type = \"deviance\") is not defined for a trunc()ed ",
             "response: the unit deviance compares against the untruncated ",
             "family, not the likelihood the model was fitted with. Use ",
             "type = \"osa\", which builds its CDF on [lb, ub]", call. = FALSE)
  }
  if (!is.null(aterms[["cens"]]) && any(aterms[["cens"]] != 0)) {
    frm_stop("residuals(type = \"deviance\") is not defined on a cens()ed ",
             "response: a censored row observes an event, not a value, so it ",
             "has no unit deviance. Use type = \"osa\"", call. = FALSE)
  }
  d <- dev(y, dpars, aterms)
  w <- aterms[["weights"]] %||% 1
  # rounding can push an exactly saturated row a few ulps below zero
  sign(y - response_mean(fam, dpars, aterms)) * sqrt(pmax(w * d, 0))
}

#' Per-observation subset of dpar / aterm vectors, for resampling the
#' rows a rejection step has not accepted yet.
#'
#' @noRd
subset_obs <- function(x, idx, n) {
  lapply(x, function(v) {
    if (is.numeric(v) && length(v) %in% c(1L, n)) {
      rep(v, length.out = n)[idx]
    } else {
      v
    }
  })
}

#' One simulated response vector, respecting `trunc()` bounds by
#' rejection: out-of-bounds draws are redrawn until every row is inside
#' its own interval. Bounds that exclude nearly all the family's mass
#' never converge, so the iteration cap reports the acceptance rate
#' instead of spinning.
#'
#' @noRd
sim_response <- function(fam, dpars, aterms, n, max_iter = 100L,
                         extra = NULL) {
  # families whose draws need parameters outside the dpar system
  # (ordinal thresholds) declare a fourth argument; the rest keep the
  # three-argument contract
  fam_sim <- if (length(formals(fam[["sim"]])) >= 4L) {
    function(dp, av, nn) fam[["sim"]](dp, av, nn, extra)
  } else {
    fam[["sim"]]
  }
  y <- fam_sim(dpars, aterms, n)
  tb <- trunc_bounds(aterms, n)
  if (is.null(tb)) return(y)
  drawn <- n
  bad <- which(y < tb$lb | y > tb$ub)
  it <- 0L
  while (length(bad)) {
    it <- it + 1L
    if (it > max_iter) {
      frm_stop("trunc(): rejection sampling did not fill ", length(bad),
               " of ", n, " rows in ", max_iter, " passes (acceptance rate ",
               format((n - length(bad)) / drawn, digits = 2),
               "). The bounds exclude nearly all of the fitted ",
               "distribution's mass.", call. = FALSE)
    }
    drawn <- drawn + length(bad)
    yb <- fam_sim(subset_obs(dpars, bad, n), subset_obs(aterms, bad, n),
                  length(bad))
    y[bad] <- yb
    bad <- bad[yb < tb$lb[bad] | yb > tb$ub[bad]]
  }
  y
}

# --- the structured simulator contract --------------------------------
#
# `fam$sim(dpars, aterms, n)` is stateless and ROWWISE: it sees one
# column of numbers per distributional parameter and draws each row on
# its own. Three model classes cannot be written that way.
#
#   - family-level EXTRAS. An lca() draw needs the item-profile
#     parameters, which live outside the dpar system. That case already
#     has an answer: a four-argument `sim` receives the extra vector.
#   - GROUP or TIME structure. A group-level mixture draws one class per
#     GROUP, and an hmm() draw walks a Markov chain per SEQUENCE. Both
#     structures are resolved at frame-assembly time and live on the
#     frame, as the structured family's `blocks` entry, not on the
#     family object.
#   - CROSS-ROW dependence. An autocor residual is one multivariate draw
#     per group, so no per-row simulator exists at all.
#
# `fam$sim_ctx(ctx)` is the extension: one function per family, taking a
# CONTEXT that carries the fit-like object as well as the evaluated
# dpars. The context is built the same way from all three entry points -
# `simulate()` on a fit, `posterior_predict()` on one draw (where
# `draws_fit_at()` makes a real fit out of it), and `frm_simulate()` on
# the de novo shim - so a structured family has exactly ONE
# implementation and every entry point reaches it.
#
# `[[ ]]`, never `$`, on the context and on frame fields: `$` partial
# matching once let a `ctx$mix` read return the group structure instead.

#' The simulation context: everything a simulator may read, assembled
#' identically from a fit, from one posterior draw, and from the de novo
#' shim. `fit` is any object carrying `spec`, `frame` and `estimates`.
#'
#' @noRd
sim_context <- function(fit, rspec, dpars, aterms = NULL, n = NULL,
                        extra = NULL, max_iter = NULL) {
  resp <- rspec[["resp_name"]]
  frame <- fit[["frame"]]
  list(fit = fit,
       family = rspec[["family"]],
       rspec = rspec,
       resp = resp,
       dpars = dpars,
       aterms = aterms %||% frame[["aterm_values"]][[resp]] %||% list(),
       n = n %||% frame[["n_obs"]],
       # a multivariate fit's extras are namespaced by response; the
       # simulator reads its own block under its family's names
       extra = resp_extras(frame, extra %||% fit_extras(fit), resp),
       autocor = frame[["autocor"]][[resp]],
       # the rejection-sampling limit a trunc()ed response needs, which
       # predict() exposes as brms's `ntrys`; absent means the
       # simulator's own default
       max_iter = max_iter,
       block = frame_block_of(frame, resp))
}

#' The family's structured simulator, from its structure or, for a
#' family that has one without a structure, from its own slot.
#'
#' @noRd
fam_sim_ctx <- function(fam) {
  fam_structure(fam)[["sim_ctx"]] %||% fam[["sim_ctx"]]
}

#' Whether a family can simulate at all, by either half of the contract.
#'
#' @noRd
sim_can <- function(fam) {
  !is.null(fam[["sim"]]) || !is.null(fam_sim_ctx(fam))
}

#' The family's own explanation for having no simulator, appended by
#' whichever entry point refused. Empty for a family that simply has not
#' been given one yet.
#'
#' @noRd
sim_note <- function(fam) {
  note <- fam[["sim_refusal"]]
  if (is.null(note)) "" else paste0(". ", note)
}

#' Whether the draws of this context come whole rather than row by row,
#' which is what makes `trunc()` rejection and `newdata` inapplicable.
#'
#' @noRd
sim_is_structured <- function(ctx) {
  !is.null(ctx[["autocor"]]) || !is.null(fam_sim_ctx(ctx[["family"]]))
}

#' One simulated response for a context: the single implementation every
#' entry point calls. Cross-row structure is resolved first (a residual
#' correlation belongs to the frame, not to the family), then the
#' family's structured simulator, then the rowwise contract.
#'
#' @noRd
sim_draw <- function(ctx) {
  ac <- ctx[["autocor"]]
  sf <- fam_sim_ctx(ctx[["family"]])
  if (is.null(ac) && is.null(sf)) {
    return(sim_response(ctx[["family"]], ctx[["dpars"]], ctx[["aterms"]],
                        ctx[["n"]],
                        max_iter = ctx[["max_iter"]] %||% 100L,
                        extra = ctx[["extra"]]))
  }
  # drawn one within-group position at a time, each row by the family's
  # own rowwise simulator, so trunc() rejection applies row by row
  if (autocor_is_cond(ac)) return(sim_autocor_cond(ctx, ac))
  if (!is.null(trunc_bounds(ctx[["aterms"]], ctx[["n"]]))) {
    frm_stop("trunc() cannot be combined with a structured draw (here: '",
             ctx[["family"]][["family"]], "'): a hidden state sequence, a ",
             "group-level latent class and a correlated residual are each ",
             "drawn whole, so a row outside its bounds cannot be redrawn ",
             "on its own and the rejection step has nothing to resample",
             call. = FALSE,
             package = frm_family_package(ctx[["family"]]))
  }
  if (!is.null(ac)) return(sim_autocor_rows(ctx, ac))
  sf(ctx)
}

#' The autocor branch: one multivariate residual draw per group added to
#' the mean predictor, rather than n independent family draws.
#'
#' @noRd
sim_autocor_rows <- function(ctx, ac) {
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  th <- ctx[["fit"]][["estimates"]][["thetaac"]]
  R <- autocor_cor(th[ac[["theta_idx"]]], ac)
  rep(dp[["mu"]], length.out = n) +
    autocor_draw_resid(ac, R, rep(dp[["sigma"]], length.out = n), n,
                       nu = if (isTRUE(ac[["student"]])) dp[["nu"]][1])
}

#' The `cov = FALSE` branch: a fresh replicate of brms's recursion.
#'
#' A replicate has no observed past, so the lagged errors are the ones
#' it draws itself: at each within-group position the rows present get
#' `mu` plus the ARMA term of their own earlier draws, are drawn by the
#' family's rowwise simulator, and leave `err = y - mu - MA` behind for
#' the next position. The first row of a group has no lagged term, which
#' is exactly the conditioning the likelihood makes.
#'
#' @noRd
sim_autocor_cond <- function(ctx, ac) {
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  av <- ctx[["aterms"]]
  th <- ctx[["fit"]][["estimates"]][["thetaac"]][ac[["theta_idx"]]]
  cf <- autocor_cond_coefs(th, ac)
  mu <- rep_len(as.numeric(dp[["mu"]]), n)
  y <- numeric(n)
  err <- vector("list", length(ac[["pos_rows"]]))
  for (t in seq_along(ac[["pos_rows"]])) {
    rows <- ac[["pos_rows"]][[t]]
    k <- seq_along(rows)
    sma <- 0
    for (i in seq_len(min(length(cf$ma), t - 1L))) {
      sma <- sma + cf$ma[i] * err[[t - i]][k]
    }
    sar <- 0
    for (i in seq_len(min(length(cf$ar), t - 1L))) {
      sar <- sar + cf$ar[i] * err[[t - i]][k]
    }
    dpt <- subset_obs(dp, rows, n)
    dpt[["mu"]] <- mu[rows] + sma + sar
    yt <- sim_response(ctx[["family"]], dpt, subset_obs(av, rows, n),
                       length(rows), max_iter = ctx[["max_iter"]] %||% 100L,
                       extra = ctx[["extra"]])
    y[rows] <- yt
    err[[t]] <- yt - mu[rows] - sma
  }
  y
}

#' Gaussian family, dpars `mu` and `sigma`. A known `se()` term enters
#' the residual SD, which is what makes meta-analysis work. It supplies
#' a CDF, so `cens()` and `trunc()` apply to it.
#'
#' @noRd
fam_gaussian <- function(link = "identity", link_sigma = "log") {
  lk_sigma <- dpar_link(link_sigma, "sigma", "gaussian", dpar_links_positive)
  frmtmb_family(
    "gaussian",
    accepts_aterms = c("weights", "cens", "trunc", "se", "mi"),
    se_dpar = "sigma",
    dpars = c("mu", "sigma"),
    links = list(mu = mu_link(link, "gaussian"), sigma = lk_sigma),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], resid_sd(dpars[["sigma"]], aterms),
                log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::pnorm((q - dpars[["mu"]]) / resid_sd(dpars[["sigma"]], aterms))
    },
    lccdf = function(q, dpars, aterms) {
      RTMB::pnorm((q - dpars[["mu"]]) / resid_sd(dpars[["sigma"]], aterms),
                  lower.tail = FALSE, log.p = TRUE)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) resid_sd(dpars[["sigma"]], aterms)^2,
      # sigma is the dispersion, so it divides out of the unit deviance.
      # se() breaks that: the residual sd is row-specific, and a raw
      # squared residual would then compare rows of different precision
      # on one scale. The known variance enters exactly as a glm prior
      # weight sigma^2 / s_i^2, which is 1 without se(); se() alone maps
      # sigma out at 1, so the row weight is the usual 1 / se_i^2.
      dev_fn = function(y, dpars, aterms) {
        s <- resid_sd(dpars[["sigma"]], aterms)
        (y - dpars[["mu"]])^2 * (dpars[["sigma"]] / s)^2
      },
      # mu + sigma * (phi(a) - phi(b)) / (Phi(b) - Phi(a))
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        s <- resid_sd(dpars[["sigma"]], aterms)
        a <- (lb - dpars[["mu"]]) / s
        b <- (ub - dpars[["mu"]]) / s
        dpars[["mu"]] + s * (stats::dnorm(a) - stats::dnorm(b)) /
          pnorm_diff(a, b)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnorm(n, dpars[["mu"]], resid_sd(dpars[["sigma"]], aterms))
    }
  )
}

#' brms's `rate(denom)` on a count family: the exposure multiplies the
#' mean, `mu * denom`. Without the term this returns `dpars$mu` itself,
#' so a model without it is unchanged.
#'
#' Under a log link the product is taken on the linear-predictor scale,
#' `exp(eta + log(denom))`, which is brms's `eta + log_denom` and is
#' the same arithmetic as `offset(log(denom))`; the numeric post-fit
#' paths, where no linear predictor rides along, multiply.
#'
#' @noRd
rate_mu <- function(dpars, link, aterms) {
  d <- aterms[["rate"]]
  if (is.null(d)) return(dpars[["mu"]])
  lm <- robust_logmu(dpars, link)
  if (!is.null(lm)) return(exp(lm + log(d)))
  dpars[["mu"]] * d
}

#' The mean of a count family under `rate()`, for the post-fit paths.
#' A mean function spelled this way still counts as "the mean is mu"
#' (`mean_is_mu()`), because it is exactly mu when the term is absent;
#' `has_rate()` routes a response that carries it through the mean.
#'
#' @noRd
rate_mean <- function(dpars, aterms) {
  d <- aterms[["rate"]]
  if (is.null(d)) dpars[["mu"]] else dpars[["mu"]] * d
}

#' A count family's shape under `rate()`: brms multiplies it by the
#' exposure too (`shape .* denom`), so the variance becomes
#' `mu d (1 + mu / shape)`.
#'
#' @noRd
rate_shape <- function(shape, aterms) {
  d <- aterms[["rate"]]
  if (is.null(d)) shape else shape * d
}

#' Poisson family, single dpar `mu` (the mean). It supplies a CDF and a
#' truncated mean, so `cens()` and `trunc()` apply to it. `rate()`
#' multiplies the mean by an exposure, as in brms.
#'
#' @noRd
fam_poisson <- function(link = "log") {
  lk <- mu_link(link, "poisson")
  frmtmb_family(
    "poisson",
    accepts_aterms = c("weights", "cens", "trunc", "rate"),
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, rate_mu(dpars, lk, aterms), log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::ppois(q, rate_mu(dpars, lk, aterms))
    },
    valid_y = count_y("poisson"),
    init_dpars = list(mu = function(y, aterms) {
      mean(y / (aterms[["rate"]] %||% 1)) + 0.1
    }),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) rate_mean(dpars, aterms),
      var_fn = function(dpars, aterms) rate_mean(dpars, aterms),
      dev_fn = function(y, dpars, aterms) {
        mu <- rate_mean(dpars, aterms)
        2 * (ylogy_mu(y, mu) - (y - mu))
      },
      # sum_{lb}^{ub} y dpois(y) = mu * (F(ub-1) - F(lb-2)), over the
      # same F(ub) - F(lb-1) normalizer the likelihood uses: the
      # inclusive lower bound keeps its own mass (brms#1903)
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        mu <- rate_mean(dpars, aterms)
        mu * (stats::ppois(ub - 1, mu) - stats::ppois(lb - 2, mu)) /
          (stats::ppois(ub, mu) - stats::ppois(lb - 1, mu))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rpois(n, rate_mean(dpars, aterms))
    }
  )
}

#' Binomial family, single dpar `mu` (the success probability). The
#' number of trials comes from the `trials()` addition term and defaults
#' to one, so the response is a count out of `trials`.
#'
#' On the tape it uses `RTMB::dbinom_robust()`, TMB's log-odds
#' parameterization of the same density, which is what glmmTMB fits
#' with. `dbinom(y, n, plogis(eta))` is -Inf from `|eta| = 37` and
#' already wrong in the second decimal at 30, so a separated predictor
#' takes the optimizer into a region where it has no gradient at all.
#'
#' @noRd
fam_binomial <- function(link = "logit") {
  lk <- mu_link(link, "binomial")
  frmtmb_family(
    "binomial",
    accepts_aterms = c("weights", "trials"),
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]] %||% 1
      lo <- robust_logit(dpars, lk)
      if (is.null(lo)) {
        return(RTMB::dbinom(y, size, dpars[["mu"]], log = TRUE))
      }
      RTMB::dbinom_robust(y, size, lo, log = TRUE)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]] %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        frm_stop("binomial: response must be integer counts in [0, trials]",
                 call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        size <- aterms[["trials"]] %||% 1
        p <- mean(y / size)
        min(max(p, 0.02), 0.98)
      }
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars,
                         aterms) dpars[["mu"]] * (aterms[["trials"]] %||% 1),
      var_fn = function(dpars, aterms) {
        (aterms[["trials"]] %||% 1) * dpars[["mu"]] * (1 - dpars[["mu"]])
      },
      dev_fn = function(y, dpars, aterms) {
        binomial_deviance(y, dpars[["mu"]], aterms[["trials"]] %||% 1)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, aterms[["trials"]] %||% 1, dpars[["mu"]])
    }
  )
}

#' Gamma family in the mean parameterization, dpars `mu` and `shape`.
#' The scale is `mu / shape`, so `mu` stays the mean at any shape.
#'
#' @noRd
fam_Gamma <- function(link = "log", link_shape = "log") {
  lk_shape <- dpar_link(link_shape, "shape", "Gamma", dpar_links_positive)
  frmtmb_family(
    "Gamma",
    accepts_aterms = "weights",
    dpars = c("mu", "shape"),
    links = list(mu = mu_link(link, "Gamma"), shape = lk_shape),
    lpdf = function(y, dpars, aterms) {
      RTMB::dgamma(y, shape = dpars[["shape"]],
                                       scale = dpars[["mu"]] / dpars[["shape"]],
                   log = TRUE)
    },
    valid_y = positive_y("Gamma"),
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      shape = function(y, aterms) {
        max(mean(y)^2 / stats::var(y), 0.1)
      }
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) dpars[["mu"]]^2 / dpars[["shape"]],
      # 1 / shape is the dispersion and divides out
      dev_fn = function(y, dpars, aterms) gamma_deviance(y, dpars[["mu"]])
    ),
    sim = function(dpars, aterms, n) {
      stats::rgamma(n, shape = dpars[["shape"]],
                                      scale = dpars[["mu"]] / dpars[["shape"]])
    }
  )
}

#' Lognormal family, dpars `mu` and `sigma`. As in brms, they are the
#' mean and SD on the LOG scale, which is why `mu` takes the identity
#' link. It supplies a CDF and a truncated mean.
#'
#' @noRd
fam_lognormal <- function(link = "identity", link_sigma = "log") {
  lk_sigma <- dpar_link(link_sigma, "sigma", "lognormal", dpar_links_positive)
  frmtmb_family(
    "lognormal",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = c("mu", "sigma"),
    links = list(mu = mu_link(link, "lognormal"), sigma = lk_sigma),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(log(y), dpars[["mu"]], dpars[["sigma"]], log = TRUE) - log(y)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::pnorm((log(q) - dpars[["mu"]]) / dpars[["sigma"]])
    },
    lccdf = function(q, dpars, aterms) {
      RTMB::pnorm((log(q) - dpars[["mu"]]) / dpars[["sigma"]],
                  lower.tail = FALSE, log.p = TRUE)
    },
    valid_y = positive_y("lognormal"),
    init_dpars = list(
      mu = function(y, aterms) mean(log(y)),
      sigma = function(y, aterms) stats::sd(log(y))
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        exp(dpars[["mu"]] + dpars[["sigma"]]^2 / 2)
      },
      var_fn = function(dpars, aterms) {
        (exp(dpars[["sigma"]]^2) - 1) * exp(2 * dpars[["mu"]] +
          dpars[["sigma"]]^2)
      },
      # E[Y] * (Phi(b - sigma) - Phi(a - sigma)) / (Phi(b) - Phi(a)),
      # a, b the log-scale standardized bounds
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        sg <- dpars[["sigma"]]
        a <- (log(pmax(lb, 0)) - dpars[["mu"]]) / sg
        b <- (log(ub) - dpars[["mu"]]) / sg
        exp(dpars[["mu"]] + sg^2 / 2) * pnorm_diff(a - sg, b - sg) /
          pnorm_diff(a, b)
      }
    ),
    sim = function(dpars, aterms,
                   n) stats::rlnorm(n, dpars[["mu"]], dpars[["sigma"]])
  )
}

#' Binet's `mu(x)`: what `log Gamma(x)` has left over its Stirling
#' part, `(x - 1/2) log x - x + log(2 pi) / 2`.
#'
#' The series is asymptotic, so it has a smallest usable argument
#' rather than a largest, and that argument falls fast with the number
#' of terms. Smallest `x` whose truncation is 1e-16 or better, against
#' a 400-bit reference: 28 for four terms, 17 for five, 12 for six, 10
#' for seven, 8 for eight. `lgamma_shift_diff()` shifts by 12, so it
#' never calls this below `x = 12`, where these six terms are 5.8e-17
#' out.
#'
#' The coefficients are `B_{2n} / (2n(2n - 1))`, and each one is
#' confirmed by the truncation it predicts: the remainder after `k`
#' terms has the sign of, and is bounded by, the first omitted term,
#' and the measured error over that term is -0.96 at every `k` from 4
#' to 7 and every `x` from 8 to 28.
#'
#' Six terms with a shift of 12 replaced four with a shift of 25. The
#' pair is no less accurate anywhere measured, over 800 `(nu, z)`
#' points in six bands with neither ordering dominating, and it sits 4
#' times further from its own truncation floor at the smallest
#' argument it can be handed (5.8e-17 at `x = 12` against 2.2e-16 at
#' `x = 25`, where four terms were already inside where they reach
#' 1e-16). What the change buys is cost: 12 recurrence nodes rather
#' than 25 halves the objective when `nu` carries a linear predictor.
#' `dev/reviews/2026-09-08-remlopt.md` has the sweep.
#'
#' @noRd
lgamma_binet <- function(x) {
  ix2 <- 1 / (x * x)
  (1 / 12 - ix2 * (1 / 360 - ix2 * (1 / 1260 - ix2 * (1 / 1680 -
    ix2 * (1 / 1188 - ix2 * (691 / 360360)))))) / x
}

#' `log Gamma(a + s) - log Gamma(a)`, formed so that nothing large
#' cancels.
#'
#' Subtracting the two `lgamma()` values directly loses the answer once
#' `a` is large: at `a = 1e10` each one is about 2.3e11, where a double
#' is spaced 3e-5 apart, and the difference they have to produce is
#' 11.5. Measured against a 300-bit reference the naive difference is
#' wrong by 1.4e-5 at `a = 1e10` and by 57 at `a = 1e50`. That is what
#' turns the flat large-`nu` tail of a student likelihood into noise:
#' the optimizer stops at whichever sign change of a noise-dominated
#' gradient it reaches first, so a permuted row order lands it
#' somewhere else. This form holds 1.3e-14 relative over 369 `(a, s)`
#' points, `a` from 1e-300 to 1e300 crossed with `s` from 1/2 to 100,
#' which covers both the density's `s = 1/2` and the `s = k/2` the
#' multivariate-t blocks ask for.
#'
#' Two steps. Push `a` up with `Gamma(x + 1) = x Gamma(x)`, one exact
#' `log1p` per step, until Binet's remainder series holds. Then take
#' the Stirling difference in a shape whose big terms are never formed:
#' `(a - 1/2) log1p(s/a) + s log(a + s) - s` is the whole of it, and
#' every piece stays the size of the answer.
#'
#' Branchless on purpose, for two reasons, neither of which is that a
#' branch is impossible. RTMB 1.9 exports no `CondExp*`, and a
#' comparison on an advector ERRORS under the default
#' `TapeConfig(comparison = "forbid")` ("Comparison is generally unsafe
#' for AD types"); under `TapeConfig(comparison = "tape")` it tapes
#' correctly, verified on an indicator that switched at the right
#' argument on re-evaluation, with jacobians 1 and 2 either side. But
#' `TapeConfig` is process-global, so a package
#' must not flip it under the user. And the one conditional that needs
#' no setting, a multiplicative blend of the naive and the asymptotic
#' form, evaluates both arms: `lgamma(a)` overflows to `Inf` at
#' `a = 2.533e305`, which is `nu = 5.07e305` and `log(nu - 1) = 703.9`,
#' inside the 709.78 the `logm1` link reaches, so the blend returns
#' `NaN` at a reachable `nu` where this form still gives the gaussian
#' limit to 2.0e-14.
#'
#' `m` is 12 for every `a`, and it cannot be shortened by reading one.
#' The shift needed falls with `a`, but `a` is `nu / 2` and `nu` is a
#' fitted parameter that the optimizer moves across the whole range
#' after taping. The only facts fixed at tape time are the link's, and
#' `logm1` bounds `a > 0.5`, which buys ONE step. The lever that does
#' work is the series length, which depends on no parameter at all:
#' `m` is whatever `lgamma_binet()` needs, so six terms buy 12 where
#' four needed 25.
#'
#' @noRd
lgamma_shift_diff <- function(a, s) {
  # steps of the recurrence that put any positive `a` above the point
  # where lgamma_binet() is accurate. m and the series length are one
  # choice: see lgamma_binet() for the table that pairs them.
  m <- 12L
  acc <- 0
  for (j in seq_len(m) - 1L) acc <- acc + log1p(s / (a + j))
  b <- a + m
  (b - 0.5) * log1p(s / b) + s * log(b + s) - s +
    lgamma_binet(b + s) - lgamma_binet(b) - acc
}

#' Log density of a standard student-t, stable for every `nu`.
#'
#' `RTMB::dt()` sends a double straight to `stats::dt()`, which is
#' accurate, and an AD number to a tape that forms
#' `lgamma((nu + 1) / 2) - lgamma(nu / 2)` as written. Every fit runs on
#' the tape, so the objective the optimizer sees is the inaccurate one;
#' `lgamma_shift_diff()` carries the measurement. `log(nu) + log(pi)`
#' rather than `log(nu * pi)` because the product overflows near the top
#' of the double range while the density there is still defined.
#'
#' @noRd
dt_stable <- function(z, nu) {
  lgamma_shift_diff(nu / 2, 0.5) - (log(nu) + log(pi)) / 2 -
    (nu + 1) / 2 * log1p(z^2 / nu)
}

#' Student-t family, dpars `mu`, `sigma` and `nu`. The `logm1` link on
#' `nu` holds the degrees of freedom above one. A known `se()` term
#' enters the scale, as in the gaussian family.
#'
#' @noRd
fam_student <- function(link = "identity", link_sigma = "log",
                        link_nu = "logm1") {
  lk_sigma <- dpar_link(link_sigma, "sigma", "student", dpar_links_positive)
  lk_nu <- dpar_link(link_nu, "nu", "student", dpar_links_above1)
  frmtmb_family(
    "student",
    accepts_aterms = c("weights", "se", "mi"),
    se_dpar = "sigma",
    dpars = c("mu", "sigma", "nu"),
    links = list(mu = mu_link(link, "student"), sigma = lk_sigma, nu = lk_nu),
    lpdf = function(y, dpars, aterms) {
      sd_t <- resid_sd(dpars[["sigma"]], aterms)
      dt_stable((y - dpars[["mu"]]) / sd_t, dpars[["nu"]]) - log(sd_t)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y),
      nu = function(y, aterms) 4
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) {
        sd_t <- resid_sd(dpars[["sigma"]], aterms)
        ifelse(dpars[["nu"]] > 2, sd_t^2 * dpars[["nu"]] / (dpars[["nu"]] - 2),
               NA_real_)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars[["mu"]] + resid_sd(dpars[["sigma"]], aterms) * stats::rt(n,
        dpars[["nu"]])
    }
  )
}

#' Negative binomial with quadratic variance (brms negbinomial, glmmTMB
#' nbinom2), dpars `mu` and `shape`. The variance is
#' `mu + mu^2 / shape`.
#'
#' @noRd
fam_negbinomial <- function(link = "log", link_shape = "log") {
  lk_shape <- dpar_link(link_shape, "shape", "negbinomial", dpar_links_positive)
  lk <- mu_link(link, "negbinomial")
  frmtmb_family(
    "negbinomial",
    accepts_aterms = c("weights", "rate"),
    dpars = c("mu", "shape"),
    links = list(mu = lk, shape = lk_shape),
    lpdf = function(y, dpars, aterms) {
      # dnbinom2 forms p = mu / var, which is 0 / 0 = NaN once exp(eta)
      # has underflowed; dnbinom_robust takes log(mu) and
      # log(var - mu) = 2 log(mu) - log(shape) and never divides.
      lmu <- robust_logmu(dpars, lk)
      d <- aterms[["rate"]]
      if (is.null(lmu)) {
        mu <- rate_mu(dpars, lk, aterms)
        return(RTMB::dnbinom2(y, mu,
                              mu + mu^2 / rate_shape(dpars[["shape"]],
                                                     aterms),
                              log = TRUE))
      }
      lsh <- dpar_log(dpars, "shape", lk_shape)
      # rate(): brms's neg_binomial_2_log(eta + log d, shape * d)
      if (!is.null(d)) {
        lmu <- lmu + log(d)
        lsh <- lsh + log(d)
      }
      RTMB::dnbinom_robust(y, lmu, 2 * lmu - lsh, log = TRUE)
    },
    valid_y = count_y("negbinomial"),
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      shape = function(y, aterms) {
        m <- mean(y)
        v <- stats::var(y)
        if (v > m) max(m^2 / (v - m), 0.1) else 10
      }
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) rate_mean(dpars, aterms),
      var_fn = function(dpars, aterms) {
        mu <- rate_mean(dpars, aterms)
        mu + mu^2 / rate_shape(dpars[["shape"]], aterms)
      },
      dev_fn = function(y, dpars, aterms) {
        nbinom_deviance(y, rate_mean(dpars, aterms),
                        rate_shape(dpars[["shape"]], aterms))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnbinom(n, size = rate_shape(dpars[["shape"]], aterms),
                     mu = rate_mean(dpars, aterms))
    }
  )
}

#' Negative binomial with linear variance (glmmTMB nbinom1), dpars `mu`
#' and `phi`. The variance is `mu * (1 + phi)`, so `phi` is a
#' quasi-Poisson style dispersion.
#'
#' @noRd
fam_nbinom1 <- function(link = "log", link_phi = "log") {
  lk_phi <- dpar_link(link_phi, "phi", "nbinom1", dpar_links_positive)
  lk <- mu_link(link, "nbinom1")
  frmtmb_family(
    "nbinom1",
    accepts_aterms = "weights",
    dpars = c("mu", "phi"),
    links = list(mu = lk, phi = lk_phi),
    lpdf = function(y, dpars, aterms) {
      # var - mu = mu * phi, so log(var - mu) = log(mu) + log(phi)
      lmu <- robust_logmu(dpars, lk)
      if (is.null(lmu)) {
        return(RTMB::dnbinom2(y, dpars[["mu"]],
               dpars[["mu"]] * (1 + dpars[["phi"]]),
                              log = TRUE))
      }
      RTMB::dnbinom_robust(y, lmu, lmu + dpar_log(dpars, "phi", lk_phi),
                           log = TRUE)
    },
    valid_y = count_y("nbinom1"),
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      phi = function(y, aterms) {
        m <- mean(y)
        v <- stats::var(y)
        max(v / max(m, 0.1) - 1, 0.1)
      }
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) dpars[["mu"]] * (1 + dpars[["phi"]]),
      # glmmTMB's convention: the negative-binomial unit deviance with
      # the size held at the FITTED row's mu / phi. Letting the size
      # follow the saturated mean instead is not a deviance at all - the
      # nbinom1 log-likelihood in mu is not maximized at mu = y once the
      # size moves with it, and the difference goes negative.
      dev_fn = function(y, dpars, aterms) {
        nbinom_deviance(y, dpars[["mu"]], dpars[["mu"]] / dpars[["phi"]])
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnbinom(n,
                     size = dpars[["mu"]] / dpars[["phi"]], mu = dpars[["mu"]])
    }
  )
}

#' Beta family in the mean-precision parameterization, dpars `mu` and
#' `phi`. The two shapes are `mu * phi` and `(1 - mu) * phi`, and the
#' response must lie strictly inside `(0, 1)`.
#'
#' @noRd
fam_beta <- function(link = "logit", link_phi = "log") {
  lk_phi <- dpar_link(link_phi, "phi", "Beta", dpar_links_positive)
  lk <- mu_link(link, "beta")
  frmtmb_family(
    "beta",
    accepts_aterms = "weights",
    dpars = c("mu", "phi"),
    links = list(mu = lk, phi = lk_phi),
    lpdf = function(y, dpars, aterms) {
      # the SECOND shape is (1 - mu) * phi, and 1 - plogis(eta) is
      # exactly 0 past eta = 37: dbeta at shape 0 is -Inf. Taking the
      # pair off the log-odds keeps both shapes strictly positive.
      mp <- dpar_complement(dpars, "mu", lk)
      RTMB::dbeta(y, mp$p * dpars[["phi"]], mp$q * dpars[["phi"]], log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y <= 0) || any(y >= 1)) {
        frm_stop("beta: response must lie strictly in (0, 1)", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      phi = function(y, aterms) {
        m <- mean(y)
        max(m * (1 - m) / stats::var(y) - 1, 0.5)
      }
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) {
        dpars[["mu"]] * (1 - dpars[["mu"]]) / (1 + dpars[["phi"]])
      },
      # 2 (ll(y; y, phi) - ll(y; mu, phi)); the terms in y alone cancel
      # (the betareg deviance-residual definition). glmmTMB returns NA
      # for beta, so there is nothing to match there.
      dev_fn = function(y, dpars, aterms) {
        mu <- dpars[["mu"]]
        ph <- dpars[["phi"]]
        2 * (lgamma(mu * ph) + lgamma((1 - mu) * ph) -
               lgamma(y * ph) - lgamma((1 - y) * ph) +
               (y - mu) * ph * log(y / (1 - y)))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbeta(n, dpars[["mu"]] * dpars[["phi"]],
        (1 - dpars[["mu"]]) * dpars[["phi"]])
    }
  )
}

#' Tweedie family for a non-negative response with a mass at zero,
#' dpars `mu`, `phi` and `power`. The `power12` link holds `power`
#' between 1 and 2, the compound Poisson-gamma range.
#'
#' @noRd
fam_tweedie <- function(link = "log", link_phi = "log") {
  lk_phi <- dpar_link(link_phi, "phi", "tweedie", dpar_links_positive)
  frmtmb_family(
    "tweedie",
    accepts_aterms = "weights",
    dpars = c("mu", "phi", "power"),
    links = list(mu = mu_link(link, "tweedie"),
                 phi = lk_phi, power = "power12"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dtweedie(y, dpars[["mu"]], dpars[["phi"]], dpars[["power"]],
                                       log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y < 0)) {
        frm_stop("tweedie: response must be non-negative", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      phi = function(y, aterms) 1,
      power = function(y, aterms) 1.5
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars,
                        aterms) dpars[["phi"]] * dpars[["mu"]]^dpars[["power"]],
      # the standard Tweedie unit deviance (phi is the dispersion and
      # divides out); at y = 0 the first two terms vanish
      dev_fn = function(y, dpars, aterms) {
        mu <- dpars[["mu"]]
        p <- dpars[["power"]]
        2 * (y^(2 - p) / ((1 - p) * (2 - p)) -
               y * mu^(1 - p) / (1 - p) +
               mu^(2 - p) / (2 - p))
      }
    ),
    # The generative process, not the density inverted. The power12
    # link guarantees 1 < p < 2, and on that range the Tweedie IS a
    # compound Poisson sum of gamma variables, so a draw is a Poisson
    # count of gamma jumps. That is also where its point mass at zero
    # comes from, and it comes out on its own: a row that draws no
    # jumps has a gamma shape of exactly zero, which rgamma() returns
    # as an exact zero rather than something merely small.
    sim = function(dpars, aterms, n) {
      mu <- rep(dpars[["mu"]], length.out = n)
      phi <- rep(dpars[["phi"]], length.out = n)
      p <- rep(dpars[["power"]], length.out = n)
      # the reparameterization that makes the sum have mean mu and
      # variance phi * mu^p
      jumps <- stats::rpois(n, mu^(2 - p) / (phi * (2 - p)))
      stats::rgamma(n, shape = jumps * (2 - p) / (p - 1),
                    scale = phi * (p - 1) * mu^(p - 1))
    }
  )
}

#' The COM-Poisson log rate that puts the distribution's MEAN at `mu`,
#' at dispersion `nu`, on the support `0..ymax`.
#'
#' The mean parameterization is implicit: the rate is whatever makes
#' `sum(y w_y) / sum(w_y)` equal `mu`, for the unnormalized weights
#' `w_y = lambda^y / (y!)^nu`. The mean is strictly increasing in the
#' rate, so bisection converges from any bracket that contains the
#' answer, and the bracket is grown rather than guessed: `nu` well away
#' from 1 moves the rate roughly to `nu * log(mu)`, which leaves a fixed
#' bracket wrong for a large mean at a large dispersion.
#'
#' @noRd
compois_loglambda <- function(mu, nu, ymax) {
  y <- seq.int(0L, ymax)
  lgy <- nu * lgamma(y + 1)
  excess <- function(ll) {
    lw <- y * ll - lgy
    w <- exp(lw - max(lw))
    sum(y * w) / sum(w) - mu
  }
  lo <- -8
  while (excess(lo) > 0 && lo > -800) lo <- lo * 2
  hi <- 8
  while (excess(hi) < 0 && hi < 800) hi <- hi * 2
  # the bracket stops halving once it is narrower than the double it is
  # converging to, which matters because each step costs one pass over
  # the whole support and an overdispersed nu makes that support long
  for (i in seq_len(200L)) {
    if (hi - lo < 1e-12) break
    md <- (lo + hi) / 2
    if (excess(md) > 0) hi <- md else lo <- md
  }
  (lo + hi) / 2
}

#' The COM-Poisson probability vector at `(mu, nu)`, built from the
#' distribution's own weights and truncated where they have decayed to
#' nothing. The support is unbounded, so `ymax` grows until the top
#' weight is `e^-40` below the mode rather than being fixed at a size
#' that a small `nu` (a long right tail) would overrun.
#'
#' @noRd
compois_probs <- function(mu, nu) {
  ymax <- max(64L, as.integer(ceiling(4 * mu / min(nu, 1) + 32)))
  repeat {
    ll <- compois_loglambda(mu, nu, ymax)
    y <- seq.int(0L, ymax)
    lw <- y * ll - nu * lgamma(y + 1)
    lw <- lw - max(lw)
    if (lw[ymax + 1L] < -40 || ymax >= 262144L) break
    ymax <- ymax * 2L
  }
  w <- exp(lw)
  list(y = y, p = w / sum(w))
}

#' Conway-Maxwell-Poisson family for counts that are under- or
#' over-dispersed, dpars `mu` (the mean) and `nu` (the dispersion,
#' with `nu = 1` giving the Poisson).
#'
#' @noRd
fam_compois <- function(link = "log", link_nu = "log") {
  lk_nu <- dpar_link(link_nu, "nu", "compois", dpar_links_positive)
  frmtmb_family(
    "compois",
    accepts_aterms = "weights",
    dpars = c("mu", "nu"),
    links = list(mu = mu_link(link, "compois"), nu = lk_nu),
    lpdf = function(y, dpars, aterms) {
      RTMB::dcompois2(y, dpars[["mu"]], dpars[["nu"]], log = TRUE)
    },
    valid_y = count_y("compois"),
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      nu = function(y, aterms) 1
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]]
    ),
    # The COM-Poisson has no constructive generative process the way a
    # Poisson or a gamma does: every published sampler works from the
    # unnormalized weights. This one does the same, building the
    # distribution from its defining weights and its own solve for the
    # rate that hits the requested MEAN, which is the convention the
    # density's mean parameterization claims and therefore the thing
    # worth checking independently.
    #
    # Rows sharing a parameter pair share one solve, because the solve
    # is the expensive part. That pays on a factor design, where a
    # handful of pairs cover every row. A continuous predictor gives
    # every row its own pair and the cache never hits: n = 1000 costs
    # about 1 s, linearly in n, against milliseconds for poisson() on
    # the same design. Solving on a grid of log(mu) per distinct nu and
    # interpolating would fix that, since the mean is smooth and
    # strictly monotone in the rate at fixed nu.
    sim = function(dpars, aterms, n) {
      mu <- rep(dpars[["mu"]], length.out = n)
      nu <- rep(dpars[["nu"]], length.out = n)
      out <- integer(n)
      key <- paste(mu, nu)
      for (k in unique(key)) {
        idx <- which(key == k)
        g <- compois_probs(mu[idx[1L]], nu[idx[1L]])
        out[idx] <- sample(g[["y"]], length(idx), replace = TRUE,
                           prob = g[["p"]])
      }
      out
    }
  )
}

#' Zero-inflated Poisson, dpars `mu` (the Poisson mean) and `zi` (the
#' probability of a structural zero).
#'
#' @noRd
fam_zi_poisson <- function(link = "log", link_zi = "logit") {
  lk_zi <- dpar_link(link_zi, "zi", "zero_inflated_poisson", dpar_links_unit)
  frmtmb_family(
    "zero_inflated_poisson",
    accepts_aterms = "weights",
    dpars = c("mu", "zi"),
    links = list(mu = mu_link(link, "zero_inflated_poisson"), zi = lk_zi),
    lpdf = function(y, dpars, aterms) {
      # y == 0 is data, so the mixture stays branch-free in parameters.
      # The whole mixture runs in log space: log(1 - zi) is -Inf once
      # the zi predictor separates, and the Poisson's own log P(0) is
      # -mu exactly, with no exp() to underflow.
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "zi", lk_zi)
      i0 * RTMB::logspace_add(g$l, g$l1m - dpars[["mu"]]) +
        (1 - i0) * (g$l1m + RTMB::dpois(y, dpars[["mu"]], log = TRUE))
    },
    valid_y = count_y("zero_inflated_poisson"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars[["zi"]]) * dpars[["mu"]],
      var_fn = function(dpars, aterms) {
        (1 - dpars[["zi"]]) * dpars[["mu"]] * (1 + dpars[["zi"]] *
          dpars[["mu"]])
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars[["zi"]]) * stats::rpois(n, dpars[["mu"]])
    }
  )
}

#' Zero-inflated negative binomial (nbinom2 variance), dpars `mu`,
#' `shape` and `zi` (the probability of a structural zero).
#'
#' @noRd
fam_zi_negbinomial <- function(link = "log", link_shape = "log",
                               link_zi = "logit") {
  lk_shape <- dpar_link(
    link_shape, "shape", "zero_inflated_negbinomial", dpar_links_positive)
  lk_zi <- dpar_link(
    link_zi, "zi", "zero_inflated_negbinomial", dpar_links_unit)
  lk <- mu_link(link, "zero_inflated_negbinomial")
  frmtmb_family(
    "zero_inflated_negbinomial",
    accepts_aterms = "weights",
    dpars = c("mu", "shape", "zi"),
    links = list(mu = lk, shape = lk_shape, zi = lk_zi),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "zi", lk_zi)
      lmu <- robust_logmu(dpars, lk)
      if (is.null(lmu)) {
        lp0 <- dpars[["shape"]] * (log(dpars[["shape"]]) -
                                log(dpars[["shape"]] + dpars[["mu"]]))
        base <- RTMB::dnbinom2(y, dpars[["mu"]],
                               dpars[["mu"]] + dpars[["mu"]]^2 /
                                 dpars[["shape"]],
                               log = TRUE)
      } else {
        # log P(0) = shape * log(shape / (shape + mu)), which is
        # -shape * log(1 + mu / shape) with the ratio taken in logs
        lsh <- dpar_log(dpars, "shape", lk_shape)
        lp0 <- -dpars[["shape"]] *
          RTMB::logspace_add(0 * lmu, lmu - lsh)
        base <- RTMB::dnbinom_robust(y, lmu, 2 * lmu - lsh, log = TRUE)
      }
      i0 * RTMB::logspace_add(g$l, g$l1m + lp0) + (1 - i0) * (g$l1m + base)
    },
    valid_y = count_y("zero_inflated_negbinomial"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      shape = function(y, aterms) 1,
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars[["zi"]]) * dpars[["mu"]]
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars[["zi"]]) *
        stats::rnbinom(n, size = dpars[["shape"]], mu = dpars[["mu"]])
    }
  )
}

#' Hurdle Poisson, dpars `mu` and `hu`. `hu` is the probability of a
#' zero and the positive part is a zero-truncated Poisson, so zeros come
#' only from the hurdle.
#'
#' @noRd
fam_hurdle_poisson <- function(link = "log", link_hu = "logit") {
  lk_hu <- dpar_link(link_hu, "hu", "hurdle_poisson", dpar_links_unit)
  frmtmb_family(
    "hurdle_poisson",
    accepts_aterms = "weights",
    dpars = c("mu", "hu"),
    links = list(mu = mu_link(link, "hurdle_poisson"), hu = lk_hu),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "hu", lk_hu)
      # nonzero part is a zero-truncated poisson. Its normalizer
      # log(1 - exp(-mu)) cancels to log(0) once mu underflows below
      # the double epsilon; expm1 keeps it, and the truth there is
      # simply log(mu).
      lztrunc <- log(-expm1(-dpars[["mu"]]))
      i0 * g$l +
        (1 - i0) * (g$l1m + RTMB::dpois(y, dpars[["mu"]], log = TRUE) -
                      lztrunc)
    },
    valid_y = count_y("hurdle_poisson"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars[["hu"]]) * dpars[["mu"]] / (1 - exp(-dpars[["mu"]]))
      }
    ),
    # The generative process the family's name states: clear the hurdle
    # with probability 1 - hu, and if you do, draw a count that cannot
    # be zero.
    #
    # The positive part is drawn by inverse transform on the Poisson
    # CDF ABOVE its own zero, not by drawing Poisson variates until one
    # comes out positive. Both are correct; only this one has a bounded
    # cost, and the rejection loop degenerates exactly where hurdle
    # models are used, at a small mu, where almost every plain draw is
    # the zero the hurdle already accounts for.
    sim = function(dpars, aterms, n) {
      mu <- rep(dpars[["mu"]], length.out = n)
      hu <- rep(dpars[["hu"]], length.out = n)
      p0 <- exp(-mu)
      # The positive part's support starts at 1, so clamp it there.
      # Two things push qpois() under that floor. The pmin() keeps the
      # probability strictly below one, and at a tiny mu the double
      # grid just under one is too coarse to resolve
      # p0 + u (1 - p0), because 1 - p0 is about mu: below mu of about
      # 1e-14 the smallest representable draws round down to zero,
      # which is the hurdle's own atom rather than the positive part's.
      u <- p0 + stats::runif(n) * (1 - p0)
      (1 - stats::rbinom(n, 1L, hu)) *
        pmax(stats::qpois(pmin(u, 1 - .Machine$double.eps), mu), 1L)
    }
  )
}

#' Hurdle negative binomial (nbinom2 variance), dpars `mu`, `shape` and
#' `hu`. As in `fam_hurdle_poisson()`, `hu` is the probability of a zero
#' and the positive part is the negative binomial truncated at zero.
#'
#' @noRd
fam_hurdle_negbinomial <- function(link = "log", link_shape = "log",
                                   link_hu = "logit") {
  lk_shape <- dpar_link(
    link_shape, "shape", "hurdle_negbinomial", dpar_links_positive)
  lk_hu <- dpar_link(link_hu, "hu", "hurdle_negbinomial", dpar_links_unit)
  lk <- mu_link(link, "hurdle_negbinomial")
  frmtmb_family(
    "hurdle_negbinomial",
    accepts_aterms = "weights",
    dpars = c("mu", "shape", "hu"),
    links = list(mu = lk, shape = lk_shape, hu = lk_hu),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "hu", lk_hu)
      # One expression for every mean link. Off the log link, log(mu)
      # is taken here rather than handing dnbinom2() the variance: that
      # forms var - mu by subtraction, which at mu = exp(-25) keeps
      # about five digits of mu^2 / shape (dev/fams-validate.R, 1)
      lmu <- robust_logmu(dpars, lk) %||% log(dpars[["mu"]])
      lsh <- dpar_log(dpars, "shape", lk_shape)
      base <- RTMB::dnbinom_robust(y, lmu, 2 * lmu - lsh, log = TRUE)
      # log P(0) = -shape * log(1 + mu / shape), with the ratio in logs
      lp0 <- -dpars[["shape"]] * RTMB::logspace_add(0 * lmu, lmu - lsh)
      # The truncation normalizer log(1 - P(0)). P(0) tends to one as mu
      # falls, where 1 - exp(lp0) cancels; logspace_sub() switches to
      # expm1() there. brms writes log1m((shape / (mu + shape))^shape),
      # the same quantity with the cancellation left in.
      i0 * g$l +
        (1 - i0) * (g$l1m + base - RTMB::logspace_sub(0 * lp0, lp0))
    },
    valid_y = count_y("hurdle_negbinomial"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      shape = function(y, aterms) 1,
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "discrete",
    post = list(
      # brms's posterior_epred_hurdle_negbinomial(), with 1 - P(0)
      # formed without the subtraction, as in the density
      mean_fn = function(dpars, aterms) {
        (1 - dpars[["hu"]]) * dpars[["mu"]] /
          nb_mass_above_zero(dpars[["mu"]], dpars[["shape"]])
      },
      var_fn = function(dpars, aterms) {
        mu <- dpars[["mu"]]
        q <- (1 - dpars[["hu"]]) / nb_mass_above_zero(mu, dpars[["shape"]])
        # E[Y^2] of the untruncated NB is mu + mu^2 (1 + 1 / shape)
        q * (mu + mu^2 * (1 + 1 / dpars[["shape"]])) - (q * mu)^2
      }
    ),
    # Inverse transform on the NB CDF above its own zero, for the reason
    # fam_hurdle_poisson() gives: a rejection loop degenerates at the
    # small mu where hurdle models are used. brms's own
    # posterior_predict_hurdle_negbinomial() draws rnbinom(mu - t) + 1,
    # which is not the zero-truncated NB, so frmtmb does not copy it.
    sim = function(dpars, aterms, n) {
      mu <- rep(dpars[["mu"]], length.out = n)
      shape <- rep(dpars[["shape"]], length.out = n)
      hu <- rep(dpars[["hu"]], length.out = n)
      p0 <- 1 - nb_mass_above_zero(mu, shape)
      u <- p0 + stats::runif(n) * (1 - p0)
      (1 - stats::rbinom(n, 1L, hu)) *
        pmax(stats::qnbinom(pmin(u, 1 - .Machine$double.eps),
                            size = shape, mu = mu), 1L)
    }
  )
}

#' `1 - P(Y = 0)` of a negative binomial with mean `mu` and size `shape`,
#' off the tape. Formed without subtracting from one, so a hurdle mean
#' at a small `mu` keeps its digits.
#'
#' @noRd
nb_mass_above_zero <- function(mu, shape) -expm1(-shape * log1p(mu / shape))

#' Bernoulli family for a 0/1 response, single dpar `mu` (the success
#' probability). Robust in the linear predictor, as [fam_binomial()]
#' describes; separation is this family's normal failure mode, so the
#' saturating form is exactly the wrong one here.
#'
#' @noRd
fam_bernoulli <- function(link = "logit") {
  lk <- mu_link(link, "bernoulli")
  frmtmb_family(
    "bernoulli",
    accepts_aterms = "weights",
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      lo <- robust_logit(dpars, lk)
      if (is.null(lo)) return(RTMB::dbinom(y, 1, dpars[["mu"]], log = TRUE))
      RTMB::dbinom_robust(y, 1, lo, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (!all(y %in% c(0, 1))) {
        frm_stop("bernoulli: response must be 0/1", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) min(max(mean(y), 0.02), 0.98)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) dpars[["mu"]] * (1 - dpars[["mu"]]),
      dev_fn = function(y, dpars, aterms) {
        binomial_deviance(y, dpars[["mu"]], 1)
      }
    ),
    sim = function(dpars, aterms, n) stats::rbinom(n, 1, dpars[["mu"]])
  )
}

#' Geometric family, single dpar `mu` (the mean). It is the negative
#' binomial with the shape held at one. Under `rate(d)` brms multiplies
#' that shape by `d` as well as the mean, so the density becomes the
#' negative binomial of size `d`.
#'
#' @noRd
fam_geometric <- function(link = "log") {
  lk <- mu_link(link, "geometric")
  frmtmb_family(
    "geometric",
    accepts_aterms = c("weights", "rate"),
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      # var - mu = mu^2 (the negative binomial at shape 1)
      lmu <- robust_logmu(dpars, lk)
      d <- aterms[["rate"]]
      if (!is.null(d)) {
        # brms's neg_binomial_2_log(eta + log d, 1 .* d)
        if (is.null(lmu)) {
          mu <- rate_mu(dpars, lk, aterms)
          return(RTMB::dnbinom2(y, mu, mu + mu^2 / d, log = TRUE))
        }
        lmu <- lmu + log(d)
        return(RTMB::dnbinom_robust(y, lmu, 2 * lmu - log(d), log = TRUE))
      }
      if (is.null(lmu)) {
        return(RTMB::dnbinom2(y, dpars[["mu"]],
               dpars[["mu"]] * (1 + dpars[["mu"]]),
                              log = TRUE))
      }
      RTMB::dnbinom_robust(y, lmu, 2 * lmu, log = TRUE)
    },
    valid_y = count_y("geometric"),
    init_dpars = list(mu = function(y, aterms) mean(y) + 0.1),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) rate_mean(dpars, aterms),
      var_fn = function(dpars, aterms) {
        mu <- rate_mean(dpars, aterms)
        mu * (1 + mu / rate_shape(1, aterms))
      },
      # negbinomial with the shape fixed at 1
      dev_fn = function(y, dpars, aterms) {
        nbinom_deviance(y, rate_mean(dpars, aterms), rate_shape(1, aterms))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnbinom(n, size = rate_shape(1, aterms),
                     mu = rate_mean(dpars, aterms))
    }
  )
}

#' Exponential family in the mean parameterization, single dpar `mu`
#' (the mean, that is one over the rate). It supplies a CDF and a
#' truncated mean.
#'
#' @noRd
fam_exponential <- function(link = "log") {
  frmtmb_family(
    "exponential",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = "mu",
    links = list(mu = mu_link(link, "exponential")),
    lpdf = function(y, dpars, aterms) {
      -log(dpars[["mu"]]) - y / dpars[["mu"]]
    },
    lcdf = function(q, dpars, aterms) {
      1 - exp(-q / dpars[["mu"]])
    },
    # log S = -q / mu in closed form: no complement is ever formed
    lccdf = function(q, dpars, aterms) {
      -q / dpars[["mu"]]
    },
    valid_y = positive_y("exponential"),
    init_dpars = list(mu = function(y, aterms) mean(y)),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) dpars[["mu"]]^2,
      # Gamma with the shape fixed at 1
      dev_fn = function(y, dpars, aterms) gamma_deviance(y, dpars[["mu"]]),
      # int_lb^ub y f(y) dy = (lb + mu) e^{-lb/mu} - (ub + mu) e^{-ub/mu}
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        mu <- dpars[["mu"]]
        lo <- pmax(lb, 0)
        el <- exp(-lo / mu)
        eu <- exp(-ub / mu)
        # Inf * 0 at an absent upper bound; the term is zero there
        hi_term <- ifelse(is.finite(ub), (ub + mu) * eu, 0)
        ((lo + mu) * el - hi_term) / (el - eu)
      }
    ),
    sim = function(dpars, aterms, n) stats::rexp(n, 1 / dpars[["mu"]])
  )
}

#' Weibull family, dpars `mu` and `shape`. As in brms, `mu` is the mean
#' and the scale is `mu / gamma(1 + 1 / shape)`. It supplies a CDF and a
#' truncated mean.
#'
#' @noRd
fam_weibull <- function(link = "log", link_shape = "log") {
  lk_shape <- dpar_link(link_shape, "shape", "weibull", dpar_links_positive)
  frmtmb_family(
    "weibull",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = c("mu", "shape"),
    links = list(mu = mu_link(link, "weibull"), shape = lk_shape),
    lpdf = function(y, dpars, aterms) {
      # brms parameterization: mu is the mean, scale = mu/gamma(1+1/k)
      sc <- dpars[["mu"]] / exp(lgamma(1 + 1 / dpars[["shape"]]))
      RTMB::dweibull(y, shape = dpars[["shape"]], scale = sc, log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      sc <- dpars[["mu"]] / exp(lgamma(1 + 1 / dpars[["shape"]]))
      1 - exp(-(q / sc)^dpars[["shape"]])
    },
    lccdf = function(q, dpars, aterms) {
      sc <- dpars[["mu"]] / exp(lgamma(1 + 1 / dpars[["shape"]]))
      -(q / sc)^dpars[["shape"]]
    },
    valid_y = positive_y("weibull"),
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      shape = function(y, aterms) 1.2
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) {
        g1 <- exp(lgamma(1 + 1 / dpars[["shape"]]))
        g2 <- exp(lgamma(1 + 2 / dpars[["shape"]]))
        dpars[["mu"]]^2 * (g2 / g1^2 - 1)
      },
      # int_lb^ub y f(y) dy = scale * gamma(1 + 1/k) * (P(1 + 1/k, zu) -
      # P(1 + 1/k, zl)) with z = (y/scale)^k; scale * gamma(1 + 1/k) = mu
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        k <- dpars[["shape"]]
        sc <- dpars[["mu"]] / exp(lgamma(1 + 1 / k))
        zl <- (pmax(lb, 0) / sc)^k
        zu <- (ub / sc)^k
        dpars[["mu"]] * (stats::pgamma(zu, 1 + 1 / k) -
                      stats::pgamma(zl, 1 + 1 / k)) /
          (exp(-zl) - exp(-zu))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rweibull(n, shape = dpars[["shape"]],
                      scale = dpars[["mu"]] / gamma(1 + 1 / dpars[["shape"]]))
    }
  )
}

#' Shifted lognormal for response times, dpars `mu`, `sigma` (both on
#' the log scale) and `ndt`, the non-decision time the distribution
#' starts at.
#'
#' @noRd
fam_shifted_lognormal <- function(link = "identity", link_sigma = "log",
                                  link_ndt = "log") {
  lk_sigma <- dpar_link(
    link_sigma, "sigma", "shifted_lognormal", dpar_links_positive)
  lk_ndt <- dpar_link(link_ndt, "ndt", "shifted_lognormal", dpar_links_positive)
  frmtmb_family(
    "shifted_lognormal",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma", "ndt"),
    links = list(mu = mu_link(link, "shifted_lognormal"),
                 sigma = lk_sigma, ndt = lk_ndt),
    lpdf = function(y, dpars, aterms) {
      # y <= ndt gives NaN, which the optimizer treats as a rejected
      # step; the ndt init keeps the start feasible
      RTMB::dnorm(log(y - dpars[["ndt"]]), dpars[["mu"]], dpars[["sigma"]],
                                                      log = TRUE) -
        log(y - dpars[["ndt"]])
    },
    valid_y = positive_y("shifted_lognormal"),
    init_dpars = list(
      mu = function(y, aterms) mean(log(y - min(y) / 2)),
      sigma = function(y, aterms) stats::sd(log(y - min(y) / 2)),
      ndt = function(y, aterms) min(y) / 2
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        dpars[["ndt"]] + exp(dpars[["mu"]] + dpars[["sigma"]]^2 / 2)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars[["ndt"]] + stats::rlnorm(n, dpars[["mu"]], dpars[["sigma"]])
    }
  )
}

#' Hurdle gamma for a non-negative response, dpars `mu`, `shape` and
#' `hu`. `hu` is the probability of an exact zero and the positive part
#' is a mean-parameterized gamma.
#'
#' @noRd
fam_hurdle_gamma <- function(link = "log", link_shape = "log",
                             link_hu = "logit") {
  lk_shape <- dpar_link(
    link_shape, "shape", "hurdle_gamma", dpar_links_positive)
  lk_hu <- dpar_link(link_hu, "hu", "hurdle_gamma", dpar_links_unit)
  frmtmb_family(
    "hurdle_gamma",
    accepts_aterms = "weights",
    dpars = c("mu", "shape", "hu"),
    links = list(mu = mu_link(link, "hurdle_gamma"),
                 shape = lk_shape, hu = lk_hu),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      yp <- y + i0   # dodge dgamma(0) = -Inf; the term carries weight 0
      g <- dpar_log_complement(dpars, "hu", lk_hu)
      i0 * g$l +
        (1 - i0) * (g$l1m +
                      RTMB::dgamma(yp, shape = dpars[["shape"]],
                                   scale = dpars[["mu"]] / dpars[["shape"]],
                                   log = TRUE))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0)) {
        frm_stop("hurdle_gamma: response must be non-negative", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      shape = function(y, aterms) 1,
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars[["hu"]]) * dpars[["mu"]]
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars[["hu"]]) *
        stats::rgamma(n, shape = dpars[["shape"]],
                      scale = dpars[["mu"]] / dpars[["shape"]])
    }
  )
}

#' Hurdle lognormal for a non-negative response, dpars `mu`, `sigma`
#' (both on the log scale) and `hu`, the probability of an exact zero.
#'
#' @noRd
fam_hurdle_lognormal <- function(link = "identity", link_sigma = "log",
                                 link_hu = "logit") {
  lk_sigma <- dpar_link(
    link_sigma, "sigma", "hurdle_lognormal", dpar_links_positive)
  lk_hu <- dpar_link(link_hu, "hu", "hurdle_lognormal", dpar_links_unit)
  frmtmb_family(
    "hurdle_lognormal",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma", "hu"),
    links = list(mu = mu_link(link, "hurdle_lognormal"),
                 sigma = lk_sigma, hu = lk_hu),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      yp <- y + i0
      g <- dpar_log_complement(dpars, "hu", lk_hu)
      i0 * g$l +
        (1 - i0) * (g$l1m +
                      RTMB::dnorm(log(yp), dpars[["mu"]], dpars[["sigma"]],
                                  log = TRUE) - log(yp))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0)) {
        frm_stop("hurdle_lognormal: response must be non-negative",
                 call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        if (any(y > 0)) mean(log(y[y > 0])) else 0
      },
      sigma = function(y, aterms) {
        if (sum(y > 0) > 1) stats::sd(log(y[y > 0])) else 1
      },
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars[["hu"]]) * exp(dpars[["mu"]] + dpars[["sigma"]]^2 / 2)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars[["hu"]]) *
        stats::rlnorm(n, dpars[["mu"]], dpars[["sigma"]])
    }
  )
}

#' Zero-inflated binomial, dpars `mu` (the success probability) and `zi`
#' (the probability of a structural zero). Trials come from the
#' `trials()` addition term.
#'
#' @noRd
fam_zi_binomial <- function(link = "logit", link_zi = "logit") {
  lk_zi <- dpar_link(link_zi, "zi", "zero_inflated_binomial", dpar_links_unit)
  lk <- mu_link(link, "zero_inflated_binomial")
  frmtmb_family(
    "zero_inflated_binomial",
    accepts_aterms = c("weights", "trials"),
    dpars = c("mu", "zi"),
    links = list(mu = lk, zi = lk_zi),
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]] %||% 1
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "zi", lk_zi)
      lo <- robust_logit(dpars, lk)
      if (is.null(lo)) {
        lp0 <- size * log(1 - dpars[["mu"]])
        base <- RTMB::dbinom(y, size, dpars[["mu"]], log = TRUE)
      } else {
        # log P(0) = size * log(1 - mu), both terms off the log-odds
        lp0 <- size * log1m_inv_logit(lo)
        base <- RTMB::dbinom_robust(y, size, lo, log = TRUE)
      }
      i0 * RTMB::logspace_add(g$l, g$l1m + lp0) + (1 - i0) * (g$l1m + base)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]] %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        frm_stop("zero_inflated_binomial: response must be integer counts ",
                 "in [0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        size <- aterms[["trials"]] %||% 1
        min(max(mean(y / size), 0.05), 0.95)
      },
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars[["zi"]]) * dpars[["mu"]] * (aterms[["trials"]] %||% 1)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars[["zi"]]) *
        stats::rbinom(n, aterms[["trials"]] %||% 1, dpars[["mu"]])
    }
  )
}

#' Zero-inflated beta for a response in `[0, 1)`, dpars `mu`, `phi` and
#' `zi`, the probability of an exact zero.
#'
#' @noRd
fam_zi_beta <- function(link = "logit", link_phi = "log", link_zi = "logit") {
  lk_phi <- dpar_link(
    link_phi, "phi", "zero_inflated_beta", dpar_links_positive)
  lk_zi <- dpar_link(link_zi, "zi", "zero_inflated_beta", dpar_links_unit)
  lk <- mu_link(link, "zero_inflated_beta")
  frmtmb_family(
    "zero_inflated_beta",
    accepts_aterms = "weights",
    dpars = c("mu", "phi", "zi"),
    links = list(mu = lk, phi = lk_phi, zi = lk_zi),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      ya <- y + i0 * 0.5   # dodge dbeta(0) = -Inf; term carries weight 0
      g <- dpar_log_complement(dpars, "zi", lk_zi)
      mp <- dpar_complement(dpars, "mu", lk)
      i0 * g$l +
        (1 - i0) * (g$l1m +
                      RTMB::dbeta(ya, mp$p * dpars[["phi"]],
                                  mp$q * dpars[["phi"]], log = TRUE))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y >= 1)) {
        frm_stop("zero_inflated_beta: response must be in [0, 1)",
                 call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        if (any(y > 0)) min(max(mean(y[y > 0]), 0.05), 0.95) else 0.5
      },
      phi = function(y, aterms) 5,
      zi = function(y, aterms) min(max(mean(y == 0), 0.05), 0.9)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars[["zi"]]) * dpars[["mu"]]
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars[["zi"]]) *
        stats::rbeta(n, dpars[["mu"]] * dpars[["phi"]],
          (1 - dpars[["mu"]]) * dpars[["phi"]])
    }
  )
}

#' Zero-one-inflated beta for a response in `[0, 1]`, dpars `mu`, `phi`,
#' `zoi` and `coi`. `zoi` is the probability of an exact 0 or 1, and
#' `coi` is the probability that such a value is 1.
#'
#' @noRd
fam_zoi_beta <- function(link = "logit", link_phi = "log",
                         link_zoi = "logit", link_coi = "logit") {
  nm <- "zero_one_inflated_beta"
  lk_phi <- dpar_link(link_phi, "phi", nm, dpar_links_positive)
  lk_zoi <- dpar_link(link_zoi, "zoi", nm, dpar_links_unit)
  lk_coi <- dpar_link(link_coi, "coi", nm, dpar_links_unit)
  lk <- mu_link(link, nm)
  frmtmb_family(
    nm,
    accepts_aterms = "weights",
    dpars = c("mu", "phi", "zoi", "coi"),
    links = list(mu = lk, phi = lk_phi, zoi = lk_zoi, coi = lk_coi),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      i1 <- as.numeric(y == 1)
      ib <- i0 + i1
      # dodge dbeta() at 0 and 1, where it is -Inf; the term carries
      # weight 0 on those rows
      ya <- y + 0.5 * (i0 - i1)
      gz <- dpar_log_complement(dpars, "zoi", lk_zoi)
      gc <- dpar_log_complement(dpars, "coi", lk_coi)
      mp <- dpar_complement(dpars, "mu", lk)
      ib * gz$l + i0 * gc$l1m + i1 * gc$l +
        (1 - ib) * (gz$l1m +
                      RTMB::dbeta(ya, mp$p * dpars[["phi"]],
                                  mp$q * dpars[["phi"]], log = TRUE))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y > 1)) {
        frm_stop(nm, ": response must be in [0, 1]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        yi <- y[y > 0 & y < 1]
        if (length(yi)) min(max(mean(yi), 0.05), 0.95) else 0.5
      },
      phi = function(y, aterms) 5,
      zoi = function(y, aterms) min(max(mean(y == 0 | y == 1), 0.05), 0.9),
      coi = function(y, aterms) {
        yb <- y[y == 0 | y == 1]
        if (length(yb)) min(max(mean(yb), 0.05), 0.95) else 0.5
      }
    ),
    type = "continuous",
    post = list(
      # brms's posterior_epred_zero_one_inflated_beta()
      mean_fn = function(dpars, aterms) {
        dpars[["zoi"]] * dpars[["coi"]] +
          (1 - dpars[["zoi"]]) * dpars[["mu"]]
      },
      var_fn = function(dpars, aterms) {
        zoi <- dpars[["zoi"]]
        mu <- dpars[["mu"]]
        m <- zoi * dpars[["coi"]] + (1 - zoi) * mu
        # E[Y^2]: the atom at 1 gives zoi * coi, the beta part its
        # variance mu (1 - mu) / (1 + phi) plus mu^2
        zoi * dpars[["coi"]] +
          (1 - zoi) * (mu * (1 - mu) / (1 + dpars[["phi"]]) + mu^2) - m^2
      },
      fit_check = zoi_beta_fit_check
    ),
    sim = function(dpars, aterms, n) {
      mu <- rep(dpars[["mu"]], length.out = n)
      phi <- rep(dpars[["phi"]], length.out = n)
      edge <- stats::rbinom(n, 1L, dpars[["zoi"]])
      one <- stats::rbinom(n, 1L, dpars[["coi"]])
      ifelse(edge == 1L, one, stats::rbeta(n, mu * phi, (1 - mu) * phi))
    }
  )
}

#' The zero-one-inflated beta's fit-end check: whether the data can
#' identify `coi` at all.
#'
#' `coi` is scored on the rows at exactly 0 or 1 and on nothing else.
#' With no such row its likelihood is flat, the Hessian is singular, and
#' every standard error of the fit comes back `NaN`. `frm()` itself
#' said nothing; `vcov()` and `summary()` warn later, naming `zoi` and
#' `coi` as flat but not why. With rows at one end only, `coi` runs to
#' that boundary without any warning. brms fits both through its
#' `beta(1, 1)` prior on `coi`; maximum likelihood cannot, so the fit
#' says which case it is in and what to use instead. A `coi` held at a
#' constant, `bf(coi = 0.5)`, has nothing to identify.
#'
#' @noRd
zoi_beta_fit_check <- function(fit, resp) {
  lp <- fit$frame[["linpreds"]][[linpred_key(resp, "coi")]]
  if (is.null(lp) || !is.null(lp[["constant"]])) return(invisible(NULL))
  y <- fit$frame[["y"]][[resp]]
  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  if (n0 > 0L && n1 > 0L) return(invisible(NULL))
  what <- if (n0 + n1 == 0L) {
    paste0("has no exact 0 or 1, so zoi goes to 0, coi has no data, the ",
           "Hessian is singular and no standard error of the fit is ",
           "usable. Use Beta(), or hold coi at a value with ",
           "bf(coi = 0.5)")
  } else {
    paste0("has exact ", if (n0 > 0L) "0s but no 1" else "1s but no 0",
           ", so coi goes to ", if (n0 > 0L) "0" else "1", " and its ",
           "standard error is not usable. Use zero_inflated_beta() on ",
           if (n0 > 0L) "y" else "1 - y", ", or hold coi at a value ",
           "with bf(coi = )")
  }
  frm_warning("zero_one_inflated_beta: the response ",
              if (length(fit$spec$responses) > 1L) paste0(resp, " "),
              what, call. = FALSE)
  invisible(NULL)
}

#' Extended-support beta (Kosmidis and Zeileis 2024), brms's `xbeta`,
#' dpars `mu`, `phi` and `kappa`, for a response in `[0, 1]`.
#'
#' A latent `Z ~ Beta(mu phi, (1 - mu) phi)` is stretched to
#' `(1 + 2 kappa) Z - kappa` and censored at 0 and 1, so the ends carry
#' the latent mass beyond them: `P(Y = 0) = I_q(a, b)` and
#' `P(Y = 1) = I_q(b, a)` with `q = kappa / (1 + 2 kappa)`, the second
#' because `1 - (1 + kappa) / (1 + 2 kappa)` is that same `q`. `mu` is
#' the latent mean, not the mean of `Y`; `post$mean_fn` is the mean of
#' `Y`, brms's `posterior_epred_xbeta()`.
#'
#' @noRd
fam_xbeta <- function(link = "logit", link_phi = "log", link_kappa = "log") {
  nm <- "xbeta"
  lk_phi <- dpar_link(link_phi, "phi", nm, dpar_links_positive)
  lk_kappa <- dpar_link(link_kappa, "kappa", nm, dpar_links_positive)
  lk <- mu_link(link, nm)
  frmtmb_family(
    nm,
    accepts_aterms = "weights",
    dpars = c("mu", "phi", "kappa"),
    links = list(mu = lk, phi = lk_phi, kappa = lk_kappa),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y <= 0)
      i1 <- as.numeric(y >= 1)
      ib <- i0 + i1
      phi <- dpars[["phi"]]
      kap <- dpars[["kappa"]]
      mp <- dpar_complement(dpars, "mu", lk)
      a <- mp$p * phi
      b <- mp$q * phi
      d <- 1 + 2 * kap
      # a boundary row is moved to the middle for the interior term, which
      # carries weight 0 there, so that the density is not read at an edge
      z <- (y * (1 - ib) + 0.5 * ib + kap) / d
      # The beta log density written out rather than RTMB::dbeta(): with
      # its first argument on the tape (z moves with kappa) that returns
      # a NaN gradient once the shapes pass about 1e3, and every xbeta fit
      # that reached phi 1e4 stopped at "NA/NaN gradient evaluation"
      # (dev/fams2-p1-dbeta.R, dev/reviews/2026-09-29-fams2.md, B1)
      out <- (1 - ib) * ((a - 1) * log(z) + (b - 1) * log1p(-z) -
                           lbeta_ad(a, b) - log(d))
      # The incomplete beta is taken on the boundary rows alone, at about
      # 15 us a row in a gradient sweep; log_ibeta_half() says why it is
      # not RTMB::pbeta().
      w0 <- which(i0 == 1)
      w1 <- which(i1 == 1)
      if (length(w0) || length(w1)) {
        # scalars broadcast to one value per row before rows are picked
        q <- kap / d + 0 * y
        a <- a + 0 * y
        b <- b + 0 * y
        if (length(w0)) out[w0] <- log_ibeta_half(q[w0], a[w0], b[w0])
        if (length(w1)) out[w1] <- log_ibeta_half(q[w1], b[w1], a[w1])
      }
      out
    },
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y > 1)) {
        frm_stop(nm, ": response must be in [0, 1], where 0 and 1 are ",
                 "the latent mass beyond each end", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) min(max(mean(y), 0.05), 0.95),
      phi = function(y, aterms) 5,
      kappa = function(y, aterms) 0.1
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        xbeta_moments(dpars[["mu"]], dpars[["phi"]], dpars[["kappa"]])$m1
      },
      var_fn = function(dpars, aterms) {
        m <- xbeta_moments(dpars[["mu"]], dpars[["phi"]], dpars[["kappa"]])
        m$m2 - m$m1^2
      },
      fit_check = xbeta_fit_check
    ),
    sim = function(dpars, aterms, n) {
      mu <- rep(dpars[["mu"]], length.out = n)
      phi <- rep(dpars[["phi"]], length.out = n)
      kap <- rep(dpars[["kappa"]], length.out = n)
      z <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
      pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
    }
  )
}

#' The extended-support beta's fit-end check: is `kappa` placed?
#'
#' Without a row at exactly 0 or 1 only the interior shape places
#' `kappa`, and it can go either way. It can run to zero, where the
#' model is `Beta()` (dev/fams2-xbeta-noends.R: log kappa -22.9 with a
#' standard error of 1.1e4, and `Beta()`'s log-likelihood to the digit),
#' or it can grow with `phi` along a ridge and beat `Beta()`: the
#' reviewer's A1 data, drawn with kappa 1 and phi 200, reach their
#' optimum near kappa 47 and phi 2e5, 6.04 log-likelihood units above
#' `Beta()` (dev/fams2-rev-a1.R). So the check reads where the fit
#' landed, not only the data. With no row at 0 or 1 it warns when
#' `kappa` is below 1e-6 at every row, where the ends carry no mass, or
#' when a coefficient of a log-link `kappa` has a standard error above
#' 10 (a 95% interval wider than a factor of e^39) or none at all.
#'
#' The fitted value alone missed two cases (the re-check's n3,
#' dev/fams2-rev2-guards.txt): a `kappa` that stopped at 3.4e-5, still
#' on its way to 0, with a standard error of 96, and `kappa ~ x`, whose
#' largest row reached 0.024 while the smallest reached 1.3e-10, with
#' standard errors of 22 and 9.4. The standard error catches both. The
#' eight fits there whose `kappa` the interior shape really placed had
#' standard errors of 1.3 to 6.9. A rule on the SMALLEST row kappa was
#' tried and dropped: `kappa ~ x` with a real slope (log kappa = -1 + 3x,
#' x in (-4.7, 0)) runs some rows to 2e-10 while its standard errors are
#' 0.26 and 1.38, 28 log-likelihood units above `Beta()`
#' (dev/fams2-rev3-kappa.txt, seeds 402 and 403).
#' The standard errors cost a `sdreport()`, which the fit keeps for
#' `summary()`, and are read only when no row is at 0 or 1. brms fits
#' such data through its `gamma(0.01, 0.01)` prior on `kappa`. One end
#' in the data is enough to place `kappa`, and a `kappa` held at a
#' constant has nothing to place.
#'
#' @noRd
xbeta_fit_check <- function(fit, resp) {
  lp <- fit$frame[["linpreds"]][[linpred_key(resp, "kappa")]]
  if (is.null(lp) || !is.null(lp[["constant"]])) return(invisible(NULL))
  y <- fit$frame[["y"]][[resp]]
  if (any(y <= 0 | y >= 1)) return(invisible(NULL))
  kap <- eval_dpars(fit)[[resp]][["kappa"]]
  why <- if (max(kap) < 1e-6) {
    paste0("kappa ran to 0 (at most ", signif(max(kap), 2), "), where ",
           "the model is Beta()")
  } else {
    fam <- fit$spec$responses[[resp]]$family
    se <- if (identical(fam[["links"]][["kappa"]][["name"]], "log")) {
      xbeta_kappa_se(fit, lp)
    }
    bad <- which(!is.finite(se) | se > 10)
    if (length(bad)) {
      j <- bad[which.max(replace(se[bad], !is.finite(se[bad]), Inf))]
      paste0("the standard error of ", names(se)[j], " is ",
             if (is.finite(se[j])) signif(se[j], 2) else "not finite",
             " on the log scale")
    }
  }
  if (is.null(why)) return(invisible(NULL))
  frm_warning("xbeta: the response ",
              if (length(fit$spec$responses) > 1L) paste0(resp, " "),
              "has no exact 0 or 1, so only the interior shape places ",
              "kappa, and this fit does not place it: ", why, ". Its ",
              "standard error is not usable. Compare with Beta(), or ",
              "hold kappa at a value with bf(kappa = )", call. = FALSE)
  invisible(NULL)
}

#' The standard errors of `kappa`'s coefficients, on the link scale,
#' named as `confint()` names them; NULL when the fit has no
#' `cov.fixed` block for them. A negative variance comes back `NaN`.
#'
#' @noRd
xbeta_kappa_se <- function(fit, lp) {
  cp <- lp[["par"]]
  V <- tryCatch(sdr_of(fit)$cov.fixed, error = function(e) NULL)
  if (is.null(V) || !cp %in% rownames(V)) return(NULL)
  pos <- seq_along(fit$frame[["par_template"]][[cp]])
  if (cp == "betad") pos <- setdiff(pos, fit$frame[["betad_fixed_idx"]])
  rows <- which(rownames(V) == cp)[match(lp[["idx"]], pos)]
  rows <- rows[!is.na(rows)]
  if (!length(rows)) return(NULL)
  v <- diag(V)[rows]
  se <- sqrt(abs(v))
  se[!is.finite(v) | v < 0] <- NaN
  onm <- outer_par_names(fit)
  names(se) <- if (length(onm) == nrow(V)) {
    onm[rows]
  } else {
    paste0("kappa[", seq_along(rows), "]")
  }
  se
}

#' `log I_x(a, b)`, the log regularized incomplete beta, for `x` below
#' one half, on the tape or off it.
#'
#' `RTMB::pbeta()` is exact in its value and first two derivatives, and
#' its THIRD derivatives come back `NaN` at ordinary points: 60 of 3840
#' random points an xbeta fit visits, most where one shape is small and
#' `I` is near one (dev/fams2-pbeta-nan-rate.R). The Laplace
#' approximation needs third derivatives, so an xbeta fit with a random
#' effect and rows at 0 or 1 stopped at "NA/NaN gradient evaluation" on
#' its first try.
#'
#' Two methods, blended so that the value and its first two derivatives
#' join everywhere:
#'
#' - The continued fraction of Numerical Recipes' `betacf()`, to a fixed
#'   101 terms and evaluated from the bottom up, which is arithmetic and
#'   differentiates to any order. It converges for `x` below
#'   `m = (a + 1) / (a + b + 2)`, and its complement `1 - I_{1-x}(b, a)`
#'   above. The two are blended over `[m + sd / 2, m + 3 sd / 2]` (`sd`
#'   the beta's standard deviation), where both hold, each read at `x`
#'   clamped to that band, so the one weighted zero cannot put a `NaN`
#'   through its weight. Top down (Lentz's order) the first partial
#'   denominator vanishes at `x = (a + 1) / (a + b)`, inside the band:
#'   evaluated that way it gave `NaN` at `x = 0.041`, `a = 40`,
#'   `b = 960`. The fraction needs about `sqrt(max(a, b))` steps within
#'   a few `sd` of `m` (the reviewer's Rmpfr run: over 1000 at shape
#'   1e6), so on its own it was wrong there at large shapes: 2.5e-2 in
#'   the log value at a shape sum of 1e5, and a jump in the gradient
#'   across `x = m` (dev/reviews/2026-09-29-fams2.md, B2).
#' - `RTMB::pbeta()` within 6 `sd` of `m` once `a b / (a + b)` passes 150,
#'   the region where the fraction is short. There the probability is
#'   moderate, so `log()` of it neither underflows nor loses digits,
#'   and its third derivatives were finite at all of 20000 random points
#'   (dev/fams2-p1-pthird.R). They are all `NaN` at the mean itself,
#'   which `log_pbeta_ad()` steps around. Where it carries no weight it
#'   is read at shapes raised to 150 and `x` held within 6.5 `sd` of that
#'   mean, a point of the same region, for the same reason.
#'
#' The blends are quintic smoothsteps, in `|x - m| / sd` between 5 and 6
#' and in `log(a b / (a + b))` between `log(150)` and `log(450)`, where
#' both methods hold. Measured against `stats::pbeta(log.p = TRUE)`,
#' relative and floored at one (the reviewer found it within 1e-13 of
#' Rmpfr at the worst points): on the reviewer's grid of 4867 points,
#' shapes 1e-3 to 1e7 in half decades with `x` at `m`, near it and at
#' every blend edge (dev/fams2-rev2-ibeta.R, output
#' dev/fams2-p2-rev2-ibeta.txt), the worst is 3.1e-14 while both shapes
#' are below 1e4, and by the larger shape 4.9e-14 to 1e5, 1.7e-13 to 1e6
#' and 3.0e-13 to 1e7. Against Rmpfr the bound at 1e7 is 7.2e-13, at
#' (1e7, 1e7) just below the tie of `log_pbeta_ad()`, where
#' `stats::pbeta()` itself is 2.9e-13 out (the final check's
#' dev/fams2-rev3-clamp.txt). At 3167 random points (dev/fams2-p1-sweep.R,
#' output dev/fams2-p2-sweep-pkg.txt) the worst is 5.5e-13 and 3.3e-13
#' at `x == m` exactly; the gradient is within 6.7e-14 of
#' `RTMB::pbeta()`'s exact one where the two are compared. First to
#' third derivatives are finite at every point of both. The value moves
#' across each blend edge as its gradient predicts, to 1.6e-13
#' (dev/fams2-p1-joins.R, output dev/fams2-p2-joins.txt). A gradient
#' sweep costs 14.5 us a row, 1.04 times round 1's fraction and 9.1 times
#' `log(RTMB::pbeta())` (dev/fams2-rev2-mean.R, section 4, output
#' dev/fams2-p2-rev2-mean.txt).
#'
#' @noRd
log_ibeta_half <- function(x, a, b, N = 50L) {
  s <- a + b
  m <- (a + 1) / (s + 2)
  sd <- sqrt(a * b / (s * s * (s + 1)))
  # The fraction and its complement, blended over [m + sd / 2,
  # m + 3 sd / 2]. The band sits above m because there the fraction still
  # holds (it reads 1e-11 at worst at m + 1.5 sd) and below m + sd / 2 the
  # complement does not: it lost 2.7e-4 of the log value at m - 2 sd with
  # a = 0.3, b = 7, and 1.9e-10 at m - sd with a = 60, b = 8e5
  # (dev/fams2-p1-sides.R). The complement's prefactor takes log(1 - y)
  # and log(y) of y = 1 - xc from xc itself, which the rounded 1 - xc
  # cost 1e-11 of at a = 0.3, b = 8e5, and its odd terms take xc too
  # (log_ibeta_cf()).
  sh <- ibeta_shape(a, b)
  wd <- ad_smooth01((m + 1.5 * sd - x) / sd)
  xd <- ad_min(x, m + 1.5 * sd)
  ld <- log_ibeta_cf(xd, a, b, N, log_ibeta_pre(xd, sh))
  xc <- ad_max(x, m + 0.5 * sd)
  # the prefactor of I_{1 - xc}(b, a) is the one of I_xc(a, b)
  lc <- log_ibeta_cf(1 - xc, b, a, N, log_ibeta_pre(xc, sh), xc = xc)
  # I_{1-x}(b, a) is at most one half on its own side of the band; past
  # it the weight is zero and the cap only keeps log1p() finite. The cap
  # side has to come out exactly: cap - pos(cap - lc) is cap there, where
  # the other spelling rounded to 0 and gave log(0)
  cap <- log1p(-2^-53)
  lc <- cap - ad_pos(cap - lc)
  lcf <- wd * ld + (1 - wd) * log1p(-exp(lc))
  # within 6 sd of the mean at large shapes, RTMB::pbeta()
  k <- abs(x - m) / sd
  # a b / (a + b) lies between min(a, b) / 2 and min(a, b) and, unlike
  # min(), is smooth where a = b
  ws <- ad_smooth01((log(a * b / s) - log(150)) / (log(450) - log(150)))
  wp <- ws * ad_smooth01(6 - k)
  a2 <- ad_max(a, 150)
  b2 <- ad_max(b, 150)
  s2 <- a2 + b2
  m2 <- a2 / s2
  sd2 <- sqrt(a2 * b2 / (s2 * s2 * (s2 + 1)))
  x2 <- ad_max(ad_min(x, m2 + 6.5 * sd2), m2 - 6.5 * sd2)
  (1 - wp) * lcf + wp * log_pbeta_ad(x2, a2, b2)
}

#' `log(RTMB::pbeta(x, a, b))` with finite derivatives at the mean.
#'
#' `RTMB::pbeta()`'s gradient, Hessian and third derivatives are all
#' `NaN` at `x == a / (a + b)` to within an ulp, and finite from a
#' relative distance of 1e-15 (dev/fams2-p2-tie.R; the re-check's n1).
#' Below the mean this uses `I_x(a, b) = I_x(a + 1, b) +
#' x^a (1 - x)^b / (a B(a, b))`, whose `pbeta()` has ITS tie at
#' `(a + 1) / (a + b + 1)`, above the mean by `b / ((a + b) (a + b + 1))`.
#' The two forms are blended over the middle half of the gap between the
#' two ties, each read at `x` clamped out of its own tie, so the one
#' weighted zero never reaches a `NaN`. Both are exact, so the blend
#' costs nothing in the value.
#'
#' It does cost the higher derivatives inside the blend: the weight's
#' k-th derivative scales as `gap^-k`, with the gap about `1 / (a + b)`,
#' and multiplies the rounding difference of the two forms. The third
#' derivatives there are off by up to 0.12 relative at (1e6, 3e6) and
#' 7.7e-5 at (3e4, 7e4), the second by up to 1.3e-6, where plain
#' `log(RTMB::pbeta())` is good to 3e-9 (the final check's
#' dev/fams2-rev3-near.txt). The blend cannot leave the gap between the
#' two ties, so widening it within the gap gains at most a factor of 4.
#' The band is about `1 / sqrt(a)` sd wide, few rows land in it, and
#' values and first derivatives are unaffected.
#'
#' @noRd
log_pbeta_ad <- function(x, a, b) {
  s <- a + b
  t1 <- a / s
  gap <- b / (s * (s + 1))
  lo <- t1 + 0.25 * gap
  w <- ad_smooth01((x - lo) / (0.5 * gap))
  xa <- ad_max(x, lo)
  xb <- ad_min(x, lo + 0.5 * gap)
  la <- log(RTMB::pbeta(xa, a, b))
  lb <- log(RTMB::pbeta(xb, a + 1, b) +
              exp(log_ibeta_pre(xb, ibeta_shape(a, b, large = TRUE)) -
                    log(a)))
  w * la + (1 - w) * lb
}

#' Tape-safe `max(u, 0)`, `min(x, m)` and `max(x, m)`. The last two
#' return `x` bit for bit on its own side, because `(x - m) + (m - x)` is
#' exactly zero; the symmetric `0.5 (x + m - |m - x|)` rounds `x + m`
#' and cost `x` eleven digits at `x = 6.5e-6` against `m = 0.69`.
#'
#' @noRd
ad_pos <- function(u) 0.5 * (u + abs(u))

#' Tape-safe `min(x, m)`, which is `x` bit for bit when `x <= m`.
#'
#' @noRd
ad_min <- function(x, m) x - ad_pos(x - m)

#' Tape-safe `max(x, m)`, which is `x` bit for bit when `x >= m`.
#'
#' @noRd
ad_max <- function(x, m) x + ad_pos(m - x)

#' 0 below 0, 1 above 1, and the quintic smoothstep `10 t^3 - 15 t^4 +
#' 6 t^5` between, whose first and second derivatives vanish at both
#' ends: a blend weight whose use keeps a value and its first two
#' derivatives continuous.
#'
#' @noRd
ad_smooth01 <- function(t) {
  t <- ad_max(ad_min(t, 1), 0)
  t * t * t * (10 - 15 * t + 6 * t * t)
}

#' `log B(a, b)` for the tape, formed as `lgamma(a) - (lgamma(a + b) -
#' lgamma(b))` with the difference from `lgamma_shift_diff()` based at
#' the LARGER shape, so a small shape beside a large one does not cancel
#' two large `lgamma()` values. `RTMB::lbeta()` did, and cost the
#' gradient of `log_ibeta_half()` seven digits at `a = 0.01`, `b = 5e6`
#' (dev/fams2-p1-sweep.R).
#'
#' The two orderings are blended over `|b - a| < 0.1 (a + b)`, where both
#' are as good, rather than switched through `min()` and `max()`: RTMB's
#' `abs()` has derivative one at zero, so at `a == b` exactly both of
#' those followed `b` alone and the derivative in `a` came out 0.
#'
#' @noRd
lbeta_ad <- function(a, b) {
  w <- ad_smooth01((b - a) / (0.2 * (a + b)) + 0.5)
  w * (lgamma(a) - lgamma_shift_diff(b, a)) +
    (1 - w) * (lgamma(b) - lgamma_shift_diff(a, b))
}

#' `log I_x(a, b)` from the continued fraction alone, valid for
#' `x < (a + 1) / (a + b + 2)`. See `log_ibeta_half()`. `lpre` is
#' `log_ibeta_pre(x, ibeta_shape(a, b))`. `xc`, when given, is `1 - x` held more
#' exactly than `x` holds it, and the odd terms are then formed from it.
#'
#' Each odd level is `1 + e_odd(k) / t`, where `e_odd(k)` is near -1
#' once `a` is large, so `1 + e_odd(k)` cancels. From `x` that loses
#' about `a` ulps; from `xc` it is
#' `((a + k) (2k + 1 - b + (a + b + k) xc) + k (k + 1)) /
#' ((a + 2k) (a + 2k + 1))`, whose large terms are gone. It matters for
#' the complement at `1 - xc` with a large first shape and a small
#' `xc`: 1.1e-11 of the log value at `a = 1`, `b = 1e7` one `sd` above
#' the mean, formed from `x` (dev/fams2-p2-rev2-ibeta.txt).
#'
#' @noRd
log_ibeta_cf <- function(x, a, b, N, lpre, xc = NULL) {
  qab <- a + b
  qap <- a + 1
  qam <- a - 1
  # Numerical Recipes' betacf() fraction, 1 / (1 + e1 / (1 + e2 / ...)),
  # to 2 N + 1 terms, evaluated from the bottom up. Lentz's top-down
  # evaluation of the same approximant forms 1 / (1 + e1) on the way,
  # which is infinite at x = (a + 1) / (a + b); bottom up only the
  # approximant's own poles are left. Each odd level 1 + e_odd / t is
  # formed as (p_odd + r) / t with p_odd = 1 + e_odd and r = t - 1,
  # both small where they cancel
  p_odd <- if (is.null(xc)) {
    function(k) 1 - (a + k) * (qab + k) * x / ((a + 2 * k) * (qap + 2 * k))
  } else {
    function(k) {
      ((a + k) * (2 * k + 1 - b + (qab + k) * xc) + k * (k + 1)) /
        ((a + 2 * k) * (qap + 2 * k))
    }
  }
  e_even <- function(k) k * (b - k) * x / ((qam + 2 * k) * (a + 2 * k))
  t <- p_odd(N)
  for (k in rev(seq_len(N))) {
    r <- e_even(k) / t
    t <- (p_odd(k - 1) + r) / (1 + r)
  }
  lpre - log(a) - log(abs(t))
}

#' `log(x^a (1 - x)^b / B(a, b))`, the prefactor of the continued
#' fraction, which is symmetric in `(x, a)` and `(1 - x, b)`. `sh` is
#' `ibeta_shape(a, b)`. `lx` and `l1mx` are `log(x)` and `log(1 - x)`,
#' for a caller that holds `1 - x` more exactly than `x`.
#'
#' Written as `a log x + b log(1 - x) - log B(a, b)` its three terms are
#' each about `a + b` in size and cancel down to a few units near the
#' mean, so at shapes of 1e7 the rounding left 6.5e-10 of the log
#' incomplete beta (dev/reviews/2026-09-29-fams2.md, re-check n2). Once
#' both shapes pass 12, where `lgamma_binet()` is exact, it is formed
#' around the mean `mu = a / (a + b)` instead:
#' `a log(x / mu) + b log((1 - x) / (1 - mu)) + C(a, b)`, with
#' `C = log(a b / (2 pi (a + b))) / 2 + mu(a + b) - mu(a) - mu(b)` from
#' Stirling's formula (`mu()` Binet's remainder), exact, and the two
#' logs from `log1p()` of the relative distance to the mean. An error in
#' the rounded `mu` drops out to first order, because the sum is
#' stationary in `mu` there. The two forms are blended in
#' `a b / (a + b)` over 15 to 30, where both hold; the old one is
#' accurate while a shape is small, because `b log(1 - x)` is then about
#' `a`. `log1p()` of a relative distance near -1 loses what `x` held, so
#' from a relative distance of 1/4 to 1/2 each log is blended into a
#' difference of logs, which is by then far from cancelling.
#'
#' @noRd
log_ibeta_pre <- function(x, sh, lx = log(x), l1mx = log1p(-x)) {
  a <- sh$a
  b <- sh$b
  d <- x - sh$mua
  u <- d / sh$mua
  v <- -d / sh$mub
  wu <- ad_smooth01(4 * abs(u) - 1)
  wv <- ad_smooth01(4 * abs(v) - 1)
  # the argument is held off -1 where its weight is zero, so log1p()
  # stays finite there
  lu <- (1 - wu) * log1p(ad_max(u, -0.75)) + wu * (lx - log(sh$mua))
  lv <- (1 - wv) * log1p(ad_max(v, -0.75)) + wv * (l1mx - log(sh$mub))
  near <- a * lu + b * lv + sh$cst
  if (is.null(sh$lb)) return(near)
  (1 - sh$w) * (a * lx + b * l1mx - sh$lb) + sh$w * near
}

#' The parts of `log_ibeta_pre()` that depend on the shapes alone,
#' formed once for the several points one `log_ibeta_half()` reads.
#' `large = TRUE` promises `a b / (a + b)` of 30 or more, where the
#' form around the mean has all the weight, and skips the other.
#'
#' @noRd
ibeta_shape <- function(a, b, large = FALSE) {
  s <- a + b
  a2 <- ad_max(a, 12)
  b2 <- ad_max(b, 12)
  s2 <- a2 + b2
  list(a = a, b = b, mua = a / s, mub = b / s,
       cst = 0.5 * log(a2 * b2 / (2 * pi * s2)) + lgamma_binet(s2) -
         lgamma_binet(a2) - lgamma_binet(b2),
       lb = if (!large) lbeta_ad(a, b),
       w = if (!large) {
         ad_smooth01((log(a * b / s) - log(15)) / (log(30) - log(15)))
       })
}

#' The first two moments of an extended-support beta response, off the
#' tape.
#'
#' `m1` is brms's `posterior_epred_xbeta()` term for term. With
#' `Y = min(max(d Z - kappa, 0), 1)`, `d = 1 + 2 kappa`,
#' `E[Y] = P(Z > q1) + d E[Z; q0 < Z < q1] - kappa P(q0 < Z < q1)`, and
#' `E[Z^k; .]` is a beta moment times the incomplete beta with the first
#' shape raised by `k`.
#'
#' @noRd
xbeta_moments <- function(mu, phi, kappa) {
  a <- mu * phi
  b <- (1 - mu) * phi
  d <- 1 + 2 * kappa
  q0 <- kappa / d
  q1 <- (1 + kappa) / d
  # P(q0 < Z < q1) and the truncated first and second moments of Z
  p <- stats::pbeta(q1, a, b) - stats::pbeta(q0, a, b)
  e1 <- mu * (stats::pbeta(q1, a + 1, b) - stats::pbeta(q0, a + 1, b))
  e2 <- mu * (a + 1) / (phi + 1) *
    (stats::pbeta(q1, a + 2, b) - stats::pbeta(q0, a + 2, b))
  # P(Z > q1) = I_q0(b, a), without the subtraction from one
  top <- stats::pbeta(q0, b, a)
  list(m1 = top + d * e1 - kappa * p,
       m2 = top + d^2 * e2 - 2 * d * kappa * e1 + kappa^2 * p)
}

#' Asymmetric Laplace family, dpars `mu`, `sigma` and `quantile`. At a
#' fixed `quantile` the maximum-likelihood fit gives the quantile
#' regression point estimates.
#'
#' @noRd
fam_asym_laplace <- function(link = "identity", link_sigma = "log",
                             link_quantile = "logit") {
  lk_sigma <- dpar_link(
    link_sigma, "sigma", "asym_laplace", dpar_links_positive)
  lk_quantile <- dpar_link(
    link_quantile, "quantile", "asym_laplace", dpar_links_unit)
  frmtmb_family(
    "asym_laplace",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma", "quantile"),
    links = list(mu = mu_link(link, "asym_laplace"),
                 sigma = lk_sigma, quantile = lk_quantile),
    lpdf = function(y, dpars, aterms) {
      # rho_p(u) = 0.5 * (|u| + (2p - 1) u); ML at fixed quantile p
      # reproduces quantile-regression point estimates
      p <- dpars[["quantile"]]
      u <- (y - dpars[["mu"]]) / dpars[["sigma"]]
      # log(p) + log(1 - p) is the normalizer; an extreme quantile makes
      # one of them -Inf through the logit round trip
      gq <- dpar_log_complement(dpars, "quantile", lk_quantile)
      gq$l + gq$l1m - log(dpars[["sigma"]]) -
        0.5 * (abs(u) + (2 * p - 1) * u)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y),
      sigma = function(y, aterms) stats::sd(y) / 2,
      quantile = function(y, aterms) 0.5
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        dpars[["mu"]] + dpars[["sigma"]] * (1 - 2 * dpars[["quantile"]]) /
          (dpars[["quantile"]] * (1 - dpars[["quantile"]]))
      }
    ),
    sim = function(dpars, aterms, n) {
      p <- dpars[["quantile"]]
      dpars[["mu"]] + dpars[["sigma"]] *
        (stats::rexp(n) / p - stats::rexp(n) / (1 - p))
    }
  )
}

#' Zero-inflated asymmetric Laplace (brms spelling): a point mass at
#' zero mixed with the continuous ALD, so the mixture density is
#' well-defined without the +i0 dodge the discrete zi families need.
#'
#' @noRd
fam_zi_asym_laplace <- function(link = "identity", link_sigma = "log",
                                link_quantile = "logit", link_zi = "logit") {
  lk_sigma <- dpar_link(
    link_sigma, "sigma", "zero_inflated_asym_laplace", dpar_links_positive)
  lk_quantile <- dpar_link(
    link_quantile, "quantile", "zero_inflated_asym_laplace", dpar_links_unit)
  lk_zi <- dpar_link(
    link_zi, "zi", "zero_inflated_asym_laplace", dpar_links_unit)
  frmtmb_family(
    "zero_inflated_asym_laplace",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma", "quantile", "zi"),
    links = list(mu = mu_link(link, "zero_inflated_asym_laplace"),
                 sigma = lk_sigma, quantile = lk_quantile,
                 zi = lk_zi),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      p <- dpars[["quantile"]]
      u <- (y - dpars[["mu"]]) / dpars[["sigma"]]
      gq <- dpar_log_complement(dpars, "quantile", lk_quantile)
      ald <- gq$l + gq$l1m - log(dpars[["sigma"]]) -
        0.5 * (abs(u) + (2 * p - 1) * u)
      g <- dpar_log_complement(dpars, "zi", lk_zi)
      i0 * g$l + (1 - i0) * (g$l1m + ald)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y[y != 0]),
      sigma = function(y, aterms) {
        s <- stats::sd(y[y != 0]) / 2
        if (is.finite(s) && s > 0) s else 1
      },
      quantile = function(y, aterms) 0.5,
      zi = function(y, aterms) min(max(mean(y == 0), 0.05), 0.9)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars[["zi"]]) *
          (dpars[["mu"]] + dpars[["sigma"]] * (1 - 2 * dpars[["quantile"]]) /
             (dpars[["quantile"]] * (1 - dpars[["quantile"]])))
      }
    ),
    sim = function(dpars, aterms, n) {
      p <- dpars[["quantile"]]
      stats::rbinom(n, 1, 1 - dpars[["zi"]]) *
        (dpars[["mu"]] + dpars[["sigma"]] *
           (stats::rexp(n) / p - stats::rexp(n) / (1 - p)))
    }
  )
}

#' Huber's rho at tuning constant `k`, branch-free:
#' `rho(u) = (u^2 - max(|u| - k, 0)^2) / 2`. Inside the kink the second
#' term is zero and this is the gaussian `u^2 / 2`; outside it the
#' quadratics cancel and leave `k|u| - k^2/2`, the Laplace tail. `max(x,
#' 0)` is written `(x + |x|) / 2`, so nothing branches on a parameter
#' value - the switch between the two regimes has to happen ON the tape,
#' once per observation per gradient evaluation, and an `if` there would
#' bake one regime into the tape for good.
#'
#' @noRd
huber_rho <- function(u, k) {
  a <- abs(u)
  e <- a - k
  0.5 * (a * a - 0.25 * (e + abs(e))^2)
}

#' The normalizing constant of `exp(-rho_k(u))`:
#' `sqrt(2 pi) (2 Phi(k) - 1) + (2 / k) exp(-k^2 / 2)`, the gaussian
#' piece over `[-k, k]` plus the two exponential tails. Written through
#' `dnorm(k)` because `exp(-k^2 / 2) = sqrt(2 pi) phi(k)`; the tests
#' check it against a numeric integral at 1e-10.
#'
#' @noRd
huber_norm <- function(k) {
  sqrt(2 * pi) * (2 * stats::pnorm(k) - 1 + 2 * stats::dnorm(k) / k)
}

#' `E[U^2]` under the standard (sigma = 1) Huber density, for the
#' pearson-residual variance. Both integrals are closed form: the
#' truncated gaussian second moment plus the two exponential tails'.
#'
#' @noRd
huber_var_u <- function(k) {
  ph <- stats::dnorm(k)
  c0 <- 2 * stats::pnorm(k) - 1
  (c0 + 2 * ph * (2 / k + 2 / k^3)) / (c0 + 2 * ph / k)
}

#' Draws from the standard Huber density by inverting its CDF. This runs
#' off the tape, so the three pieces (exponential lower tail, gaussian
#' middle, exponential upper tail) branch freely.
#'
#' @noRd
rhuber_u <- function(n, k) {
  z <- huber_norm(k)
  p0 <- exp(-k^2 / 2) / (k * z)            # F(-k), and 1 - F(k)
  p <- stats::runif(n)
  u <- numeric(n)
  lo <- p < p0
  hi <- p > 1 - p0
  mid <- !lo & !hi
  u[lo] <- (log(p[lo] * k * z) - k^2 / 2) / k
  u[hi] <- (k^2 / 2 - log((1 - p[hi]) * k * z)) / k
  u[mid] <- stats::qnorm(stats::pnorm(-k) +
                           (p[mid] - p0) * z / sqrt(2 * pi))
  u
}

#' Huber's least-favorable density, dpars `mu` and `sigma`. `k` is a
#' fixed tuning constant of the FAMILY, not a dpar: it says where the
#' analyst draws the line between "residual" and "outlier", which is a
#' choice, not a quantity the data identifies. `MASS::rlm()` treats it
#' the same way. Estimating it would let the likelihood buy fit by
#' widening the gaussian core, which is the opposite of what the family
#' is for.
#'
#' @noRd
fam_huber <- function(link = "identity", k = 1.345, link_sigma = "log") {
  lk_sigma <- dpar_link(link_sigma, "sigma", "huber", dpar_links_positive)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    frm_stop("huber(k =): the tuning constant must be one finite positive ",
             "number; 1.345 (the default, and MASS::rlm()'s) gives 95% ",
             "efficiency against a gaussian", call. = FALSE)
  }
  lognorm <- log(huber_norm(k))
  varu <- huber_var_u(k)
  frmtmb_family(
    "huber",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma"),
    links = list(mu = mu_link(link, "huber"), sigma = lk_sigma),
    lpdf = function(y, dpars, aterms) {
      s <- resid_sd(dpars[["sigma"]], aterms)
      -log(s) - lognorm - huber_rho((y - dpars[["mu"]]) / s, k)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y),
      sigma = function(y, aterms) {
        s <- stats::mad(y)
        if (is.finite(s) && s > 0) s else max(stats::sd(y), 1e-3)
      }
    ),
    type = "continuous",
    post = list(
      # symmetric about mu, so the mean is mu
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) {
        resid_sd(dpars[["sigma"]], aterms)^2 * varu
      },
      # 2 (ll at mu = y - ll at mu), the scale held fixed: 2 sigma^2
      # rho_k(u), which collapses to the gaussian (y - mu)^2 as k grows
      dev_fn = function(y, dpars, aterms) {
        s <- resid_sd(dpars[["sigma"]], aterms)
        2 * s^2 * huber_rho((y - dpars[["mu"]]) / s, k)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars[["mu"]] + resid_sd(dpars[["sigma"]], aterms) * rhuber_u(n, k)
    }
  )
}

# --- RTMBdist-backed families ---

#' Beta-binomial family in the mean-precision parameterization, dpars
#' `mu` and `phi`. Trials come from the `trials()` addition term.
#'
#' @noRd
fam_beta_binomial <- function(link = "logit", link_phi = "log") {
  lk_phi <- dpar_link(link_phi, "phi", "beta_binomial", dpar_links_positive)
  lk <- mu_link(link, "beta_binomial")
  frmtmb_family(
    "beta_binomial",
    accepts_aterms = c("weights", "trials"),
    dpars = c("mu", "phi"),
    links = list(mu = lk, phi = lk_phi),
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]] %||% 1
      # RTMBdist has no log-odds form, so the robustness has to go into
      # the shapes: a zero second shape gives lgamma(0) = NaN
      mp <- dpar_complement(dpars, "mu", lk)
      RTMBdist::dbetabinom(y, size, mp$p * dpars[["phi"]],
                mp$q * dpars[["phi"]],
                           log = TRUE)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]] %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        frm_stop("beta_binomial: response must be integer counts in ",
                 "[0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        p <- mean(y / (aterms[["trials"]] %||% 1))
        min(max(p, 0.02), 0.98)
      },
      phi = function(y, aterms) 5
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        dpars[["mu"]] * (aterms[["trials"]] %||% 1)
      },
      # without it pearson residuals were refused for want of a variance
      var_fn = function(dpars, aterms) {
        beta_binomial_var(dpars[["mu"]], dpars[["phi"]],
                          aterms[["trials"]] %||% 1)
      }
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rbetabinom(n, aterms[["trials"]] %||% 1,
                           dpars[["mu"]] * dpars[["phi"]],
                           (1 - dpars[["mu"]]) * dpars[["phi"]])
    }
  )
}

#' The variance of a beta-binomial count with mean `size * mu` and
#' precision `phi`, off the tape.
#'
#' @noRd
beta_binomial_var <- function(mu, phi, size) {
  size * mu * (1 - mu) * (phi + size) / (phi + 1)
}

#' Zero-inflated beta-binomial, dpars `mu`, `phi` and `zi`, as brms
#' 2.23.0 defines it: `P(Y = 0) = zi + (1 - zi) BB(0)` and
#' `P(Y = y) = (1 - zi) BB(y)` above zero, where `BB` is
#' `fam_beta_binomial()`'s density with trials from `trials()`.
#'
#' @noRd
fam_zi_beta_binomial <- function(link = "logit", link_phi = "log",
                                 link_zi = "logit") {
  nm <- "zero_inflated_beta_binomial"
  lk_phi <- dpar_link(link_phi, "phi", nm, dpar_links_positive)
  lk_zi <- dpar_link(link_zi, "zi", nm, dpar_links_unit)
  lk <- mu_link(link, nm)
  frmtmb_family(
    nm,
    accepts_aterms = c("weights", "trials"),
    dpars = c("mu", "phi", "zi"),
    links = list(mu = lk, phi = lk_phi, zi = lk_zi),
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]] %||% 1
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "zi", lk_zi)
      # the shapes off the log-odds, for the reason fam_beta_binomial()
      # gives; at y = 0 `base` is the beta-binomial's own P(0)
      mp <- dpar_complement(dpars, "mu", lk)
      base <- RTMBdist::dbetabinom(y, size, mp$p * dpars[["phi"]],
                                   mp$q * dpars[["phi"]], log = TRUE)
      i0 * RTMB::logspace_add(g$l, g$l1m + base) + (1 - i0) * (g$l1m + base)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]] %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        frm_stop(nm, ": response must be integer counts in [0, trials]",
                 call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        p <- mean(y / (aterms[["trials"]] %||% 1))
        min(max(p, 0.02), 0.98)
      },
      phi = function(y, aterms) 5,
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      # brms's posterior_epred_zero_inflated_beta_binomial()
      mean_fn = function(dpars, aterms) {
        (1 - dpars[["zi"]]) * dpars[["mu"]] * (aterms[["trials"]] %||% 1)
      },
      var_fn = function(dpars, aterms) {
        size <- aterms[["trials"]] %||% 1
        mu <- dpars[["mu"]]
        q <- 1 - dpars[["zi"]]
        # E[Y^2] = (1 - zi) (Var BB + E[BB]^2)
        q * (beta_binomial_var(mu, dpars[["phi"]], size) + (size * mu)^2) -
          (q * size * mu)^2
      }
    ),
    sim = function(dpars, aterms, n) {
      size <- aterms[["trials"]] %||% 1
      (1 - stats::rbinom(n, 1L, dpars[["zi"]])) *
        RTMBdist::rbetabinom(n, size, dpars[["mu"]] * dpars[["phi"]],
                             (1 - dpars[["mu"]]) * dpars[["phi"]])
    }
  )
}

#' Skew-normal family in the mean parameterization, dpars `mu` (the
#' mean), `sigma` and `alpha` (the skewness). The `alpha` start avoids
#' zero, which is a stationary point of the likelihood, and the family
#' declares that point so a fit that lands on it is sent out from both
#' sides.
#'
#' @noRd
fam_skew_normal <- function(link = "identity", link_sigma = "log",
                            link_alpha = "identity") {
  lk_sigma <- dpar_link(link_sigma, "sigma", "skew_normal", dpar_links_positive)
  lk_alpha <- dpar_link(link_alpha, "alpha", "skew_normal", dpar_links_signed)
  frmtmb_family(
    "skew_normal",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma", "alpha"),
    links = list(mu = mu_link(link, "skew_normal"),
                 sigma = lk_sigma, alpha = lk_alpha),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dskewnorm2(y, dpars[["mu"]], dpars[["sigma"]], dpars[["alpha"]],
                           log = TRUE)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      # sigma here is the CONDITIONAL standard deviation, so the
      # marginal sd(y) overstates it by whatever the mu predictor
      # explains: 2.0x on the design in dev/skewinit-findings.md. With
      # no mu covariate `resid` IS the response, so this is bitwise the
      # number sd(y) gave before.
      sigma = function(y, aterms, resid) stats::sd(resid),
      # The third argument asks make_start() for the response with its
      # mu predictor taken out. alpha describes the RESIDUAL skew, and
      # a covariate can give the raw response the opposite sign: the
      # start then lands on the wrong side of alpha = 0 and the
      # optimizer stops there. Measured in dev/skewinit-findings.md.
      alpha = function(y, aterms, resid) {
        s <- stats::sd(resid)
        m3 <- if (is.finite(s) && s > 0) {
          mean((resid - mean(resid))^3) / s^3
        } else 0
        # a zero third moment would start AT the stationary point, which
        # is the one value from which no optimizer can leave
        if (!is.finite(m3) || m3 == 0) m3 <- 1
        2 * sign(m3) + 0.5 * m3
      }
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      # alpha = 0 is a stationary point of the skew-normal likelihood at
      # EVERY sample, and its information is singular there, so the
      # gradient and the curvature both vanish and no convergence test
      # can see it. An optimizer that stops there has to be sent out
      # from both sides to find out whether it is the maximum.
      stationary = list(alpha = list(at = 0, tol = 0.05, from = c(2, -2)))
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rskewnorm2(n, dpars[["mu"]], dpars[["sigma"]], dpars[["alpha"]])
    }
  )
}

#' Inverse gaussian family, dpars `mu` (the mean) and `shape`. It
#' supplies a CDF and a truncated mean.
#'
#' @noRd
fam_inverse_gaussian <- function(link = "1/mu^2", link_shape = "log") {
  lk_shape <- dpar_link(
    link_shape, "shape", "inverse.gaussian", dpar_links_positive)
  frmtmb_family(
    "inverse.gaussian",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = c("mu", "shape"),
    links = list(mu = mu_link(link, "inverse.gaussian"), shape = lk_shape),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dinvgauss(y, mean = dpars[["mu"]], shape = dpars[["shape"]],
                          log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMBdist::pinvgauss(q, mean = dpars[["mu"]], shape = dpars[["shape"]])
    },
    valid_y = positive_y("inverse.gaussian"),
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      shape = function(y, aterms) mean(y)^3 / max(stats::var(y), 1e-8)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) dpars[["mu"]]^3 / dpars[["shape"]],
      # 1 / shape is the dispersion and divides out
      dev_fn = function(y, dpars, aterms) {
        (y - dpars[["mu"]])^2 / (y * dpars[["mu"]]^2)
      },
      # Partial moment (Jorgensen): E[Y 1{Y <= x}] =
      # mu (Phi(z1) - e^{2 lambda / mu} Phi(-z2)), with the CDF the same
      # pair added instead of subtracted; the exponential factor
      # overflows on its own, so it rides in log space
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        mu <- dpars[["mu"]]
        lam <- dpars[["shape"]]
        part <- function(x) {
          r <- sqrt(lam / x)
          mu * (stats::pnorm(r * (x / mu - 1)) -
                  exp(2 * lam / mu +
                        stats::pnorm(-r * (x / mu + 1), log.p = TRUE)))
        }
        cdf <- function(x) {
          r <- sqrt(lam / x)
          stats::pnorm(r * (x / mu - 1)) +
            exp(2 * lam / mu +
                  stats::pnorm(-r * (x / mu + 1), log.p = TRUE))
        }
        lo <- pmax(lb, 0)
        pl <- ifelse(lo > 0, part(lo), 0)
        fl <- ifelse(lo > 0, cdf(lo), 0)
        pu <- ifelse(is.finite(ub), part(ub), mu)
        fu <- ifelse(is.finite(ub), cdf(ub), 1)
        (pu - pl) / (fu - fl)
      }
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rinvgauss(n, mean = dpars[["mu"]], shape = dpars[["shape"]])
    }
  )
}

#' Ex-gaussian family, dpars `mu`, `sigma` and `beta`.
#' brms parameterization: `mu` is the DISTRIBUTION mean, `beta` the
#' scale of the exponential component (gaussian component sits at
#' `mu - beta`).
#'
#' @noRd
fam_exgaussian <- function(link = "identity", link_sigma = "log",
                           link_beta = "log") {
  lk_sigma <- dpar_link(link_sigma, "sigma", "exgaussian", dpar_links_positive)
  lk_beta <- dpar_link(link_beta, "beta", "exgaussian", dpar_links_positive)
  frmtmb_family(
    "exgaussian",
    accepts_aterms = "weights",
    dpars = c("mu", "sigma", "beta"),
    links = list(mu = mu_link(link, "exgaussian"),
                 sigma = lk_sigma, beta = lk_beta),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dexgauss(y, dpars[["mu"]] - dpars[["beta"]], dpars[["sigma"]],
                         1 / dpars[["beta"]], log = TRUE)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y) / 2,
      beta = function(y, aterms) stats::sd(y) / 2
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) dpars[["sigma"]]^2 + dpars[["beta"]]^2
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rexgauss(n, dpars[["mu"]] - dpars[["beta"]], dpars[["sigma"]],
                         1 / dpars[["beta"]])
    }
  )
}

# Exact 1{y == k} for k = 1..K built from arithmetic alone. RTMB
# advectors carry no comparison operators, and oneStepPredict re-tapes
# the objective with the response promoted to a parameter, which is what
# breaks the indexing forms below. At integer y this Lagrange basis is
# exact in floating point - off the diagonal one factor is exactly 0, on
# it every factor is exactly 1 - so the taped path and the data path
# agree bit for bit.
#' During oneStepPredict RTMB hands the response to the lpdf as an "osa"
#' object: the taped value in `@x` and the per-row data-term indicator
#' in `@keep`. RTMB's own densities apply the indicator through
#' dGenericOSA, so a hand-written lpdf has to do it itself or every
#' observation stays switched on and the one-step sequence collapses.
#' This returns the pair, or `NULL` when the response is plain data.
#'
#' @noRd
osa_unwrap <- function(y) {
  if (!methods::is(y, "osa")) return(NULL)
  keep <- y@keep
  if (ncol(keep) != 1L) {
    frm_stop("osa_method = \"cdf\" is not supported for this family",
             call. = FALSE)
  }
  list(y = y@x, keep = keep[, 1])
}

#' The `1{y == k}` indicators for `k = 1..K`, one vector per category,
#' from the arithmetic Lagrange basis described above. It works with the
#' response on the AD tape, where comparison operators are not available.
#'
#' @noRd
ord_cat_sel <- function(y, K) {
  lapply(seq_len(K), function(k) {
    s <- 1
    for (j in seq_len(K)) {
      if (j != k) s <- s * ((y - j) / (k - j))
    }
    s
  })
}

#' `1{j < y}` for `j = 1..K-1`, from the same basis.
#'
#' @noRd
ord_cat_below <- function(sel, K) {
  lapply(seq_len(K - 1L), function(j) {
    b <- sel[[j + 1L]]
    if (j + 1L < K) for (k in (j + 2L):K) b <- b + sel[[k]]
    b
  })
}

#' Cumulative ordinal: response in `1..K` (or an ordered factor). The
#' linear predictor has no intercept; K-1 ordered thresholds take its
#' place, parameterized as (tau_1, log increments) in `extra_pars` under
#' `threshold = "flexible"`. The distribution function is read at
#' `disc * (tau_k - eta)`, with `disc` held at one unless the formula
#' models it, as in brms.
#'
#' @noRd
fam_cumulative <- function(link = "logit", link_disc = "log",
                           threshold = "flexible") {
  threshold <- ord_threshold_arg(threshold, "cumulative")
  # brms allows softit here and not for the sequential pair
  lk <- ord_link(link, "cumulative")
  lk_disc <- dpar_link(link_disc, "disc", "cumulative", dpar_links_positive)
  Fcdf <- lk$linkinv
  make_lpdf <- function(tmap) {
    force(tmap)
    function(y, dpars, aterms, extra) {
      tau <- tmap(extra$tau_raw)
      eta <- dpars[["mu"]]
      disc <- dpars[["disc"]] %||% 1
      K1 <- length(tau)
      K <- K1 + 1L
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        # OSA re-tape: pick the category probability arithmetically
        sel <- ord_cat_sel(ov$y, K)
        Fk <- lapply(seq_len(K1), function(k) Fcdf(disc * (tau[k] - eta)))
        dens <- sel[[1]] * Fk[[1]]
        if (K1 > 1) {
          for (k in 2:K1) dens <- dens + sel[[k]] * (Fk[[k]] - Fk[[k - 1]])
        }
        return(log(dens + sel[[K]] * (1 - Fk[[K1]])) * ov$keep)
      }
      ord_cumulative_logpmf(y, eta, tau, lk, disc)
    }
  }
  fam <- frmtmb_family(
    "cumulative",
    accepts_aterms = c("weights", "thres"),
    family_finalize = thres_finalizer("cumulative", ordered = TRUE,
                                      link = lk),
    dpars = c("mu", "disc"),
    links = list(mu = "identity", disc = lk_disc),
    lpdf = make_lpdf(ord_tau_from_raw_ad),
    valid_y = ord_valid_y("cumulative"),
    type = "ordinal",
    extra_pars = function(y, aterms) {
      ord_tau_init(y, ordered = TRUE, link = lk)
    },
    sim = ord_sim("cumulative", link = lk,
                  tmap = function(raw) ord_tau_from_raw(raw, TRUE)),
    post = list(ord_thresholds = ord_threshold_map(TRUE),
                ord_thresholds_raw = ord_threshold_raw_map(TRUE),
                fit_check = ord_fit_check),
    drop_intercept = TRUE
  )
  ord_family_tail(fam, lk, threshold, make_lpdf,
                  function(tmap) ord_sim("cumulative", link = lk,
                                         tmap = tmap))
}

#' The fields every ordinal constructor sets after `frmtmb_family()`:
#' `disc` held at one unless the formula models it, the threshold
#' structure, and the two factories that `thres_finalizer()` calls to
#' rebuild the density and the simulator once the threshold count is
#' known. A structure other than `"flexible"` needs that count, and the
#' count is a fact of the response.
#'
#' `disc` held at one is also a HIDDEN default: brms shows it nowhere
#' unless the formula models it (no Links entry, no coefficient), so
#' the output a fit prints leaves it out (`lp_hidden_fixed()`) while
#' the objective keeps the mapped coefficient.
#'
#' @noRd
ord_family_tail <- function(fam, lk, threshold, make_lpdf, make_sim) {
  fam[["fixed_dpars"]] <- list(disc = 1)
  fam[["hidden_fixed_dpars"]] <- "disc"
  fam[["threshold"]] <- threshold
  fam[["ord_lpdf_make"]] <- make_lpdf
  fam[["ord_sim_make"]] <- make_sim
  ord_tag_link(fam, lk)
}

#' The code of an ordinal family's first category: 1, or 0 for a family
#' with brms's `extra_cat` special, whose response has a category 0
#' below the ordinal ones (`hurdle_cumulative()`). An ordinal fit's
#' category probabilities, simulated codes and category means read their
#' codes from here rather than from `seq_len(K)`.
#'
#' @noRd
ord_code0 <- function(fam) if (isTRUE(fam[["extra_cat"]])) 0L else 1L

#' The ordinal families whose thresholds are held ordered, as
#' `(tau_1, log increments)`: brms declares them `ordered` for these two
#' (`brms:::has_ordered_thres()`). `cs()` is refused for both, as its
#' offsets could make a difference of their probabilities negative.
#'
#' @noRd
ord_ordered_families <- c("cumulative", "hurdle_cumulative")

#' The cumulative log-probability of categories `y` in `1..K`, with the
#' distribution function read at `disc * (tau - eta)`: the data path of
#' the densities of `fam_cumulative()` and `fam_hurdle_cumulative()`.
#'
#' @noRd
ord_cumulative_logpmf <- function(y, eta, tau, lk, disc = 1) {
  Fcdf <- lk$linkinv
  # Any link carrying `logit_eta` has an exact log-space difference
  # (see below), because that field turns its CDF into a logistic one.
  # cauchit is the ordinal link that does not, and does not need one:
  # its tails are polynomial, so the plain difference never saturates.
  q <- lk[["logit_eta"]]
  K1 <- length(tau)
  K <- K1 + 1L
  iK <- as.numeric(y == K)
  i1 <- as.numeric(y == 1)
  if (is.null(q)) {
    up <- Fcdf(disc * (tau[pmin(y, K1)] - eta)) * (1 - iK) + iK
    lo <- Fcdf(disc * (tau[pmax(y - 1, 1)] - eta)) * (1 - i1)
    return(log(up - lo))
  }
  # log(F(a) - F(b)) entirely in log space. `q` is the log odds of F,
  # so F is the logistic of it exactly, and for a logistic
  # F(a) - F(b) = (e^-qb - e^-qa) / ((1 + e^-qa)(1 + e^-qb)). The
  # difference of two saturated CDFs loses every digit from |eta| = 20
  # and is exactly 0 by 40; this form has no such point. For the logit
  # q is the identity and the two arguments of logspace_sub differ by
  # disc (tau_k - tau_{k-1}), which does not move with eta; for the
  # others the gap does move with eta, so the conditioning is merely
  # good rather than fixed.
  out <- i1 * log_inv_logit(q(disc * (tau[1] - eta))) +
    iK * log1m_inv_logit(q(disc * (tau[K1] - eta)))
  if (K1 >= 2L) {
    # data-only clamp into the interior categories, so the masked rows
    # still evaluate a legal (strictly ordered) threshold pair
    ym <- pmin(pmax(y, 2L), K1)
    a <- q(disc * (tau[ym] - eta))
    b <- q(disc * (tau[ym - 1L] - eta))
    out <- out + (1 - i1 - iK) *
      (RTMB::logspace_sub(-b, -a) + log_inv_logit(a) + log_inv_logit(b))
  }
  out
}

#' Hurdle cumulative (brms 2.23.0 `hurdle_cumulative`), dpars `mu`, `hu`
#' and `disc`. The response is `0..K`: `P(Y = 0) = hu`, and
#' `P(Y = k) = (1 - hu) P_cum(k)` for `k` in `1..K`, where `P_cum` is
#' `fam_cumulative()`'s model with the distribution function read at
#' `disc * (tau_k - eta)`. `disc` is held at one unless the formula
#' models it, as in brms. An ordered-factor response takes its first
#' level as the category 0.
#'
#' @noRd
fam_hurdle_cumulative <- function(link = "logit", link_hu = "logit",
                                  link_disc = "log", threshold = "flexible") {
  threshold <- ord_threshold_arg(threshold, "hurdle_cumulative")
  nm <- "hurdle_cumulative"
  lk <- ord_link(link, nm)
  lk_hu <- dpar_link(link_hu, "hu", nm, dpar_links_unit)
  lk_disc <- dpar_link(link_disc, "disc", nm, dpar_links_positive)
  Fcdf <- lk$linkinv
  make_lpdf <- function(tmap) {
    force(tmap)
    function(y, dpars, aterms, extra) {
      tau <- tmap(extra$tau_raw)
      i0 <- as.numeric(y == 0)
      g <- dpar_log_complement(dpars, "hu", lk_hu)
      # a zero row reads category 1 in the ordinal term, which carries
      # weight 0 there
      base <- ord_cumulative_logpmf(pmax(y, 1), dpars[["mu"]], tau, lk,
                                    dpars[["disc"]] %||% 1)
      i0 * g$l + (1 - i0) * (g$l1m + base)
    }
  }
  make_sim <- function(tmap) {
    force(tmap)
    function(dpars, aterms, n, extra) {
      tau <- tmap(extra$tau_raw)
      P <- hurdle_cum_probs(rep(dpars[["mu"]], length.out = n),
                            rep(dpars[["disc"]] %||% 1, length.out = n),
                            tau, Fcdf)
      K <- ncol(P)
      cp <- t(apply(P, 1L, cumsum))
      if (n == 1L) cp <- matrix(cp, 1L, K)
      cat_ <- pmin(1L + rowSums(cp < stats::runif(n)), K)
      hu <- rep(dpars[["hu"]], length.out = n)
      ifelse(stats::runif(n) < hu, 0L, cat_)
    }
  }
  fam <- frmtmb_family(
    nm,
    accepts_aterms = c("weights", "thres"),
    family_finalize = hurdle_thres_finalizer(lk),
    dpars = c("mu", "hu", "disc"),
    links = list(mu = "identity", hu = lk_hu, disc = lk_disc),
    lpdf = make_lpdf(ord_tau_from_raw_ad),
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y != round(y))) {
        frm_stop("Family '", nm, "' requires either non-negative ",
                 "integers or ordered factors as responses: 0 for the ",
                 "hurdle and 1..K for the ordinal categories",
                 call. = FALSE)
      }
      if (max(y) < 2) {
        frm_stop("Could not extract the number of thresholds. Family '",
                 nm, "' needs at least 2 ordinal categories above the ",
                 "hurdle, and the response reaches ", max(y),
                 call. = FALSE)
      }
    },
    init_dpars = list(
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "ordinal",
    extra_pars = function(y, aterms) {
      ord_tau_init(y[y > 0], ordered = TRUE, link = lk, K = max(y))
    },
    sim = make_sim(function(raw) ord_tau_from_raw(raw, TRUE)),
    post = list(ord_thresholds = ord_threshold_map(TRUE),
                ord_thresholds_raw = ord_threshold_raw_map(TRUE),
                fit_check = ord_fit_check),
    drop_intercept = TRUE
  )
  fam[["extra_cat"]] <- TRUE
  ord_family_tail(fam, lk, threshold, make_lpdf, make_sim)
}

#' `ord_tau_from_raw()` for an ordered threshold vector, written so that
#' it tapes: the same map `fam_cumulative()`'s density forms inline.
#'
#' @noRd
ord_tau_from_raw_ad <- function(raw) {
  "[<-" <- RTMB::ADoverload("[<-")
  K1 <- length(raw)
  tau <- rep(raw[1], K1)
  if (K1 > 1) for (k in 2:K1) tau[k] <- tau[k - 1] + exp(raw[k])
  tau
}

#' The `n x K` probabilities of the ordinal categories `1..K` of a
#' cumulative model with discrimination `disc`, in plain doubles.
#'
#' @noRd
hurdle_cum_probs <- function(eta, disc, tau, Fcdf) {
  n <- length(eta)
  M <- disc * (matrix(tau, n, length(tau), byrow = TRUE) - eta)
  Fm <- cbind(0, Fcdf(M), 1)
  Fm[, -1L, drop = FALSE] - Fm[, -ncol(Fm), drop = FALSE]
}

#' `hurdle_cumulative()`'s `family_finalize()` slot: `thres(x = )` and
#' the threshold structure are resolved by the ordinal families' own
#' finalizer on the rows above the hurdle, and `thres(gr = )` is
#' refused.
#'
#' The grouped densities `thres_finalizer()` builds read `mu` alone and
#' have no hurdle, so accepting `gr = ` would silently fit the model
#' without it.
#'
#' @noRd
hurdle_thres_finalizer <- function(lk) {
  inner <- thres_finalizer("cumulative", ordered = TRUE, link = lk)
  function(fam, y, aterms) {
    if (!is.null(aterms[["thres_gr"]])) {
      frm_stop("hurdle_cumulative() takes thres(x = ) but not ",
               "thres(gr = ): grouped thresholds are not implemented ",
               "for the hurdle family. brms fits them", call. = FALSE)
    }
    if (is.null(aterms[["thres"]]) &&
          identical(fam[["threshold"]] %||% "flexible", "flexible")) {
      return(fam)
    }
    pos <- y > 0
    ax <- aterms
    if (length(ax[["thres"]]) == length(y)) ax[["thres"]] <- ax[["thres"]][pos]
    out <- inner(fam, y[pos], ax)
    # the ordinal finalizer's start values count the categories over
    # every row it is handed, and a hurdle zero is not one of them
    ep <- out[["extra_pars"]]
    out[["extra_pars"]] <- function(y, aterms) ep(y[y > 0], aterms)
    out
  }
}

#' An ordinal family's fit-end check: the ordinal finalizer's
#' unplaced-threshold check when `thres(x = )` asked for one, and an
#' intercept in `disc` that nothing identifies.
#'
#' `disc * (tau - eta)` is unchanged when `disc` is divided by a
#' constant and the thresholds and `mu`'s coefficients multiplied by it,
#' so the likelihood cannot place an intercept in `disc`: a fit of
#' `disc ~ 1 + z` finished silently with standard errors of 12 to 67
#' (dev/fams2-rev-disc1.txt). brms accepts the same formula without a
#' word and identifies it by its default `normal(0, 1)` prior on the
#' intercept (dev/fams2-p1-brmsdisc.txt), which frmtmb fits too when
#' that prior is set, so this warns rather than refuses, and is silent
#' when a prior holds the intercept. A prior on the thresholds or on
#' `mu`'s coefficients pins the scale through that prior, and then it
#' says so rather than that no standard error is usable. The same
#' holds for every ordinal family, because each reads its distribution
#' function at `disc` times a difference of thresholds and `mu`.
#'
#' @noRd
ord_fit_check <- function(fit, resp) {
  fam <- fit$spec$responses[[resp]]$family
  if (length(fam[["thres"]][["unident"]])) thres_fit_check(fit, resp)
  lp <- fit$frame[["linpreds"]][[linpred_key(resp, "disc")]]
  if (is.null(lp) || !is.null(lp[["constant"]])) return(invisible(NULL))
  j <- which(colnames(lp[["X"]]) == "(Intercept)")
  # a design without an intercept column can still span one: the cell
  # means of `disc ~ 0 + h` add up to it, and then a common shift of
  # those coefficients is the same scale the likelihood cannot place
  if (!length(j)) {
    X <- as.matrix(lp[["X"]][, seq_len(lp[["n_param_cols"]] %||%
                                          ncol(lp[["X"]])), drop = FALSE])
    if (!ncol(X)) return(invisible(NULL))
    one <- rep(1, nrow(X))
    res <- qr.resid(qr(X), one)
    if (sqrt(sum(res^2)) > sqrt(.Machine$double.eps) * sqrt(nrow(X))) {
      return(invisible(NULL))
    }
    j <- seq_len(ncol(X))
    spans <- TRUE
  } else {
    spans <- FALSE
  }
  ent <- if (!is.null(fit$prior)) {
    resolve_prior_input(list(frame = fit$frame, spec = fit$spec),
                        fit$prior)$entries
  }
  held <- any(vapply(ent, function(e) {
    identical(e$comp, lp[["par"]]) && any(lp[["idx"]][j] %in% e$idx)
  }, NA))
  if (held) return(invisible(NULL))
  # a prior on the thresholds or on mu's coefficients pins the common
  # scale, and with it the intercept, if only through that prior: the
  # re-check's fit with a prior on class Intercept alone had standard
  # errors of at most 1.71 (dev/fams2-rev2-guards.txt, section 4)
  lpm <- fit$frame[["linpreds"]][[linpred_key(resp, "mu")]]
  tau_nm <- extra_tpl_name(fit$frame, resp, "tau_raw")
  pinned <- any(vapply(ent, function(e) {
    identical(e$comp, tau_nm) ||
      (!is.null(lpm) && identical(e$comp, lpm[["par"]]) &&
         any(lpm[["idx"]] %in% e$idx))
  }, NA))
  frm_warning(fam[["family"]], ": disc has an intercept",
              if (spans) {
                paste0(" (its columns, ",
                       paste(colnames(lp[["X"]])[j], collapse = ", "),
                       ", add up to one)")
              },
              ", which the ",
              "likelihood cannot tell apart from the scale of the ",
              "thresholds, so ",
              if (pinned) {
                paste0("only the priors on the thresholds or on the ",
                       "coefficients of mu place it. ")
              } else {
                paste0("it, the thresholds and the coefficients of mu ",
                       "have no usable standard error. ")
              },
              if (spans) {
                paste0("Leave out one of those columns, or hold them with ",
                       "a prior on class = \"b\", dpar = \"disc\"")
              } else {
                paste0("Write disc ~ 0 + ..., or hold the intercept with a ",
                       "prior, as brms does with its default ",
                       "set_prior(\"normal(0, 1)\", class = \"Intercept\", ",
                       "dpar = \"disc\")")
              }, call. = FALSE)
  invisible(NULL)
}

#' Shared scaffolding for the sequential ordinal families: an
#' n x (K-1) matrix of `(tau_j - eta_i)` or `(eta_i - tau_j)`, and
#' data-only indicator matrices selecting the observed category
#' (branch-free).
#'
#' @noRd
ord_indicators <- function(y, K1) {
  n <- length(y)
  jj <- rep(seq_len(K1), each = n)
  yy <- rep(y, K1)
  list(
    sel = matrix(as.numeric(yy == jj), n, K1),   # j == y (y <= K-1)
    below = matrix(as.numeric(jj < yy), n, K1)   # j < y
  )
}

#' The n x (K-1) matrix of `tau_j - eta_i` the ordinal lpdfs work over.
#' It broadcasts by matrix multiplication because `rep()` would strip the
#' advector class.
#'
#' @noRd
ord_eta_mat <- function(eta, tau, n, K1) {
  # broadcast via matmul: rep() strips the advector class
  TM <- RTMB::matrix(1, n, 1) %*% RTMB::matrix(tau, 1, K1)
  TM - eta   # column-wise recycling: row i is tau_j - eta_i
}

#' Sequential (sratio/cratio) log-density with the response on the tape.
#' Column-at-a-time so nothing indexes an advector and nothing needs the
#' n x (K-1) matrices, which the data path builds for speed.
#'
#' @noRd
ord_seq_lpdf_ad <- function(y, eta, tau, K1, cs, Fcdf, stopping,
                            disc = 1) {
  sel <- ord_cat_sel(y, K1 + 1L)
  below <- ord_cat_below(sel, K1 + 1L)
  out <- 0
  for (j in seq_len(K1)) {
    Mj <- tau[j] - eta
    if (!is.null(cs)) Mj <- Mj - cs[, j]
    Mj <- disc * Mj
    Pj <- if (stopping) Fcdf(Mj) else 1 - Fcdf(-Mj)
    out <- out + sel[[j]] * log(Pj) + below[[j]] * log(1 - Pj)
  }
  out
}

#' Make the response validator the ordinal families share: the response
#' must be integers `1..K`, and at least two categories must be
#' observed.
#'
#' @noRd
ord_valid_y <- function(name) {
  force(name)
  # brms's sentence opens each refusal, so a message matched on in a
  # ported script reads the same
  function(y, aterms) {
    if (any(y < 1) || any(y != round(y))) {
      frm_stop("Family '", name, "' requires either positive integers or ",
               "ordered factors as responses: integer codes 1..K",
               call. = FALSE)
    }
    if (length(unique(y)) < 2) {
      frm_stop("Family '", name, "' needs at least 2 observed response ",
               "categories", call. = FALSE)
    }
  }
}

#' Starting values for the K-1 ordinal thresholds, from the observed
#' cumulative category frequencies. With `ordered = TRUE` they come back
#' in the (first threshold, log increments) parameterization. `K` is
#' larger than `max(y)` when `thres(x = )` names categories nobody chose.
#'
#' @noRd
ord_tau_init <- function(y, ordered = TRUE, link = "logit", K = max(y)) {
  p <- cumsum(tabulate(y, K) / length(y))[-K]
  p <- pmin(pmax(p, 0.01), 0.99)
  # the thresholds live on the link's own scale, so the observed
  # cumulative proportions are mapped through that link and not always
  # through qlogis: a probit threshold is about 0.6 of a logit one, and
  # starting a cloglog fit on logit thresholds starts it skewed
  tau0 <- get_link(link)$linkfun(p)
  if (!ordered) return(list(tau_raw = tau0))
  incr <- pmax(diff(tau0), 0.05)
  list(tau_raw = c(tau0[1], log(incr)))
}

# --- ordinal simulators ----------------------------------------------
# The lpdfs are taped, so they answer "how likely is this category"; a
# simulator needs the whole category distribution instead. These build
# it in plain doubles, one branch per family, matching each lpdf term
# for term.

#' The map from the internal threshold vector to the thresholds
#' themselves, as a family declares it in `post$ord_thresholds`.
#'
#' `cumulative()` estimates `(tau_1, log increments)` so that its
#' thresholds cannot cross, and the other three ordinal families
#' estimate the thresholds directly, as brms declares them. `variables()`
#' and `hypothesis()` report the THRESHOLDS, brms's `b_Intercept[k]`, so
#' they need the map; asking the family for it keeps an ordinal family
#' shipped by another package working, which a list of four names here
#' would not.
#'
#' @noRd
ord_threshold_map <- function(ordered) {
  function(raw) ord_tau_from_raw(raw, ordered)
}

#' The inverse of `ord_threshold_map()`, from the thresholds back to the
#' internal vector, as a family declares it in `post$ord_thresholds_raw`.
#'
#' frmtmb.sample stores a draw's thresholds under brms's names and on
#' brms's scale, `b_Intercept[k]`, and hands the internal vector back to
#' the model at every draw. A family that declares only the forward map
#' keeps its internal names in the draws, because a draws column must
#' not carry a name that its values do not have.
#'
#' @noRd
ord_threshold_raw_map <- function(ordered) {
  function(tau) ord_raw_from_tau(tau, ordered)
}

#' The inverse of `ord_tau_from_raw()`.
#'
#' @noRd
ord_raw_from_tau <- function(tau, ordered) {
  if (!ordered || length(tau) < 2L) return(tau)
  c(tau[1L], log(diff(tau)))
}

#' cumulative stores ordered thresholds as (tau_1, log increments);
#' sratio, cratio and acat store them raw. This returns the thresholds
#' themselves.
#'
#' @noRd
ord_tau_from_raw <- function(raw, ordered) {
  if (!ordered || length(raw) < 2L) return(raw)
  c(raw[1L], raw[1L] + cumsum(exp(raw[-1L])))
}

#' Resolve an ordinal family's `link` against the distribution
#' functions its thresholds can be read through.
#'
#' An ordinal `link` is not a link on a mean: it names the cumulative
#' distribution function applied to `tau_j - eta`, so only a link whose
#' INVERSE maps onto the unit interval can serve. The set is brms's own,
#' from `brms_mu_links`: brms allows softit for cumulative and refuses it
#' for the sequential pair, and that difference is kept rather than
#' tidied away, because the point of the roster is that a brms model
#' ports.
#'
#' @noRd
ord_link <- function(link, family, choices = brms_mu_links[[family]]) {
  if (!is.character(link) || length(link) != 1L || is.na(link) ||
      !link %in% choices) {
    head <- if (is.character(link) && length(link) == 1L && !is.na(link)) {
      paste0("'", link, "' is not a supported link for family '", family,
             "'. ")
    } else {
      paste0(family, "(link =) takes a single link name, not ",
             arg_desc(link), ". ")
    }
    frm_stop(head, "Supported links are: ", link_set_text(choices),
             ". An ordinal link names the distribution function the ",
             "thresholds are read through, so it has to map onto (0, 1). ",
             "See ?`frmtmb-links`", call. = FALSE)
  }
  get_link(link)
}

#' Check an ordinal constructor's `threshold` argument, brms's
#' `"flexible"`, `"equidistant"` or `"sum_to_zero"`, and return it. The
#' structure itself is built by `thres_finalizer()`, because every
#' structure but `"flexible"` needs the threshold count, which is a fact
#' of the response.
#'
#' @noRd
ord_threshold_arg <- function(threshold, family) {
  choices <- c("flexible", "equidistant", "sum_to_zero")
  if (!is.character(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !threshold %in% choices) {
    frm_stop(family, "(threshold =) takes one of ",
             paste0("'", choices, "'", collapse = ", "), ", not ",
             arg_desc(threshold), call. = FALSE)
  }
  threshold
}

#' The CDF an ordinal family reads its thresholds through: the inverse
#' of its link. One function serves the taped density and the plain
#' numeric simulators alike, because every inverse link in the registry
#' is written over arithmetic that RTMB overloads and that base R
#' evaluates unchanged on doubles.
#'
#' @noRd
ord_cdf <- function(link) get_link(link)$linkinv

#' Resolve `acat()`'s link against brms's roster for the family: the
#' logit, which reads the categories in brms's log-linear form, and
#' every other unit-interval link, which reads them in brms's second
#' form (`acat_general_E()`).
#'
#' @noRd
acat_link <- function(link) {
  ord_link(link, "acat")
}

#' The unnormalized log category probabilities of brms's `acat` off the
#' logit, `brms:::inv_link_acat()`: category `k` of `K` is
#' `prod_{j < k} F(x_j) * prod_{j >= k} (1 - F(x_j))`, with
#' `x_j = disc * (eta + cs_j - tau_j)`. Written as the running sum
#' `E_1 = sum_j log(1 - F(x_j))`, `E_k = E_{k-1} + log F(x_{k-1}) -
#' log(1 - F(x_{k-1}))`, one threshold at a time, so it tapes and works
#' on plain doubles alike. Under the logit `log F - log(1 - F)` is `x`
#' itself, which is why the log-linear form is the same model there.
#'
#' `xs` is the list of the `K - 1` columns `x_j`, `live` (optional) the
#' list of data masks that switch a threshold position off for a row
#' whose own vector is shorter (grouped thresholds). Log space through
#' the link's log-odds form where it has one, and the plain logarithms
#' otherwise.
#'
#' @noRd
acat_general_E <- function(xs, lk, live = NULL) {
  lg <- ord_log_cdf_pair(lk)
  if (is.null(lg)) {
    Fcdf <- lk$linkinv
    lg <- list(lF = function(v) log(Fcdf(v)),
               l1mF = function(v) log(1 - Fcdf(v)))
  }
  K1 <- length(xs)
  lF <- lapply(xs, lg$lF)
  l1 <- lapply(xs, lg$l1mF)
  if (!is.null(live)) {
    lF <- Map(`*`, lF, live)
    l1 <- Map(`*`, l1, live)
  }
  E <- vector("list", K1 + 1L)
  E1 <- l1[[1L]]
  if (K1 > 1L) for (j in 2:K1) E1 <- E1 + l1[[j]]
  E[[1L]] <- E1
  for (k in 2:(K1 + 1L)) {
    E[[k]] <- E[[k - 1L]] + lF[[k - 1L]] - l1[[k - 1L]]
  }
  E
}

#' Record the distribution function an ordinal family reads its
#' thresholds through, so `summary()` can name the scale its
#' coefficients are on.
#'
#' It cannot go in `links`. That slot is what the objective applies to
#' a linear predictor, and an ordinal `mu` genuinely has an identity
#' link: the CDF is applied to `tau_j - eta`, never to `eta` alone.
#'
#' @noRd
ord_tag_link <- function(fam, lk) {
  fam[["ord_link"]] <- lk
  fam
}

#' n x K matrix of category probabilities, one branch per ordinal
#' family, in plain doubles.
#'
#' @noRd
ord_cat_probs <- function(family, eta, tau, cs, link, disc = 1) {
  n <- length(eta)
  K1 <- length(tau)
  K <- K1 + 1L
  if (identical(family, "acat") &&
        !identical(if (is.list(link)) link$name else link, "logit")) {
    xs <- lapply(seq_len(K1), function(j) {
      v <- eta - tau[j]
      if (!is.null(cs)) v <- v + cs[, j]
      disc * v
    })
    E <- do.call(cbind, acat_general_E(xs, get_link(link)))
    ex <- exp(E - apply(E, 1L, max))
    return(ex / rowSums(ex))
  }
  if (identical(family, "acat")) {
    # P(y=r) proportional to exp((r-1) eta - cumsum(tau)[r]); the row
    # maximum comes out before exp() so a wide eta cannot overflow
    ct0 <- c(0, cumsum(tau))
    E <- outer(eta, seq_len(K) - 1) - matrix(ct0, n, K, byrow = TRUE)
    if (!is.null(cs)) {
      acc <- rep(0, n)
      for (r in seq.int(2L, K)) {
        acc <- acc + cs[, r - 1L]
        E[, r] <- E[, r] + acc
      }
    }
    # brms reads acat at disc * (eta - tau_j), so disc scales the
    # whole exponent of every category
    E <- disc * E
    ex <- exp(E - apply(E, 1L, max))
    return(ex / rowSums(ex))
  }
  Fcdf <- ord_cdf(link)
  M <- matrix(tau, n, K1, byrow = TRUE) - eta   # tau_j - eta_i
  if (!is.null(cs)) M <- M - cs
  M <- disc * M
  if (identical(family, "cumulative")) {
    # P(y=k) = F(tau_k - eta) - F(tau_{k-1} - eta), with the two
    # boundary values pinned at 0 and 1: a column-wise difference of the
    # K+1 cumulative probabilities
    Fm <- cbind(0, Fcdf(M), 1)
    return(Fm[, -1L, drop = FALSE] - Fm[, -ncol(Fm), drop = FALSE])
  }
  # sequential families: h_j is the chance of stopping at category j
  # given survival past j-1, on each family's own scale
  h <- if (identical(family, "sratio")) Fcdf(M) else 1 - Fcdf(-M)
  P <- matrix(0, n, K)
  surv <- rep(1, n)
  for (j in seq_len(K1)) {
    P[, j] <- surv * h[, j]
    surv <- surv * (1 - h[, j])
  }
  P[, K] <- surv
  P
}

#' Make the simulator for an ordinal family. It draws one category per
#' row by inverse-CDF sampling of the full category distribution, which
#' the taped log-density does not give. `tmap` maps the internal
#' threshold vector to the thresholds, in plain doubles.
#'
#' @noRd
ord_sim <- function(family, link, tmap) {
  force(tmap)
  function(dpars, aterms, n, extra) {
    tau <- tmap(extra$tau_raw)
    P <- ord_cat_probs(family, rep(dpars[["mu"]], length.out = n), tau,
                       dpars[[".cs"]], link,
                       rep(dpars[["disc"]] %||% 1, length.out = n))
    K <- ncol(P)
    # inverse-CDF sampling, one uniform per row
    cp <- t(apply(P, 1L, cumsum))
    if (n == 1L) cp <- matrix(cp, 1L, K)
    pmin(1L + rowSums(cp < stats::runif(n)), K)
  }
}

#' The log hazard pair an ordinal family's robust branch needs, or
#' `NULL` when the link has no exact log-odds form.
#'
#' `logit_eta` turns the link's own CDF into a logistic one: if `q` is
#' the log odds of `F(x)` then `F(x)` is `plogis(q(x))` exactly, so
#' `log F` and `log(1 - F)` both come out of `logspace_add()` and
#' neither saturates. For the logit `q` is the identity and this
#' reduces to the pair the sequential families already used.
#'
#' @noRd
ord_log_cdf_pair <- function(lk) {
  q <- lk[["logit_eta"]]
  if (is.null(q)) return(NULL)
  list(lF = function(x) log_inv_logit(q(x)),
       l1mF = function(x) log1m_inv_logit(q(x)))
}

#' The sequential families' log-density from log-space hazards. `lstop`
#' gives `log P(stop at j)` and `lgo` gives `log P(continue past j)`,
#' both as functions of the `tau_j - eta` matrix column; sratio and
#' cratio differ only in which way round they are. One column at a time
#' rather than one matrix operation, because `logspace_add()` works on
#' the flat vector and would drop the dimensions the indicator matrices
#' need.
#'
#' @noRd
ord_log_hazard_sum <- function(M, ind, K1, lstop, lgo) {
  out <- 0
  for (j in seq_len(K1)) {
    Mj <- M[, j]
    out <- out + ind$sel[, j] * lstop(Mj) + ind$below[, j] * lgo(Mj)
  }
  out
}

#' Stopping ratio (brms sratio): `P(y=k) = F(tau_k - eta) *
#' prod_{j<k} (1 - F(tau_j - eta))`; unordered thresholds, like cratio.
#' Every distribution function is read at `disc` times its argument.
#'
#' Each category probability is a product of hazards, so it is positive
#' whatever order the thresholds take, and brms 2.23.0 declares sratio's
#' thresholds as a plain `vector` (`brms:::has_ordered_thres()` is FALSE
#' here). Holding them ordered put the estimate on the ordering boundary
#' wherever the unconstrained optimum has crossing thresholds, away from
#' brms's mode (dev/sratio-findings.md).
#'
#' @noRd
fam_sratio <- function(link = "logit", link_disc = "log",
                       threshold = "flexible") {
  fam_sequential("sratio", link, link_disc, threshold)
}

#' Continuation ratio (brms cratio): `P(y=k) = (1 - F(eta - tau_k)) *
#' prod_{j<k} F(eta - tau_j)`; unordered thresholds. Every distribution
#' function is read at `disc` times its argument.
#'
#' @noRd
fam_cratio <- function(link = "logit", link_disc = "log",
                       threshold = "flexible") {
  fam_sequential("cratio", link, link_disc, threshold)
}

#' The sequential pair, which differ only in which way round the
#' distribution function is read: sratio stops at `F(M)` and cratio at
#' `1 - F(-M)`, with `M = disc * (tau_j - eta - cs_j)`.
#'
#' @noRd
fam_sequential <- function(name, link, link_disc, threshold) {
  stopping <- identical(name, "sratio")
  threshold <- ord_threshold_arg(threshold, name)
  lk <- ord_link(link, name)
  lk_disc <- dpar_link(link_disc, "disc", name, dpar_links_positive)
  Fcdf <- lk$linkinv
  lgm <- ord_log_cdf_pair(lk)
  # sratio's h_j = F(M) is read at M itself. cratio reads its CDF at -M,
  # and the pair has to be read there too. The logit-only version relied
  # on 1 - F(-x) = F(x) to stay at M instead, which is true of the
  # logistic, the normal, the Cauchy and the cubic of probit_approx, and
  # false of the cloglog. Composing with the negation costs nothing and
  # holds for an asymmetric CDF.
  lg <- if (is.null(lgm) || stopping) lgm else {
    list(lF = function(x) lgm$l1mF(-x), l1mF = function(x) lgm$lF(-x))
  }
  make_lpdf <- function(tmap) {
    force(tmap)
    function(y, dpars, aterms, extra) {
      tau <- tmap(extra$tau_raw)
      K1 <- length(tau)
      n <- length(y)
      disc <- dpars[["disc"]] %||% 1
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        return(ord_seq_lpdf_ad(ov$y, dpars[["mu"]], tau, K1,
                               dpars[[".cs"]], Fcdf, stopping = stopping,
                               disc = disc) * ov$keep)
      }
      M <- ord_eta_mat(dpars[["mu"]], tau, n, K1)   # tau_j - eta
      if (!is.null(dpars[[".cs"]])) M <- M - dpars[[".cs"]]
      # column-wise recycling: row i is scaled by disc_i
      M <- M * disc
      ind <- ord_indicators(y, K1)
      if (!is.null(lg)) {
        # log F and log(1 - F) straight out of logspace_add: the naive
        # pair is -Inf on whichever side the CDF saturated, which for
        # sratio is the negative eta tail and for cratio its mirror
        return(ord_log_hazard_sum(M, ind, K1, lg$lF, lg$l1mF))
      }
      ones <- rep(1, K1)   # rowSums strips the advector class
      if (stopping) {
        P <- Fcdf(M)
        return(as.vector((log(P) * ind$sel) %*% ones) +
                 as.vector((log(1 - P) * ind$below) %*% ones))
      }
      P <- Fcdf(-M)                            # F(disc (eta + cs - tau))
      as.vector((log(1 - P) * ind$sel) %*% ones) +
        as.vector((log(P) * ind$below) %*% ones)
    }
  }
  make_sim <- function(tmap) ord_sim(name, link = lk, tmap = tmap)
  fam <- frmtmb_family(
    name,
    accepts_aterms = c("weights", "thres"),
    family_finalize = thres_finalizer(name, ordered = FALSE, link = lk),
    dpars = c("mu", "disc"),
    links = list(mu = "identity", disc = lk_disc),
    lpdf = make_lpdf(identity),
    valid_y = ord_valid_y(name),
    type = "ordinal",
    extra_pars = function(y, aterms) {
      ord_tau_init(y, ordered = FALSE, link = lk)
    },
    sim = make_sim(identity),
    post = list(ord_thresholds = ord_threshold_map(FALSE),
                ord_thresholds_raw = ord_threshold_raw_map(FALSE),
                fit_check = ord_fit_check),
    drop_intercept = TRUE
  )
  ord_family_tail(fam, lk, threshold, make_lpdf, make_sim)
}

#' Adjacent category (brms acat): under the logit link `P(y=k)`
#' proportional to `exp(disc * sum_{j<k} (eta - tau_j))`; under any other
#' link brms's second form (`acat_general_E()`). Unordered thresholds.
#'
#' @noRd
fam_acat <- function(link = "logit", link_disc = "log",
                     threshold = "flexible") {
  threshold <- ord_threshold_arg(threshold, "acat")
  lk <- acat_link(link)
  lk_disc <- dpar_link(link_disc, "disc", "acat", dpar_links_positive)
  make_lpdf <- function(tmap) {
    force(tmap)
    if (!identical(lk$name, "logit")) return(acat_general_lpdf(tmap, lk))
    function(y, dpars, aterms, extra) {
      tau <- tmap(extra$tau_raw)
      K <- length(tau) + 1L
      n <- length(y)
      eta <- dpars[["mu"]]
      disc <- dpars[["disc"]] %||% 1
      "c" <- RTMB::ADoverload("c")
      ct0 <- c(0, cumsum(tau))                 # length K
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        sel <- ord_cat_sel(ov$y, K)
        acc <- 0 * eta
        num <- 0
        den <- NULL
        for (r in seq_len(K)) {
          Er <- (r - 1) * eta - ct0[r]
          if (!is.null(dpars[[".cs"]]) && r >= 2L) {
            acc <- acc + dpars[[".cs"]][, r - 1L]
            Er <- Er + acc
          }
          Er <- disc * Er
          num <- num + sel[[r]] * Er
          # logsumexp fold: the top category's exponent is (K-1) * eta,
          # so a plain sum of exp() overflows at eta = 709 / (K - 1)
          den <- if (is.null(den)) Er else RTMB::logspace_add(den, Er)
        }
        return((num - den) * ov$keep)
      }
      # E[i, r] = (r-1) * eta_i - cumsum tau, r = 1..K; broadcast by
      # matmul (rep() strips the advector class)
      Rm <- matrix(rep(seq_len(K) - 1L, each = n), n, K)
      E <- Rm * eta -
        RTMB::matrix(1, n, 1) %*% RTMB::matrix(ct0, 1, K)
      if (!is.null(dpars[[".cs"]])) {
        "[<-" <- RTMB::ADoverload("[<-")
        CS <- dpars[[".cs"]]
        acc <- 0 * eta
        for (r in seq.int(2L, K)) {
          acc <- acc + CS[, r - 1L]
          E[, r] <- E[, r] + acc
        }
      }
      # brms: cumulative_sum(disc * (mu - thres)), so disc scales every
      # category's exponent; column-wise recycling scales row i by disc_i
      E <- E * disc
      jj <- rep(seq_len(K), each = n)
      S <- matrix(as.numeric(rep(y, K) == jj), n, K)
      ones <- rep(1, K)   # rowSums strips the advector class
      den <- E[, 1L]
      for (r in seq_len(K)[-1L]) den <- RTMB::logspace_add(den, E[, r])
      as.vector((E * S) %*% ones) - den
    }
  }
  make_sim <- function(tmap) ord_sim("acat", link = lk, tmap = tmap)
  fam <- frmtmb_family(
    "acat",
    accepts_aterms = c("weights", "thres"),
    family_finalize = thres_finalizer("acat", ordered = FALSE, link = lk),
    dpars = c("mu", "disc"),
    links = list(mu = "identity", disc = lk_disc),
    lpdf = make_lpdf(identity),
    valid_y = ord_valid_y("acat"),
    type = "ordinal",
    extra_pars = function(y, aterms) ord_tau_init(y, ordered = FALSE),
    sim = make_sim(identity),
    post = list(ord_thresholds = ord_threshold_map(FALSE),
                ord_thresholds_raw = ord_threshold_raw_map(FALSE),
                fit_check = ord_fit_check),
    drop_intercept = TRUE
  )
  ord_family_tail(fam, lk, threshold, make_lpdf, make_sim)
}

#' The log-density of `acat()` off the logit, from `acat_general_E()`,
#' on the data path and on the one-step path alike: the category is
#' picked by indicators rather than by indexing, and the normalizer is
#' a log-space fold over the categories.
#'
#' @noRd
acat_general_lpdf <- function(tmap, lk) {
  force(tmap)
  function(y, dpars, aterms, extra) {
    tau <- tmap(extra$tau_raw)
    K1 <- length(tau)
    eta <- dpars[["mu"]]
    disc <- dpars[["disc"]] %||% 1
    cs <- dpars[[".cs"]]
    xs <- lapply(seq_len(K1), function(j) {
      v <- eta - tau[j]
      if (!is.null(cs)) v <- v + cs[, j]
      disc * v
    })
    E <- acat_general_E(xs, lk)
    ov <- osa_unwrap(y)
    sel <- if (is.null(ov)) {
      lapply(seq_len(K1 + 1L), function(k) as.numeric(y == k))
    } else {
      ord_cat_sel(ov$y, K1 + 1L)
    }
    num <- 0
    den <- NULL
    for (k in seq_len(K1 + 1L)) {
      num <- num + sel[[k]] * E[[k]]
      den <- if (is.null(den)) E[[k]] else RTMB::logspace_add(den, E[[k]])
    }
    out <- num - den
    if (is.null(ov)) out else out * ov$keep
  }
}

#' Where a mixture component's mean starts. A component whose mean
#' lives on (0, 1) takes its start from the proportion y / trials, not
#' from y: the quantile of a count lies outside the logit's range, so
#' every component fell back to the link origin and
#' mixture(beta_binomial, beta_binomial) sat at two identical
#' components with no warning (the Bayesian Cognitive Modeling port
#' found it on the malingering data, logLik -64.26 against -58.71 from
#' a separated start). The quantile is clamped the way the component
#' families clamp their own start, because eight of twenty-two
#' respondents scoring 45 of 45 put the two-thirds quantile at exactly
#' 1, where the logit is infinite.
#'
#' @noRd
mixture_mu_bounded <- function(comp) {
  lk <- comp$links[["mu"]]
  nm <- if (is.list(lk)) lk[["name"]] else lk
  isTRUE(nm %in% c("logit", "probit", "probit_approx", "cauchit",
                   "cloglog", "softit"))
}

#' Refuse the component sets brms refuses, in brms's order and words.
#'
#' A mixture density is a weighted sum of component densities, and that
#' sum is a density only when every component is one with respect to the
#' same measure. A real-support density and an integer-support mass
#' function are not, so `mixture(lognormal, exgaussian, poisson())`
#' returned a family whose likelihood added a probability to a density
#' (dev/famlink-findings.md, item 5). A continuous hurdle and a
#' zero-inflated beta each put an atom at zero beside a density, which is
#' the same defect inside one component, and a categorical or
#' multinomial response is not a number the other components could
#' share. The barred families are `brms_no_mixture`, generated from
#' `brms:::no_mixture()` rather than listed here by hand.
#'
#' The type of a frmtmb family is its support: `continuous` is brms's
#' `real`, and `discrete`, `ordinal` and `categorical` are all `int`.
#'
#' @noRd
mixture_check_components <- function(comps) {
  types <- vapply(comps, function(cp) cp[["type"]] %||% "continuous", "")
  fams <- vapply(comps, function(cp) cp[["family"]], "")
  listed <- paste0(fams, " (", types, ")", collapse = ", ")
  if (any(types == "continuous") && any(types != "continuous")) {
    frm_stop("Cannot mix families with real and integer support. The ",
             "components are ", listed, ", and a weighted sum of a density ",
             "and a probability mass is not a likelihood", call. = FALSE)
  }
  ord <- types == "ordinal"
  if (any(ord) && !all(ord)) {
    frm_stop("Cannot mix ordinal and non-ordinal families. The components ",
             "are ", listed, call. = FALSE)
  }
  barred <- types == "categorical" | fams %in% brms_no_mixture
  if (any(barred)) {
    frm_stop("Some of the families are not allowed in mixture models: ",
             paste(unique(fams[barred]), collapse = ", "), ". A categorical ",
             "or multinomial response is not a value the other components ",
             "can share, and a continuous hurdle or a zero-inflated or ",
             "zero-one-inflated beta puts a probability mass at zero ",
             "(or one) beside a density",
             call. = FALSE)
  }
  invisible(NULL)
}

mixture_mu_start <- function(y, aterms, p, bounded) {
  if (!bounded) return(stats::quantile(y, p, names = FALSE))
  size <- aterms[["trials"]] %||% 1
  q <- stats::quantile(y / size, p, names = FALSE)
  min(max(q, 0.02), 0.98)
}

#' Finite mixture families
#'
#' `mixture(fam1, fam2, ...)` builds a K-component mixture: each
#' component keeps its own distributional parameters, suffixed by the
#' component index (`mu1`, `sigma1`, `mu2`, ...), and the mixing
#' proportions come from `theta1 ... theta{K-1}` (multinomial logit
#' against the last component, each with its own linear predictor, so
#' mixing weights may depend on covariates). The main model formula
#' applies to every component mean; override per component with
#' `bf(y ~ x, mu2 ~ 1)`.
#'
#' A formula for a mixing weight follows brms's rule. Write one for all
#' of `theta1 ... thetaK` but one, and the one left out is the reference
#' component, whose linear predictor is 0: `bf(y ~ 1, theta2 ~ x)` on
#' two components makes component 1 the reference, so `theta2`'s
#' coefficients are the log odds of component 2 against component 1.
#' A formula for fewer of them is refused, as brms refuses it.
#'
#' A mixing weight's RESPONSE scale is the softmax over the component
#' predictors, so `frm_linpred(type = "response", dpar = "theta1")` is a
#' probability while `type = "link"` stays the predictor the density
#' works on. Under `se.fit = TRUE` that probability's standard error is
#' the delta method through its OWN predictor, `p (1 - p)` times the
#' predictor's standard error. For two components that is exact. For
#' three or more the softmax also moves with the other components'
#' predictors, and those terms are dropped, so the standard error is
#' CONSERVATIVE: measured 5.5% to 26.1% wider than the joint delta
#' method on a three-component fit, never narrower.
#'
#' The likelihood is a parameter-branch-free logsumexp, so Laplace
#' machinery is untouched; the usual finite-mixture ML caveats apply
#' instead: the likelihood is invariant to component relabeling, so the
#' component means are initialized on spread-out response quantiles,
#' and multimodality is real (compare starts, or order the intercepts
#' with bounds: `set_prior("", class = "Intercept", dpar = ..., lb = )`
#' per component). Component families with extra parameters
#' (ordinal) are not supported.
#'
#' With `groups = ~g` the mixture moves to the group level (latent
#' classes): every observation of a group shares one class draw, and
#' the marginal likelihood sums the class assignment per group.
#' Continuous random effects, smooths, and gp() terms are allowed in
#' the component formulas - the class sum happens conditional on the
#' latent effects, so one Laplace approximation integrates them
#' (growth-mixture models). Random effects written in a component
#' formula are class-specific by construction; the Laplace
#' approximation of the class-mixture integrand is not exact even for
#' gaussian responses (a fraction of a log-likelihood unit in typical
#' well-separated problems). `quadrature = TRUE` makes the integral
#' numerically exact when the per-group integrand is univariate (one
#' scalar random intercept, in one class); with class-specific
#' intercepts in several classes the coordinates couple and quadrature
#' remains approximate - use `frmtmb.sample::check_laplace()` to judge.
#' Mixing-weight predictors are evaluated at each group's first row
#' (use group-constant covariates). [mixture_probs()] returns the
#' posterior class probabilities per group (or per observation for
#' ordinary mixtures), conditional on the random-effect modes.
#'
#' @param ... Two or more component families.
#' @param groups Optional one-sided formula naming the latent-class
#'   grouping factor.
#' @param order brms's argument. `NULL`, `"none"` or `FALSE` leaves the
#'   components unordered, which is what a fit here does. `"mu"` or
#'   `TRUE`, brms's ordering of the `mu` intercepts, is refused, and any
#'   other value is refused as invalid, in brms's words.
#' @return A `frmtmb_family`.
#' @examples
#' # two well-separated gaussian components
#' set.seed(3)
#' dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)),
#'                  x = rnorm(160))
#' fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
#' # one mu and sigma per component, plus the mixing weight theta1
#' fixef(fit)
#' # posterior class probability per observation
#' head(mixture_probs(fit))
#'
#' # the mixing weight can take its own predictor
#' frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()), data = dd)
#'
#' \donttest{
#' # latent classes: every observation of a group shares one class
#' set.seed(4)
#' n_g <- 40
#' cls <- rep(c(1, 2), each = n_g / 2)
#' dg <- data.frame(g = factor(rep(seq_len(n_g), each = 5)))
#' dg$y <- rnorm(nrow(dg), c(0, 4)[cls[as.integer(dg$g)]], 1)
#' fg <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
#'           data = dg)
#' head(mixture_probs(fg))   # one row per group, not per observation
#' }
#' @export
mixture <- function(..., groups = NULL, order = NULL) {
  named <- ...names()
  named <- named[!is.na(named) & nzchar(named)]
  if (length(named)) {
    # a name in the dots is an argument meant for mixture() itself, and
    # read as a component it came back as "not a supported family"
    frm_stop("mixture() has no argument ", paste0("`", named, "`",
                                                    collapse = ", "),
             ". It takes the component families, unnamed, `groups` and ",
             "`order`. brms's `nmix` is not supported: repeat a ",
             "component to use it twice", call. = FALSE)
  }
  mixture_check_order(order)
  comps <- lapply(list(...), as_frmtmb_family)
  K <- length(comps)
  if (K < 2L) {
    frm_stop("Expecting at least 2 mixture components. mixture() needs two ",
             "or more component families", call. = FALSE)
  }
  mixture_check_components(comps)
  for (cp in comps) {
    if (!is.null(cp$extra_pars) || isTRUE(cp$drop_intercept)) {
      frm_stop("mixture() does not support component family '", cp$family,
               "'", call. = FALSE)
    }
    if (!"mu" %in% cp$dpars) {
      frm_stop("mixture() components need a 'mu' parameter", call. = FALSE)
    }
  }
  mixture_build(comps, groups, ref = K)
}

#' brms's `mixture(order = )`, validated as brms validates it.
#'
#' brms orders the components by their `mu` intercepts under
#' `order = "mu"` (or `TRUE`), which identifies them for a sampler, and
#' leaves them alone under `"none"` (or `FALSE`). A fit here never
#' constrains them, so `"none"` is what it does and is accepted. `"mu"`
#' is refused rather than ignored: the caller asked for an ordering the
#' fit would not deliver. Any other value is refused in brms's words.
#'
#' @noRd
mixture_check_order <- function(order) {
  if (is.null(order)) return(invisible(NULL))
  if (length(order) != 1L) {
    frm_stop("Argument 'order' must be of length 1.", call. = FALSE)
  }
  if (is.character(order)) {
    if (!order %in% c("none", "mu")) {
      frm_stop("Argument 'order' is invalid. Valid options are: none, mu",
               call. = FALSE)
    }
  } else if (!is.logical(order) || is.na(order)) {
    frm_stop("Argument 'order' is invalid. Valid options are: none, mu, ",
             "or TRUE and FALSE for them", call. = FALSE)
  } else {
    order <- if (order) "mu" else "none"
  }
  if (identical(order, "mu")) {
    frm_stop("mixture(order = \"mu\") is not supported: this package does ",
             "not constrain the order of the components' mu intercepts. ",
             "The likelihood is the same under any labeling, so a fit ",
             "needs no ordering; to fix one, bound the intercepts with ",
             "set_prior(\"\", class = \"Intercept\", dpar = \"mu2\", lb = ) ",
             "and so on per component. order = \"none\" is accepted",
             call. = FALSE)
  }
  invisible(NULL)
}

#' The mixture family a formula's written thetas call for.
#'
#' brms predicts the mixing weights only when a formula is written for
#' all of them but one, and the one left unwritten is the reference
#' whose linear predictor is 0 (`stan_mixture()`'s `missing_id`,
#' measured on brms 2.23.0 in dev/correct-log/brms-theta.txt). So
#' `theta2 ~ x` on two components makes component 1 the reference, and
#' `theta1 ~ x` on three is refused. Without a written theta the family
#' is returned unchanged, against its last component, where brms's
#' simplex `theta1 ... thetaK` is read off it (`brms_coef_table()`).
#'
#' @noRd
mixture_theta_reference <- function(fam, written) {
  rebuild <- fam[["mix_rebuild"]]
  if (!is.function(rebuild)) return(fam)
  K <- fam[["mix"]][["K"]]
  th <- grep("^theta[0-9]+$", written, value = TRUE)
  if (!length(th)) return(fam)
  ids <- unique(as.integer(sub("^theta", "", th)))
  # a theta past K is not this family's, and the dpar check names it
  if (any(ids < 1L | ids > K)) return(fam)
  if (length(ids) != K - 1L) {
    frm_stop("Can only predict all but one mixing proportion: a mixture of ",
             K, " components takes a formula for ", K - 1L, " of theta1 ",
             "... theta", K, ", and the one left out is the reference ",
             "whose linear predictor is 0, as in brms. This formula has ",
             paste(sort(th), collapse = ", "), call. = FALSE)
  }
  ref <- setdiff(seq_len(K), ids)
  if (identical(ref, fam[["mix"]][["ref"]])) return(fam)
  rebuild(ref)
}

#' A mixture family whose mixing weights are the multinomial logit
#' against component `ref`: `theta<k>` for every `k` but `ref`, and
#' `ref`'s own log ratio held at 0.
#'
#' `mixture()` builds it against the last component. A formula written
#' for every theta but one moves the reference to that one, which is
#' brms's rule (`mixture_theta_reference()`), and the family is rebuilt
#' through `mix_rebuild` with the same components.
#'
#' @noRd
mixture_build <- function(comps, groups, ref) {
  K <- length(comps)
  dpars <- character(0)
  links <- list()
  for (k in seq_len(K)) {
    for (dp in comps[[k]]$dpars) {
      nm <- paste0(dp, k)
      dpars <- c(dpars, nm)
      links[[nm]] <- comps[[k]]$links[[dp]]
    }
  }
  theta_ids <- setdiff(seq_len(K), ref)
  for (k in theta_ids) {
    nm <- paste0("theta", k)
    dpars <- c(dpars, nm)
    links[[nm]] <- "identity"
  }

  comp_dpars <- function(dpars_all, k) {
    stats::setNames(
      lapply(comps[[k]]$dpars, function(dp) dpars_all[[paste0(dp, k)]]),
      comps[[k]]$dpars
    )
  }
  # log mixing weights: multinomial logit against the reference
  # component, one list entry per component in component order
  log_pi <- function(dpars_all) {
    Ts <- vector("list", K)
    for (k in theta_ids) Ts[[k]] <- dpars_all[[paste0("theta", k)]]
    Ts[[ref]] <- 0 * dpars_all$mu1
    lse <- Ts[[1L]]
    for (k in seq.int(2L, K)) lse <- RTMB::logspace_add(lse, Ts[[k]])
    lapply(Ts, function(t_) t_ - lse)
  }

  init <- list()
  stat <- list()
  for (k in seq_len(K)) {
    kk <- k
    init[[paste0("mu", k)]] <- local({
      k_ <- kk
      bounded <- mixture_mu_bounded(comps[[k_]])
      function(y, aterms) {
        mixture_mu_start(y, aterms, k_ / (K + 1), bounded)
      }
    })
    st_k <- comps[[k]][["post"]][["stationary"]]
    for (dp in setdiff(comps[[k]]$dpars, "mu")) {
      fn <- comps[[k]]$init_dpars[[dp]]
      if (!is.null(fn)) init[[paste0(dp, k)]] <- fn
      # A component's stationary point is still one of the mixture's,
      # in that component's own dpar, and the mixture renames the dpar
      # by appending the component index. Without this,
      # mixture(skew_normal(), ...) keeps the alpha = 0 stall that the
      # component declares an escape for.
      if (is.list(st_k) && is.list(st_k[[dp]])) {
        stat[[paste0(dp, k)]] <- st_k[[dp]]
      }
    }
  }

  types <- unique(vapply(comps, `[[`, "", "type"))
  # A mixture reads a per-row datum only through its components, and a
  # term reaches the density if ANY of them reads it, so the allow-list
  # is their union. One component that declares nothing (a custom
  # family) leaves the mixture undeclared too: a union with "every
  # term" is every term.
  comp_acc <- lapply(comps, accepted_aterm_names)
  fam <- frmtmb_family(
    paste0("mixture(", paste(vapply(comps, `[[`, "", "family"),
                             collapse = ", "), ")"),
    # `se` is the one exception to the union. A component DOES read it
    # (through resid_sd), but reading it is only half of what se()
    # means: the other half is that the residual scale it replaces is
    # mapped out, and that step names the dpar `sigma`, which a
    # mixture's `sigma1` and `sigma2` are not. Measured before this
    # line existed: mixture(gaussian, gaussian) with se() fitted, both
    # component sigmas sat at their shared starting value 0.9794869
    # having never moved, and every standard error was NaN. Claiming
    # the term here would be claiming a capability the mixture has half
    # of, so it is refused by name instead.
    accepts_aterms = if (!any(vapply(comp_acc, is.null, NA))) {
      setdiff(unique(unlist(comp_acc, use.names = FALSE)), "se")
    },
    # Exclusivity is the INTERSECTION where the allow-list is the union,
    # and the asymmetry is not an oversight: a term reaches the density
    # if any component reads it, but two spellings are interchangeable
    # only if every component treats them so. One component that reads
    # both is a component the second spelling is data for.
    exclusive_aterms = mixture_exclusive_aterms(comps),
    dpars = dpars,
    links = links,
    lpdf = function(y, dpars, aterms) {
      lp <- log_pi(dpars)
      ll <- NULL
      for (k in seq_len(K)) {
        llk <- comps[[k]]$lpdf(y, comp_dpars(dpars, k), aterms) + lp[[k]]
        ll <- if (is.null(ll)) llk else RTMB::logspace_add(ll, llk)
      }
      ll
    },
    valid_y = function(y, aterms) {
      for (cp in comps) {
        if (!is.null(cp$valid_y)) cp$valid_y(y, aterms)
      }
    },
    init_dpars = init,
    type = if (length(types) == 1L) types else "continuous",
    post = c(if (length(stat)) list(stationary = stat), list(
      mean_fn = function(dpars, aterms) {
        lp <- log_pi(dpars)
        out <- 0
        for (k in seq_len(K)) {
          mk <- comps[[k]]$post$mean_fn
          if (is.null(mk)) return(NULL)
          out <- out + exp(lp[[k]]) * mk(comp_dpars(dpars, k), aterms)
        }
        out
      },
      # A mixing weight's RESPONSE scale is the softmax over the
      # component predictors, not the predictor itself. theta's link is
      # identity because that is the scale the multinomial logit inside
      # the density works on, and reporting it as "the response scale"
      # handed back a number that was not a probability and, on this
      # family, was not bounded by one either. The density is untouched:
      # every consumer that feeds dpar values to lpdf(), sim() or
      # mean_fn() reads them through dpars_natural().
      dpar_response = list(
        dpars = paste0("theta", theta_ids),
        value = function(dpars, dnm) {
          exp(log_pi(dpars)[[as.integer(sub("^theta", "", dnm))]])
        },
        deriv = function(dpars, dnm) {
          p <- exp(log_pi(dpars)[[as.integer(sub("^theta", "", dnm))]])
          p * (1 - p)
        }
      )
    )),
    sim = function(dpars, aterms, n) {
      lp <- log_pi(dpars)
      P <- vapply(lp, function(l) rep(exp(l), length.out = n),
                  numeric(n))
      ks <- vapply(seq_len(n), function(i) {
        sample.int(K, 1L, prob = P[i, ])
      }, integer(1))
      out <- numeric(n)
      for (k in seq_len(K)) {
        sk <- comps[[k]]$sim
        if (is.null(sk)) frm_stop("Mixture component ", k, " ('",
                                  comps[[k]]$family,
                                  "') has no simulator", call. = FALSE)
        idx <- which(ks == k)
        if (length(idx)) {
          dk <- lapply(comp_dpars(dpars, k), function(v) {
            rep(v, length.out = n)[idx]
          })
          out[idx] <- sk(dk, aterms, length(idx))
        }
      }
      out
    },
    primary_dpars = paste0("mu", seq_len(K))
  )
  # read by the frame checks brms applies through a mixture's components
  # (trials(), the bernoulli message)
  fam[["component_families"]] <- vapply(comps, `[[`, "", "family")
  # internals for the objective's group-level branch and for
  # mixture_probs(): per-component log-densities and log mixing weights
  fam[["mix"]] <- list(
    K = K,
    comp_lpdf = function(y, dpars, aterms, k) {
      comps[[k]]$lpdf(y, comp_dpars(dpars, k), aterms)
    },
    comp_dpars = comp_dpars,
    comp_sim = function(dpars_k, aterms, n, k) {
      sk <- comps[[k]]$sim
      if (is.null(sk)) {
        frm_stop("Mixture component ", k, " ('", comps[[k]]$family,
                 "') has no simulator for class-wise draws", call. = FALSE)
      }
      sk(dpars_k, aterms, n)
    },
    log_pi = log_pi,
    ref = ref
  )
  # the same components against another reference, for a formula that
  # leaves a different theta unwritten (mixture_theta_reference())
  fam[["mix_rebuild"]] <- function(r) mixture_build(comps, groups, r)
  if (is.null(groups)) {
    # A rowwise mixture still refuses the two fitting options that
    # expand about a single inner mode, and that refusal is a property
    # of the family, so it travels with it rather than sitting in a
    # named gate the core has to run.
    fam[["structure"]] <- frmtmb_structure(
      latent_probs = function(fit, block) mixture_posterior(fit),
      supports = structure_supports_all(reml = FALSE, profile = FALSE),
      refusals = mixture_multimodal_refusals("a mixture() family")
    )
  } else {
    if (!inherits(groups, "formula") || length(groups) != 2L) {
      frm_stop("groups must be a one-sided formula: groups = ~g",
               call. = FALSE)
    }
    fam[["mix_groups"]] <- groups
    # A class belongs to the GROUP, so neither the likelihood nor the
    # draw is rowwise: the group's per-observation log-densities sum
    # BEFORE the logsumexp over classes, and the rowwise `sim` above
    # would give one group two classes. That is what a structure is.
    fam[["structure"]] <- mixture_structure(fam[["mix"]])
  }
  fam
}

#' The two fitting options every mixture-type family refuses, in its own
#' name.
#'
#' A mixture likelihood is invariant to permuting its components, so the
#' location coefficients enter a multimodal objective. Both `REML` and
#' `profile = TRUE` integrate those coefficients out with a Laplace
#' approximation about a single inner mode, which is not defined here:
#' the inner Newton solve walks between the component modes and the fit
#' either dies at "NA/NaN gradient evaluation" or reports an optimum
#' with a gradient near 1e9. Quadrature is unaffected, because it
#' marginalizes the random effects, not the coefficients.
#'
#' `what` is the family as the sentence should name it. These used to be
#' one message listing every mixture-type family, raised from one gate
#' in fit.R that knew all three by name; each family states its own now.
#'
#' @rdname frmtmb-extension-api
#' @export
mixture_multimodal_refusals <- function(what) {
  list(
    reml = paste0(
      "REML = TRUE cannot be combined with ", what, ": the mixture ",
      "likelihood is multimodal in the fixed effects REML integrates ",
      "out, so the restricted likelihood is not defined. Use ",
      "REML = FALSE"),
    profile = paste0(
      "frmtmb_control(profile = TRUE) cannot be combined with ", what,
      ": profiling moves the fixed effects into the inner Laplace ",
      "problem, and the mixture likelihood is multimodal in them. Use ",
      "profile = FALSE")
  )
}

#' The structured-family protocol for a group-level (latent-class)
#' mixture: `mixture(groups = ~g)`.
#'
#' This is the protocol's reference implementation, and it is short on
#' purpose. Only the likelihood, the grouping data it reads and the draw
#' are not rowwise. Everything else about a mixture IS rowwise - its
#' per-row mean is the mixing-weighted mean of the component means, one
#' row at a time - so the four capability flags that depend on that mean
#' are the ones it opts into, `fitted_mean` stays `NULL`, and the core
#' keeps using the family's own `post$mean_fn`.
#'
#' @noRd
mixture_structure <- function(mx) {
  force(mx)
  frmtmb_structure(
    frame_vars = function(fam) list(fam[["mix_groups"]][[2L]]),
    check_spec = mixture_check_spec,
    frame_block = mixture_frame_block,
    unit = "a group-level mixture (mixture(groups = ))",
    loglik = function(y, dpars, aterms, weights, block, extra) {
      lps_pi <- mx[["log_pi"]](dpars)
      total <- NULL
      for (k in seq_len(mx[["K"]])) {
        ll_k <- mx[["comp_lpdf"]](y, dpars, aterms, k)
        g_k <- as.vector(block[["Gt"]] %*% (weights * ll_k)) +
          lps_pi[[k]][block[["first"]]]
        total <- if (is.null(total)) g_k else {
          RTMB::logspace_add(total, g_k)
        }
      }
      sum(total)
    },
    latent_probs = function(fit, block) mixture_posterior(fit),
    sim_ctx = mixture_sim_groups,
    # The four the doc names are the four that follow from the per-row
    # mean being rowwise. Most of the rest are TRUE not because a
    # group-level mixture can do them but because something ELSE
    # already refuses them in words of its own, and a flag here would
    # only shadow the better message: the missing CDF states cens(),
    # trunc() and mi(). `osa` is the one refusal with no better home:
    # the group branch of the objective registers no observation
    # vector, so there is nothing to step through. `cluster_robust` is
    # allowed and then checked, because the sandwich needs every
    # mixture group to sit inside one cluster rather than needing no
    # groups at all.
    # multivariate = FALSE states what mixture_check_spec() enforces
    # upstream (a latent class belongs to the group; a second response
    # would need a joint class process): the flag must not promise what
    # a future consumer of it would then wrongly allow
    supports = structure_supports_all(reml = FALSE, profile = FALSE,
                                      osa = FALSE, multivariate = FALSE),
    refusals = mixture_multimodal_refusals("a mixture() family")
  )
}

#' The model shape a group-level mixture cannot carry: a latent class
#' belongs to the GROUP, and a second response would need a joint class
#' process across responses.
#'
#' `cens()`, `trunc()` and `mi()` are refused too, but not here. A
#' mixture density carries no CDF, so the generic addition-term guard
#' refuses the first two by naming the missing CDF, and `mi()` by naming
#' the families it needs - and those are the messages a user actually
#' meets, because a `check_spec` runs BEFORE those guards and a
#' duplicate refusal here would only shadow the more specific one.
#'
#' @noRd
mixture_check_spec <- function(resp, spec, av) {
  if (length(spec$responses) > 1L || spec$rescor) {
    frm_stop("Group-level mixtures support univariate models",
             call. = FALSE)
  }
  invisible(NULL)
}

#' The group incidence matrix and the row indices a group-level mixture
#' folds its per-observation densities over. All data, resolved once.
#'
#' @noRd
mixture_frame_block <- function(resp, spec, av, mf, y, n) {
  gv <- factor(eval(resp$family[["mix_groups"]][[2L]], mf,
                    resp$formula_env))
  G <- Matrix::sparseMatrix(i = seq_len(n), j = as.integer(gv),
                            x = 1, dims = c(n, nlevels(gv)))
  list(G = G, Gt = Matrix::t(G), first = match(levels(gv), gv),
       gindex = as.integer(gv), levels = levels(gv))
}

#' A group-level mixture draw: one class per group from the group's own
#' mixing weights, then each row from its group's component. The rowwise
#' `sim` cannot express this - it would resample the class per row and
#' wash the grouping out - which is why `groups =` installs this
#' instead.
#'
#' @noRd
mixture_sim_groups <- function(ctx) {
  mg <- ctx[["block"]]
  fam <- ctx[["family"]]
  mx <- fam[["mix"]]
  if (is.null(mg)) {
    frm_stop("mixture(groups =) simulated without the group structure: the ",
             "frame carries no latent-class grouping for response '",
             ctx[["resp"]], "'. Rebuild the model frame from the same ",
             "formula and data", call. = FALSE)
  }
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  av <- ctx[["aterms"]]
  K <- mx[["K"]]
  lps <- mx[["log_pi"]](dp)
  Pg <- vapply(lps, function(l) {
    exp(rep(l, length.out = n)[mg[["first"]]])
  }, numeric(length(mg[["first"]])))
  kg <- vapply(seq_len(nrow(Pg)), function(g_) {
    sample.int(K, 1L, prob = Pg[g_, ])
  }, integer(1))
  kk <- kg[mg[["gindex"]]]
  ys <- numeric(n)
  for (k in seq_len(K)) {
    idx <- which(kk == k)
    if (!length(idx)) next
    dk <- lapply(mx[["comp_dpars"]](dp, k), function(v) {
      rep(v, length.out = n)[idx]
    })
    ys[idx] <- mx[["comp_sim"]](dk, av, length(idx), k)
  }
  ys
}

#' Posterior class probabilities of a mixture fit
#'
#' For an ordinary [mixture()] or [mixture_mvn()] fit, one row per
#' observation; for a group-level mixture (`groups = ~g`), one row per
#' group.
#'
#' @param fit A `frmtmb_fit` with a mixture family.
#' @return A matrix of class probabilities (rows sum to one).
#' @examples
#' set.seed(3)
#' dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)))
#' fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
#'
#' p <- mixture_probs(fit)
#' head(p)
#' rowSums(p)[1:3]                    # rows sum to one
#'
#' # the hard assignment, and how well it recovers the truth
#' cl <- max.col(p)
#' table(cl, truth = rep(1:2, each = 80))
#' @export
mixture_probs <- function(fit) {
  rspec <- single_response(fit, "mixture_probs()")
  if (is.null(rspec$family[["mix"]])) {
    frm_stop("mixture_probs() needs a mixture() family fit", call. = FALSE)
  }
  latent_probs(fit)
}

#' The posterior class probabilities of any family implementing the
#' `fam$mix` component interface, which is `mixture()`, `mixture_mvn()`
#' and `lca()`. This is the structures' `latent_probs` slot for all
#' three, and `mixture_probs()` and `lca_probs()` reach it through the
#' `latent_probs()` generic.
#'
#' @rdname frmtmb-extension-api
#' @export
mixture_posterior <- function(fit) {
  rspec <- single_response(fit, "mixture_probs()")
  fam <- rspec$family
  dp <- eval_dpars(fit)[[rspec$resp_name]]
  av <- fit$frame[["aterm_values"]][[rspec$resp_name]]
  yv <- fit$frame[["y"]][[rspec$resp_name]]
  lps_pi <- fam[["mix"]]$log_pi(dp)
  K <- fam[["mix"]]$K
  # the group incidence of a group-level mixture, or NULL for a rowwise
  # one, where the posterior is per observation
  mg <- frame_block_of(fit$frame, rspec$resp_name)
  w <- av[["weights"]] %||% 1
  # matrix responses (mixture_mvn) have one density per ROW, so NROW,
  # not length; their class covariances live in the extra parameters
  ex <- if (!is.null(fam[["extra_pars"]])) {
    fit$estimates[fit$frame[["extra_names"]]]
  }
  M <- vapply(seq_len(K), function(k) {
    ll_k <- if (is.null(ex)) {
      fam[["mix"]]$comp_lpdf(yv, dp, av, k)
    } else {
      fam[["mix"]]$comp_lpdf(yv, dp, av, k, ex)
    }
    if (!is.null(mg)) {
      as.vector(Matrix::t(mg$G) %*% (w * ll_k)) + lps_pi[[k]][mg$first]
    } else {
      ll_k + rep(lps_pi[[k]], length.out = length(ll_k))
    }
  }, numeric(if (!is.null(mg)) length(mg$first) else NROW(yv)))
  P <- exp(M - apply(M, 1, max))
  P <- P / rowSums(P)
  rownames(P) <- if (!is.null(mg)) mg$levels
  colnames(P) <- paste0("class", seq_len(K))
  P
}

# mclust's covariance taxonomy for mixture_mvn(). mclust writes
# Sigma_k = lambda_k * D_k * A_k * D_k' (volume * orientation * shape);
# each letter of the model name says whether volume, shape and
# orientation are Equal across classes or Vary. The subset here is the
# one whose orientation is either the identity (spherical and diagonal
# models) or completely shared/free, which is exactly the subset a
# log-SD / scaled-Cholesky parameterization expresses without
# constrained eigenvector machinery.
mvn_cov_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
                    "EEE", "VVV")

#' The extras a covariance model needs and how they assemble class k's
#' covariance. `pars` is the (name, length) template used for the error
#' messages and the documentation; `init` turns the response into the
#' starting values; `sigma` runs on the tape, so it must stay AD-safe.
#'
#' @noRd
mvn_cov_spec <- function(model, K, D) {
  if (!is.character(model) || length(model) != 1L ||
        !model %in% mvn_cov_models) {
    frm_stop("mixture_mvn(): unknown covariance model '",
             paste(model, collapse = ", "), "'. Supported models: ",
             paste(mvn_cov_models, collapse = ", "), call. = FALSE)
  }
  us_len <- as.integer(D + D * (D - 1L) / 2L)
  # per-column response SDs; the spherical and volume-shape models
  # collapse them to their log-scale mean. unname(): response column
  # names would otherwise leak into the parameter template and give
  # confint()/frm_sample() ragged parameter labels.
  base_ls <- function(y) unname(log(pmax(apply(y, 2, stats::sd), 1e-3)))
  # diagonal covariance from log-SDs (length D, or length 1 recycled)
  diag_S <- function(ls) {
    "[<-" <- RTMB::ADoverload("[<-")
    S <- RTMB::matrix(0, D, D)
    one <- length(ls) == 1L
    for (j in seq_len(D)) S[j, j] <- exp(2 * ls[if (one) 1L else j])
    S
  }
  # volume-shape diagonal: the shape's last log entry is minus the sum
  # of the free ones, which is what pins det(A) = 1 and keeps the
  # volume identified separately from the shape
  volshape_S <- function(vol, sh) {
    "[<-" <- RTMB::ADoverload("[<-")
    S <- RTMB::matrix(0, D, D)
    a_last <- -sum(sh)
    for (j in seq_len(D)) {
      S[j, j] <- exp(2 * (vol[1] + (if (j < D) sh[j] else a_last)))
    }
    S
  }
  cls <- function(prefix) paste0(prefix, seq_len(K))
  spread <- function(nms, v) {
    stats::setNames(rep(list(v), length(nms)), nms)
  }
  free_shape <- function(y) {
    ls <- base_ls(y)
    (ls - mean(ls))[seq_len(D - 1L)]
  }
  switch(
    model,
    EII = list(
      init = function(y) list(sigmaraw = mean(base_ls(y))),
      sigma = function(extra, k) diag_S(extra[["sigmaraw"]])
    ),
    VII = list(
      init = function(y) spread(cls("sigmaraw"), mean(base_ls(y))),
      sigma = function(extra, k) diag_S(extra[[paste0("sigmaraw", k)]])
    ),
    EEI = list(
      init = function(y) list(sigmaraw = base_ls(y)),
      sigma = function(extra, k) diag_S(extra[["sigmaraw"]])
    ),
    VEI = list(
      init = function(y) {
        c(spread(cls("sigmavol"), mean(base_ls(y))),
          list(sigmashape = free_shape(y)))
      },
      sigma = function(extra, k) {
        volshape_S(extra[[paste0("sigmavol", k)]], extra[["sigmashape"]])
      }
    ),
    EVI = list(
      init = function(y) {
        c(list(sigmavol = mean(base_ls(y))),
          spread(cls("sigmashape"), free_shape(y)))
      },
      sigma = function(extra, k) {
        volshape_S(extra[["sigmavol"]], extra[[paste0("sigmashape", k)]])
      }
    ),
    VVI = list(
      init = function(y) spread(cls("sigmaraw"), base_ls(y)),
      sigma = function(extra, k) diag_S(extra[[paste0("sigmaraw", k)]])
    ),
    EEE = list(
      init = function(y) {
        list(sigmaraw = c(base_ls(y), numeric(us_len - D)))
      },
      sigma = function(extra, k) us_sigma(extra[["sigmaraw"]], D)
    ),
    VVV = list(
      init = function(y) {
        spread(cls("sigmaraw"), c(base_ls(y), numeric(us_len - D)))
      },
      sigma = function(extra, k) us_sigma(extra[[paste0("sigmaraw", k)]], D)
    )
  )
}

#' Multivariate gaussian mixture family
#'
#' `mixture_mvn(K, D)` does model-based clustering of an n x D matrix
#' response (mclust-style): K classes, each with its own D-dimensional
#' mean and a D x D covariance from mclust's model taxonomy. Every
#' class mean is a full
#' linear predictor - the main model formula applies to all of them -
#' so cluster means may depend on covariates, which mclust cannot do.
#' The location dpars are named `mu<k>d<j>` (class k, response column
#' j) and are individually overridable, e.g. `bf(Y ~ x, mu2d1 ~ 1)`
#' (all except the first, `mu1d1`). Mixing weights are `theta1 ...
#' theta{K-1}`, multinomial logit against class K, each with its own
#' linear predictor - so gating on covariates works like [mixture()].
#'
#' Class covariances are family-level extra parameters, covariate-free,
#' and their structure follows `model`, mclust's volume-shape-orientation
#' taxonomy for `Sigma_k = lambda_k * D_k * A_k * D_k'`:
#'
#' \tabular{lll}{
#'   `EII` \tab spherical, equal volume
#'     \tab `sigmaraw`, one log-SD \cr
#'   `VII` \tab spherical, varying volume
#'     \tab `sigmaraw<k>`, one log-SD each \cr
#'   `EEI` \tab diagonal, equal volume and shape
#'     \tab `sigmaraw`, D log-SDs \cr
#'   `VEI` \tab diagonal, varying volume, equal shape
#'     \tab `sigmavol<k>` plus `sigmashape` (D - 1) \cr
#'   `EVI` \tab diagonal, equal volume, varying shape
#'     \tab `sigmavol` plus `sigmashape<k>` (D - 1 each) \cr
#'   `VVI` \tab diagonal, free
#'     \tab `sigmaraw<k>`, D log-SDs each \cr
#'   `EEE` \tab one shared full covariance
#'     \tab `sigmaraw`, a `us` block \cr
#'   `VVV` \tab free full covariance per class (default)
#'     \tab `sigmaraw<k>`, one `us` block each
#' }
#'
#' A `us` block is D log-SDs then the scaled-Cholesky correlation
#' entries, as in [frm()]'s `us()` covariance structure. The
#' `sigmashape` vectors hold the first D - 1 log-shape entries; the last
#' is minus their sum, which is what fixes `det(A_k) = 1`. Log-SDs start
#' at the per-column response SDs and correlations at zero; class means
#' start on spread-out per-column response quantiles to break the label
#' symmetry. The usual finite-mixture ML caveats apply: the
#' likelihood is invariant to
#' relabeling and can be multimodal (compare starts via
#' [frm_allfit()]). [mixture_probs()] returns posterior class
#' probabilities per row; [fitted()] returns the n x D mixture-mean
#' matrix. Covariances take no linear predictor (no covariance
#' regression), and the models with a class-varying eigenvector basis
#' (`EEV`, `VEV`, `EVE`, `VEE`, `VVE`, `EVV`) are not available.
#' `cens()`/`trunc()`, [mvbf()], and `simulate()` are not supported.
#'
#' @param K Number of mixture classes (at least 2).
#' @param D Number of response columns (at least 2; for `D = 1` use
#'   `mixture(gaussian(), ...)`).
#' @param model Covariance model name from mclust's vocabulary: one of
#'   `"EII"`, `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`,
#'   `"VVV"` (the default, a free covariance per class).
#' @return A `frmtmb_family`.
#' @examples
#' set.seed(1)
#' Y <- rbind(matrix(rnorm(60, 0), ncol = 2),
#'            matrix(rnorm(60, 4), ncol = 2))
#' dd <- data.frame(row = seq_len(nrow(Y)))
#' dd$Y <- Y
#' fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
#' fixef(fit)
#' head(mixture_probs(fit))
#' # a shared spherical covariance (mclust's EII, k-means-like)
#' frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = "EII"), data = dd)
#' @export
mixture_mvn <- function(K, D, model = "VVV") {
  if (missing(K) || missing(D) || K < 2 || D < 2) {
    frm_stop("mixture_mvn() needs K >= 2 classes and D >= 2 response ",
             "columns (for D = 1 use mixture(gaussian(), ...))",
             call. = FALSE)
  }
  K <- as.integer(K)
  D <- as.integer(D)
  cspec <- mvn_cov_spec(model, K, D)
  mu_names <- paste0("mu", rep(seq_len(K), each = D),
                     "d", rep(seq_len(D), K))
  dpars <- c(mu_names, paste0("theta", seq_len(K - 1L)))
  links <- stats::setNames(rep(list("identity"), length(dpars)), dpars)

  # n x D class-mean matrix from the class's D location dpars
  class_mean <- function(dpars_all, k, n) {
    "[<-" <- RTMB::ADoverload("[<-")
    M <- RTMB::matrix(0, n, D)
    for (j in seq_len(D)) {
      M[, j] <- dpars_all[[paste0("mu", k, "d", j)]]
    }
    M
  }
  # log mixing weights: multinomial logit, last class reference
  log_pi <- function(dpars_all) {
    Ts <- lapply(seq_len(K - 1L), function(k) {
      dpars_all[[paste0("theta", k)]]
    })
    Ts[[K]] <- 0 * dpars_all[[mu_names[1]]]
    lse <- Ts[[1L]]
    for (k in seq.int(2L, K)) lse <- RTMB::logspace_add(lse, Ts[[k]])
    lapply(Ts, function(t_) t_ - lse)
  }
  # per-row class log-density; extra carries the raw covariance
  # parameters, whose layout the covariance model decides
  comp_lpdf <- function(y, dpars_all, aterms, k, extra) {
    Sk <- cspec$sigma(extra, k)
    R <- y - class_mean(dpars_all, k, nrow(y))
    RTMB::dmvnorm(R, 0, Sk, log = TRUE)
  }

  init <- list()
  for (k in seq_len(K)) {
    for (j in seq_len(D)) {
      init[[paste0("mu", k, "d", j)]] <- local({
        k_ <- k
        j_ <- j
        function(y, aterms) {
          stats::quantile(y[, j_], k_ / (K + 1), names = FALSE)
        }
      })
    }
  }

  fam <- frmtmb_family(
    paste0("mixture_mvn(K = ", K, ", D = ", D, ", model = \"",
           model, "\")"),
    accepts_aterms = "weights",
    dpars = dpars,
    links = links,
    lpdf = function(y, dpars, aterms, extra) {
      lp <- log_pi(dpars)
      ll <- NULL
      for (k in seq_len(K)) {
        llk <- comp_lpdf(y, dpars, aterms, k, extra) + lp[[k]]
        ll <- if (is.null(ll)) llk else RTMB::logspace_add(ll, llk)
      }
      ll
    },
    valid_y = function(y, aterms) {
      if (!is.matrix(y) || ncol(y) != D) {
        frm_stop("mixture_mvn(K = ", K, ", D = ", D, "): response must be ",
                 "an n x ", D, " numeric matrix", call. = FALSE)
      }
    },
    init_dpars = init,
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        lp <- log_pi(dpars)
        n <- length(dpars[[mu_names[1]]])
        out <- 0
        for (k in seq_len(K)) {
          Mk <- vapply(seq_len(D), function(j) {
            dpars[[paste0("mu", k, "d", j)]]
          }, numeric(n))
          out <- out + exp(lp[[k]]) * Mk
        }
        out
      },
      # the same softmax reporting scale mixture() gives its mixing
      # weights, for the same reason
      dpar_response = list(
        dpars = paste0("theta", seq_len(K - 1L)),
        value = function(dpars, dnm) {
          exp(log_pi(dpars)[[as.integer(sub("^theta", "", dnm))]])
        },
        deriv = function(dpars, dnm) {
          p <- exp(log_pi(dpars)[[as.integer(sub("^theta", "", dnm))]])
          p * (1 - p)
        }
      )
    ),
    extra_pars = function(y, aterms) {
      # identical covariance starts across classes; the quantile-spread
      # mean inits break the label symmetry
      cspec$init(y)
    },
    primary_dpars = mu_names
  )
  # internals for mixture_probs(): same shape as mixture()'s, with the
  # extra (covariance) parameters as a fifth comp_lpdf argument
  fam[["mix"]] <- list(
    K = K,
    D = D,
    cov_model = model,
    comp_lpdf = comp_lpdf,
    comp_dpars = function(dpars_all, k) {
      stats::setNames(
        lapply(seq_len(D), function(j) {
          dpars_all[[paste0("mu", k, "d", j)]]
        }),
        paste0("mud", seq_len(D))
      )
    },
    # class k's covariance from the fit's extra parameters; the same
    # assembler the objective uses, so it also runs off the tape
    sigma = cspec$sigma,
    log_pi = log_pi
  )
  # the class covariances are family-level EXTRAS, not dpars, so the
  # rowwise contract cannot see them; the structured one can
  fam[["sim_ctx"]] <- mvn_sim_rows
  # the likelihood is rowwise (one D-variate density per row), so the
  # structure carries nothing but the two multimodality refusals
  fam[["structure"]] <- frmtmb_structure(
    latent_probs = function(fit, block) mixture_posterior(fit),
    supports = structure_supports_all(reml = FALSE, profile = FALSE),
    refusals = mixture_multimodal_refusals("a mixture_mvn() family")
  )
  fam
}

#' One `mixture_mvn()` draw: a class per row from its mixing weights,
#' then a D-variate normal about that class's mean with that class's
#' covariance, assembled from the fit's extra parameters by the same
#' `sigma()` the likelihood uses.
#'
#' @noRd
mvn_sim_rows <- function(ctx) {
  mx <- ctx[["family"]][["mix"]]
  K <- mx[["K"]]
  D <- mx[["D"]]
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  ex <- ctx[["extra"]]
  P <- vapply(mx[["log_pi"]](dp), function(l) {
    rep(exp(l), length.out = n)
  }, numeric(n))
  ks <- vapply(seq_len(n), function(i) {
    sample.int(K, 1L, prob = P[i, ])
  }, integer(1))
  out <- matrix(NA_real_, n, D)
  for (k in seq_len(K)) {
    idx <- which(ks == k)
    if (!length(idx)) next
    # upper Cholesky: Z %*% R has covariance R'R = Sigma_k
    R <- chol(as.matrix(mx[["sigma"]](ex, k)))
    M <- matrix(NA_real_, length(idx), D)
    for (j in seq_len(D)) {
      M[, j] <- rep(dp[[paste0("mu", k, "d", j)]], length.out = n)[idx]
    }
    Z <- matrix(stats::rnorm(length(idx) * D), length(idx), D)
    out[idx, ] <- M + Z %*% R
  }
  out
}

#' von Mises family for a circular response in `(-pi, pi]`, dpars `mu`
#' (the mean direction, through the tan-half link) and `kappa` (the
#' concentration, log link). brms's parameterization exactly.
#'
#' The density needs `log I0(kappa)`, the modified Bessel function of the
#' first kind of order zero. `base::besselI()` is not AD-safe, but RTMB
#' carries its own `besselI` method for advectors, and `RTMBdist::dvm()`
#' already uses it in the exponentially scaled form
#' (`log I0(k) = log besselI(k, 0, expon.scaled = TRUE) + k`), which is
#' what keeps a large concentration from overflowing. So the whole
#' density, and its exact derivative, come off the tape with no series
#' approximation of our own.
#'
#' @noRd
fam_von_mises <- function(link = "tan_half", link_kappa = "log") {
  lk_kappa <- dpar_link(link_kappa, "kappa", "von_mises", dpar_links_positive)
  frmtmb_family(
    "von_mises",
    accepts_aterms = "weights",
    dpars = c("mu", "kappa"),
    links = list(mu = mu_link(link, "von_mises"), kappa = lk_kappa),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dvm(y, dpars[["mu"]], dpars[["kappa"]], log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y < -pi) || any(y > pi)) {
        frm_stop("von_mises: response must be angles in radians on ",
                 "(-pi, pi]; wrap the response first, for example ",
                 "atan2(sin(y), cos(y))", call. = FALSE)
      }
    },
    init_dpars = list(
      # the resultant vector's direction and length: the circular
      # analogues of the mean and (through Fisher's approximation of
      # A(kappa) = I1/I0) the precision
      mu = function(y, aterms) atan2(mean(sin(y)), mean(cos(y))),
      kappa = function(y, aterms) {
        rbar <- sqrt(mean(sin(y))^2 + mean(cos(y))^2)
        rbar <- min(max(rbar, 0.02), 0.99)
        if (rbar < 0.53) {
          2 * rbar + rbar^3 + 5 * rbar^5 / 6
        } else if (rbar < 0.85) {
          -0.4 + 1.39 * rbar + 0.43 / (1 - rbar)
        } else {
          1 / (rbar^3 - 4 * rbar^2 + 3 * rbar)
        }
      }
    ),
    type = "continuous",
    post = list(
      # the mean DIRECTION, which is what brms's posterior_epred()
      # reports for this family; a circular response has no mean in the
      # arithmetic sense
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      # the circular variance 1 - A(kappa), on (0, 1)
      var_fn = function(dpars, aterms) {
        k <- dpars[["kappa"]]
        1 - besselI(k, 1, expon.scaled = TRUE) /
          besselI(k, 0, expon.scaled = TRUE)
      }
    ),
    sim = function(dpars, aterms, n) {
      rvon_mises(n, rep(dpars[["mu"]], length.out = n),
                 rep(dpars[["kappa"]], length.out = n))
    }
  )
}

#' von Mises draws with a per-row mean direction and concentration, by
#' Best and Fisher's (1979) wrapped-Cauchy rejection scheme. `RTMBdist`
#' ships `rvm()`, but it takes scalar parameters only, and a
#' distributional fit (or `dharma_residuals()`, which draws hundreds of
#' replicates) needs one draw per row of a varying `kappa`. Rejection
#' runs vectorized over the rows still outstanding.
#'
#' @noRd
rvon_mises <- function(n, mu, kappa) {
  mu <- rep(mu, length.out = n)
  kappa <- rep(kappa, length.out = n)
  out <- numeric(n)
  # kappa near zero is the uniform distribution on the circle, and the
  # rejection constants below divide by it
  flat <- kappa < 1e-8
  if (any(flat)) out[flat] <- stats::runif(sum(flat), -pi, pi)
  todo <- which(!flat)
  if (!length(todo)) return(out)
  k <- kappa[todo]
  a <- 1 + sqrt(1 + 4 * k^2)
  b <- (a - sqrt(2 * a)) / (2 * k)
  r <- (1 + b^2) / (2 * b)
  left <- seq_along(todo)
  it <- 0L
  while (length(left)) {
    it <- it + 1L
    if (it > 1000L) {
      frm_stop("von Mises simulation did not converge for ", length(left),
               " of ", n, " rows", call. = FALSE)
    }
    m <- length(left)
    z <- cos(pi * stats::runif(m))
    f <- (1 + r[left] * z) / (r[left] + z)
    c_ <- k[left] * (r[left] - f)
    u2 <- stats::runif(m)
    ok <- c_ * (2 - c_) - u2 > 0 | log(c_ / u2) + 1 - c_ >= 0
    if (any(ok)) {
      idx <- left[ok]
      sgn <- sign(stats::runif(sum(ok)) - 0.5)
      out[todo[idx]] <- mu[todo[idx]] + sgn * acos(pmin(pmax(f[ok], -1), 1))
    }
    left <- left[!ok]
  }
  # back onto (-pi, pi], the support the family declares
  (out + pi) %% (2 * pi) - pi
}

# --- categorical (nominal) responses ---------------------------------

#' Multinomial-logit log-density over a CATEGORY INDEX response `1..K`.
#' `eta` is the list of K-1 non-reference linear predictors in category
#' order, category 1 being the reference with `eta_1 = 0`; the shared
#' code is factored out so the taped path (`ord_cat_sel`, which needs no
#' comparison operators) and the data path agree term for term.
#'
#' @noRd
cat_logit_lpdf <- function(y, eta, K, sel = NULL) {
  if (is.null(sel)) sel <- lapply(seq_len(K), function(k) {
    as.numeric(y == k)
  })
  # the denominator accumulates in LOG space from the reference
  # category's implicit zero, so a wide predictor cannot overflow it;
  # summing exp(eta_j) first loses the model at eta = 709
  lden <- 0 * eta[[1L]]
  num <- 0
  for (j in seq_len(K - 1L)) {
    lden <- RTMB::logspace_add(lden, eta[[j]])
    num <- num + sel[[j + 1L]] * eta[[j]]
  }
  num - lden
}

#' Categorical (nominal) family over a `1..K` category index, category 1
#' the reference. One linear predictor per non-reference category, all
#' receiving the main model formula unless overridden.
#'
#' `dpar_names` is the brms spelling `mu<Level>` when the levels are
#' known and the positional `mu2 ... muK` when only `K` is.
#'
#' @noRd
fam_categorical_impl <- function(dpar_names, levels = NULL,
                                 link = "logit") {
  categorical_link_check(link)
  K <- length(dpar_names) + 1L
  fam <- frmtmb_family(
    "categorical",
    accepts_aterms = "weights",
    dpars = dpar_names,
    links = stats::setNames(rep(list("identity"), K - 1L), dpar_names),
    lpdf = function(y, dpars, aterms) {
      eta <- lapply(dpar_names, function(nm) dpars[[nm]])
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        # the response is on the tape: pick the category arithmetically
        return(cat_logit_lpdf(NULL, eta, K,
                              sel = ord_cat_sel(ov$y, K)) * ov$keep)
      }
      cat_logit_lpdf(y, eta, K)
    },
    valid_y = function(y, aterms) {
      if (any(y < 1) || any(y > K) || any(y != round(y))) {
        frm_stop("categorical: response must be a factor with ", K,
                 " levels, or integer category codes 1..", K, call. = FALSE)
      }
      if (length(unique(y)) < 2) {
        frm_stop("categorical: the response takes only one value; there is ",
                 "nothing to model", call. = FALSE)
      }
    },
    type = "categorical",
    sim = function(dpars, aterms, n) {
      eta <- lapply(dpar_names, function(nm) {
        rep(dpars[[nm]], length.out = n)
      })
      E <- cbind(0, do.call(cbind, eta))
      P <- exp(E - apply(E, 1L, max))
      P <- P / rowSums(P)
      cp <- t(apply(P, 1L, cumsum))
      if (n == 1L) cp <- matrix(cp, 1L, K)
      pmin(1L + rowSums(cp < stats::runif(n)), K)
    },
    primary_dpars = dpar_names
  )
  fam[["cat_levels"]] <- levels
  fam[["cat_K"]] <- K
  fam
}

#' Refuse any categorical link but the logit, in brms's words.
#'
#' @noRd
categorical_link_check <- function(link) {
  if (identical(link, "logit")) return(invisible(NULL))
  shown <- if (is.character(link) && length(link) == 1L) link else {
    arg_desc(link)
  }
  frm_stop("'", shown, "' is not a supported link for family ",
           "'categorical'. Supported links are: 'logit'. categorical() ",
           "takes the 'logit' link only, the multinomial logit its ",
           "category probabilities are defined by", call. = FALSE)
}

#' Placeholder returned by a bare `categorical()`: the category count is
#' a property of the data, and the formula grammar is resolved before any
#' data arrives, so `frm()` swaps this for the real family through
#' `resolve_deferred_families()` once it can see the response. It carries
#' one `mu` so a spec built from it is still well formed; anything that
#' reaches the likelihood with it in place says so.
#'
#' @noRd
fam_categorical_deferred <- function(link = "logit") {
  # checked here and not only once the levels are known: the levels come
  # from the data, and a data problem would otherwise mask the bad link
  categorical_link_check(link)
  fam <- frmtmb_family(
    "categorical",
    accepts_aterms = "weights",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms) {
      frm_stop("categorical(): the response categories were never resolved. ",
               "They are read from the data by frm(); on another entry point ",
               "name them, categorical(levels = c(\"a\", \"b\", \"c\")) or ",
               "categorical(K = 3)", call. = FALSE)
    },
    type = "categorical"
  )
  fam[["defer"]] <- function(formula, data) {
    lv <- categorical_levels(formula, data)
    fam_categorical_impl(paste0("mu", lv[-1L]), levels = lv, link = link)
  }
  fam
}

#' The response's category levels, for the deferred `categorical()`. The
#' response expression is evaluated against the data with the addition
#' terms stripped; a character response becomes a factor here, and says
#' which order it took.
#'
#' @noRd
categorical_levels <- function(formula, data) {
  ri <- parse_response(formula)
  y <- tryCatch(eval(ri$resp, data, environment(formula) %||% globalenv()),
                error = function(e) NULL)
  if (is.null(y)) {
    frm_stop("categorical(): the response '", deparse1(ri$resp),
             "' could not be evaluated on the data to find its categories. ",
             "Name them instead: categorical(levels = c(\"a\", \"b\"))",
             call. = FALSE)
  }
  lv <- categorical_y_levels(y, deparse1(ri$resp))
  if (length(lv) < 2L) {
    frm_stop("categorical(): the response '", deparse1(ri$resp),
             "' has fewer than two categories", call. = FALSE)
  }
  lv
}

#' The category labels of a categorical response, in the order that fixes
#' the reference category and the dpar names. A character vector is
#' coerced with a message naming the order it took, because that order is
#' the model; a factor keeps its own levels; numeric codes are read in
#' numeric order, as brms reads them.
#'
#' @noRd
categorical_y_levels <- function(y, label) {
  if (is.character(y) || is.logical(y)) {
    lv <- levels(factor(y))
    frm_message("Categorical response '", label, "' is a ",
                if (is.logical(y)) "logical" else "character",
                " vector; it is read as a factor with levels ",
                paste(lv, collapse = ", "), " and '", lv[1L],
                "' as the reference category. Set the order with factor() ",
                "if that is not what you want.")
    return(lv)
  }
  if (is.factor(y)) return(levels(y))
  # Numeric codes are categories too, in numeric order, as brms reads
  # them: 11:20 is ten categories with 11 the reference. This returned
  # NULL once, and the caller then claimed "fewer than two categories"
  # of a response that had ten.
  if (is.numeric(y) && !is.matrix(y)) return(levels(factor(y)))
  NULL
}

#' Swap every deferred family of a bform for the concrete one its data
#' implies. Called by `frm()` before parsing, so a per-category dpar
#' formula (`bf(y ~ x, mub ~ z)`) reaches the parser with the dpar
#' already in the family's vocabulary.
#'
#' @noRd
resolve_deferred_families <- function(bform, data) {
  if (is.null(data)) return(bform)
  one <- function(f) {
    if (is.null(f$family) || is.null(f$family[["defer"]])) return(f)
    f$family <- f$family[["defer"]](f$formula, data)
    f
  }
  if (inherits(bform, "frmtmb_mvformula")) {
    bform$forms <- lapply(bform$forms, one)
  } else {
    bform <- one(bform)
  }
  bform
}

#' Carry a frame's finalized responses back into the spec beside it.
#'
#' A family with a `family_finalize` slot is not itself until the
#' response has been seen: assembly is where it derives its links and
#' any slot it computes from the data. Every stage that keeps the spec
#' separately from the frame has to take the finalized copy, or it
#' reads the family as written rather than the one the likelihood
#' scores. Only responses whose family declares the slot are replaced;
#' the rest of the spec stays as parsed, which is what every stage has
#' always read. A bernoulli response's coding of its two values
#' (`bin_levels`) is carried the same way, because refits and newdata
#' are coded against it.
#'
#' @noRd
carry_finalized_responses <- function(spec, frame) {
  for (rn_ in names(spec$responses)) {
    fr_ <- frame[["spec"]]$responses[[rn_]]
    if (is.null(spec$responses[[rn_]]$family[["family_finalize"]])) {
      if (!is.null(fr_$family[["bin_levels"]])) {
        spec$responses[[rn_]]$family[["bin_levels"]] <-
          fr_$family[["bin_levels"]]
      }
      next
    }
    spec$responses[[rn_]] <- fr_
  }
  spec
}

# --- Cox proportional hazards ----------------------------------------

#' B-spline basis of a given order (degree + 1) over a full knot vector,
#' by the Cox-de Boor recursion. Written out rather than taken from
#' `splines::splineDesign()` so the baseline hazard needs no dependency
#' outside the ones the package already carries; it agrees with
#' `splines2` to machine precision, which the tests assert.
#'
#' @noRd
bspline_basis <- function(x, knots, order) {
  n <- length(x)
  nb1 <- length(knots) - 1L
  # order 1: the indicator of [t_j, t_{j+1}); the last non-degenerate
  # interval closes on the right so the upper boundary is covered
  B <- matrix(0, n, nb1)
  last <- max(which(knots < knots[length(knots)]))
  for (j in seq_len(nb1)) {
    B[, j] <- if (j == last) {
      as.numeric(x >= knots[j] & x <= knots[j + 1L])
    } else {
      as.numeric(x >= knots[j] & x < knots[j + 1L])
    }
  }
  if (order == 1L) return(B[, seq_len(length(knots) - order), drop = FALSE])
  for (k in 2:order) {
    Bn <- matrix(0, n, length(knots) - k)
    for (j in seq_len(ncol(Bn))) {
      d1 <- knots[j + k - 1L] - knots[j]
      d2 <- knots[j + k] - knots[j + 1L]
      # a repeated knot gives a zero-width span; its term drops out
      t1 <- if (d1 > 0) (x - knots[j]) / d1 * B[, j] else 0
      t2 <- if (d2 > 0) (knots[j + k] - x) / d2 * B[, j + 1L] else 0
      Bn[, j] <- t1 + t2
    }
    B <- Bn
  }
  B
}

#' The knot placement brms uses for the Cox baseline: internal knots on
#' equally spaced response quantiles, boundary knots just outside the
#' observed range so every observation sits strictly inside the basis
#' (`splines2`'s default through brms's `bhaz_basis_matrix()`).
#'
#' @noRd
bhaz_spec <- function(y, df, degree, intercept) {
  n_int <- df - degree - as.integer(intercept)
  if (n_int < 0L) {
    frm_stop("cox(): df must be at least degree + 1 (", degree + 1L,
             ") with an intercept in the baseline basis", call. = FALSE)
  }
  rng <- range(y)
  d <- rng[2L] - rng[1L]
  boundary <- c(max(rng[1L] - d / 50, 0), rng[2L] + d / 50)
  internal <- if (n_int > 0L) {
    p <- seq(0, 1, length.out = n_int + 2L)[-c(1L, n_int + 2L)]
    unname(stats::quantile(y, probs = p))
  } else {
    numeric(0)
  }
  list(internal = internal, boundary = boundary, degree = degree,
       intercept = intercept, df = df)
}

#' The full knot vector of a baseline-hazard spec: boundary knots
#' repeated to the spline's order at each end.
#'
#' @noRd
bhaz_knots <- function(sp, extra = 0L) {
  k <- sp$degree + 1L + extra
  c(rep(sp$boundary[1L], k), sp$internal, rep(sp$boundary[2L], k))
}

#' The M-spline basis of the baseline hazard: the B-spline basis scaled
#' so each function integrates to one over its support, which is what
#' makes a simplex of coefficients a density-like hazard.
#'
#' @noRd
mspline_design <- function(x, sp) {
  k <- sp$degree + 1L
  tt <- bhaz_knots(sp)
  N <- bspline_basis(x, tt, k)
  M <- N
  for (j in seq_len(ncol(N))) M[, j] <- N[, j] * k / (tt[j + k] - tt[j])
  if (!sp$intercept) M <- M[, -1L, drop = FALSE]
  M
}

#' The I-spline basis: `I_j(x)` is the integral of `M_j` from the lower
#' boundary knot, so it is the CUMULATIVE baseline hazard's basis and it
#' is monotone by construction. It equals the reverse cumulative sum of
#' the order-(degree + 2) B-spline basis on the once-more-repeated
#' boundary knots (Ramsay 1988); the tests check it against both
#' `splines2::iSpline()` and a numeric integral of `mspline_design()`.
#'
#' @noRd
ispline_design <- function(x, sp) {
  N2 <- bspline_basis(x, bhaz_knots(sp, extra = 1L), sp$degree + 2L)
  R <- t(apply(N2, 1L, function(v) rev(cumsum(rev(v)))))
  if (length(x) == 1L) R <- matrix(R, 1L, ncol(N2))
  I <- R[, -1L, drop = FALSE]
  if (!sp$intercept) I <- I[, -1L, drop = FALSE]
  I
}

#' Cox proportional hazards with a flexible (M-spline) baseline hazard,
#' brms's construction. The baseline hazard is `Zbhaz %*% s` and the
#' cumulative baseline hazard `Zcbhaz %*% s` over the SAME simplex `s`,
#' with `Zcbhaz` the I-spline integral of the M-spline `Zbhaz`. The
#' simplex is what identifies the baseline against the intercept, and it
#' rides in the parameter vector as `K - 1` unconstrained values through
#' a softmax with the first component pinned at zero.
#'
#' @noRd
fam_cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE) {
  df <- as.integer(df)
  degree <- as.integer(degree)
  intercept <- isTRUE(intercept)
  # the simplex from its K-1 free coordinates; softmax with the first
  # component pinned at zero covers the whole open simplex
  sbhaz <- function(raw) {
    "c" <- RTMB::ADoverload("c")   # base c() strips the advector class
    e <- exp(c(0, raw))
    e / sum(e)
  }
  fam <- frmtmb_family(
    "cox",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = "mu",
    links = list(mu = mu_link(link, "cox")),
    lpdf = function(y, dpars, aterms, extra) {
      s <- sbhaz(extra$sbhaz_raw)
      bhaz <- as.vector(aterms[["Zbhaz"]] %*% s)
      cbhaz <- as.vector(aterms[["Zcbhaz"]] %*% s)
      # log h(t) + log S(t), with h(t) = h0(t) * mu and mu = exp(eta)
      log(bhaz) + log(dpars[["mu"]]) - cbhaz * dpars[["mu"]]
    },
    lcdf = function(q, dpars, aterms, extra) {
      # F(t) = 1 - exp(-H0(t) * mu). Every bound a censored or truncated
      # row asks about is plain data, so its I-spline row block was
      # built once at frame assembly and is looked up by value here.
      cbhaz <- as.vector(cox_cbhaz_design(q, aterms) %*%
                           sbhaz(extra$sbhaz_raw))
      1 - exp(-cbhaz * dpars[["mu"]])
    },
    lccdf = function(q, dpars, aterms, extra) {
      # log S = -H0(t) * mu, the cumulative hazard itself
      cbhaz <- as.vector(cox_cbhaz_design(q, aterms) %*%
                           sbhaz(extra$sbhaz_raw))
      -cbhaz * dpars[["mu"]]
    },
    valid_y = function(y, aterms) {
      if (any(y <= 0)) {
        frm_stop("cox: the response is a survival time and must be ",
                 "strictly positive", call. = FALSE)
      }
      if (length(unique(y)) < df) {
        frm_stop("cox: fewer distinct event times (", length(unique(y)),
                 ") than baseline basis functions (", df,
                 "); lower df in cox(df = )", call. = FALSE)
      }
    },
    init_dpars = list(
      # the baseline integrates to one over the observed window, so the
      # hazard ratio starts at the crude event rate
      mu = function(y, aterms) {
        ev <- if (is.null(aterms[["cens"]])) 1 else mean(aterms[["cens"]] == 0)
        max(ev, 0.05) / mean(y)
      }
    ),
    type = "continuous",
    post = list(
      # mu is the hazard RATIO, not a mean, and reporting it as one is
      # exactly the mistake this family invites. brms refuses the same
      # question (posterior_epred has no cox method).
      mean_fn = function(dpars, aterms) {
        frm_stop("cox: a survival time has no mean on the response scale ",
                 "here - the model fits a hazard, and the mean survival ",
                 "time would be an integral over the baseline that the ",
                 "censored rows do not identify. frm_linpred(type = \"link\") ",
                 "gives the log hazard ratio and cox_baseline() the fitted ",
                 "baseline weights", call. = FALSE)
      }
    ),
    extra_pars = function(y, aterms) {
      list(sbhaz_raw = rep(0, ncol(aterms[["Zbhaz"]]) - 1L))
    },
    # a deliberate omission, not a gap: every entry point repeats this
    # reason after its own refusal
    sim_refusal = paste0(
      "Drawing a survival time means inverting the cumulative baseline ",
      "hazard, and cox() carries no quantile function for it: the ",
      "I-spline baseline is identified only on the observed time ",
      "window, so a draw beyond the last event time has no defined ",
      "distribution"
    )
  )
  # the baseline bases are data, not parameters: they are built once
  # from the observed times and ride with the response's addition terms
  fam[["aterm_data"]] <- function(y, aterms) {
    sp <- bhaz_spec(y, df, degree, intercept)
    out <- list(Zbhaz = mspline_design(y, sp), bhaz_spec = sp,
                cox_cb = list(cox_cb_entry(y, sp)))
    for (nm in c("cens_y2", "trunc_lb", "trunc_ub")) {
      if (!is.null(aterms[[nm]])) {
        out$cox_cb[[length(out$cox_cb) + 1L]] <-
          cox_cb_entry(aterms[[nm]], sp)
      }
    }
    out$Zcbhaz <- out$cox_cb[[1L]]$Z
    out
  }
  fam[["cox_sbhaz"]] <- sbhaz
  fam
}

#' One cached (bound, I-spline design) pair. A bound outside the
#' boundary knots is clamped first: below the lower knot the cumulative
#' baseline hazard is zero and above the upper one it is its total, and
#' the raw B-spline recursion returns an all-zero row for both, which
#' would read as "no risk accumulated" at an upper bound.
#'
#' @noRd
cox_cb_entry <- function(q, sp) {
  qq <- pmin(pmax(as.numeric(q), sp$boundary[1L]), sp$boundary[2L])
  list(q = as.numeric(q), Z = ispline_design(qq, sp))
}

#' The cumulative-baseline-hazard design for a bound the Cox likelihood
#' asks about. Frame assembly cached one per bound; anything else (a
#' post-fit query) is built on the spot.
#'
#' @noRd
cox_cbhaz_design <- function(q, aterms) {
  qn <- as.numeric(q)
  for (e in aterms[["cox_cb"]] %||% list()) {
    if (identical(e$q, qn)) return(e$Z)
  }
  cox_cb_entry(qn, aterms[["bhaz_spec"]])$Z
}

#' Call a family's log-CDF, passing the family-level extra parameters to
#' the families that declare them. The Cox baseline lives there, so its
#' survivor function needs them; every other CDF keeps the three-argument
#' contract.
#'
#' @noRd
fam_lcdf <- function(fam, q, dpars, aterms, extra) {
  if (length(formals(fam[["lcdf"]])) >= 4L) {
    fam[["lcdf"]](q, dpars, aterms, extra)
  } else {
    fam[["lcdf"]](q, dpars, aterms)
  }
}

#' Call a family's log survivor function, with the same extra-parameter
#' arity shim `fam_lcdf()` carries.
#'
#' @noRd
fam_lccdf <- function(fam, q, dpars, aterms, extra) {
  if (length(formals(fam[["lccdf"]])) >= 4L) {
    fam[["lccdf"]](q, dpars, aterms, extra)
  } else {
    fam[["lccdf"]](q, dpars, aterms)
  }
}

#' Whether a family can score a RIGHT-censored row exactly, which it can
#' when it declares `lccdf`.
#'
#' @noRd
has_lccdf <- function(fam) !is.null(fam[["lccdf"]])

#' The fitted baseline-hazard simplex of a `cox()` fit.
#'
#' @param fit A `frmtmb_fit` with a [cox()] family.
#' @return The `Kbhaz` M-spline weights, summing to one.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200))
#' dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
#' fit <- frm(bf(time ~ x), family = cox(), data = dd)
#' cox_baseline(fit)
#' @export
cox_baseline <- function(fit) {
  rspec <- single_response(fit, "cox_baseline()")
  fam <- rspec$family
  if (is.null(fam[["cox_sbhaz"]])) {
    frm_stop("cox_baseline() needs a cox() family fit", call. = FALSE)
  }
  s <- fam[["cox_sbhaz"]](fit$estimates[["sbhaz_raw"]])
  stats::setNames(s, paste0("s", seq_along(s)))
}

#' Matrix-response multinomial: y is an n x K count matrix, category 1
#' is the reference. One linear predictor per non-reference category
#' (`mu2`, ..., `muK`), all receiving the main model formula unless
#' overridden.
#'
#' @noRd
fam_multinomial <- function(K) {
  if (missing(K) || K < 2) {
    frm_stop("multinomial() needs the number of categories, e.g. ",
             "multinomial(K = 3)", call. = FALSE)
  }
  dpn <- paste0("mu", seq_len(K)[-1])
  frmtmb_family(
    "multinomial",
    accepts_aterms = c("weights", "trials"),
    dpars = dpn,
    links = stats::setNames(rep(list("identity"), K - 1L), dpn),
    lpdf = function(y, dpars, aterms) {
      # log-space denominator, as in cat_logit_lpdf(): a count response
      # multiplies it by the row total, so an overflow there is fatal
      lden <- 0 * dpars[[dpn[1L]]]
      for (k in dpn) lden <- RTMB::logspace_add(lden, dpars[[k]])
      nvec <- rowSums(y)
      ll <- -nvec * lden
      for (j in seq_along(dpn)) {
        ll <- ll + y[, j + 1L] * dpars[[dpn[j]]]
      }
      # multinomial coefficient: data-only
      ll + lgamma(nvec + 1) - rowSums(lgamma(y + 1))
    },
    valid_y = function(y, aterms) {
      if (!is.matrix(y) || ncol(y) != K) {
        frm_stop("multinomial(K = ", K, "): response must be an n x ", K,
                 " count matrix", call. = FALSE)
      }
      if (any(y < 0) || any(y != round(y))) {
        frm_stop("multinomial: response must be non-negative integer counts",
                 call. = FALSE)
      }
      tr <- aterms[["trials"]]
      if (!is.null(tr) && any(rowSums(y) != tr)) {
        frm_stop("Number of trials does not match the number of events. ",
                 "multinomial: rowSums(response) must equal trials()",
                 call. = FALSE)
      }
    },
    type = "discrete",
    sim = function(dpars, aterms, n) {
      size <- aterms[["trials"]]
      if (is.null(size)) {
        frm_stop("simulate(): a multinomial fit needs trials() to know how ",
                 "many draws each row gets", call. = FALSE)
      }
      denom <- 1
      for (k in dpn) denom <- denom + exp(dpars[[k]])
      P <- matrix(0, n, K)
      P[, 1L] <- rep(1 / denom, length.out = n)
      for (j in seq_along(dpn)) {
        P[, j + 1L] <- rep(exp(dpars[[dpn[j]]]) / denom, length.out = n)
      }
      size <- rep(size, length.out = n)
      out <- matrix(0L, n, K)
      for (i in seq_len(n)) {
        out[i, ] <- stats::rmultinom(1L, size[i], P[i, ])
      }
      out
    },
    primary_dpars = dpn
  )
}

family_registry <- list(
  gaussian    = fam_gaussian,
  poisson     = fam_poisson,
  binomial    = fam_binomial,
  Gamma       = fam_Gamma,
  lognormal   = fam_lognormal,
  student     = fam_student,
  negbinomial = fam_negbinomial,
  nbinom2     = fam_negbinomial,
  nbinom1     = fam_nbinom1,
  beta        = fam_beta,
  Beta        = fam_beta,
  tweedie     = fam_tweedie,
  compois     = fam_compois,
  zero_inflated_poisson     = fam_zi_poisson,
  zero_inflated_negbinomial = fam_zi_negbinomial,
  hurdle_poisson            = fam_hurdle_poisson,
  hurdle_negbinomial        = fam_hurdle_negbinomial,
  multinomial               = fam_multinomial,
  cumulative                = fam_cumulative,
  hurdle_cumulative         = fam_hurdle_cumulative,
  beta_binomial             = fam_beta_binomial,
  zero_inflated_beta_binomial = fam_zi_beta_binomial,
  skew_normal               = fam_skew_normal,
  inverse.gaussian          = fam_inverse_gaussian,
  exgaussian                = fam_exgaussian,
  bernoulli                 = fam_bernoulli,
  geometric                 = fam_geometric,
  exponential               = fam_exponential,
  weibull                   = fam_weibull,
  shifted_lognormal         = fam_shifted_lognormal,
  hurdle_gamma              = fam_hurdle_gamma,
  hurdle_lognormal          = fam_hurdle_lognormal,
  zero_inflated_binomial    = fam_zi_binomial,
  zero_inflated_beta        = fam_zi_beta,
  zero_one_inflated_beta    = fam_zoi_beta,
  xbeta                     = fam_xbeta,
  asym_laplace              = fam_asym_laplace,
  zero_inflated_asym_laplace = fam_zi_asym_laplace,
  huber                     = fam_huber,
  sratio                    = fam_sratio,
  cratio                    = fam_cratio,
  acat                      = fam_acat,
  von_mises                 = fam_von_mises,
  cox                       = fam_cox,
  categorical               = fam_categorical_deferred
)

# An unknown value errors and names the supported families.

#' Find a family constructor by name, under brms's spelling rules.
#'
#' brms lower-cases the name and reads `normal` as `gaussian`, `zi_` as
#' `zero_inflated_` and `hu_` as `hurdle_`, so `brmsfamily("zi_poisson")`
#' and `family = "Gamma"` both work there. The exact registry name is
#' tried first, so frmtmb's own spellings keep resolving as they did.
#'
#' @noRd
family_ctor <- function(name) {
  ctor <- family_registry[[name]]
  if (!is.null(ctor)) return(ctor)
  low <- tolower(name)
  low <- sub("^normal$", "gaussian", low)
  low <- sub("^zi_", "zero_inflated_", low)
  low <- sub("^hu_", "hurdle_", low)
  # brms's name for the family frmtmb calls compois
  low <- sub("^com_poisson$", "compois", low)
  keys <- names(family_registry)
  hit <- keys[tolower(keys) == low]
  if (length(hit)) return(family_registry[[hit[1L]]])
  frm_stop(name, " is not a supported family. Supported families are: ",
           paste(unique(names(family_registry)), collapse = ", "),
           call. = FALSE)
}

#' Build a family from its name and links, the shared body of
#' [brmsfamily()] and a character `family =`.
#'
#' `link` is already a value here: the callers that take an unquoted link
#' have resolved it. `NULL` leaves the constructor's default.
#'
#' @noRd
family_from_name <- function(family, link = NULL, args = list(),
                             caller = "brmsfamily") {
  ctor <- family_ctor(family)
  if (length(args) &&
      (is.null(names(args)) || !all(nzchar(names(args))))) {
    frm_stop(caller, "(): every argument after `family` and `link` has to ",
             "be named, as in link_sigma = \"softplus\"", call. = FALSE)
  }
  unknown <- setdiff(names(args), setdiff(names(formals(ctor)), "link"))
  if (length(unknown)) {
    frm_stop(caller, "(\"", family, "\") has no argument ",
             paste0("`", unknown, "`", collapse = ", "), ". It takes ",
             paste0("`", setdiff(names(formals(ctor)), "..."), "`",
                    collapse = ", "), ". See ?`frmtmb-links`", call. = FALSE)
  }
  if (!is.null(link)) args <- c(list(link = link), args)
  do.call(ctor, args)
}

#' Build a family by name, with a link for any of its parameters
#'
#' `brmsfamily()` is brms's generic family constructor. It names the
#' family with a string, takes the link for the mean as the second
#' argument, and takes a `link_<dpar>` argument for each of the
#' family's other distributional parameters.
#'
#' The family constructors take the same arguments, so
#' `student(link_sigma = "softplus")` does not need it. They exist
#' for the families frmtmb has no constructor of its own for.
#' `gaussian()`, `poisson()`, `binomial()` and `Gamma()` come from
#' 'stats', which knows nothing of a `link_sigma` and refuses several
#' links brms allows for the mean. frmtmb does NOT shadow them: a
#' `gaussian()` that answered a `frmtmb_family` would break every
#' `glm()` call in an attached session.
#'
#' The name follows brms's spelling rules: case does not matter,
#' `"normal"` is `"gaussian"`, `"com_poisson"` is `"compois"`, and the
#' prefixes `"zi_"` and `"hu_"` stand for `"zero_inflated_"` and
#' `"hurdle_"`. The link may be unquoted, as in
#' `brmsfamily("gaussian", inverse)`.
#'
#' The family and its link as one character vector,
#' `c("weibull", "log")`, is the form brms accepts for `family =` in a
#' model call, and [frm()] accepts it there. `brmsfamily()` refuses it,
#' as brms does.
#'
#' @section Differences from brms:
#' An argument the family does not have is refused by name. brms
#' ignores it, so `brmsfamily("poisson", link_sigma = "log")` builds a
#' poisson family there. `threshold` is an argument of the ordinal
#' families alone, and `refcat` of none.
#'
#' @param family Family name, as a single string.
#' @param link The link for the mean, quoted or not. `NULL` gives the
#'   family's default. See [frmtmb-links] for the set each family takes.
#' @param ... Any `link_<dpar>` argument the named family takes, and its
#'   other constructor arguments. Every one has to be named.
#' @return A `frmtmb_family` object, ready for [frm()].
#' @seealso [frmtmb-links] for the links and which parameter takes
#'   which, [frmtmb-families] for the constructors.
#' @examples
#' # sigma on a softplus: gaussian() from 'stats' has no link_sigma
#' brmsfamily("gaussian", link_sigma = "softplus")
#'
#' # a mean link stats::poisson() refuses, and brms's short names
#' brmsfamily("poisson", softplus)$link
#' brmsfamily("zi_poisson")$link_zi
#'
#' set.seed(2)
#' d <- data.frame(x = rnorm(60))
#' d$y <- rnorm(60, 1 + d$x, exp(0.3 * d$x))
#' frm(bf(y ~ x, sigma ~ x),
#'     family = brmsfamily("gaussian", link_sigma = "softplus"), data = d)
#' @export
brmsfamily <- function(family, link = NULL, ...) {
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    frm_stop("brmsfamily(): `family` names a family with a single string, ",
             "as in brmsfamily(\"gaussian\", link_sigma = \"softplus\"); not ",
             arg_desc(family), ". frm(family =) also takes ",
             "c(family, link)", call. = FALSE)
  }
  link <- link_arg_value(substitute(link), link, names(frmtmb_links), NULL)
  family_from_name(family, link, list(...), caller = "brmsfamily")
}

#' Split `family = c(name, link)`, brms's one-vector spelling, and
#' refuse every other shape of a family name.
#'
#' brms reads the second element as the link and ignores any third; a
#' third element is refused here, because it means the call says
#' something that would be dropped.
#'
#' @noRd
family_link_pair <- function(family, caller) {
  if (!is.character(family) || !length(family) %in% 1:2 ||
      anyNA(family)) {
    frm_stop(caller, "(): `family` names a family, as \"weibull\" or as ",
             "c(\"weibull\", \"log\") with its link; not ", arg_desc(family),
             call. = FALSE)
  }
  list(family = family[1L],
       link = if (length(family) == 2L) family[2L])
}

#' @rdname frmtmb-extension-api
#' @export
as_frmtmb_family <- function(x) {
  if (inherits(x, "frmtmb_family")) return(x)
  # A bare constructor (`family = cumulative`, `family = gaussian`) is
  # legal in brms and in stats::glm, so call it with its defaults. The
  # frmtmb constructors return a `frmtmb_family` rather than a
  # `stats::family`, so the result is re-checked here.
  if (is.function(x)) {
    x <- x()
    if (inherits(x, "frmtmb_family")) return(x)
  }
  if (inherits(x, "family")) {
    # [[ ]]: a stats family is a plain list, and `$` on it partial-matches
    ctor <- family_ctor(x[["family"]])
    # a brms family object also carries link_<dpar>; carry the ones the
    # constructor takes rather than silently dropping them
    extra <- intersect(grep("^link_", names(x), value = TRUE),
                       names(formals(ctor)))
    return(do.call(ctor, c(list(link = x[["link"]]), unclass(x)[extra])))
  }
  if (is.character(x)) {
    fl <- family_link_pair(x, "frm")
    return(family_from_name(fl$family, fl$link, caller = "frm"))
  }
  frm_stop("Cannot interpret `family` of class ",
           paste(class(x), collapse = "/"),
           ": pass a family constructor or its name as the `family` ",
           "argument of frm(), or attach it with `bf(...) + gaussian()`",
           call. = FALSE)
}

#' Additional response families
#'
#' Family constructors without a [stats::family] equivalent, following
#' brms naming. `gaussian()`, `poisson()`, `binomial()`, and `Gamma()`
#' from 'stats' are accepted directly by [frm()] and `+`.
#'
#' An ordinal family (`cumulative()`, `sratio()`, `cratio()`, `acat()`)
#' takes the response's level order as the category order. Supply an
#' ordered factor, or integer codes `1..K`. An unordered factor is
#' refused, as brms refuses it, because its level order is alphabetical
#' unless someone set it, and that order is the model.
#'
#' A `bernoulli()` response can hold any two values, which are coded 0
#' and 1 by level order as brms codes them: the levels of a factor, the
#' sorted values of a number or a character vector (so `-2` and `-1`
#' become 0 and 1, and so do `1` and `2`), and `FALSE` and `TRUE` of a
#' logical. A number with one value is the 1 of 0 and 1 unless the value
#' is 0. The fit keeps the coding, and codes newdata's response with it.
#' A third value is refused with brms's message. Two values that both
#' lie strictly between 0 and 1, such as 0.1 and 0.9, are coded too, as
#' brms codes them, with a warning: they are usually proportions, which
#' `binomial()` with `trials()` or `Beta()` fits. [simulate()] and the
#' draws return the codes 0 and 1, as brms's `posterior_predict()` does.
#'
#' A response with only two outcomes gets brms's message suggesting
#' `bernoulli()`: an ordinal or categorical response with two
#' categories, and a `binomial()`, `beta_binomial()` or
#' `zero_inflated_binomial()` response whose trials are all one.
#'
#' Every constructor has brms's fields: `$link` is the name of the link
#' for the mean, and `$link_<dpar>` the link of each other parameter,
#' as in `beta_binomial()$link_phi`. See [frmtmb_family()].
#'
#' @section Hurdle and zero-one-inflated responses:
#' A hurdle family gives every zero to `hu`, the hurdle probability, and
#' truncates its count or continuous part at zero. For
#' `hurdle_negbinomial()`, `P(Y = 0) = hu` and a positive count `y` has
#' probability `(1 - hu) NB(y | mu, shape) / (1 - NB(0 | mu, shape))`,
#' so `mu` is the mean of the untruncated negative binomial and not the
#' mean of the positive counts.
#'
#' `zero_one_inflated_beta()` takes a response in `[0, 1]`. `zoi` is the
#' probability of an exact 0 or 1, `coi` is the probability that such a
#' value is 1, and a value strictly between 0 and 1 follows a beta
#' distribution with mean `mu` and precision `phi`. The expected
#' response, which `fitted()` returns, is `zoi * coi + (1 - zoi) * mu`.
#' `coi` is scored only on the values at exactly 0 or 1. A response with
#' 0s but no 1 (or the reverse) sends `coi` to that boundary, and one
#' with neither gives `coi` no data at all, so that no standard error of
#' the fit is usable; the fit warns in both cases. brms fits such data
#' through its prior on `coi`. Here, hold `coi` at a value with
#' `bf(coi = 0.5)`.
#'
#' `zero_inflated_beta_binomial()` is `beta_binomial()` with a
#' zero-inflation probability `zi`: `P(Y = 0) = zi + (1 - zi) BB(0)`,
#' and a count `y > 0` has probability `(1 - zi) BB(y)`, with the trials
#' from `trials()`. `fitted()` returns `(1 - zi) * mu * trials`.
#'
#' `hurdle_cumulative()` is brms's hurdle over an ordinal response. The
#' response is `0..K`: 0 is the hurdle, with probability `hu`, and the
#' categories `1..K` follow `cumulative()` with probabilities scaled by
#' `1 - hu`. An ordered factor takes its FIRST level as the category 0.
#' `fitted()` returns the `K + 1` category probabilities, 0 first. The
#' expected category that `conditional_effects()` shows by default is
#' scored by these codes, so the hurdle scores 0; brms scores each
#' column by its position and reads one higher. `disc` and `threshold`
#' work as for the other ordinal families (see Ordinal thresholds and
#' discrimination). `thres(x = )` works; `thres(gr = )` and `cs()` are
#' refused.
#'
#' @section Ordinal thresholds and discrimination:
#' Every ordinal family (`cumulative()`, `sratio()`, `cratio()`,
#' `acat()` and `hurdle_cumulative()`) has brms's discrimination
#' parameter `disc`. The distribution function is read at
#' `disc * (tau_k - eta)` (at `disc * (eta - tau_k)` for `cratio()` and
#' `acat()`), with `disc` held at 1 unless the formula models it, as in
#' brms: `bf(y ~ x, disc ~ z)`. Its link is `link_disc`, `"log"` by
#' default. The likelihood cannot tell an intercept in `disc` apart from
#' the scale of the thresholds, so write `disc ~ 0 + z`. brms accepts
#' `disc ~ 1 + z` and holds the intercept with its default
#' `normal(0, 1)` prior; frmtmb fits it too, and warns unless a prior
#' holds the intercept, as it does for cell-means columns that add up to
#' one (`disc ~ 0 + h`). Held at 1, `disc` shows nowhere in the output,
#' as in brms: not on the `Links:` line, not as a fixed parameter, and
#' not in `fixef(flatten = TRUE)` or `coef()`.
#'
#' `threshold` selects brms's threshold structure. `"flexible"`
#' estimates each threshold. `"equidistant"` estimates the first
#' threshold and the distance `delta` between neighbors, so that
#' `tau_k = tau_1 + (k - 1) * delta`; `variables()` and `summary()` name
#' it `delta` as brms does, and `set_prior(class = "delta")` gives it a
#' prior, while class `"Intercept"` addresses the first threshold only.
#' `delta` is positive for `cumulative()` and `hurdle_cumulative()`,
#' whose thresholds are ordered. `"sum_to_zero"` holds the thresholds
#' to a sum of zero and does not center the design, as in brms. brms
#' declares one parameter per threshold and relies on its prior for the
#' common location, which the likelihood does not see; frmtmb estimates
#' the `K - 2` free directions instead, so class `"Intercept"` has
#' nothing to address there and is refused. `thres(gr = )` gives each
#' level its own first threshold and `delta`, or its own zero sum.
#' `"equidistant"` needs two thresholds per vector, since `delta` has
#' nothing to measure on one. `"sum_to_zero"` holds a single threshold
#' at zero, a vector with no parameter, which brms runs under
#' `thres(gr = )` and frmtmb also fits without it.
#'
#' @section Extended-support beta:
#' `xbeta()` is the extended-support beta of Kosmidis and Zeileis
#' (2024), brms's `xbeta`, for a response in `[0, 1]` with exact 0s and
#' 1s. A latent beta variable with mean `mu` and precision `phi` is
#' stretched by the exceedance `kappa` to `(1 + 2 kappa) Z - kappa` and
#' censored at 0 and 1, so the ends carry the latent mass beyond them.
#' `mu` is the mean of the latent variable, not of the response;
#' `fitted()` returns the mean of the response, brms's
#' `posterior_epred()`. The mass at each end is an incomplete beta
#' function, which is taken mostly from a continued fraction rather
#' than from `RTMB::pbeta()`, whose third derivatives are not finite
#' everywhere and which therefore failed the Laplace approximation of a
#' model with random effects.
#'
#' `kappa` is placed mostly by the rows at exactly 0 or 1. With none,
#' only the shape of the interior places it, and it can run to 0,
#' where the model is `Beta()`, or grow without bound. When it runs to
#' 0, or a coefficient of a log-link `kappa` has a standard
#' error above 10 or none, the fit warns; brms fits such data through its
#' `gamma(0.01, 0.01)` prior on `kappa`. Compare the fit with `Beta()`,
#' or hold `kappa` at a value with `bf(kappa = )`.
#'
#' @section Categorical (nominal) responses:
#' `categorical()` fits a multinomial logit to an unordered factor. The
#' FIRST level is the reference category, its linear predictor is held
#' at zero, and each remaining level gets its own predictor named
#' `mu<Level>`, as in brms. The main model formula applies to every one
#' of them, and any single category is overridden by naming it:
#' `bf(y ~ x, mub ~ z)` gives category `"b"` its own predictor and
#' leaves the rest on `~ x`.
#'
#' The categories come from the data, so `frm()` reads them off the
#' response before it parses the formula. Constructing the family away
#' from a data set (to inspect it, or to reach the parser through
#' another entry point) needs them stated: `categorical(levels = c("a",
#' "b", "c"))`, or `categorical(K = 3)` for a response already coded
#' `1..K`, which names the dpars `mu2 ... muK`. A character or logical
#' response is coerced to a factor with a message naming the level
#' order, because that order is the model.
#'
#' `fitted()` and `frm_linpred(type = "response")` return the `n x K` matrix
#' of category probabilities, columns named by the response's own
#' levels and rows summing to one - the same convention the ordinal
#' families follow. `frm_linpred(type = "link")` and `frm_linpred(dpar =)` give
#' the per-category latent predictors, which is where `se.fit` lives.
#' `simulate()` draws factor levels. The same likelihood is available on
#' a count-matrix response as `multinomial(K)`, and a one-hot matrix
#' gives an identical log-likelihood.
#'
#' @section Circular responses:
#' `von_mises()` models an angle in radians on `(-pi, pi]`. Its `mu` is
#' the mean direction and takes the `tan_half` link, which maps the
#' whole line onto that interval; `kappa` is the concentration and takes
#' a log link, with `kappa = 0` the uniform distribution on the circle.
#' Both are brms's choices. `fitted()` and `frm_linpred(type = "response")`
#' report the mean direction. The normalizing constant needs
#' `log I0(kappa)`, which RTMB differentiates exactly through its own
#' `besselI` method, so nothing here is a series approximation.
#' Residuals are differences of angles and are NOT wrapped, so read
#' `residuals()` on a von Mises fit with that in mind.
#'
#' @section Cox proportional hazards:
#' `cox()` is the flexible-parametric proportional hazards model brms
#' fits: the baseline hazard is an M-spline in time,
#' `h0(t) = sum_j s_j M_j(t)`, and the cumulative baseline hazard is the
#' I-spline integral of the same basis over the same weights, which
#' makes it monotone by construction. The weights `s` form a simplex -
#' that is what identifies the baseline against the intercept - and are
#' estimated as `sbhaz_raw`, their `Kbhaz - 1` free softmax
#' coordinates; [cox_baseline()] returns the simplex itself. The
#' default basis is brms's: `df = 5` cubic M-splines with an intercept,
#' internal knots on response quantiles and boundary knots just outside
#' the observed range.
#'
#' The hazard is `h0(t) exp(eta)`, so a coefficient is a log hazard
#' ratio, exactly as in `survival::coxph()`. Right, left, and interval
#' censoring come through the ordinary `cens()` addition term: an event
#' contributes the density and a censored observation the survivor
#' function, which is what this family's log-density and log-CDF are.
#' Random effects are the point: `time | cens(c) ~ x + (1 | g)` is a
#' frailty model, and the Laplace approximation integrates the frailties
#' out.
#'
#' The baseline is semiparametric only in spirit - it has `df`
#' parameters, not one per event time - so coefficients agree with
#' `coxph()` closely rather than exactly. A survival response has no
#' mean, so `fitted()` and `frm_linpred(type = "response")` are refused;
#' `frm_linpred(type = "link")` gives the log hazard ratio. `simulate()` is
#' not available.
#'
#' Maximum likelihood often puts one or more baseline weights ON the
#' simplex boundary, at exactly zero. Their softmax coordinates then run
#' off to minus infinity along a flat ridge, the Hessian is singular in
#' those directions, and the optimizer reports singular convergence even
#' though the gradient is zero and the regression coefficients are at
#' their optimum - [diagnose()] names `sbhaz_raw` as the culprit. This
#' is what an unpenalized flexible baseline does; brms does not meet it
#' because its Dirichlet prior keeps the weights interior. Lower `df`
#' until the baseline is one the data supports, and read
#' [cox_baseline()] to see which weights collapsed.
#'
#' @section Quantile regression inference:
#' `asym_laplace()` and `zero_inflated_asym_laplace()` fit quantile
#' regression through a WORKING likelihood: at a fixed `quantile` the
#' point estimates are consistent quantile estimates (they match
#' `quantreg::rq()`), but the asymmetric Laplace is not the data's
#' true density, so Wald standard errors and `confint()` intervals
#' computed from it are not calibrated. This is a property of the
#' asymmetric-Laplace approach, shared by brms. Use
#' [frm_bootstrap()] for intervals you can defend. The check
#' function's kink can also produce a benign false-convergence
#' warning near the optimum; `frm_allfit()` confirms the fit when in
#' doubt.
#'
#' @section Degrees of freedom that run off:
#' `student()` estimates `nu` on the `logm1` link, which holds it above
#' one but puts no ceiling on it. A student-t reaches `gaussian()` only
#' in the limit `nu -> Inf`, so on data with no heavy tails the
#' likelihood keeps rising as `nu` grows and the maximum is never
#' attained. The fit then reports whatever `nu` the optimizer last
#' reached, with a standard error to match: values of `1e9` and above
#' are ordinary, and two runs of the same model on the same data in a
#' different row order can differ by orders of magnitude in `nu` while
#' agreeing on every coefficient to the last bit.
#'
#' That is the model saying the data shows no heavy tails, not a
#' failure to converge. [diagnose()] names it under "Distributional
#' parameter at the end of its link". Read the fit as the gaussian one
#' it has become. If you want a number you can report, refit with
#' `gaussian()`, or hold `nu` somewhere finite with a prior (see
#' [set_prior()]); brms does the same with its default
#' `gamma(2, 0.1)`.
#'
#' `nu` runs off the OTHER end too, and there the reading is opposite.
#' On tails heavier than any identified `nu` can hold, Cauchy data for
#' instance, `nu` goes down to one instead: `log(nu - 1)` around -20,
#' a standard error in the thousands, and a natural-scale value that
#' prints as `1`. [diagnose()] names that under the same heading. The
#' remedy is not `gaussian()`, which is the worst fit available for
#' such data. There the data is the message; a prior is still the way
#' to hold `nu` finite.
#'
#' The likelihood itself stays accurate the whole way. The log density
#' is formed so that `log Gamma((nu + 1) / 2) - log Gamma(nu / 2)`
#' never cancels, which holds it to 7e-15 of a 300-bit reference for
#' every `nu` up to `1e50`.
#'
#' @section Robust regression:
#' `huber()` fits Huber's least-favorable distribution: gaussian within
#' `k` residual standard deviations of `mu` and Laplace outside, so a
#' far-out point pulls on the fit with a bounded influence instead of
#' its squared distance. It is a proper normalized density, not a
#' penalty, so this is ordinary maximum likelihood and `logLik()`,
#' `AIC()` and the likelihood-ratio machinery all mean what they say.
#'
#' `k` is a fixed constant of the family, `huber(k = 1.345)`, not a
#' distributional parameter. It states where the analyst draws the line
#' between a residual and an outlier, which is a modelling choice rather
#' than something the data identifies; `MASS::rlm()` treats it the same
#' way. Estimating it would let the likelihood buy fit by widening the
#' gaussian core, which is the opposite of the point. As `k` grows the
#' family collapses to `gaussian()`.
#'
#' Point estimates track `MASS::rlm(psi = psi.huber)` closely but not
#' exactly, and the difference is the scale: `rlm()` fixes the scale at
#' a MAD-type estimate and iterates the location, while `huber()`
#' estimates `sigma` by maximum likelihood jointly with `mu`. The
#' coefficients agree to about `1e-2` on well-behaved data; the
#' `sigma` estimates need not.
#'
#' The working-likelihood caveat that applies to `asym_laplace()`
#' applies here for the same reason. If the data are not actually
#' Huber-distributed - and the family is chosen precisely because the
#' error distribution is unknown - then the model is misspecified, the
#' information matrix is not the variance of the score, and Wald
#' standard errors and `confint()` intervals from it are not
#' calibrated. The point estimates stay consistent for the location.
#' Use [frm_bootstrap()] for intervals you can defend.
#' `cens()` and `trunc()` are unavailable: the CDF is a three-piece
#' function of a parameter-dependent residual and has no branch-free
#' form for the tape.
#'
#' `rho` has a kink at `|u| = k`, so the objective is only piecewise
#' smooth and the optimizer often stops with a maximum absolute
#' gradient around `1e-4` and the accompanying false-convergence
#' warning. That is the kink, not a bad fit: the same thing happens to
#' `asym_laplace()`. The estimates satisfy Huber's own estimating
#' equations, `X' psi(u) = 0` with `psi(u) = min(max(u, -k), k)`, to
#' the same order. [frm_allfit()] confirms the fit when in doubt.
#'
#' @section Links:
#' Every constructor takes `link` for the mean and a `link_<dpar>` for
#' each of its other distributional parameters, following brms:
#' `student(link_sigma = "softplus")`, `zero_inflated_poisson(link_zi =
#' "identity")`. [frmtmb-links] lists the whole roster, what each link
#' maps, which families take it for the mean, and which set each
#' parameter admits.
#'
#' The four families 'stats' owns, `gaussian()`, `poisson()`,
#' `binomial()` and `Gamma()`, have no frmtmb constructor to carry
#' these. Reach their links through [brmsfamily()]:
#' `brmsfamily("gaussian", link_sigma = "softplus")`.
#'
#' The link for the mean may be unquoted, `student(identity)`, and is
#' refused by name when brms does not allow it for the family:
#' `bernoulli("sqrt")` is an error.
#'
#' An ordinal family's `link` is not a link on a mean. It names the
#' distribution function the thresholds are read through, so
#' all five ordinal families take `logit`, `probit`, `probit_approx`,
#' `cloglog` and `cauchit` (and `cumulative()`, `hurdle_cumulative()`
#' and `acat()` also take `softit`), and refuse anything else. Off the
#' logit, `acat()` reads its categories in brms's second form, a product
#' of distribution functions times a reversed product of survivals (see
#' [frmtmb-links]).
#'
#' `cumulative()` and `hurdle_cumulative()` keep their thresholds
#' increasing, because their category probabilities are differences of
#' the distribution function. The
#' thresholds of `sratio()`, `cratio()` and `acat()` are unconstrained,
#' as in brms: their category probabilities are positive for any
#' thresholds, so a fit may have two of them cross.
#'
#' @param link Link for `mu`. See [frmtmb-links].
#' @param link_sigma,link_shape,link_phi,link_kappa,link_ndt,link_beta,link_disc
#'   Link for a strictly positive parameter: one of `"log"` (the
#'   default), `"identity"`, `"softplus"` or `"squareplus"`.
#' @param threshold For the ordinal families: the threshold structure,
#'   as in brms. `"flexible"` (the default) estimates every threshold.
#'   `"equidistant"` estimates the first threshold and `delta`, the
#'   distance between neighboring thresholds. `"sum_to_zero"` holds
#'   the thresholds to a sum of zero. See Ordinal thresholds and
#'   discrimination.
#' @param link_nu Link for `nu`. `student()`'s degrees of freedom take
#'   `"logm1"` (the default) or `"identity"`, which keeps them above
#'   one; `compois()`'s dispersion is an ordinary positive parameter
#'   and takes the positive set.
#' @param link_zi,link_hu,link_quantile,link_zoi,link_coi Link for a
#'   parameter on the unit interval: `"logit"` (the default) or
#'   `"identity"`. In `zero_one_inflated_beta()`, `zoi` is the
#'   probability of an exact 0 or 1, and `coi` is the probability that
#'   such a value is 1.
#' @param link_alpha Link for `skew_normal()`'s skewness, which is
#'   signed: `"identity"` (the default), `"log"`, `"softplus"` or
#'   `"squareplus"`.
#' @return A `frmtmb_family` object.
#' @examples
#' set.seed(4)
#' n <- 120
#' dd <- data.frame(x = rnorm(n))
#'
#' # heavier tails than gaussian(), with an estimated df
#' dd$y <- 1 + 0.8 * dd$x + rt(n, df = 4)
#' fixef(frm(bf(y ~ x) + student(), data = dd))
#'
#' # counts with more spread than poisson() allows
#' dd$cnt <- rnbinom(n, mu = exp(0.5 + 0.4 * dd$x), size = 2)
#' fit <- frm(bf(cnt ~ x) + negbinomial(), data = dd)
#' fixef_by_dpar(fit)$mu
#'
#' # a zero-inflated count: the zi dpar gets its own predictor
#' dd$zi <- ifelse(runif(n) < 0.3, 0, dd$cnt)
#' frm(bf(zi ~ x, zi ~ 1) + zero_inflated_poisson(), data = dd)
#'
#' # an ordered response: level order is the category order
#' dd$grade <- cut(1 + 0.8 * dd$x + rlogis(n), 3,
#'                 labels = c("low", "mid", "high"), ordered_result = TRUE)
#' frm(bf(grade ~ x) + cumulative(), data = dd)
#'
#' # a proportion in (0, 1)
#' dd$p <- plogis(0.2 + 0.6 * dd$x + rnorm(n, 0, 0.3))
#' frm(bf(p ~ x) + Beta(), data = dd)
#'
#' # bounded influence: a few wild points barely move the slope
#' dd$rob <- 1 + 0.8 * dd$x + rnorm(n)
#' dd$rob[1:5] <- dd$rob[1:5] + 30
#' fixef_by_dpar(frm(bf(rob ~ x), family = huber(), data = dd))$mu
#' fixef_by_dpar(frm(bf(rob ~ x), family = gaussian(), data = dd))$mu
#'
#' # an unordered factor: one predictor per non-reference category,
#' # named after the level it belongs to
#' dd$pick <- factor(sample(c("ale", "stout", "lager"), n, TRUE))
#' cat_fit <- frm(bf(pick ~ x), family = categorical(), data = dd)
#' fixef(cat_fit)                     # mulager and mustout; ale is the
#'                                    # reference
#' head(fitted(cat_fit))              # n x K category probabilities
#'
#' # one category may take its own predictor
#' dd$w <- rnorm(n)
#' frm(bf(pick ~ x, mustout ~ w), family = categorical(), data = dd)
#'
#' # an angle: mu is the mean direction, kappa the concentration
#' dd$angle <- atan2(sin(0.5 + dd$x), cos(0.5 + dd$x))
#' vm_fit <- frm(bf(angle ~ x), family = von_mises(), data = dd)
#' head(fitted(vm_fit))               # the mean direction, in radians
#'
#' # proportional hazards with a spline baseline; (1 | g) is a frailty
#' dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
#' dd$out <- rbinom(n, 1, 0.3)        # 1 = right censored
#' cox_fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
#' fixef_by_dpar(cox_fit)$mu                  # log hazard ratios
#' cox_baseline(cox_fit)              # the baseline hazard weights
#' @name frmtmb-families
NULL

#' @rdname frmtmb-families
#' @export
student <- function(link = "identity", link_sigma = "log", link_nu = "logm1") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["student"]], "identity")
  fam_student(link, link_sigma, link_nu)
}

#' @rdname frmtmb-families
#' @export
lognormal <- function(link = "identity", link_sigma = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["lognormal"]], "identity")
  fam_lognormal(link, link_sigma)
}

#' @rdname frmtmb-families
#' @export
negbinomial <- function(link = "log", link_shape = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["negbinomial"]], "log")
  fam_negbinomial(link, link_shape)
}

#' @rdname frmtmb-families
#' @export
nbinom1 <- function(link = "log", link_phi = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["nbinom1"]], "log")
  fam_nbinom1(link, link_phi)
}

#' @rdname frmtmb-families
#' @export
Beta <- function(link = "logit", link_phi = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["beta"]], "logit")
  fam_beta(link, link_phi)
}

#' @rdname frmtmb-families
#' @export
tweedie <- function(link = "log", link_phi = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["tweedie"]], "log")
  fam_tweedie(link, link_phi)
}

#' @rdname frmtmb-families
#' @export
compois <- function(link = "log", link_nu = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["compois"]], "log")
  fam_compois(link, link_nu)
}

#' @rdname frmtmb-families
#' @export
zero_inflated_poisson <- function(link = "log", link_zi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_inflated_poisson"]], "log")
  fam_zi_poisson(link, link_zi)
}

#' @rdname frmtmb-families
#' @export
zero_inflated_negbinomial <- function(link = "log", link_shape = "log",
                                      link_zi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_inflated_negbinomial"]], "log")
  fam_zi_negbinomial(link, link_shape, link_zi)
}

#' @rdname frmtmb-families
#' @export
hurdle_poisson <- function(link = "log", link_hu = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["hurdle_poisson"]], "log")
  fam_hurdle_poisson(link, link_hu)
}

#' @rdname frmtmb-families
#' @export
hurdle_negbinomial <- function(link = "log", link_shape = "log",
                               link_hu = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["hurdle_negbinomial"]], "log")
  fam_hurdle_negbinomial(link, link_shape, link_hu)
}

#' @rdname frmtmb-families
#' @param K For `multinomial()`: number of response categories (columns
#'   of the count-matrix response); category 1 is the reference.
#' @export
multinomial <- function(K) fam_multinomial(K)

#' @rdname frmtmb-families
#' @export
cumulative <- function(link = "logit", link_disc = "log",
                       threshold = "flexible") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["cumulative"]], "logit")
  fam_cumulative(link, link_disc, threshold)
}

#' @rdname frmtmb-families
#' @export
hurdle_cumulative <- function(link = "logit", link_hu = "logit",
                              link_disc = "log", threshold = "flexible") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["hurdle_cumulative"]], "logit")
  fam_hurdle_cumulative(link, link_hu, link_disc, threshold)
}

#' @rdname frmtmb-families
#' @export
beta_binomial <- function(link = "logit", link_phi = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["beta_binomial"]], "logit")
  fam_beta_binomial(link, link_phi)
}

#' @rdname frmtmb-families
#' @export
zero_inflated_beta_binomial <- function(link = "logit", link_phi = "log",
                                        link_zi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_inflated_beta_binomial"]],
                         "logit")
  fam_zi_beta_binomial(link, link_phi, link_zi)
}

#' @rdname frmtmb-families
#' @export
skew_normal <- function(link = "identity", link_sigma = "log",
                        link_alpha = "identity") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["skew_normal"]], "identity")
  fam_skew_normal(link, link_sigma, link_alpha)
}

#' @rdname frmtmb-families
#' @export
exgaussian <- function(link = "identity", link_sigma = "log",
                       link_beta = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["exgaussian"]], "identity")
  fam_exgaussian(link, link_sigma, link_beta)
}

#' @rdname frmtmb-families
#' @export
bernoulli <- function(link = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["bernoulli"]], "logit")
  fam_bernoulli(link)
}

#' @rdname frmtmb-families
#' @export
geometric <- function(link = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["geometric"]], "log")
  fam_geometric(link)
}

#' @rdname frmtmb-families
#' @export
exponential <- function(link = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["exponential"]], "log")
  fam_exponential(link)
}

#' @rdname frmtmb-families
#' @export
weibull <- function(link = "log", link_shape = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["weibull"]], "log")
  fam_weibull(link, link_shape)
}

#' @rdname frmtmb-families
#' @export
shifted_lognormal <- function(link = "identity", link_sigma = "log",
                              link_ndt = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["shifted_lognormal"]], "identity")
  fam_shifted_lognormal(link, link_sigma, link_ndt)
}

#' @rdname frmtmb-families
#' @export
hurdle_gamma <- function(link = "log", link_shape = "log", link_hu = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["hurdle_gamma"]], "log")
  fam_hurdle_gamma(link, link_shape, link_hu)
}

#' @rdname frmtmb-families
#' @export
hurdle_lognormal <- function(link = "identity", link_sigma = "log",
                             link_hu = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["hurdle_lognormal"]], "identity")
  fam_hurdle_lognormal(link, link_sigma, link_hu)
}

#' @rdname frmtmb-families
#' @export
zero_inflated_binomial <- function(link = "logit", link_zi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_inflated_binomial"]], "logit")
  fam_zi_binomial(link, link_zi)
}

#' @rdname frmtmb-families
#' @export
zero_inflated_beta <- function(link = "logit", link_phi = "log",
                               link_zi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_inflated_beta"]], "logit")
  fam_zi_beta(link, link_phi, link_zi)
}

#' @rdname frmtmb-families
#' @export
zero_one_inflated_beta <- function(link = "logit", link_phi = "log",
                                   link_zoi = "logit", link_coi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_one_inflated_beta"]], "logit")
  fam_zoi_beta(link, link_phi, link_zoi, link_coi)
}

#' @rdname frmtmb-families
#' @export
xbeta <- function(link = "logit", link_phi = "log", link_kappa = "log") {
  link <- link_arg_value(substitute(link), link, brms_mu_links[["xbeta"]],
                         "logit")
  fam_xbeta(link, link_phi, link_kappa)
}

#' @rdname frmtmb-families
#' @export
asym_laplace <- function(link = "identity", link_sigma = "log",
                         link_quantile = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["asym_laplace"]], "identity")
  fam_asym_laplace(link, link_sigma, link_quantile)
}

#' @rdname frmtmb-families
#' @export
zero_inflated_asym_laplace <- function(link = "identity", link_sigma = "log",
                                       link_quantile = "logit",
                                       link_zi = "logit") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["zero_inflated_asym_laplace"]],
                         "identity")
  fam_zi_asym_laplace(link, link_sigma, link_quantile, link_zi)
}

#' @rdname frmtmb-families
#' @param k For `huber()`: Huber's tuning constant, the residual size in
#'   units of `sigma` where the density stops being gaussian and becomes
#'   Laplace. It is FIXED, not estimated - the default 1.345 is
#'   `MASS::rlm()`'s, which gives 95% efficiency against a gaussian.
#' @export
huber <- function(link = "identity", k = 1.345, link_sigma = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["huber"]], "identity")
  fam_huber(link, k, link_sigma)
}

#' @rdname frmtmb-families
#' @export
sratio <- function(link = "logit", link_disc = "log",
                   threshold = "flexible") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["sratio"]], "logit")
  fam_sratio(link, link_disc, threshold)
}

#' @rdname frmtmb-families
#' @export
cratio <- function(link = "logit", link_disc = "log",
                   threshold = "flexible") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["cratio"]], "logit")
  fam_cratio(link, link_disc, threshold)
}

#' @rdname frmtmb-families
#' @export
acat <- function(link = "logit", link_disc = "log",
                 threshold = "flexible") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["acat"]], "logit")
  fam_acat(link, link_disc, threshold)
}

#' @rdname frmtmb-families
#' @export
von_mises <- function(link = "tan_half", link_kappa = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["von_mises"]], "tan_half")
  fam_von_mises(link, link_kappa)
}

#' @rdname frmtmb-families
#' @param levels For `categorical()`: the response's category labels, in
#'   the order that fixes the reference category (the first) and the
#'   dpar names. Only needed when the family is built away from the
#'   data; [frm()] reads them off the response.
#' @export
categorical <- function(link = "logit", levels = NULL, K = NULL) {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["categorical"]], "logit")
  categorical_link_check(link)
  if (!is.null(levels)) {
    levels <- as.character(levels)
    if (length(levels) < 2L || anyDuplicated(levels)) {
      frm_stop("categorical(levels =): needs at least two distinct category ",
               "labels", call. = FALSE)
    }
    return(fam_categorical_impl(paste0("mu", levels[-1L]), levels, link))
  }
  if (!is.null(K)) {
    if (length(K) != 1L || is.na(K) || K < 2) {
      frm_stop("categorical(K =): needs at least two categories",
               call. = FALSE)
    }
    K <- as.integer(K)
    return(fam_categorical_impl(paste0("mu", seq_len(K)[-1L]),
                                levels = NULL, link = link))
  }
  fam_categorical_deferred(link)
}

#' @rdname frmtmb-families
#' @param df For `cox()`: the number of M-spline basis functions in the
#'   baseline hazard (brms's `bhaz(df = 5)` default).
#' @param degree For `cox()`: the spline degree of that basis (cubic by
#'   default).
#' @param intercept For `cox()`: keep the basis function that is
#'   non-zero at the lower boundary knot.
#' @export
cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE) {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["cox"]], "log")
  # df and degree size the I-spline baseline basis, so a length-2 value
  # used to build a basis of one shape and record another
  check_count(df, "df", min = 1L)
  check_count(degree, "degree", min = 1L)
  check_flag(intercept, "intercept")
  fam_cox(link, df = df, degree = degree, intercept = intercept)
}

#' The links of one family, as `mu = identity; sigma = log`, in dpar
#' order.
#'
#' `summary()` and `print()` both show it. A coefficient is reported on
#' the LINK scale, and the family name alone does not say which scale
#' that is: every family in the roster now has more than one choice for
#' its mean, and most have one for each of their other parameters too.
#'
#' An ordinal family reports the distribution function its thresholds
#' are read through instead. Its `mu` really does carry an identity
#' link, so printing that would answer a question the reader did not
#' ask.
#'
#' @noRd
family_link_str <- function(fam, shown = character(0)) {
  ol <- fam[["ord_link"]]
  lk <- fam[["links"]]
  # a dpar the family holds at a default the user did not write, and
  # that brms does not name, unless a fit says it is modeled
  dp <- setdiff(fam[["dpars"]],
                setdiff(fam[["hidden_fixed_dpars"]] %||% character(0),
                        shown))
  if (!is.null(ol)) {
    # an ordinal mu's link is the identity; the cdf takes its place, and
    # a dpar beside mu (a hurdle's hu and disc) keeps its own
    rest <- setdiff(dp, "mu")
    more <- if (length(rest)) {
      paste0(rest, " = ", vapply(rest, function(d) {
        lk[[d]][["name"]] %||% "?"
      }, ""))
    }
    return(paste(c(paste0("cdf = ", ol[["name"]]), more), collapse = "; "))
  }
  if (!length(dp) || !length(lk)) return("")
  nm <- vapply(dp, function(d) {
    l <- lk[[d]]
    if (is.null(l)) NA_character_ else (l[["name"]] %||% NA_character_)
  }, character(1))
  keep <- !is.na(nm)
  if (!any(keep)) return("")
  paste(paste0(dp[keep], " = ", nm[keep]), collapse = "; ")
}

#' Print a family
#'
#' Prints the family, its link and each distributional parameter with
#' its link, in brms's layout: the first two lines are brms's, spacing
#' included, and a mixture opens with `Mixture` and names each
#' component with its own link.
#'
#' @param x A family.
#' @param links `TRUE` also prints brms's line per distributional
#'   parameter other than `mu`, "Link function of 'sigma' (if
#'   predicted): log"; a character vector prints those parameters only.
#' @param newline Print a blank line at the end, as brms does.
#' @param ... Unused; an argument given here is refused by name.
#' @return `x`, invisibly.
#' @examples
#' print(student(), links = TRUE)
#' print(mixture(gaussian(), exponential()))
#' @export
print.frmtmb_family <- function(x, links = FALSE, newline = TRUE, ...) {
  frm_check_dots(...)
  if (!(isTRUE(links) || isFALSE(links) ||
          (is.character(links) && !anyNA(links)))) {
    frm_stop("`links` must be TRUE, FALSE or the names of distributional ",
             "parameters, not ", arg_desc(links), call. = FALSE)
  }
  check_flag(newline, "newline")
  lk <- vapply(x[["links"]], function(l) {
    if (is.list(l)) l[["name"]] %||% "?" else as.character(l)
  }, "")
  comps <- x[["component_families"]]
  if (!is.null(comps)) {
    # brms's print.mixfamily: the word, then each component in turn
    cat("\nMixture\n")
    for (k in seq_along(comps)) {
      cat("\nFamily:", comps[[k]], "\n")
      cat("Link function:", lk[[paste0("mu", k)]], "\n")
    }
  } else {
    cat("\nFamily:", x[["family"]], "\n")
    cat("Link function:", x$link, "\n")
  }
  cat("Parameters: ", paste0(x[["dpars"]], " (", lk[x[["dpars"]]], ")",
                             collapse = ", "), "\n", sep = "")
  if (!isFALSE(links)) {
    dp <- setdiff(x[["dpars"]], "mu")
    if (is.character(links)) dp <- intersect(dp, links)
    for (d in dp) {
      cat("Link function of '", d, "' (if predicted): ", lk[[d]], "\n",
          sep = "")
    }
  }
  if (newline) cat("\n")
  invisible(x)
}

#' Validate `frmtmb_family(required_aterms =)`.
#'
#' Two spellings, because two things are being said. A character vector
#' is a conjunction: every value named must reach the density. A list
#' adds the disjunction: an element of length one is one required value,
#' and an element of length more than one is a set of spellings the
#' density reads ANY of. A family whose per-row datum arrives through
#' either of two addition terms needs the second form; without it the
#' refusal has to be hand-rolled inside `valid_y`, where it fires after
#' the frame is built and says nothing the framework can check.
#'
#' @noRd
check_required_aterms <- function(x) {
  ok <- function(z) {
    is.character(z) && length(z) >= 1L && !anyNA(z) && all(nzchar(z))
  }
  if (is.character(x)) {
    if (length(x) && !ok(x)) {
      frm_stop("frmtmb_family(required_aterms =) names addition-term ",
               "values, and one of the names it was given is missing or ",
               "empty.", call. = FALSE)
    }
    return(invisible(NULL))
  }
  # !is.object: a data frame is a list of character columns often
  # enough to be accepted by the test below and mean nothing here
  if (is.list(x) && !is.object(x) && all(vapply(x, ok, NA))) {
    return(invisible(NULL))
  }
  frm_stop("frmtmb_family(required_aterms =) names the addition terms the ",
           "density cannot do without: a character vector for the values it ",
           "needs ALL of, or a list whose length-one elements are required ",
           "and whose longer elements are alternatives, one of which must ",
           "be supplied. Got ", arg_desc(x), call. = FALSE)
}

#' Validate `frmtmb_family(accepts_aterms =)`.
#'
#' Membership is NOT checked against the addition-term registry: a
#' family may name a term its own package registers, and a package
#' registers from `.onLoad()` while a family object can be built in any
#' order. A name that no term ever carries costs nothing, while a name
#' checked too early would refuse a legitimate declaration.
#'
#' @noRd
check_accepts_aterms <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.character(x) || anyNA(x) || (length(x) && !all(nzchar(x)))) {
    frm_stop("frmtmb_family(accepts_aterms =) names the addition terms the ",
             "family takes, as a formula writes them and without ",
             "parentheses: c(\"weights\", \"trials\"). character(0) declares ",
             "a family that takes none, and NULL accepts every registered ",
             "term. Got ", arg_desc(x), call. = FALSE)
  }
  paren <- grepl("\\(", x, fixed = FALSE)
  if (any(paren)) {
    frm_stop("frmtmb_family(accepts_aterms =) names a term WITHOUT its ",
             "parentheses, because the name has to match the term however ",
             "many arguments it takes: write \"",
             sub("\\(.*$", "", x[paren][[1L]]),
             "\", not \"", x[paren][[1L]], "\"", call. = FALSE)
  }
  unique(x)
}

#' The addition terms a family takes, as a formula writes them: its
#' declared allow-list widened by what it requires, because a term the
#' density cannot do without is one it accepts.
#'
#' `NULL` means every term, which is what a family that declares no
#' allow-list keeps.
#'
#' @noRd
accepted_aterm_names <- function(fam) {
  acc <- fam[["accepts_aterms"]]
  if (is.null(acc)) return(NULL)
  req <- unlist(required_aterm_groups(fam[["required_aterms"]]),
                use.names = FALSE)
  sort(unique(c(acc, vapply(req %||% character(0), aterm_base, ""))))
}

#' Validate `frmtmb_family(se_dpar =)`: one dpar name, or `NA` for a
#' family whose whole scale is the known standard error.
#'
#' A location parameter is refused by name, because `se()` replaces the
#' residual SCALE and mapping out `mu` would delete the model.
#'
#' @noRd
check_se_dpar <- function(se_dpar, dpars, primary_dpars, family) {
  if (is.null(se_dpar)) return(NULL)
  if (length(se_dpar) != 1L ||
      !(is.character(se_dpar) || (is.logical(se_dpar) && is.na(se_dpar)))) {
    frm_stop("frmtmb_family(se_dpar =) names the ONE dpar that a known ",
             "standard error replaces, or is NA where it replaces none ",
             "because the whole scale is the known one; got ",
             arg_desc(se_dpar), call. = FALSE)
  }
  if (is.na(se_dpar)) return(NA_character_)
  if (!se_dpar %in% dpars) {
    frm_stop("frmtmb_family(se_dpar = \"", se_dpar, "\") names a dpar '",
             family, "' does not have (it has: ",
             paste(dpars, collapse = ", "),
             "). Name the scale se() replaces, or NA if it replaces none",
             call. = FALSE)
  }
  if (se_dpar %in% primary_dpars) {
    frm_stop("frmtmb_family(se_dpar = \"", se_dpar, "\") names a LOCATION ",
             "parameter of '", family, "'. se() carries a known standard ",
             "deviation, so it replaces a scale; mapping out the location ",
             "would leave the model with nothing to estimate", call. = FALSE)
  }
  se_dpar
}

#' The dpar a known `se()` replaces: the declaration if the family made
#' one, `NA_character_` where it declared that `se()` replaces nothing,
#' and NULL where it said nothing and has no `sigma` to fall back on.
#'
#' NULL is the case the parse-time guard refuses on, because an
#' undeclared family with a second free dpar cannot be told from one
#' whose second dpar is a shape the density reads.
#'
#' @noRd
family_se_dpar <- function(fam) {
  d <- fam[["se_dpar"]]
  if (!is.null(d)) return(d)
  if ("sigma" %in% fam[["dpars"]]) "sigma" else NULL
}

#' Whether a family EXPLICITLY declares that it reads an addition term.
#'
#' `accepted_aterm_names()` answers the allow-list's question - would
#' this term be IGNORED - and an undeclared family answers "no term is
#' ignored", which is why it returns NULL. This answers the opposite
#' question, and an undeclared family answers nothing: it is the test
#' for a term whose entire effect is inside the density, where the core
#' cannot act on the value by itself and reading it is the family's own
#' promise. `se()` is that term.
#'
#' @noRd
family_declares_aterm <- function(fam, term) {
  req <- unlist(required_aterm_groups(fam[["required_aterms"]]),
                use.names = FALSE)
  term %in% c(fam[["accepts_aterms"]] %||% character(0),
              vapply(req %||% character(0), aterm_base, ""))
}

#' Refuse an addition term the family does not declare.
#'
#' Runs at frame assembly, after every guard that can say something
#' specific about the pair, so a message here is the one nothing else
#' would have given.
#'
#' @noRd
check_accepted_aterms <- function(resp, av) {
  fam <- resp$family
  ok <- accepted_aterm_names(fam)
  if (is.null(ok)) return(invisible(NULL))
  # subset() and index() choose rows for every family and reach no
  # density, so no family declares them (R/subset.R)
  have <- setdiff(unique(c(names(resp$aterms), names(av))), row_aterms)
  bad <- setdiff(unique(vapply(have, aterm_base, "")), ok)
  if (!length(bad)) return(invisible(NULL))
  frm_stop(fam[["family"]], ": the addition term",
           if (length(bad) > 1L) "s " else " ",
           paste0("`", bad, "()`", collapse = ", "),
           if (length(bad) > 1L) " are not ones" else " is not one",
           " this family reads, so writing ",
           if (length(bad) > 1L) "them" else "it",
           " would change nothing about the fit. ",
           if (length(ok)) {
             paste0("This family takes ",
                    paste0("`", ok, "()`", collapse = ", "), ".")
           } else {
             "This family takes no addition terms."
           },
           call. = FALSE,
           package = frm_family_package(fam))
}

#' `required_aterms` as the groups frame assembly checks: each element
#' is a set of addition-term values, at least one of which must reach
#' the density.
#'
#' A character vector is a conjunction, so each of its entries becomes a
#' group of one and the argument's original meaning is unchanged.
#'
#' @noRd
required_aterm_groups <- function(x) {
  if (is.null(x) || !length(x)) return(list())
  if (is.character(x)) return(as.list(x))
  lapply(x, as.character)
}

#' `exclusive_aterms` as the groups frame assembly checks: each element
#' is a set of addition-term VALUES that mean the same thing to the
#' density, at most one of which may be supplied.
#'
#' Spelled in values rather than in terms, so that it lines up with the
#' any-of groups of `required_aterms` it is usually the other half of:
#' `dec` and `vint1` are one datum under two spellings, while `vint2`
#' beside them is a second datum.
#'
#' A bare character vector is ONE set, which is the opposite of
#' `required_aterms`'s reading of the same shape. The two arguments say
#' different things about a list of names -- everything, versus at most
#' one thing -- so a shared convention would have made one of them
#' unwritable.
#'
#' @noRd
exclusive_aterm_groups <- function(x) {
  if (is.null(x) || !length(x)) return(list())
  if (is.character(x)) return(list(as.character(x)))
  lapply(x, as.character)
}

#' Validate `frmtmb_family(exclusive_aterms =)`.
#'
#' The cross-check against `required_aterms` is the one that earns its
#' place: a family whose conjunction demands two values and whose
#' exclusivity forbids them together has declared a model nobody can
#' fit, and every fit would then fail at frame assembly with a message
#' about the user's formula rather than about the family.
#'
#' Membership in the addition-term registry is NOT checked, for the
#' reason `check_accepts_aterms()` gives: a family may name a value from
#' a term its own package registers, and registration order is not
#' fixed.
#'
#' @noRd
check_exclusive_aterms <- function(x, required = character(0)) {
  bad <- function(why = NULL, set = NULL) {
    frm_stop("frmtmb_family(exclusive_aterms =) names sets of addition-term ",
             "VALUES that say the same thing to the density, so that at ",
             "most one of each set may be supplied: a character vector for ",
             "one set, or a list of them. Each set needs at least two ",
             "distinct values. ",
             if (is.null(why)) paste0("Got ", arg_desc(x)) else {
               paste0("Set ", set, " ", why)
             },
             call. = FALSE)
  }
  if (is.null(x) || (is.list(x) && !length(x)) ||
        (is.character(x) && !length(x))) {
    return(list())
  }
  if (is.object(x) || !(is.character(x) || is.list(x))) bad()
  # The RAW elements, before exclusive_aterm_groups() coerces them.
  # Validating the coerced value cannot fail for anything as.character()
  # accepts, so list(c(1, 2)) was taken and became c("1", "2"): a set
  # matching no addition-term value, and a rule that never fires. The
  # sibling required_aterms refuses that shape, and so must this.
  raw <- if (is.character(x)) list(x) else x
  for (i in seq_along(raw)) {
    z <- raw[[i]]
    if (!is.character(z)) {
      bad(paste0("is ", arg_desc(z),
                 ", and a set is a character vector of term values"), i)
    }
    if (anyNA(z) || !all(nzchar(z))) {
      bad("has a missing or empty value in it", i)
    }
    if (length(unique(z)) < 2L) {
      bad(paste0("has ", length(unique(z)),
                 " distinct value(s), and exclusivity needs two"), i)
    }
  }
  grps <- exclusive_aterm_groups(x)
  # A length-one required group is a value the density cannot do
  # without. Two of them inside one exclusive set is unsatisfiable.
  must <- unlist(Filter(function(g) length(g) == 1L,
                        required_aterm_groups(required)),
                 use.names = FALSE)
  for (g in grps) {
    clash <- intersect(g, must)
    if (length(clash) > 1L) {
      frm_stop("frmtmb_family(exclusive_aterms =) makes ",
               paste0("`", clash, "`", collapse = " and "),
               " mutually exclusive, and required_aterms demands both of ",
               "them, so no model could satisfy the family. Make them ",
               "alternatives instead: required_aterms = list(c(",
               paste0("\"", clash, "\"", collapse = ", "), ")).",
               call. = FALSE)
    }
  }
  grps
}

#' The exclusive sets a mixture keeps: the INTERSECTION of what its
#' components declare, not the sets they agree on exactly.
#'
#' The rule is that two spellings are interchangeable for the mixture
#' only if EVERY component treats them so, and set equality is stricter
#' than that in the direction that fails open. A component declaring
#' `c("dec", "vint1", "vint2")` beside one declaring `c("dec", "vint1")`
#' does treat `dec` and `vint1` as one datum; comparing whole sets found
#' them unequal, dropped the rule, and let the mixture take both
#' spellings together, which is the defect this argument exists to
#' close.
#'
#' Intersecting against EVERY set of each later component rather than
#' the best-matching one keeps a component that splits one set in two
#' from silently costing the mixture the half it still shares. Order
#' inside a kept set comes from the first component, so the refusal
#' names the spelling that component reads.
#'
#' @noRd
mixture_exclusive_aterms <- function(comps) {
  grps <- lapply(comps, function(f) {
    exclusive_aterm_groups(f[["exclusive_aterms"]])
  })
  if (!length(grps) || any(!lengths(grps))) return(list())
  keep <- grps[[1L]]
  for (gs in grps[-1L]) {
    nxt <- list()
    for (g in keep) {
      for (s in gs) {
        hit <- g[g %in% s]
        if (length(hit) >= 2L) nxt[[length(nxt) + 1L]] <- hit
      }
    }
    keep <- nxt
    if (!length(keep)) return(list())
  }
  # A set fully contained in another says nothing the larger one does
  # not, and two components can produce the same set by two routes.
  keep <- keep[!duplicated(lapply(keep, sort))]
  Filter(function(g) {
    !any(vapply(keep, function(h) {
      length(h) > length(g) && all(g %in% h)
    }, NA))
  }, keep)
}

#' `c("a", "b", "c")` as "a, b and c", for a message that has to list
#' an arbitrary number of names and still read as a sentence.
#'
#' @noRd
comma_and <- function(x) {
  if (length(x) < 2L) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[[length(x)]])
}

#' Refuse two spellings of one datum supplied together.
#'
#' Runs immediately after the `required_aterms` check, because it is the
#' same declaration read the other way round: that one asks whether the
#' datum arrived at all, this one whether it arrived twice. Running it
#' here rather than with the allow-list keeps a family's own `valid_y`
#' from reporting on values the density was never going to read.
#'
#' @noRd
check_exclusive_aterms_supplied <- function(resp, av) {
  grps <- exclusive_aterm_groups(resp$family[["exclusive_aterms"]])
  if (!length(grps)) return(invisible(NULL))
  have <- unique(c(names(av), names(resp$aterms)))
  for (g in grps) {
    hit <- g[g %in% have]
    if (length(hit) < 2L) next
    # Group order is the declaration's own precedence, so the name kept
    # is the one the density actually reads.
    keep <- hit[[1L]]
    drop <- hit[-1L]
    frm_stop(resp$family[["family"]], ": ",
             comma_and(paste0("`", hit, "`")),
             " are spellings of the same datum for this family, and this ",
             "response supplies ",
             if (length(hit) > 2L) "all of them" else "both",
             ". The density reads `", keep, "` and would ignore ",
             comma_and(paste0("`", drop, "`")),
             ", so they could disagree with nothing to say so. Keep ",
             aterm_spelling(keep), " and drop ",
             comma_and(vapply(drop, aterm_spelling, "")),
             ".", call. = FALSE,
             package = frm_family_package(resp$family))
  }
  invisible(NULL)
}
