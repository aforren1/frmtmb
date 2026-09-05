# Helpers for the post-fit method tier (test-brms-methods.R).
#
# The log-density tier next door proves that frmtmb's objective IS the
# Stan program's log density at a point. That says nothing about what
# the two packages RETURN from fitted(), predict(), ranef() and the
# rest, and nothing in the repository compared those against brms until
# this file. See dev/brms-methods-tests.md.
#
# The mechanism is a brms fit whose draws ARE frmtmb's estimates, so
# that every brms method is evaluated at exactly the parameter vector
# frmtmb converged to and any difference in a returned value is a
# difference in the METHOD, never in the estimate. Stan's Fixed_param
# algorithm returns its initial values as the draws, and the tier's
# translator, stan_pars_from_fit(), already produces the constrained
# Stan-named list that init takes.
#
# Everything here calls the tier's helper-brms.R rather than
# reimplementing it: the flat prior set, the compile cache and the
# translator are shared, so a program compiled for the log-density tier
# is reused here for free.

# ---------------------------------------------------------------------
# The mechanism
# ---------------------------------------------------------------------

# Route (b) of two. See dev/brms-methods-tests.md for the measurement
# that chose it over route (a), brm(algorithm = "fixed_param"): route
# (a) compiles through brms's own machinery, which cannot see the
# tier's content-addressed cache, so it pays a fresh compile for a
# program the cache already holds.
#
# brms::rename_pars() is what brm() itself calls after fitting, and it
# is exported. brm(empty = TRUE) supplies the brmsfit scaffold with
# every design object built and nothing sampled, which is the same call
# the tier already makes for its ranef table. One brms internal is
# reached for, exclude_pars(), and only for fidelity: see below.
brms_fixed_fit <- function(bform, family, data, fit, ndraws = 10, ...) {
  prior <- brms_flat_prior(bform, data = data, family = family, ...)
  code <- brms::make_stancode(bform, data = data, family = family,
                              prior = prior, ...)
  sdat <- brms_standata(bform, data = data, family = family,
                        prior = prior, ...)
  scaffold <- suppressMessages(
    brms::brm(bform, data = data, family = family, prior = prior,
              empty = TRUE, ...))
  pars <- stan_pars_from_fit(fit, sdat, code, scaffold$ranef)
  mod <- brms_stan_model(code)
  # brm() passes its own exclusion set to rstan::sampling(), and
  # rename_pars() renames what survives. This is the one brms internal
  # the file reaches for, so it is guarded, and the fallback was
  # measured rather than assumed: without it rename_pars() does NOT
  # error, and on the sleepstudy shape the object carries 127 variables
  # instead of 47, the 80 extra being r_1 (36), z_1 (36), L_1 (4) and
  # Cor_1 (4) beside the renamed r_Subject[level,coef]. ranef() and
  # posterior_epred() come back identical() either way. So losing this
  # function costs fidelity to what a brm() fit looks like, and costs
  # no comparison in this file.
  exclude <- if (exists("exclude_pars", asNamespace("brms"),
                        inherits = FALSE)) {
    get("exclude_pars", asNamespace("brms"))(scaffold)
  } else {
    NULL
  }
  # warmup = 0 keeps every iteration: Fixed_param does not move the
  # parameters, so a warmup phase would only discard identical draws.
  # More than one draw is kept because several brms methods summarize
  # over draws and a single one degenerates their return shape.
  sf <- suppressMessages(
    rstan::sampling(mod, data = sdat, chains = 1, iter = ndraws,
                    warmup = 0, algorithm = "Fixed_param",
                    init = list(pars), refresh = 0, seed = 1,
                    pars = exclude, include = FALSE))
  scaffold$fit <- sf
  brms::rename_pars(scaffold)
}

# One brmsfit per (Stan program, standata) within a session. Sampling
# is under a second once the program is compiled, but brm(empty = TRUE)
# and make_stancode() are not free and a shape is used by a dozen test
# blocks.
.brms_fixed_fits <- new.env(parent = emptyenv())

# The key hashes the Stan program AND the standata, because one program
# serves every data variant of a shape (that is the compile cache's
# whole point) while a fitted object does not.
brms_fixed_key <- function(bform, family, data, ...) {
  prior <- brms_flat_prior(bform, data = data, family = family, ...)
  code <- brms::make_stancode(bform, data = data, family = family,
                              prior = prior, ...)
  sdat <- brms_standata(bform, data = data, family = family,
                        prior = prior, ...)
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(code, utils::capture.output(utils::str(sdat,
                                                      digits.d = 12))), f)
  unname(tools::md5sum(f))
}

brms_fixed_cached <- function(bform, family, data, fit, ...) {
  key <- brms_fixed_key(bform, family, data, ...)
  hit <- .brms_fixed_fits[[key]]
  if (!is.null(hit)) {
    return(hit)
  }
  out <- brms_fixed_fit(bform, family, data, fit, ...)
  .brms_fixed_fits[[key]] <- out
  out
}

# ---------------------------------------------------------------------
# The shape registry
# ---------------------------------------------------------------------
#
# One entry per row of the log-density tier's model matrix that its
# translator already handles. Each builds its own data with its own
# seed and returns the frmtmb fit beside the brms formula and family,
# so that a test block names a shape rather than repeating a fit.
#
# The data and the seeds are the tier's, deliberately: the same numbers
# under both tiers mean a divergence found here can be replayed there
# and the other way round. Row 3 is absent because a sibling lane is
# changing frmtmb's monotonic semantics this round.

# The shared data builders, so that four ordinal rows and three
# gaussian rows do not each repeat a seed.
brms_ord_data <- function() {
  set.seed(5)
  n <- 300
  do <- data.frame(x = rnorm(n))
  do$y <- ordered(cut(0.9 * do$x + rlogis(n),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  do$z <- rnorm(n)
  do
}

brms_ord_shape <- function(fam) {
  d <- brms_ord_data()
  bfam <- switch(fam,
                 cumulative = brms::cumulative(),
                 sratio = brms::sratio(),
                 cratio = brms::cratio(), acat = brms::acat())
  ffam <- switch(fam,
                 cumulative = frmtmb::cumulative(),
                 sratio = frmtmb::sratio(),
                 cratio = frmtmb::cratio(),
                 acat = frmtmb::acat())
  list(data = d, bform = brms::bf(y ~ x), family = bfam,
       fit = frm(frmtmb::bf(y ~ x) + ffam, data = d))
}

brms_cens_data <- function() {
  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
  dd$cens <- sample(c(0L, 1L), n, TRUE, prob = c(0.8, 0.2))
  dd$s <- runif(n, 0.5, 1.5)
  dd$w <- runif(n, 0.5, 2)
  dd
}

brms_zip_data <- function() {
  set.seed(13)
  n <- 300
  dz <- data.frame(x = rnorm(n))
  dz$y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.3 * dz$x)), 0L,
                 rpois(n, exp(0.6 + 0.4 * dz$x)))
  dz
}

# The roster's data, one family per call. The tier draws all four from
# ONE stream in order, so each family's y depends on the ones before
# it; that order is reproduced here rather than re-seeded per family.
brms_fam_data <- function(which) {
  set.seed(17)
  n <- 250
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.6 * d$x))
  if (identical(which, "poisson")) {
    return(d)
  }
  d$y <- rgamma(n, shape = 3, rate = 3 / exp(0.7 + 0.4 * d$x))
  if (identical(which, "gamma")) {
    return(d)
  }
  d$y <- rnbinom(n, size = 2, mu = exp(0.5 + 0.5 * d$x))
  if (identical(which, "negbinomial")) {
    return(d)
  }
  d$y <- rbinom(n, 1, plogis(0.2 + 0.8 * d$x))
  d
}

brms_methods_shapes <- list(
  # row 1: distributional gaussian
  r1 = function() {
    set.seed(11)
    n <- 150
    dd <- data.frame(x = rnorm(n), z = rnorm(n))
    dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
      rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
    list(data = dd, bform = brms::bf(y ~ x + z, sigma ~ x),
         family = gaussian(),
         fit = frm(frmtmb::bf(y ~ x + z, sigma ~ x) + gaussian(),
                   data = dd))
  },
  # row 2: monotonic effects
  r2 = function() {
    set.seed(3)
    dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
    dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
    list(data = dm, bform = brms::bf(y ~ mo(inc) + z),
         family = gaussian(),
         fit = frm(frmtmb::bf(y ~ mo(inc) + z) + gaussian(), data = dm))
  },
  # row 5: nonlinear
  r5 = function() {
    set.seed(7)
    n <- 120
    dn <- data.frame(x = runif(n, 0, 3))
    dn$y <- 2.5 * exp(-0.8 * dn$x) + rnorm(n, 0, 0.15)
    list(data = dn,
         bform = brms::bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE),
         family = gaussian(),
         fit = frm(frmtmb::bf(y ~ a * exp(-b * x), a + b ~ 1,
                              nl = TRUE) + gaussian(), data = dn))
  },
  # row 7: (1 | q | g) merged across mu and sigma
  r7 = function() {
    set.seed(43)
    ng <- 30
    nper <- 12
    n <- ng * nper
    gg <- factor(rep(seq_len(ng), each = nper))
    u <- MASS::mvrnorm(ng, c(0, 0), matrix(c(0.6, 0.25, 0.25, 0.3), 2))
    dd <- data.frame(x = rnorm(n), g = gg)
    dd$y <- 1 + 0.8 * dd$x + u[gg, 1] +
      rnorm(n, 0, exp(0.2 + 0.3 * dd$x + u[gg, 2]))
    list(data = dd,
         bform = brms::bf(y ~ x + (1 | q | g), sigma ~ x + (1 | q | g)),
         family = gaussian(),
         fit = frm(frmtmb::bf(y ~ x + (1 | q | g),
                              sigma ~ x + (1 | q | g)) + gaussian(),
                   data = dd))
  },
  # row 12: the four ordinal families, and cs()
  r12a = function() brms_ord_shape("cumulative"),
  r12b = function() brms_ord_shape("sratio"),
  r12c = function() brms_ord_shape("cratio"),
  r12d = function() brms_ord_shape("acat"),
  r12e = function() {
    d <- brms_ord_data()
    list(data = d, bform = brms::bf(y ~ x + cs(z)),
         family = brms::sratio(),
         fit = frm(frmtmb::bf(y ~ x + cs(z)) + frmtmb::sratio(), data = d))
  },
  # row 13: categorical
  r13 = function() {
    set.seed(31)
    n <- 300
    dc <- data.frame(x = rnorm(n))
    p <- cbind(1, exp(0.3 + 0.5 * dc$x), exp(-0.2 + 0.9 * dc$x))
    p <- p / rowSums(p)
    dc$y <- factor(apply(p, 1, function(pr) sample(1:3, 1, prob = pr)))
    list(data = dc, bform = brms::bf(y ~ x),
         family = brms::categorical(),
         fit = frm(frmtmb::bf(y ~ x) + frmtmb::categorical(), data = dc))
  },
  # row 14: censoring, and known standard errors
  r14a = function() {
    d <- brms_cens_data()
    list(data = d, bform = brms::bf(y | cens(cens) ~ x),
         family = gaussian(),
         fit = frm(frmtmb::bf(y | cens(cens) ~ x) + gaussian(),
                   data = d))
  },
  r14c = function() {
    d <- brms_cens_data()
    list(data = d, bform = brms::bf(y | se(s) ~ x), family = gaussian(),
         fit = frm(frmtmb::bf(y | se(s) ~ x) + gaussian(), data = d))
  },
  # row 15: binomial with trials()
  r15 = function() {
    set.seed(9)
    n <- 200
    db <- data.frame(x = rnorm(n), n = sample(3:12, n, TRUE))
    db$y <- rbinom(n, db$n, plogis(0.3 + 0.7 * db$x))
    list(data = db, bform = brms::bf(y | trials(n) ~ x),
         family = binomial(),
         fit = frm(frmtmb::bf(y | trials(n) ~ x) + binomial(),
                   data = db))
  },
  # row 16: zero-inflated poisson with a predictor on zi
  r16 = function() {
    d <- brms_zip_data()
    list(data = d, bform = brms::bf(y ~ x, zi ~ x),
         family = brms::zero_inflated_poisson(),
         fit = frm(frmtmb::bf(y ~ x, zi ~ x) +
                     frmtmb::zero_inflated_poisson(), data = d))
  },
  # row 17: two-component gaussian mixture
  r17 = function() {
    set.seed(37)
    n <- 400
    dx <- data.frame(x = rnorm(n))
    k <- rbinom(n, 1, 0.35)
    dx$y <- ifelse(k == 1, rnorm(n, 3, 1), rnorm(n, -1, 1))
    list(data = dx, bform = brms::bf(y ~ 1, theta1 ~ x),
         family = brms::mixture(gaussian(), gaussian()),
         fit = frm(frmtmb::bf(y ~ 1, theta1 ~ x) +
                     frmtmb::mixture(gaussian(), gaussian()), data = dx))
  },
  # row 20: case weights
  r20 = function() {
    d <- brms_cens_data()
    list(data = d, bform = brms::bf(y | weights(w) ~ x),
         family = gaussian(),
         fit = frm(frmtmb::bf(y | weights(w) ~ x) + gaussian(),
                   data = d))
  },
  # row 21: the family roster, one predictor each
  r21a = function() {
    d <- brms_fam_data("poisson")
    list(data = d, bform = brms::bf(y ~ x), family = poisson(),
         fit = frm(frmtmb::bf(y ~ x) + poisson(), data = d))
  },
  r21b = function() {
    d <- brms_fam_data("gamma")
    list(data = d, bform = brms::bf(y ~ x), family = Gamma(link = "log"),
         fit = frm(frmtmb::bf(y ~ x) + Gamma(link = "log"), data = d))
  },
  r21c = function() {
    d <- brms_fam_data("negbinomial")
    list(data = d, bform = brms::bf(y ~ x), family = brms::negbinomial(),
         fit = frm(frmtmb::bf(y ~ x) + frmtmb::negbinomial(), data = d))
  },
  r21d = function() {
    d <- brms_fam_data("bernoulli")
    list(data = d, bform = brms::bf(y ~ x), family = brms::bernoulli(),
         fit = frm(frmtmb::bf(y ~ x) + frmtmb::bernoulli(), data = d))
  },
  # Not a row of the log-density tier's matrix. Every shape above has
  # numeric covariates only, and conditional_effects() holds a numeric
  # covariate at its mean and a factor at something else, so without a
  # factor in the population-level formula half of that semantics is
  # unreachable. Structurally this is row 1 with an interaction, so the
  # translator needs nothing new.
  rfac = function() {
    set.seed(23)
    n <- 240
    df <- data.frame(x = rnorm(n),
                     f = factor(rep(c("a", "b", "c"), length.out = n)))
    df$y <- 0.5 + 0.8 * df$x + c(a = 0, b = 0.7, c = -0.4)[df$f] +
      0.3 * df$x * (df$f == "c") + rnorm(n, 0, 0.6)
    list(data = df, bform = brms::bf(y ~ x * f), family = gaussian(),
         fit = frm(frmtmb::bf(y ~ x * f) + gaussian(), data = df))
  },
  # the random-effect variants, where brms's draws carry frmtmb's
  # conditional modes through the non-centered z block
  rC0 = function() {
    e <- new.env()
    utils::data(sleepstudy, package = "lme4", envir = e)
    ss <- e$sleepstudy
    list(data = ss,
         bform = brms::bf(Reaction ~ Days + (Days | Subject),
                          sigma ~ Days),
         family = gaussian(),
         fit = frm(frmtmb::bf(Reaction ~ Days + (Days | Subject),
                              sigma ~ Days) + gaussian(), data = ss))
  },
  rC16 = function() {
    set.seed(13)
    n <- 300
    dz <- data.frame(x = rnorm(n), g = factor(rep(1:20, 15)))
    dz$y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.3 * dz$x)), 0L,
                   rpois(n, exp(0.6 + 0.4 * dz$x)))
    list(data = dz, bform = brms::bf(y ~ x + (1 | g), zi ~ x),
         family = brms::zero_inflated_poisson(),
         fit = frm(frmtmb::bf(y ~ x + (1 | g), zi ~ x) +
                     frmtmb::zero_inflated_poisson(), data = dz))
  }
)

# Shapes built at most once per session. A shape costs a frmtmb fit and
# a brms scaffold; the compiled program behind it is cached on disk by
# the log-density tier's brms_stan_model().
.brms_shape_cache <- new.env(parent = emptyenv())

brms_shape <- function(name) {
  hit <- .brms_shape_cache[[name]]
  if (!is.null(hit)) {
    return(hit)
  }
  builder <- brms_methods_shapes[[name]]
  if (is.null(builder)) {
    stop("no method-tier shape is registered under the name ", name)
  }
  s <- builder()
  s$name <- name
  s$brmsfit <- brms_fixed_cached(s$bform, s$family, s$data, s$fit)
  .brms_shape_cache[[name]] <- s
  s
}

# A second fit of the same shape at a different draw count. Two of the
# divergences below are Monte Carlo effects on brms's side, and the way
# to tell a Monte Carlo effect from a real one is to change the draw
# count and watch the gap move.
.brms_ndraw_fits <- new.env(parent = emptyenv())

brms_fixed_cached_n <- function(shape, ndraws) {
  key <- paste0(shape$name, "#", ndraws)
  hit <- .brms_ndraw_fits[[key]]
  if (!is.null(hit)) {
    return(hit)
  }
  out <- brms_fixed_fit(shape$bform, shape$family, shape$data,
                        shape$fit, ndraws = ndraws)
  .brms_ndraw_fits[[key]] <- out
  out
}

# ---------------------------------------------------------------------
# Reaching values frmtmb has no accessor for
# ---------------------------------------------------------------------

# frmtmb core has no log_lik() on a fit: the pointwise log density is
# frmtmb.sample's draws_row_loglik(), and core exports only the pieces
# it is built from. This composes them the same way, so the comparison
# against brms::log_lik() is available without adding an accessor under
# R/, which sibling lanes own this round. Recorded as a missing
# accessor in dev/brms-methods-tests.md.
frm_row_loglik <- function(fit, resp = NULL) {
  frame <- fit$frame
  rspecs <- fit$spec$responses
  dpv <- with_cs_offsets(fit, NULL, eval_dpars(fit))
  extra <- fit_extras(fit)
  n <- frame[["n_obs"]]
  use <- resp %||% names(rspecs)
  out <- numeric(n)
  for (r in use) {
    av <- frame[["aterm_values"]][[r]]
    fam <- rspecs[[r]][["family"]]
    yv <- frame[["y"]][[r]]
    ll <- row_lpdf(fam, yv, yv, dpv[[r]], av, extra)
    out <- out + (av[["weights"]] %||% 1) * as.numeric(ll)
  }
  out
}

# frmtmb's fixef() keyed the way brms names the same coefficients.
# frmtmb returns one named vector per linear predictor; brms flattens
# them into one matrix and puts the dpar or nlpar in FRONT of the
# coefficient name, with "Intercept" for model.matrix's "(Intercept)".
brms_flatten_fixef <- function(fit) {
  fe <- fixef(fit)
  out <- unlist(lapply(names(fe), function(dp) {
    v <- fe[[dp]]
    cn <- names(v)
    cn[cn == "(Intercept)"] <- "Intercept"
    stats::setNames(as.numeric(v),
                    if (identical(dp, "mu")) cn else paste0(dp, "_", cn))
  }))
  out
}

# frmtmb's column name for one group-level coefficient as brms's
# METHODS name it.
#
# The tier's brms_frm_coef() takes the dpar as a separate argument,
# because brms's ranef TABLE carries it in its own column. brms's
# ranef() and VarCorr() fold it into the coefficient instead, so a
# block merged across mu and sigma comes back as "Intercept" and
# "sigma_Intercept" where frmtmb has "y.mu:(Intercept)" and
# "y.sigma:(Intercept)". This splits the fold off and hands the two
# halves to the tier's mapper.
brms_re_coef_to_frm <- function(fit, cols, cn) {
  bare <- brms_coef_to_frm(cn)
  if (bare %in% cols) {
    return(bare)
  }
  for (d in names(family(fit)$links)) {
    pre <- paste0(d, "_")
    if (startsWith(cn, pre)) {
      hit <- try(brms_frm_coef(cols, substring(cn, nchar(pre) + 1L), d),
                 silent = TRUE)
      if (!inherits(hit, "try-error")) {
        return(hit)
      }
    }
  }
  brms_frm_coef(cols, cn, "")
}

# The dpars beyond mu that BOTH packages know about. frmtmb's family
# object lists every dpar the family has, including ones brms fixes
# rather than estimates (a cumulative model's disc), and asking brms
# about a parameter it did not declare is a question about brms's
# argument checking, not about the two packages agreeing.
brms_dpars_of <- function(shape) {
  dp <- setdiff(names(family(shape$fit)$links), "mu")
  # A dpar whose response scale diverges is routed out through the
  # exclusion table below rather than by a rule here, so that the guard
  # block can ask each one whether its defect is still live.
  dp <- setdiff(dp, sub("^[^:]*:", "",
                        grep(paste0("^", shape$name, ":"),
                             brms_excluded("brms_dpars_of"),
                             value = TRUE)))
  v <- brms::variables(shape$brmsfit)
  keep <- vapply(dp, function(d) {
    d %in% v || any(startsWith(v, paste0("b_", d, "_"))) ||
      any(startsWith(v, paste0("bsp_", d, "_")))
  }, logical(1))
  dp[keep]
}

# TRUE when brms declares the dpar as a bare scalar on its natural
# scale rather than building a linear predictor for it. That is the one
# case where posterior_linpred(dpar = ) and predict(type = "link",
# dpar = ) are on different scales, so it decides which comparison a
# dpar belongs in.
brms_dpar_is_scalar <- function(shape, dpar) {
  dpar %in% brms::variables(shape$brmsfit)
}

# ---------------------------------------------------------------------
# The exclusion table
# ---------------------------------------------------------------------
#
# Each agreement loop runs over a list of shapes, and a shape that
# diverges is kept OUT of that list so the loop asserts agreement where
# agreement is the claim. That routing is the one place in this tier
# where a fix could pass silently: fix the zero-inflated defect, leave
# the shapes out, and nothing fails while the tier quietly stops
# covering what it was built for.
#
# So the lists are DERIVED from one table, every row names the finding
# in dev/brms-methods-tests.md that put it there and that finding's
# class, and every row whose class is "D" carries a live-defect probe.
# The guard block in test-brms-methods.R runs those probes and fails the
# moment one starts agreeing, naming the list to edit.
#
#   D   a frmtmb defect. The exclusion is TEMPORARY, and the probe is
#       what makes it so.
#   P   a paradigm difference, right for a maximum-likelihood fit.
#   C   a design choice needing a user decision.
#
# P and C rows carry no probe: there is nothing to wait for, and a probe
# that can never flip is noise.
brms_exclusions <- function() {
  rows <- rbind(
    # conditional_effects, one-way panels
    c("brms_ce_shapes", "r2", "17", "D"),
    c("brms_ce_shapes", "r5", "18", "P"),
    c("brms_ce_shapes", "r12a", "6b", "D"),
    c("brms_ce_shapes", "r12b", "6b", "D"),
    c("brms_ce_shapes", "r12c", "6b", "D"),
    c("brms_ce_shapes", "r12d", "6b", "D"),
    c("brms_ce_shapes", "r12e", "6b", "D"),
    c("brms_ce_shapes", "r13", "6b", "D"),
    c("brms_ce_shapes", "r16", "1b", "D"),
    c("brms_ce_shapes", "rC16", "1b", "D"),
    c("brms_ce_shapes", "r17", "1d", "C"),
    # posterior_epred(dpar = ), keyed shape:dpar
    c("brms_dpars_of", "r17:theta1", "1c", "D"),
    c("brms_dpars_of", "r14c:sigma", "15", "D"),
    # posterior_linpred() against predict(type = "link")
    c("brms_linpred_shapes", "r12e", "13", "C"),
    c("brms_linpred_shapes", "r13", "13", "C"),
    c("brms_linpred_shapes", "r17", "13", "C"),
    # posterior_linpred(transform = TRUE) against the response scale
    c("brms_meanlink_shapes", "r15", "14", "P"),
    c("brms_meanlink_shapes", "r16", "14", "P"),
    c("brms_meanlink_shapes", "rC16", "14", "P"),
    c("brms_meanlink_shapes", "r12a", "14", "P"),
    c("brms_meanlink_shapes", "r12b", "14", "P"),
    c("brms_meanlink_shapes", "r12c", "14", "P"),
    c("brms_meanlink_shapes", "r12d", "14", "P"),
    c("brms_meanlink_shapes", "r12e", "14", "P"),
    c("brms_meanlink_shapes", "r13", "14", "P"),
    c("brms_meanlink_shapes", "r17", "14", "P"),
    # residuals: a nominal response has no y - E[Y] to form, and the
    # ordinal families score the categories, which brms does not, so
    # there is no shared definition to compare rather than a defect
    c("brms_resid_shapes", "r13", "omissions", "P"),
    c("brms_resid_shapes", "r12a", "omissions", "P"),
    c("brms_resid_shapes", "r12b", "omissions", "P"),
    c("brms_resid_shapes", "r12c", "omissions", "P"),
    c("brms_resid_shapes", "r12d", "omissions", "P"),
    c("brms_resid_shapes", "r12e", "omissions", "P")
  )
  data.frame(list = rows[, 1], key = rows[, 2], finding = rows[, 3],
             class = rows[, 4], stringsAsFactors = FALSE)
}

brms_excluded <- function(list_nm) {
  ex <- brms_exclusions()
  ex$key[ex$list == list_nm]
}

# TRUE when the divergence that justifies one exclusion row has gone
# away, which is the signal to put the shape back into its list. NA for
# a row with no probe.
brms_exclusion_agrees <- function(list_nm, key) {
  if (identical(list_nm, "brms_ce_shapes")) {
    return(brms_ce_agrees(brms_shape(key)))
  }
  if (identical(list_nm, "brms_dpars_of")) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1]]
    return(brms_dpar_epred_agrees(brms_shape(parts[[1]]), parts[[2]]))
  }
  NA
}

# The one-way conditional_effects comparison the loop performs, as a
# predicate rather than as assertions, so an exclusion can be asked
# whether it is still earning its place.
brms_ce_agrees <- function(shape, tol = 1e-8) {
  cb <- try(suppressWarnings(brms::conditional_effects(shape$brmsfit)),
            silent = TRUE)
  cf <- try(suppressWarnings(conditional_effects(shape$fit)),
            silent = TRUE)
  if (inherits(cb, "try-error") || inherits(cf, "try-error")) {
    return(FALSE)
  }
  if (!identical(names(cb), names(cf))) {
    return(FALSE)
  }
  for (e in names(cb)) {
    if (grepl(":", e, fixed = TRUE)) next
    if (!identical(nrow(cb[[e]]), nrow(cf[[e]]))) {
      return(FALSE)
    }
    a <- cb[[e]][["estimate__"]]
    b <- cf[[e]][["estimate__"]]
    if (max(abs(a - b) / pmax(1, pmax(abs(a), abs(b)))) > tol) {
      return(FALSE)
    }
  }
  TRUE
}

brms_dpar_epred_agrees <- function(shape, dpar, tol = 1e-8) {
  a <- try(brms::posterior_epred(shape$brmsfit, dpar = dpar),
           silent = TRUE)
  b <- try(predict(shape$fit, type = "response", dpar = dpar),
           silent = TRUE)
  if (inherits(a, "try-error") || inherits(b, "try-error")) {
    return(FALSE)
  }
  a <- as.numeric(if (length(dim(a)) == 3L) a[1, , ] else a[1, ])
  b <- as.numeric(as.matrix(b))
  if (length(a) != length(b)) {
    return(FALSE)
  }
  max(abs(a - b) / pmax(1, pmax(abs(a), abs(b)))) <= tol
}

# The shapes each comparison is defined on, derived from the table
# above so that a list and its justification cannot drift apart.

brms_resid_shapes <- function() {
  setdiff(names(brms_methods_shapes), brms_excluded("brms_resid_shapes"))
}

# The shapes with group-level effects. brms's draws carry frmtmb's
# conditional modes through the non-centered z block, so every
# CONDITIONAL quantity agrees and the marginal ones do not. Not an
# exclusion: it is a property of the models, not a deferred defect.
brms_re_shapes <- function() {
  c("r7", "rC0", "rC16")
}

# posterior_linpred(): the shapes whose mu predictor is ONE column per
# observation, so brms's draws x N array and frmtmb's N-vector are the
# same object.
brms_linpred_shapes <- function() {
  setdiff(names(brms_methods_shapes),
          brms_excluded("brms_linpred_shapes"))
}

# The shapes whose MEAN is the inverse link of the mu linear predictor,
# so brms's posterior_linpred(transform = TRUE) and frmtmb's
# predict(type = "response") are the same quantity.
brms_meanlink_shapes <- function() {
  setdiff(names(brms_methods_shapes),
          brms_excluded("brms_meanlink_shapes"))
}

# conditional_effects: the shapes whose one-way panels agree exactly.
brms_ce_shapes <- function() {
  setdiff(names(brms_methods_shapes), brms_excluded("brms_ce_shapes"))
}

# ---------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------

# The tier's claim is exactness, so the tolerance is relative and tiny.
# A value comparison that needs anything looser is a finding, not a
# reason to widen this.
expect_exact_num <- function(x, y, tol = 1e-8, label = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  testthat::expect_identical(length(x), length(y))
  scale <- pmax(1, pmax(abs(x), abs(y)))
  testthat::expect_lt(max(abs(x - y) / scale), tol,
                      label = label %||% "max relative difference")
}

# Every draw identical is the mechanism's own precondition, and a
# method comparison that ran without it would be comparing a Monte
# Carlo mean to a point estimate. Asserted per shape rather than
# assumed.
expect_draws_degenerate <- function(x, tol = 0) {
  m <- posterior::as_draws_matrix(x)
  spread <- apply(m, 2, function(v) max(v) - min(v))
  testthat::expect_lte(max(spread), tol)
}
