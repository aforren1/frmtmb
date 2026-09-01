# Grammar-aware pairwise fuzzer for the model grammar.
#
# The package's worst historical bugs were two-feature interactions
# (rescor x cens, profile x bounds, trunc x post-fit surface). A single
# example per feature never reaches them, and a full cross of the
# grammar is far too large to fit. A pairwise covering array is the
# cheap middle: every pair of feature values appears in at least one
# generated model, so any bug that needs only two features to show
# itself is reachable by construction.
#
# Nothing here runs during a normal test run. test-fuzz.R gates the
# whole tier behind FRMTMB_FUZZ=true; this file only defines functions.
#
# Entry points:
#   fuzz_plan(seed, size)   -> data frame of specs (covering array)
#   fuzz_run(plan, ...)     -> list of findings
#   fuzz_repro(spec)        -> runnable code for one spec, as text

# ---------------------------------------------------------------------
# 0. small utilities
# ---------------------------------------------------------------------

fuzz_or <- function(x, y) if (is.null(x)) y else x

# Evaluate once, keeping errors and warnings as data rather than as
# control flow: an invariant needs to know both that a fit warned and
# what it produced.
fuzz_try <- function(expr) {
  warns <- character(0)
  err <- NULL
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) { err <<- e; NULL }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  list(value = val, error = err, warnings = warns,
       ok = is.null(err))
}

# Messages that mean "an internal index or dimension went wrong", as
# opposed to "you asked for something this package refuses to do".
FUZZ_CRASH_PATTERNS <- c(
  "subscript out of bounds",
  "non-conformable",
  "NA/NaN/Inf",
  "missing value where TRUE/FALSE needed",
  "argument is of length zero",
  "undefined columns selected",
  "attempt to select (less|more) than one element",
  "invalid 'times' argument",
  "incorrect number of dimensions",
  "arguments imply differing number of rows",
  "object of type 'closure' is not subsettable",
  "invalid subscript",
  "^object '[^']*' not found$",
  "dims \\[product",
  "invalid argument to unary operator"
)

fuzz_is_crash <- function(msg) {
  any(vapply(FUZZ_CRASH_PATTERNS, function(p) grepl(p, msg), TRUE))
}

# ---------------------------------------------------------------------
# 1. feature registry (generation only)
# ---------------------------------------------------------------------
#
# Deliberately self-contained and minimal: this describes what the
# generator may emit and how to simulate a response for it, not what
# the package supports in general. A user-facing compatibility registry
# is a separate concern.
#
# Per family:
#   ctor        our family constructor, as source text
#   brms        the brms family expression, as source text (NULL = no
#               brms equivalent to compare against)
#   dpar2       second dpar the generator may give a formula to
#   reml        REML is meaningful for this family (continuous)
#   sim         simulate() is implemented
#   mean_check  mean(simulate) vs mean(predict) is the natural check
#   trials/cens/cens_chr/trunc/se  aterm applicability
#   eta0        intercept keeping the response in support at n ~ 100
#   trunc_lb    lower truncation bound used by the trunc aterm

fuzz_families <- list(
  gaussian = list(
    ctor = "gaussian()", brms = "gaussian()", dpar2 = "sigma",
    reml = TRUE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = TRUE, trunc = TRUE, se = TRUE,
    eta0 = 1, trunc_lb = -0.6
  ),
  student = list(
    ctor = "student()", brms = "brms::student()", dpar2 = "sigma",
    reml = TRUE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = FALSE, trunc = FALSE, se = TRUE,
    eta0 = 1, trunc_lb = -0.6
  ),
  lognormal = list(
    ctor = "lognormal()", brms = "brms::lognormal()", dpar2 = "sigma",
    reml = TRUE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = TRUE, trunc = TRUE, se = FALSE,
    eta0 = 0.5, trunc_lb = 0.9
  ),
  weibull = list(
    ctor = "weibull()", brms = "brms::weibull()", dpar2 = "shape",
    reml = TRUE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = TRUE, trunc = TRUE, se = FALSE,
    eta0 = 0.6, trunc_lb = 0.5
  ),
  Gamma = list(
    ctor = "Gamma(link = \"log\")", brms = "Gamma(link = \"log\")",
    dpar2 = "shape",
    reml = TRUE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = FALSE, trunc = FALSE, se = FALSE,
    eta0 = 1, trunc_lb = 0.3
  ),
  Beta = list(
    ctor = "Beta()", brms = "brms::Beta()", dpar2 = "phi",
    reml = TRUE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = FALSE, trunc = FALSE, se = FALSE,
    eta0 = 0, trunc_lb = 0.02
  ),
  poisson = list(
    ctor = "poisson()", brms = "poisson()", dpar2 = NULL,
    reml = FALSE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = FALSE, trunc = TRUE, se = FALSE,
    eta0 = 1, trunc_lb = 1
  ),
  binomial = list(
    ctor = "binomial()", brms = "binomial()", dpar2 = NULL,
    reml = FALSE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = TRUE, cens = FALSE, trunc = FALSE, se = FALSE,
    eta0 = 0, trunc_lb = 0.5
  ),
  negbinomial = list(
    ctor = "negbinomial()", brms = "brms::negbinomial()",
    dpar2 = "shape",
    reml = FALSE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = FALSE, trunc = FALSE, se = FALSE,
    eta0 = 1, trunc_lb = 1
  ),
  zero_inflated_poisson = list(
    ctor = "zero_inflated_poisson()",
    brms = "brms::zero_inflated_poisson()", dpar2 = "zi",
    reml = FALSE, sim = TRUE, mean_check = TRUE, ordinal = FALSE,
    trials = FALSE, cens = FALSE, trunc = FALSE, se = FALSE,
    eta0 = 1.2, trunc_lb = 1
  ),
  cumulative = list(
    ctor = "cumulative()", brms = "brms::cumulative()", dpar2 = NULL,
    reml = FALSE, sim = FALSE, mean_check = FALSE, ordinal = TRUE,
    trials = FALSE, cens = FALSE, trunc = FALSE, se = FALSE,
    eta0 = 0, trunc_lb = 1
  )
)

# Random-effect structures. `brms` is the same text when brms parses it
# the same way; NULL means the term uses a covstruct wrapper brms has
# no spelling for, so the brms oracle sits out.
fuzz_re <- list(
  none    = list(term = NULL,                brms = NULL,  blocks = 0),
  ri      = list(term = "(1 | g)",           brms = "(1 | g)",      blocks = 1,
                 scalar = TRUE),
  rs      = list(term = "(1 + x | g)",       brms = "(1 + x | g)",  blocks = 1),
  dbar    = list(term = "(1 + x || g)",      brms = "(1 + x || g)", blocks = 2),
  dbar_f  = list(term = "(fac || g)",        brms = "(fac || g)",   blocks = 3),
  diagbar = list(term = "diag(1 + x | g)",   brms = NULL,           blocks = 1),
  ar1     = list(term = "ar1(tt + 0 | g)",   brms = NULL,           blocks = 1),
  nested  = list(term = "(1 | ga/gb)",       brms = "(1 | ga/gb)",  blocks = 2,
                 scalar = TRUE),
  ri_fx   = list(term = "(1 | factor(gc))",  brms = "(1 | factor(gc))",
                 blocks = 1, scalar = TRUE)
)

fuzz_special <- list(
  none   = list(term = NULL,             brms = NULL,             blk = 0),
  smooth = list(term = "s(xs)",          brms = "s(xs)",          blk = 1),
  mo     = list(term = "mo(xo)",         brms = "mo(xo)",         blk = 0),
  # mo() crossed with a numeric multiplier: the interaction machinery
  # without the factor spelling, which the package refuses (probed
  # separately in fuzz_refusal_cases)
  mo_int = list(term = "mo(xo) * z",     brms = "mo(xo) * z",     blk = 0),
  gp     = list(term = "gp(xs, k = 10)", brms = "gp(xs, k = 10)", blk = 1),
  # refusal probe only (never generated): mo() crossed with a factor
  mo_fac = list(term = "mo(xo) * fac",   brms = "mo(xo) * fac",   blk = 0)
)

fuzz_aterm <- list(
  none     = list(term = NULL),
  weights  = list(term = "weights(w)"),
  trials   = list(term = "trials(nt)"),
  cens     = list(term = "cens(cc)"),
  cens_chr = list(term = "cens(ccl)"),
  trunc    = list(term = NULL),   # bound is family specific, built later
  se       = list(term = "se(sdy)"),
  # not a generated value: the unit-weights metamorphic refit reuses
  # the formula builder through this entry
  unit_w   = list(term = "weights(one)")
)

fuzz_modes <- c("ml", "reml", "quadrature", "profile", "autoscale",
                "sparse_x")
fuzz_ops <- c("predict", "simulate", "confint", "vcov")

fuzz_dims <- list(
  family  = names(fuzz_families),
  aterm   = setdiff(names(fuzz_aterm), "unit_w"),
  re      = names(fuzz_re),
  special = setdiff(names(fuzz_special), "mo_fac"),
  dpar    = c("none", "dpar_x"),
  mode    = fuzz_modes,
  op      = fuzz_ops
)

# Applicability. The generator only emits combinations the package
# CLAIMS to support; combinations it refuses by design are sampled
# separately (fuzz_refusals()) so the refusal itself gets tested.
fuzz_ok <- function(a) {
  fm <- fuzz_families[[a$family]]
  if (a$aterm == "trials"   && !fm$trials) return(FALSE)
  if (a$aterm == "cens"     && !fm$cens)   return(FALSE)
  if (a$aterm == "cens_chr" && !fm$cens)   return(FALSE)
  if (a$aterm == "trunc"    && !fm$trunc)  return(FALSE)
  if (a$aterm == "se"       && !fm$se)     return(FALSE)
  if (a$dpar == "dpar_x") {
    if (is.null(fm$dpar2)) return(FALSE)
    # se() without sigma = TRUE maps sigma out; a sigma formula would
    # put an unidentified parameter back, which is a user error, not a
    # package bug
    if (a$aterm == "se") return(FALSE)
  }
  if (a$mode == "reml" && !fm$reml) return(FALSE)
  if (a$mode == "quadrature") {
    if (!isTRUE(fuzz_re[[a$re]]$scalar)) return(FALSE)
    if (a$special != "none") return(FALSE)
  }
  if (a$op == "simulate" && !fm$sim) return(FALSE)
  TRUE
}

# ---------------------------------------------------------------------
# 2. pairwise covering plan
# ---------------------------------------------------------------------
#
# Greedy set cover over the feasible (dim=value, dim=value) pairs. The
# candidate pool is the complete set of applicable assignments, so the
# cover is exact rather than approximate in its feasibility. Once every
# pair is covered the plan is padded by greedy 3-way coverage, which
# spends the remaining budget where it can still buy something.

fuzz_all_valid <- function() {
  grid <- expand.grid(fuzz_dims, stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)
  keep <- vapply(seq_len(nrow(grid)), function(i) fuzz_ok(grid[i, ]), TRUE)
  grid[keep, , drop = FALSE]
}

# Global id for every (dimension, value): pair ids are then a pair of
# small integers and the covered set is a plain logical vector.
fuzz_value_ids <- function() {
  ids <- list()
  k <- 0L
  for (d in names(fuzz_dims)) {
    v <- fuzz_dims[[d]]
    ids[[d]] <- stats::setNames(seq_along(v) + k, v)
    k <- k + length(v)
  }
  ids
}

fuzz_combo_matrix <- function(grid, m) {
  ids <- fuzz_value_ids()
  M <- vapply(names(fuzz_dims), function(d) unname(ids[[d]][grid[[d]]]),
              integer(nrow(grid)))
  nd <- length(fuzz_dims)
  cmb <- utils::combn(nd, m)
  V <- max(unlist(ids))
  P <- matrix(0L, nrow(grid), ncol(cmb))
  for (j in seq_len(ncol(cmb))) {
    k <- cmb[, j]
    key <- M[, k[1]]
    for (i in k[-1]) key <- key * (V + 1L) + M[, i]
    P[, j] <- key
  }
  # dense-ify the keys so `covered` stays a small logical vector
  u <- sort(unique(as.vector(P)))
  matrix(match(P, u), nrow(P), ncol(P))
}

#' Build the pairwise covering plan.
#' @param seed RNG seed; the whole plan and every datum derive from it.
#' @param size Target number of specs (>= the size of the pair cover).
#' @return A data frame with one row per spec: the dimension values,
#'   a `seed`, and a `kind` ("cover", "pad", or "refusal").
fuzz_plan <- function(seed = 20260901L, size = 300L) {
  grid <- fuzz_all_valid()
  set.seed(seed)
  grid <- grid[sample.int(nrow(grid)), , drop = FALSE]
  rownames(grid) <- NULL

  P2 <- fuzz_combo_matrix(grid, 2L)
  chosen <- integer(0)
  covered <- logical(max(P2))
  repeat {
    gain <- rowSums(matrix(!covered[P2], nrow(P2)))
    b <- which.max(gain)
    if (gain[b] == 0L) break
    chosen <- c(chosen, b)
    covered[P2[b, ]] <- TRUE
  }
  n_cover <- length(chosen)

  if (size > n_cover) {
    P3 <- fuzz_combo_matrix(grid, 3L)
    cov3 <- logical(max(P3))
    cov3[as.vector(P3[chosen, , drop = FALSE])] <- TRUE
    while (length(chosen) < size) {
      gain <- rowSums(matrix(!cov3[P3], nrow(P3)))
      gain[chosen] <- -1L
      b <- which.max(gain)
      if (gain[b] <= 0L) {
        rest <- setdiff(seq_len(nrow(grid)), chosen)
        if (!length(rest)) break
        b <- rest[seq_len(min(length(rest), size - length(chosen)))]
      }
      chosen <- c(chosen, b)
      cov3[as.vector(P3[b, , drop = FALSE])] <- TRUE
    }
    chosen <- chosen[seq_len(min(length(chosen), size))]
  }

  out <- grid[chosen, , drop = FALSE]
  out$kind <- c(rep("cover", n_cover),
                rep("pad", nrow(out) - n_cover))[seq_len(nrow(out))]
  out <- rbind(out, fuzz_refusals())
  out$seed <- seed + seq_len(nrow(out)) * 977L
  out$id <- seq_len(nrow(out))
  rownames(out) <- NULL
  attr(out, "n_pairs") <- length(covered)
  attr(out, "n_pairs_covered") <- sum(covered)
  attr(out, "n_cover_rows") <- n_cover
  out
}

# Deliberately refused-by-design combinations. Each must produce a
# clean, informative error - never a crash from inside an index.
fuzz_refusal_cases <- list(
  list(name = "cens_no_cdf", expect = "CDF",
       a = list(family = "Gamma", aterm = "cens", re = "ri",
                special = "none", dpar = "none", mode = "ml",
                op = "predict")),
  list(name = "trunc_no_cdf", expect = "CDF",
       a = list(family = "Beta", aterm = "trunc", re = "none",
                special = "none", dpar = "none", mode = "ml",
                op = "predict")),
  list(name = "cens_discrete", expect = "discrete",
       a = list(family = "poisson", aterm = "cens", re = "ri",
                special = "none", dpar = "none", mode = "ml",
                op = "predict")),
  list(name = "quad_smooth", expect = "scalar random",
       a = list(family = "gaussian", aterm = "none", re = "ri",
                special = "smooth", dpar = "none", mode = "quadrature",
                op = "predict")),
  list(name = "quad_slope", expect = "scalar random",
       a = list(family = "gaussian", aterm = "none", re = "rs",
                special = "none", dpar = "none", mode = "quadrature",
                op = "predict")),
  list(name = "quad_reml", expect = "REML",
       a = list(family = "gaussian", aterm = "none", re = "ri",
                special = "none", dpar = "none", mode = "quad_reml",
                op = "predict")),
  list(name = "profile_reml", expect = "REML",
       a = list(family = "gaussian", aterm = "none", re = "ri",
                special = "none", dpar = "none", mode = "profile_reml",
                op = "predict")),
  list(name = "mo_factor", expect = "numeric multiplier",
       a = list(family = "gaussian", aterm = "none", re = "none",
                special = "mo_fac", dpar = "none", mode = "ml",
                op = "predict")),
  list(name = "bad_dpar", expect = "not available",
       a = list(family = "poisson", aterm = "none", re = "none",
                special = "none", dpar = "bad_dpar", mode = "ml",
                op = "predict")),
  list(name = "bar_crossed", expect = "term",
       a = list(family = "gaussian", aterm = "none", re = "crossed_bar",
                special = "none", dpar = "none", mode = "ml",
                op = "predict"))
)

fuzz_refusals <- function() {
  rows <- lapply(fuzz_refusal_cases, function(cs) {
    r <- as.data.frame(cs$a, stringsAsFactors = FALSE)
    r$kind <- paste0("refusal:", cs$name)
    r
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 3. data generation - simulate FROM the sampled model
# ---------------------------------------------------------------------

# Balanced group x time design so ar1() has a full time grid and the
# nested factor divides evenly. n stays in [60, 150].
fuzz_design <- function(sp) {
  rn <- (sp$seed %% 4L) + 1L
  ng <- c(12L, 15L, 20L, 24L)[rn]
  nt <- c(6L, 6L, 5L, 5L)[rn]
  n <- ng * nt
  g <- factor(rep(seq_len(ng), each = nt))
  tt <- factor(rep(seq_len(nt), ng))
  gb <- factor(rep(rep(c("i", "ii"), length.out = nt), ng))
  data.frame(
    g = g, tt = tt, ga = g, gb = gb,
    gc = as.character(g),
    fac = factor(rep(c("a", "b", "c"), length.out = n)),
    x = stats::rnorm(n), z = stats::rnorm(n),
    xs = stats::runif(n, 0, 10),
    xo = sample.int(4L, n, replace = TRUE) - 1L,
    nt_ = rep(10L, n),
    w = stats::runif(n, 0.5, 2),
    sdy = stats::runif(n, 0.2, 0.6),
    one = 1,
    stringsAsFactors = FALSE
  )
}

# Linear predictor of the sampled model, with the same structure the
# formula will fit: fixed x, the chosen special, the chosen RE.
fuzz_eta <- function(sp, d) {
  n <- nrow(d)
  eta <- fuzz_families[[sp$family]]$eta0 + 0.4 * d$x
  sm <- sp$special
  if (sm == "smooth" || sm == "gp") {
    eta <- eta + 0.5 * sin(0.6 * d$xs)
  } else if (sm == "mo") {
    eta <- eta + 0.25 * d$xo
  } else if (sm == "mo_int" || sm == "mo_fac") {
    eta <- eta + 0.25 * d$xo + 0.2 * d$z + 0.1 * d$xo * d$z
  }
  re <- sp$re
  if (re %in% c("ri", "ri_fx")) {
    eta <- eta + stats::rnorm(nlevels(d$g), 0, 0.5)[d$g]
  } else if (re %in% c("rs", "dbar", "diagbar")) {
    ng <- nlevels(d$g)
    eta <- eta + stats::rnorm(ng, 0, 0.5)[d$g] +
      stats::rnorm(ng, 0, 0.3)[d$g] * d$x
  } else if (re == "dbar_f") {
    ng <- nlevels(d$g)
    U <- matrix(stats::rnorm(ng * 3, 0, 0.4), ng, 3)
    eta <- eta + U[cbind(as.integer(d$g), as.integer(d$fac))]
  } else if (re == "ar1") {
    ng <- nlevels(d$g); nt <- nlevels(d$tt); rho <- 0.6
    U <- matrix(0, ng, nt)
    U[, 1] <- stats::rnorm(ng, 0, 0.5)
    for (t in 2:nt) {
      U[, t] <- rho * U[, t - 1] +
        sqrt(1 - rho^2) * stats::rnorm(ng, 0, 0.5)
    }
    eta <- eta + U[cbind(as.integer(d$g), as.integer(d$tt))]
  } else if (re == "nested") {
    eta <- eta + stats::rnorm(nlevels(d$ga), 0, 0.45)[d$ga] +
      stats::rnorm(nlevels(d$ga) * nlevels(d$gb), 0, 0.3)[
        as.integer(d$ga) * nlevels(d$gb) + as.integer(d$gb) -
          nlevels(d$gb)]
  }
  eta
}

# One draw of the response from the sampled model.
fuzz_draw <- function(fam, eta, d, trials) {
  n <- length(eta)
  switch(fam,
    gaussian = eta + stats::rnorm(n, 0, 1),
    student  = eta + 1.0 * stats::rt(n, 6),
    lognormal = exp(eta + stats::rnorm(n, 0, 0.5)),
    weibull  = stats::rweibull(n, shape = 2,
                               scale = exp(eta) / gamma(1.5)),
    Gamma    = stats::rgamma(n, shape = 3, rate = 3 / exp(eta)),
    Beta     = {
      mu <- stats::plogis(eta); phi <- 8
      pmin(pmax(stats::rbeta(n, mu * phi, (1 - mu) * phi), 1e-4),
           1 - 1e-4)
    },
    poisson  = stats::rpois(n, exp(eta)),
    binomial = stats::rbinom(n, trials, stats::plogis(eta)),
    negbinomial = stats::rnbinom(n, size = 2, mu = exp(eta)),
    zero_inflated_poisson =
      (1 - stats::rbinom(n, 1, 0.25)) * stats::rpois(n, exp(eta)),
    cumulative = {
      thr <- c(-1.2, 0, 1.2)
      lat <- eta + stats::rlogis(n)
      1L + as.integer(rowSums(outer(lat, thr, ">")))
    },
    stop("fuzz: no simulator for family ", fam)
  )
}

# A simulation is degenerate when the fit it feeds could not identify
# the model no matter how correct the code is. Regenerate with a bumped
# seed instead of reporting a spurious finding.
fuzz_degenerate <- function(fam, y, trials) {
  if (anyNA(y) || any(!is.finite(y))) return(TRUE)
  # a Bernoulli response has only two values by construction, so the
  # generic "needs three distinct values" rule does not apply to it
  bernoulli <- fam == "binomial" && all(trials == 1)
  if (!bernoulli && length(unique(y)) < 3) return(TRUE)
  switch(fam,
    binomial = if (bernoulli) min(mean(y == 0), mean(y == 1)) < 0.1
               else all(y == 0) || all(y == trials),
    poisson = ,
    negbinomial = ,
    zero_inflated_poisson = mean(y == 0) > 0.9 || max(y) > 5e3,
    cumulative = length(unique(y)) < 3,
    Beta = any(y <= 0) || any(y >= 1),
    Gamma = ,
    weibull = ,
    lognormal = any(y <= 0) || max(y) > 1e6,
    FALSE
  )
}

#' Build the data for one spec.
#' @param sp One row of a plan (a list or one-row data frame).
#' @return A data frame with `y` and every covariate the spec needs.
fuzz_data <- function(sp) {
  sp <- as.list(sp)
  fm <- fuzz_families[[sp$family]]
  for (attempt in 0:7) {
    set.seed(sp$seed + attempt * 7919L)
    d <- fuzz_design(sp)
    eta <- fuzz_eta(sp, d)
    trials <- if (identical(sp$aterm, "trials")) d$nt_ else 1L
    if (identical(sp$aterm, "trunc")) {
      # per-observation rejection keeps the design grid intact while
      # drawing genuinely from the truncated law
      y <- fuzz_draw(sp$family, eta, d, trials)
      lb <- fm$trunc_lb
      for (k in 1:400) {
        bad <- which(y < lb)
        if (!length(bad)) break
        y[bad] <- fuzz_draw(sp$family, eta, d, trials)[bad]
      }
      if (any(y < lb)) next
    } else {
      y <- fuzz_draw(sp$family, eta, d, trials)
    }
    if (fuzz_degenerate(sp$family, y, trials)) next
    d$nt <- trials
    d$y <- y
    if (fm$ordinal) d$y <- factor(d$y, ordered = TRUE)
    if (identical(sp$aterm, "cens") || identical(sp$aterm, "cens_chr")) {
      qs <- stats::quantile(y, c(0.15, 0.85), names = FALSE)
      cc <- ifelse(y < qs[1], -1L, ifelse(y > qs[2], 1L, 0L))
      d$y <- pmin(pmax(y, qs[1]), qs[2])
      d$cc <- cc
      d$ccl <- c("left", "none", "right")[cc + 2L]
      if (length(unique(cc)) < 3) next
    }
    return(d)
  }
  NULL
}

# ---------------------------------------------------------------------
# 4. spec -> source text
# ---------------------------------------------------------------------

fuzz_mu_formula <- function(sp) {
  fm <- fuzz_families[[sp$family]]
  at <- if (identical(sp$aterm, "trunc")) {
    paste0("trunc(lb = ", format(fm$trunc_lb), ")")
  } else {
    fuzz_aterm[[sp$aterm]]$term
  }
  lhs <- if (is.null(at)) "y" else paste("y |", at)
  rhs <- c("1", "x", fuzz_special[[sp$special]]$term)
  if (identical(sp$re, "crossed_bar")) {
    rhs <- c(rhs, NULL)
    return(paste0(lhs, " ~ ", paste(rhs, collapse = " + "),
                  " + x * (1 | g)"))
  }
  rhs <- c(rhs, fuzz_re[[sp$re]]$term)
  paste0(lhs, " ~ ", paste(rhs, collapse = " + "))
}

fuzz_dpar_formula <- function(sp) {
  if (identical(sp$dpar, "bad_dpar")) return("shape ~ 1 + z")
  if (!identical(sp$dpar, "dpar_x")) return(NULL)
  paste0(fuzz_families[[sp$family]]$dpar2, " ~ 1 + z")
}

fuzz_bf_text <- function(sp) {
  parts <- c(fuzz_mu_formula(sp), fuzz_dpar_formula(sp))
  paste0("bf(", paste(parts, collapse = ", "), ") + ",
         fuzz_families[[sp$family]]$ctor)
}

fuzz_mode_args <- function(sp) {
  m <- sp$mode
  ctl <- character(0)
  reml <- "FALSE"; quad <- "FALSE"
  if (m == "reml") reml <- "TRUE"
  if (m == "quadrature") quad <- "TRUE"
  if (m == "quad_reml") { quad <- "TRUE"; reml <- "TRUE" }
  if (m == "profile") ctl <- "profile = TRUE"
  if (m == "profile_reml") { ctl <- "profile = TRUE"; reml <- "TRUE" }
  if (m == "autoscale") ctl <- "autoscale = TRUE"
  if (m == "sparse_x") ctl <- "sparse_x = TRUE"
  list(reml = reml, quad = quad,
       control = if (length(ctl)) paste0("frmtmb_control(", ctl, ")")
                 else "frmtmb_control()")
}

fuzz_call_text <- function(sp, data = "d") {
  ma <- fuzz_mode_args(sp)
  paste0("frm(", fuzz_bf_text(sp), ",\n    data = ", data,
         ", REML = ", ma$reml, ", quadrature = ", ma$quad,
         ",\n    control = ", ma$control, ")")
}

fuzz_spec_key <- function(sp) {
  paste(vapply(c("family", "aterm", "re", "special", "dpar", "mode", "op"),
               function(d) paste0(d, "=", sp[[d]]), ""), collapse = " ")
}

#' Runnable reproduction for one spec.
#' @param sp One row of a plan.
#' @return A character vector of source lines.
fuzz_repro <- function(sp) {
  sp <- as.list(sp)
  sp <- sp[!startsWith(names(sp), ".")]
  c(paste0("# ", fuzz_spec_key(sp), "  seed=", sp$seed,
           "  kind=", sp$kind),
    "library(frmtmb); source(\"tests/testthat/helper-fuzz.R\")",
    paste0("sp  <- list(", paste(sprintf("%s = %s", names(sp),
                                         vapply(sp, function(v)
                                           if (is.character(v))
                                             paste0("\"", v, "\"")
                                           else as.character(v), "")),
                                 collapse = ", "), ")"),
    "d   <- fuzz_data(sp)",
    paste0("fit <- ", fuzz_call_text(sp)))
}

# ---------------------------------------------------------------------
# 5. invariants
# ---------------------------------------------------------------------

fuzz_finding <- function(recs, sp, invariant, class, detail,
                         values = NULL) {
  recs$add(list(id = sp$id, kind = sp$kind, seed = sp$seed,
                spec = fuzz_spec_key(sp), invariant = invariant,
                class = class, detail = detail,
                fit_warnings = fuzz_or(sp$.warned, character(0)),
                values = values,
                repro = paste(fuzz_repro(sp), collapse = "\n")))
  invisible(NULL)
}

# A tiny appendable collector; a growing list in a closure beats
# rbind-ing data frames per finding.
fuzz_recorder <- function() {
  out <- list()
  list(
    add = function(x) { out[[length(out) + 1L]] <<- x; invisible() },
    get = function() out
  )
}

fuzz_fit_one <- function(sp, d) {
  ma <- fuzz_mode_args(sp)
  ctl <- eval(parse(text = ma$control))
  bform <- eval(parse(text = fuzz_bf_text(sp)))
  fuzz_try(frm(bform, data = d, REML = ma$reml == "TRUE",
               quadrature = ma$quad == "TRUE", control = ctl))
}

# I4: the reported logLik is the objective at the optimum. Independent
# of any df convention, so REML and profile are checked the same way.
# A model whose likelihood is unbounded is a separate, harder failure
# than a bookkeeping mismatch, so it gets its own invariant name.
fuzz_inv_loglik_identity <- function(recs, sp, fit) {
  ll <- suppressWarnings(as.numeric(stats::logLik(fit)))
  if (!is.finite(ll)) {
    fuzz_finding(recs, sp, "finite_loglik", "candidate",
                 "logLik() is not finite",
                 list(logLik = ll, objective = fit$opt$objective))
  }
  fn <- fuzz_try(fit$obj$fn(fit$opt$par))
  if (!fn$ok) {
    return(fuzz_finding(recs, sp, "loglik_identity", "candidate",
                        paste("obj$fn(opt$par) errored:",
                              conditionMessage(fn$error))))
  }
  nfn <- -as.numeric(fn$value)
  if (!is.finite(ll) || !is.finite(nfn)) {
    if (!identical(ll, nfn)) {
      fuzz_finding(recs, sp, "loglik_identity", "candidate",
                   "logLik(fit) and -obj$fn(opt$par) disagree on infinity",
                   list(logLik = ll, neg_fn = nfn))
    }
    return(invisible())
  }
  # obj$fn re-solves the inner problem from wherever the tape was last
  # left, so exact equality is not owed on an ill-conditioned fit; the
  # tolerance is still far tighter than any real bookkeeping error
  dif <- abs(ll - nfn)
  if (dif > 1e-7) {
    fuzz_finding(recs, sp, "loglik_identity", "candidate",
                 "logLik(fit) != -obj$fn(opt$par)",
                 list(logLik = ll, neg_fn = nfn, diff = dif))
  }
}

# I2: silent NaN. A NaN estimate is acceptable only if the fit said so.
fuzz_inv_no_silent_na <- function(recs, sp, fit, warns) {
  est <- unlist(fit$estimates[c("beta", "betad", "theta")], use.names = FALSE)
  ll <- suppressWarnings(as.numeric(stats::logLik(fit)))
  bad <- (length(est) && any(!is.finite(est))) || !is.finite(ll)
  if (bad && !length(warns)) {
    fuzz_finding(recs, sp, "no_silent_na", "candidate",
                 "non-finite estimate or logLik with no warning",
                 list(n_nonfinite = sum(!is.finite(est)), logLik = ll))
  }
}

# I3: predict() on the training data reproduces fitted().
fuzz_inv_predict_fitted <- function(recs, sp, fit, d) {
  if (fuzz_families[[sp$family]]$ordinal) return(invisible())
  ft <- fuzz_try(stats::fitted(fit))
  pr <- fuzz_try(stats::predict(fit, newdata = d, type = "response"))
  if (!ft$ok || !pr$ok) {
    return(fuzz_finding(recs, sp, "predict_eq_fitted", "candidate",
                        paste("fitted/predict errored:",
                              conditionMessage(fuzz_or(ft$error, pr$error)))))
  }
  a <- as.numeric(ft$value); b <- as.numeric(pr$value)
  if (length(a) != length(b)) {
    return(fuzz_finding(recs, sp, "predict_eq_fitted", "candidate",
                        "length mismatch",
                        list(fitted = length(a), predict = length(b))))
  }
  dif <- max(abs(a - b))
  if (!is.finite(dif) || dif > 1e-8) {
    fuzz_finding(recs, sp, "predict_eq_fitted", "candidate",
                 "predict(newdata = training) != fitted()",
                 list(max_abs_diff = dif,
                      fitted_head = utils::head(a, 3),
                      predict_head = utils::head(b, 3)))
  }
}

# Are the two frames the same design, once the permutation is undone?
# Elementwise equality means the two fits also share a parameter
# meaning; span equality is the weaker structural statement that the
# model is the same even though the basis was rebuilt (mgcv rotates a
# smooth basis on rebuilt data, which is legitimate).
fuzz_design_match <- function(fr1, fr2, perm) {
  keys <- intersect(names(fr1$linpreds), names(fr2$linpreds))
  if (!length(keys)) return(list(exact = FALSE, span = TRUE))
  exact <- TRUE; span <- TRUE
  for (k in keys) {
    for (part in c("X", "Z")) {
      a <- fr1$linpreds[[k]][[part]]
      b <- fr2$linpreds[[k]][[part]]
      if (is.null(a) && is.null(b)) next
      if (is.null(a) || is.null(b)) { exact <- FALSE; span <- FALSE; next }
      a <- as.matrix(a)[perm, , drop = FALSE]
      b <- as.matrix(b)
      if (!identical(dim(a), dim(b))) { exact <- FALSE; span <- FALSE; next }
      if (max(abs(a - b)) > 1e-10) {
        exact <- FALSE
        if (fuzz_max_diff(fuzz_span_proj(a), fuzz_span_proj(b)) > 1e-8) {
          span <- FALSE
        }
      }
    }
  }
  list(exact = exact, span = span)
}

# I5: the likelihood does not depend on row order.
#
# Three checks, from strongest evidence to weakest. The design
# comparison is pure assembly and cannot be blamed on the optimizer.
# The shared-parameter objective comparison is exact but only means
# anything when both fits parameterize the model identically - a
# rebuilt smooth basis spans the same space in a different rotation,
# and then `beta` denotes different things in the two fits. The logLik
# comparison is what the invariant nominally says, and is the one that
# survives a reparameterization.
fuzz_inv_permutation <- function(recs, sp, fit, d) {
  set.seed(sp$seed + 31L)
  perm <- sample.int(nrow(d))
  d2 <- d[perm, , drop = FALSE]
  rownames(d2) <- NULL
  f2 <- fuzz_fit_one(sp, d2)
  if (!f2$ok) {
    return(fuzz_finding(recs, sp, "row_permutation", "candidate",
                        paste("refit on permuted rows errored:",
                              conditionMessage(f2$error))))
  }
  dm <- fuzz_try(fuzz_design_match(fit$frame, f2$value$frame, perm))
  if (isTRUE(dm$ok) && !dm$value$span) {
    return(fuzz_finding(recs, sp, "row_permutation", "candidate",
                        paste("the design spans a different space",
                              "after a row permutation")))
  }
  if (isTRUE(dm$ok) && dm$value$exact &&
      length(fit$opt$par) == length(f2$value$opt$par)) {
    par <- fit$opt$par
    s1 <- fuzz_try(as.numeric(fit$obj$fn(par)))
    s2 <- fuzz_try(as.numeric(f2$value$obj$fn(par)))
    if (s1$ok && s2$ok) {
      ds <- abs(s1$value - s2$value)
      if (!is.finite(ds) || ds > 1e-8) {
        return(fuzz_finding(recs, sp, "row_permutation", "candidate",
                            paste("objective at a shared parameter",
                                  "vector changed under a row",
                                  "permutation"),
                            list(obj = s1$value, obj_permuted = s2$value,
                                 diff = ds)))
      }
    }
  }
  l1 <- suppressWarnings(as.numeric(stats::logLik(fit)))
  l2 <- suppressWarnings(as.numeric(stats::logLik(f2$value)))
  dif <- abs(l1 - l2)
  if (!is.finite(dif) || dif > 1e-6) {
    fuzz_finding(recs, sp, "row_permutation", "candidate",
                 "logLik changed under a row permutation",
                 list(logLik = l1, logLik_permuted = l2, diff = dif))
  }
}

# I6: unit weights are the unweighted model.
fuzz_inv_unit_weights <- function(recs, sp, fit, d) {
  sp2 <- sp; sp2$aterm <- "unit_w"
  f2 <- fuzz_fit_one(sp2, d)
  if (!f2$ok) {
    return(fuzz_finding(recs, sp, "unit_weights", "candidate",
                        paste("weights(one) refit errored:",
                              conditionMessage(f2$error))))
  }
  l1 <- suppressWarnings(as.numeric(stats::logLik(fit)))
  l2 <- suppressWarnings(as.numeric(stats::logLik(f2$value)))
  dif <- abs(l1 - l2)
  if (!is.finite(dif) || dif > 1e-6) {
    fuzz_finding(recs, sp, "unit_weights", "candidate",
                 "weights = 1 differs from the unweighted fit",
                 list(logLik = l1, logLik_weighted = l2, diff = dif))
  }
}

# I7: the simulator and the fitted mean describe the same distribution.
# This is the invariant that catches a post-fit surface which ignores an
# aterm (truncation especially): the fit is right, the two summaries
# disagree.
fuzz_inv_simulate_mean <- function(recs, sp, fit, nsim = 200L) {
  fm <- fuzz_families[[sp$family]]
  if (!fm$mean_check) return(invisible())
  sm <- fuzz_try(stats::simulate(fit, nsim = nsim, seed = sp$seed + 5L))
  pm <- fuzz_try(stats::predict(fit, type = "response"))
  if (!sm$ok || !pm$ok) {
    return(fuzz_finding(recs, sp, "simulate_mean", "candidate",
                        paste("simulate/predict errored:",
                              conditionMessage(fuzz_or(sm$error, pm$error)))))
  }
  means <- colMeans(as.matrix(sm$value))
  target <- mean(as.numeric(pm$value))
  se <- stats::sd(means) / sqrt(length(means))
  if (!is.finite(se) || se <= 0) return(invisible())
  z <- abs(mean(means) - target) / se
  if (z > 4) {
    fuzz_finding(recs, sp, "simulate_mean", "candidate",
                 "mean of simulate() draws is far from mean(predict())",
                 list(sim_mean = mean(means), predict_mean = target,
                      se = se, z = z))
  }
}

# I7b: a simulator must respect the support the model was fitted under.
# Comparing simulate() with predict() cannot see an aterm that BOTH of
# them ignore, which is exactly how a truncated post-fit surface fails;
# the truncation bound is the independent reference that can.
fuzz_inv_trunc_support <- function(recs, sp, fit, nsim = 10L) {
  lb <- fuzz_families[[sp$family]]$trunc_lb
  sm <- fuzz_try(stats::simulate(fit, nsim = nsim, seed = sp$seed + 11L))
  if (!sm$ok) {
    return(fuzz_finding(recs, sp, "trunc_support", "candidate",
                        paste("simulate() errored on a truncated fit:",
                              conditionMessage(sm$error))))
  }
  v <- as.numeric(as.matrix(sm$value))
  v <- v[is.finite(v)]
  if (!length(v)) return(invisible())
  below <- sum(v < lb)
  if (below > 0) {
    fuzz_finding(recs, sp, "trunc_support", "candidate",
                 "simulate() drew values below the truncation bound",
                 list(lb = lb, n_below = below, n_draws = length(v),
                      min_draw = min(v)))
  }
}

# I8: vcov is a covariance matrix, and summary() prints.
fuzz_inv_vcov_summary <- function(recs, sp, fit) {
  V <- fuzz_try(stats::vcov(fit))
  if (!V$ok) {
    return(fuzz_finding(recs, sp, "vcov_psd", "candidate",
                        paste("vcov() errored:",
                              conditionMessage(V$error))))
  }
  M <- as.matrix(V$value)
  if (nrow(M) != ncol(M)) {
    return(fuzz_finding(recs, sp, "vcov_psd", "candidate",
                        "vcov() is not square",
                        list(dim = dim(M))))
  }
  # The reference is the optimizer's own coefficient count, not
  # length(fixef()): fixef() reports dpars held at a constant (se()
  # maps sigma out that way) and vcov() documents that it excludes
  # them. Under REML and profile the coefficients are inner
  # parameters and absent from opt$par, so there is nothing to compare.
  if (!sp$mode %in% c("reml", "profile")) {
    nc <- sum(names(fit$opt$par) %in% c("beta", "betad"))
    if (nc > 0 && nc != nrow(M)) {
      fuzz_finding(recs, sp, "vcov_dim", "candidate",
                   "vcov() dimension != estimated coefficient count",
                   list(vcov_dim = nrow(M), n_coef = nc))
    }
  }
  if (any(!is.finite(M))) {
    fuzz_finding(recs, sp, "vcov_psd", "candidate",
                 "vcov() has non-finite entries",
                 list(n_nonfinite = sum(!is.finite(M))))
  } else {
    asym <- max(abs(M - t(M)))
    if (asym > 1e-8 * max(1, max(abs(M)))) {
      fuzz_finding(recs, sp, "vcov_psd", "candidate",
                   "vcov() is not symmetric", list(asymmetry = asym))
    }
    ev <- min(eigen((M + t(M)) / 2, symmetric = TRUE,
                    only.values = TRUE)$values)
    if (ev < -1e-8 * max(1, max(abs(M)))) {
      fuzz_finding(recs, sp, "vcov_psd", "candidate",
                   "vcov() is not positive semidefinite",
                   list(min_eigenvalue = ev))
    }
  }
  s <- fuzz_try(utils::capture.output(print(summary(fit))))
  if (!s$ok) {
    fuzz_finding(recs, sp, "summary_prints", "candidate",
                 paste("summary() errored:", conditionMessage(s$error)))
  }
}

# I9: Wald intervals are finite, ordered, and named.
fuzz_inv_confint <- function(recs, sp, fit) {
  ci <- fuzz_try(stats::confint(fit, method = "wald"))
  if (!ci$ok) {
    return(fuzz_finding(recs, sp, "confint_wald", "candidate",
                        paste("confint() errored:",
                              conditionMessage(ci$error))))
  }
  M <- ci$value
  if (!nrow(M)) {
    return(fuzz_finding(recs, sp, "confint_wald", "candidate",
                        "confint() returned no rows"))
  }
  bad <- which(is.finite(M[, "lwr"]) & is.finite(M[, "upr"]) &
                 M[, "lwr"] > M[, "upr"])
  if (length(bad)) {
    fuzz_finding(recs, sp, "confint_wald", "candidate",
                 "confint() lower bound above upper bound",
                 list(rows = rownames(M)[bad]))
  }
  inside <- is.finite(M[, "est"]) &
    (M[, "est"] < M[, "lwr"] - 1e-6 | M[, "est"] > M[, "upr"] + 1e-6)
  if (any(inside, na.rm = TRUE)) {
    fuzz_finding(recs, sp, "confint_wald", "candidate",
                 "estimate outside its own interval",
                 list(rows = rownames(M)[which(inside)]))
  }
}

# ---------------------------------------------------------------------
# 6. brms differential oracle (data layer only - no Stan)
# ---------------------------------------------------------------------
#
# brms owns the meaning of this grammar, and make_standata() builds
# every design object without touching Stan, so it is a millisecond
# oracle. Estimation modes (REML, profile, quadrature, autoscale,
# sparse_x) never reach the data layer, so they are compared through
# their frame like any ML spec.

fuzz_brms_translatable <- function(sp) {
  if (!requireNamespace("brms", quietly = TRUE)) return(FALSE)
  if (grepl("^refusal", sp$kind)) return(FALSE)
  if (is.null(fuzz_families[[sp$family]]$brms)) return(FALSE)
  if (is.null(fuzz_re[[sp$re]]$term) && sp$re != "none") return(FALSE)
  if (sp$re != "none" && is.null(fuzz_re[[sp$re]]$brms)) return(FALSE)
  TRUE
}

fuzz_brms_bf_text <- function(sp) {
  fm <- fuzz_families[[sp$family]]
  at <- if (identical(sp$aterm, "trunc")) {
    paste0("trunc(lb = ", format(fm$trunc_lb), ")")
  } else {
    fuzz_aterm[[sp$aterm]]$term
  }
  lhs <- if (is.null(at)) "y" else paste("y |", at)
  rhs <- c("1", "x", fuzz_special[[sp$special]]$brms,
           if (sp$re != "none") fuzz_re[[sp$re]]$brms)
  parts <- paste0(lhs, " ~ ", paste(rhs, collapse = " + "))
  if (identical(sp$dpar, "dpar_x")) {
    parts <- c(parts, paste0(fm$dpar2, " ~ 1 + z"))
  }
  paste0("brms::bf(", paste(parts, collapse = ", "), ")")
}

fuzz_span_proj <- function(M) {
  M <- as.matrix(M)
  q <- qr(M)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  Q %*% t(Q)
}

fuzz_max_diff <- function(a, b) {
  a <- unname(as.matrix(a)); b <- unname(as.matrix(b))
  if (!identical(dim(a), dim(b))) return(Inf)
  max(abs(a - b))
}

# Compare our assembled frame against brms's standata for the same
# model. Only value-level conventions the two packages actually share
# are asserted; naming and storage differences are documented above
# each block.
fuzz_inv_brms <- function(recs, sp, d, frame, fit_ok) {
  fm <- fuzz_families[[sp$family]]
  bt <- fuzz_brms_bf_text(sp)
  sdat <- fuzz_try(suppressMessages(suppressWarnings(
    brms::make_standata(eval(parse(text = bt)), data = d,
                        family = eval(parse(text = fm$brms))))))
  if (!sdat$ok) {
    if (fit_ok) {
      fuzz_finding(recs, sp, "brms_translation", "grammar_divergence",
                   paste0("frmtmb accepted the model, brms rejected it: ",
                          conditionMessage(sdat$error)),
                   list(brms_formula = bt))
    }
    return(invisible())
  }
  if (!fit_ok) {
    return(fuzz_finding(recs, sp, "brms_translation",
                        "grammar_divergence",
                        "brms accepted the model, frmtmb did not",
                        list(brms_formula = bt)))
  }
  sd <- sdat$value
  key <- "y.mu"
  lp <- frame$linpreds[[key]]
  if (is.null(lp)) return(invisible())

  # --- population-level design -------------------------------------
  # brms keeps mo() coefficients out of X (they live in bsp); we carry
  # a zero placeholder column so the beta bookkeeping stays uniform.
  X <- as.matrix(lp$X)
  if (length(lp$mo)) {
    drop_cols <- vapply(lp$mo, function(m) m$col, 1L)
    X <- X[, -drop_cols, drop = FALSE]
  }
  if (sp$special == "smooth") {
    # brms calls smoothCon(diagonal.penalty = TRUE) and ships the
    # smooth's fixed part separately as Xs: same model, rotated basis
    if (!is.null(sd$Xs) &&
        fuzz_max_diff(fuzz_span_proj(X),
                      fuzz_span_proj(cbind(sd$X, sd$Xs))) > 1e-8) {
      fuzz_finding(recs, sp, "brms_X", "candidate",
                   "population-level column span differs from brms",
                   list(ncol_ours = ncol(X),
                        ncol_brms = ncol(sd$X) + ncol(sd$Xs)))
    }
  } else {
    dif <- fuzz_max_diff(X, sd$X)
    if (dif > 1e-10) {
      fuzz_finding(recs, sp, "brms_X", "candidate",
                   "population-level design differs from brms",
                   list(max_abs_diff = dif, dim_ours = dim(X),
                        dim_brms = dim(as.matrix(sd$X))))
    }
  }

  # --- dpar design --------------------------------------------------
  if (identical(sp$dpar, "dpar_x")) {
    kd <- paste0("y.", fm$dpar2)
    Xb <- sd[[paste0("X_", fm$dpar2)]]
    if (!is.null(frame$linpreds[[kd]]) && !is.null(Xb)) {
      dif <- fuzz_max_diff(frame$linpreds[[kd]]$X, Xb)
      if (dif > 1e-10) {
        fuzz_finding(recs, sp, "brms_X_dpar", "candidate",
                     paste0("dpar design (", fm$dpar2,
                            ") differs from brms"),
                     list(max_abs_diff = dif))
      }
    }
  }

  # --- addition terms ------------------------------------------------
  av <- frame$aterm_values$y
  cmp <- function(ours, theirs, nm) {
    if (is.null(ours) && is.null(theirs)) return(invisible())
    if (is.null(ours) || is.null(theirs)) {
      return(fuzz_finding(recs, sp, "brms_aterm", "candidate",
                          paste0(nm, ": present in one package only"),
                          list(ours = !is.null(ours),
                               brms = !is.null(theirs))))
    }
    # brms ships a per-observation vector where a constant bound is
    # enough for us; recycling makes the values comparable
    ours <- as.numeric(ours); theirs <- as.numeric(theirs)
    n <- max(length(ours), length(theirs))
    dif <- fuzz_max_diff(rep_len(ours, n), rep_len(theirs, n))
    if (dif > 1e-10) {
      fuzz_finding(recs, sp, "brms_aterm", "candidate",
                   paste0(nm, " differs from brms"),
                   list(max_abs_diff = dif))
    }
  }
  if (sp$aterm == "trials")  cmp(av$trials, sd$trials, "trials")
  if (sp$aterm == "weights") cmp(av$weights, sd$weights, "weights")
  if (sp$aterm %in% c("cens", "cens_chr")) cmp(av$cens, sd$cens, "cens")
  if (sp$aterm == "trunc")   cmp(av$trunc_lb, sd$lb, "trunc lb")
  if (sp$aterm == "se")      cmp(av$se, sd$se, "se")

  # --- random-effect blocks ------------------------------------------
  # DIVERGENCE (convention, never asserted): brms numbers grouping
  # factors and keeps `||` terms in ONE block with zero correlations,
  # where lme4-style expansion (ours) splits them into independent
  # blocks. Block SHAPES therefore differ legitimately. What cannot
  # differ is how many random coefficients the model has in total, and
  # over how many levels - so that is what is compared. s()/gp() add
  # blocks brms represents outside its grouping machinery, so this runs
  # only when the special contributes no block.
  if (sp$re != "none" && fuzz_special[[sp$special]]$blk == 0) {
    ours_n <- sum(vapply(frame$re_blocks,
                         function(b) b$n_levels * b$dim, 1))
    ours_lv <- sort(unique(vapply(frame$re_blocks,
                                  function(b) as.integer(b$n_levels), 1L)))
    idx <- 1L; theirs_n <- 0; theirs_lv <- integer(0)
    while (!is.null(sd[[paste0("N_", idx)]])) {
      N <- as.integer(sd[[paste0("N_", idx)]])
      M <- as.integer(sd[[paste0("M_", idx)]])
      theirs_n <- theirs_n + N * M
      theirs_lv <- c(theirs_lv, N)
      idx <- idx + 1L
    }
    theirs_lv <- sort(unique(theirs_lv))
    if (!isTRUE(all.equal(ours_n, theirs_n)) ||
        !identical(ours_lv, theirs_lv)) {
      fuzz_finding(recs, sp, "brms_re_blocks", "candidate",
                   "random-effect size differs from brms",
                   list(n_coef_ours = ours_n, n_coef_brms = theirs_n,
                        levels_ours = ours_lv, levels_brms = theirs_lv))
    }
  }

  # --- specials ------------------------------------------------------
  if (sp$special %in% c("mo", "mo_int") && length(lp$mo)) {
    if (!is.null(sd$Xmo_1)) {
      dif <- fuzz_max_diff(lp$mo[[1]]$codes, sd$Xmo_1)
      if (dif > 1e-10) {
        fuzz_finding(recs, sp, "brms_mo", "candidate",
                     "mo() codes differ from brms Xmo_1",
                     list(max_abs_diff = dif))
      }
      if (!identical(as.integer(lp$mo[[1]]$D), as.integer(sd$Jmo[1]))) {
        fuzz_finding(recs, sp, "brms_mo", "candidate",
                     "mo() category count differs from brms Jmo",
                     list(ours = lp$mo[[1]]$D, brms = as.integer(sd$Jmo[1])))
      }
    }
  }
  # lp$Z spans every random column of the model, so a special's own
  # columns have to be selected through its block before the bases can
  # be compared with brms's per-term matrices
  block_cols <- function(ids) {
    unlist(lapply(ids, function(i) frame$re_blocks[[i]]$c_idx))
  }
  if (sp$special == "gp" && length(lp$gps) && !is.null(sd$Lgp_1)) {
    gi <- lp$gps[[1]]
    if (fuzz_max_diff(as.numeric(gi$L), as.numeric(sd$Lgp_1)) > 1e-10) {
      fuzz_finding(recs, sp, "brms_gp", "candidate",
                   "HSGP boundary L differs from brms Lgp_1",
                   list(ours = as.numeric(gi$L),
                        brms = as.numeric(sd$Lgp_1)))
    }
    if (!is.null(sd$Xgp_1) && !is.null(sd$Jgp_1) && !is.null(lp$Z) &&
        !is.null(gi$block_id)) {
      Zg <- as.matrix(lp$Z[, block_cols(gi$block_id), drop = FALSE])
      Zb <- sd$Xgp_1[sd$Jgp_1, , drop = FALSE]
      if (fuzz_max_diff(Zg, Zb) > 1e-10) {
        fuzz_finding(recs, sp, "brms_gp", "candidate",
                     "HSGP basis differs from brms Xgp_1[Jgp_1, ]",
                     list(dim_ours = dim(Zg), dim_brms = dim(Zb),
                          max_abs_diff = fuzz_max_diff(Zg, Zb)))
      }
    }
  }
  if (sp$special == "smooth" && length(lp$smooths) && !is.null(sd$knots_1)) {
    sm <- lp$smooths[[1]]
    if (!identical(as.integer(sm$nr), as.integer(sd$knots_1))) {
      fuzz_finding(recs, sp, "brms_smooth", "candidate",
                   "smooth wiggly dimension differs from brms knots_1",
                   list(ours = as.integer(sm$nr),
                        brms = as.integer(sd$knots_1)))
    }
    # DIVERGENCE (convention): brms passes diagonal.penalty = TRUE to
    # smoothCon, so the wiggly bases are a reparameterization of each
    # other. Only the span is comparable.
    if (!is.null(sd$Zs_1_1) && !is.null(lp$Z) && !is.null(sm$block_ids)) {
      Zs <- as.matrix(lp$Z[, block_cols(sm$block_ids), drop = FALSE])
      if (ncol(Zs) != ncol(sd$Zs_1_1)) {
        fuzz_finding(recs, sp, "brms_smooth", "candidate",
                     "smooth wiggly basis has a different width than brms",
                     list(ours = ncol(Zs), brms = ncol(sd$Zs_1_1)))
      } else if (fuzz_max_diff(fuzz_span_proj(Zs),
                               fuzz_span_proj(sd$Zs_1_1)) > 1e-8) {
        fuzz_finding(recs, sp, "brms_smooth", "candidate",
                     "smooth wiggly basis spans a different space than brms")
      }
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------
# 7. refusal probes
# ---------------------------------------------------------------------

fuzz_check_refusal <- function(recs, sp, res) {
  cs <- Filter(function(k) paste0("refusal:", k$name) == sp$kind,
               fuzz_refusal_cases)[[1]]
  if (res$ok) {
    return(fuzz_finding(recs, sp, "refusal_is_error", "candidate",
                        "a refused-by-design combination was accepted",
                        list(case = cs$name)))
  }
  msg <- conditionMessage(res$error)
  if (fuzz_is_crash(msg)) {
    return(fuzz_finding(recs, sp, "refusal_is_clean", "candidate",
                        "refusal surfaced as an internal crash",
                        list(case = cs$name, message = msg)))
  }
  if (!grepl(cs$expect, msg, ignore.case = TRUE)) {
    fuzz_finding(recs, sp, "refusal_message", "candidate",
                 "refusal message does not name the reason",
                 list(case = cs$name, expected = cs$expect,
                      message = msg))
  }
}

# ---------------------------------------------------------------------
# 8. runner
# ---------------------------------------------------------------------

#' Run a fuzz plan.
#' @param plan A data frame from [fuzz_plan()].
#' @param brms_oracle Compare designs against brms::make_standata().
#' @param nsim Draws for the simulate-mean invariant.
#' @param progress Print one line per spec.
#' @return A list with `findings` (list of records) and `timing`.
fuzz_run <- function(plan, brms_oracle = TRUE, nsim = 200L,
                     progress = FALSE) {
  recs <- fuzz_recorder()
  use_brms <- brms_oracle && requireNamespace("brms", quietly = TRUE)
  t_start <- proc.time()[["elapsed"]]
  n_fit <- 0L
  for (i in seq_len(nrow(plan))) {
    sp <- as.list(plan[i, ])
    t0 <- proc.time()[["elapsed"]]
    dg <- fuzz_try(fuzz_data(sp))
    d <- dg$value
    if (is.null(d)) {
      fuzz_finding(recs, sp, "data_generation", "generator",
                   if (dg$ok) "no non-degenerate dataset after 8 seeds"
                   else paste("data generation errored:",
                              conditionMessage(dg$error)))
      next
    }
    res <- fuzz_fit_one(sp, d)
    n_fit <- n_fit + 1L
    # every finding for this spec carries whatever the fit itself said,
    # so a "singular vcov" record can be read against a convergence
    # warning that already fired
    sp$.warned <- utils::head(res$warnings, 3)

    if (grepl("^refusal", sp$kind)) {
      fuzz_check_refusal(recs, sp, res)
      if (progress) {
        cat(sprintf("[%3d/%3d] %-6s %s\n", i, nrow(plan), "refuse",
                    fuzz_spec_key(sp)))
      }
      next
    }

    if (!res$ok) {
      msg <- conditionMessage(res$error)
      fuzz_finding(recs, sp,
                   if (fuzz_is_crash(msg)) "fit_crash" else "fit_error",
                   "candidate",
                   paste("frm() failed:", msg))
    } else {
      fit <- res$value
      # an invariant that blows up is a harness defect, not a package
      # one: record it as such and keep the run going
      guard <- function(nm, f, ...) {
        g <- fuzz_try(f(recs, sp, ...))
        if (!g$ok) {
          fuzz_finding(recs, sp, nm, "generator",
                       paste("invariant errored:",
                             conditionMessage(g$error)))
        }
      }
      guard("loglik_identity", fuzz_inv_loglik_identity, fit)
      guard("no_silent_na", fuzz_inv_no_silent_na, fit, res$warnings)
      guard("predict_eq_fitted", fuzz_inv_predict_fitted, fit, d)
      guard("row_permutation", fuzz_inv_permutation, fit, d)
      if (sp$aterm == "none") {
        guard("unit_weights", fuzz_inv_unit_weights, fit, d)
      }
      if (sp$aterm == "trunc" && fuzz_families[[sp$family]]$sim) {
        guard("trunc_support", fuzz_inv_trunc_support, fit)
      }
      if (sp$op == "simulate") {
        guard("simulate_mean", fuzz_inv_simulate_mean, fit, nsim)
      }
      if (sp$op == "vcov")    guard("vcov_psd", fuzz_inv_vcov_summary, fit)
      if (sp$op == "confint") guard("confint_wald", fuzz_inv_confint, fit)
    }

    if (use_brms && fuzz_brms_translatable(sp)) {
      # grammar acceptance is a parse/assembly question, so the oracle
      # compares against dry_run = "frame", never against whether the
      # optimizer succeeded: an estimation failure is not a refusal
      fr <- fuzz_try(frm(eval(parse(text = fuzz_bf_text(sp))), data = d,
                         dry_run = "frame"))
      bo <- fuzz_try(fuzz_inv_brms(recs, sp, d, fr$value, fr$ok))
      if (!bo$ok) {
        fuzz_finding(recs, sp, "brms_oracle", "generator",
                     paste("brms oracle errored:",
                           conditionMessage(bo$error)))
      }
    }
    if (progress) {
      cat(sprintf("[%3d/%3d] %5.1fs %s\n", i, nrow(plan),
                  proc.time()[["elapsed"]] - t0, fuzz_spec_key(sp)))
      utils::flush.console()
    }
  }
  list(findings = recs$get(),
       elapsed = proc.time()[["elapsed"]] - t_start,
       n_specs = nrow(plan), n_fits = n_fit)
}

# ---------------------------------------------------------------------
# 9. triage
# ---------------------------------------------------------------------
#
# Findings that reproduce a defect already fixed on an unmerged sibling
# branch are marked KNOWN-PENDING so they do not compete for attention
# with anything new. Matching is on the spec plus the invariant. A rule
# that could swallow an unrelated failure of the same shape also names
# the symptom, because a wrong KNOWN-PENDING hides a real defect while
# a wrong REAL-NEW only costs a second look.

# quadrature is implicated in defects of its own (see
# dev/fuzz-findings.md), so no in-flight rule may claim a quadrature
# fit: those findings must stand on their own.
fuzz_not_quad <- function(f) !grepl("mode=quadrature", f$spec)

FUZZ_KNOWN_PENDING <- list(
  list(id = "trunc-postfit",
       why = "truncation ignored by fitted/predict/simulate/residuals",
       match = function(f) grepl("aterm=trunc", f$spec) &&
         (f$invariant %in% c("trunc_support", "simulate_mean") ||
            (f$invariant == "predict_eq_fitted" && fuzz_not_quad(f)))),
  list(id = "hsgp-brms-scaling",
       why = "HSGP basis scaling matches brms only after the pending merge",
       match = function(f) f$invariant == "brms_gp"),
  list(id = "double-bar-factor",
       why = "(f || g) builds a correlated block instead of independent ones",
       match = function(f) grepl("re=dbar_f", f$spec) &&
         f$invariant == "brms_re_blocks"),
  list(id = "re-factor-call",
       why = "(1 | factor(x)) errors instead of grouping by the factor",
       match = function(f) grepl("re=ri_fx", f$spec) &&
         ((f$invariant %in% c("fit_error", "fit_crash") &&
             grepl("grouping factor", f$detail)) ||
            f$invariant == "brms_translation")),
  list(id = "bar-crossed-star",
       why = "bar terms crossed with * are accepted instead of refused",
       match = function(f) identical(f$kind, "refusal:bar_crossed")),
  list(id = "cens-character-codes",
       why = "character cens() codes fail",
       match = function(f) grepl("aterm=cens_chr", f$spec) &&
         f$invariant %in% c("fit_error", "fit_crash") &&
         grepl("cens", f$detail)),
  list(id = "mo-factor-interaction",
       why = "mo() crossed with a factor is refused rather than fitted",
       match = function(f) identical(f$kind, "refusal:mo_factor"))
)

# Invariants that only mean anything about a fit that actually reached
# an optimum. A boundary or non-converged fit can violate them without
# any code being wrong, so they are reported under "unconverged" and do
# not compete with real defects. Everything else (assembly identities,
# design agreement, refusals, crashes) is judged regardless.
# loglik_identity belongs here because obj$fn re-solves the inner
# problem from wherever the tape was last left, which only drifts when
# that inner problem is ill-conditioned. finite_loglik deliberately
# does NOT: reporting an infinite logLik is a defect whatever the
# optimizer said about its gradient.
FUZZ_CONVERGENCE_SENSITIVE <- c(
  "row_permutation", "vcov_psd", "vcov_dim", "confint_wald",
  "simulate_mean", "predict_eq_fitted", "unit_weights",
  "loglik_identity"
)

FUZZ_NONCONVERGENCE <- paste(
  c("did not report convergence", "maximum absolute gradient",
    "NA/NaN function evaluation", "singular", "not positive definite"),
  collapse = "|")

fuzz_triage <- function(findings) {
  lapply(findings, function(f) {
    if (identical(f$class, "generator")) return(f)
    hit <- Filter(function(k) isTRUE(k$match(f)), FUZZ_KNOWN_PENDING)
    if (length(hit)) {
      f$class <- "known_pending"
      f$known_id <- hit[[1]]$id
      f$known_why <- hit[[1]]$why
      return(f)
    }
    if (identical(f$class, "candidate")) {
      shaky <- length(f$fit_warnings) &&
        any(grepl(FUZZ_NONCONVERGENCE, f$fit_warnings))
      f$class <- if (shaky && f$invariant %in% FUZZ_CONVERGENCE_SENSITIVE)
        "unconverged" else "real_new"
    }
    f
  })
}

fuzz_summary <- function(res) {
  tri <- fuzz_triage(res$findings)
  cls <- vapply(tri, function(f) f$class, "")
  list(triaged = tri, counts = table(cls),
       elapsed = res$elapsed, n_specs = res$n_specs)
}

# Compact, greppable rendering of one finding.
fuzz_format <- function(f) {
  vals <- if (is.null(f$values)) "" else
    paste0("    values: ",
           paste(sprintf("%s=%s", names(f$values),
                         vapply(f$values, function(v)
                           paste(format(utils::head(v, 4), digits = 8),
                                 collapse = ","), "")),
                 collapse = "; "), "\n")
  wr <- if (!length(f$fit_warnings)) "" else
    paste0("    fit warned: ", paste(f$fit_warnings, collapse = " | "), "\n")
  paste0("[", toupper(gsub("_", "-", f$class)), "] ", f$invariant, "\n",
         "    spec: ", f$spec, " (seed ", f$seed, ", ", f$kind, ")\n",
         "    ", f$detail, "\n", vals, wr)
}
