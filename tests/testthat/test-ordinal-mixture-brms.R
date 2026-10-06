# The log density of ordinal mixtures and of hurdle_cumulative() with
# thres(gr = ) and cs() against brms 2.23.0's compiled Stan program,
# at matched parameter values: frmtmb's estimates and the estimates
# moved by a fixed pattern, each translated to brms's parameters. Every
# prior is flat, so brms's log density with adjust_transform = FALSE is
# the log likelihood, up to the constant Gamma(K) of the dirichlet(1)
# brms keeps on a simplex theta. Gated: each program compiles once, then
# comes from FRMTMB_STAN_CACHE. dev/ordmix-lpcheck.R runs the same
# translation over more shapes and records the ulp counts.

omb_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < 0.2, 0L, d$y)
  d
}

# the names declared in the Stan parameters block, in order
omb_par_names <- function(code) {
  lines <- sub("//.*$", "", strsplit(code, "\n", fixed = TRUE)[[1]])
  start <- grep("^parameters[[:space:]]*[{]", lines)
  out <- character()
  depth <- 0L
  for (i in seq(start[[1]], length(lines))) {
    ln <- lines[[i]]
    depth <- depth + lengths(regmatches(ln, gregexpr("[{]", ln))) -
      lengths(regmatches(ln, gregexpr("[}]", ln)))
    m <- regmatches(ln, regexpr(
      "([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*;[[:space:]]*$", ln))
    if (length(m)) out <- c(out, trimws(sub(";.*$", "", m)))
    if (i > start[[1]] && depth <= 0L) break
  }
  out
}

# frmtmb's parameter list `est` as brms's constrained parameters
omb_translate <- function(fit, est, code, sdat) {
  fam <- family(fit)
  mx <- fam[["mix"]][["ord"]]
  is_mix <- !is.null(mx)
  rn <- names(fit$spec$responses)[1L]
  lpk <- function(dp) {
    Filter(function(l) identical(l$dpar, dp), fit$frame$linpreds)[[1L]]
  }
  coef_of <- function(dp) {
    lp <- lpk(dp)
    v <- est[[lp$par]][lp$idx]
    names(v) <- colnames(lp$X)[seq_along(v)]
    v
  }
  dsx <- function(dp) if (dp == "mu" && !is_mix) "" else paste0("_", dp)
  bX <- function(dp) {
    X <- sdat[[paste0("X", dsx(dp))]]
    cn <- colnames(X)
    X <- as.matrix(X)
    colnames(X) <- cn
    X
  }
  centered <- function(dp) !is.null(sdat[[paste0("Kc", dsx(dp))]])
  shift_of <- function(dp) {
    if (!centered(dp)) return(0)
    X <- bX(dp)
    b <- coef_of(dp)
    cn <- setdiff(colnames(X), "Intercept")
    sum(colMeans(X)[cn] * unname(b[cn]))
  }
  comp_tau <- function(k) {
    if (!is_mix) return(ord_threshold_values(fam, est[["tau_raw"]]))
    ord_threshold_values(mx$views[[k]], est[[mx$tau_names[k]]])
  }
  nth <- fam[["thres"]][["nthres"]] %||% length(comp_tau(1L))
  lay <- thres_layout(nth)
  slice <- function(tau, gi) {
    if (length(nth) == 1L) tau else tau[lay$start[gi]:lay$end[gi]]
  }
  kg <- function(nm, pre) {
    rest <- sub(paste0("^", pre), "", nm)
    k <- if (grepl("^_mu[0-9]+", rest)) {
      as.integer(sub("^_mu([0-9]+).*$", "\\1", rest))
    } else 1L
    rest <- sub("^_mu[0-9]+", "", rest)
    gi <- if (grepl("^_[0-9]+$", rest)) as.integer(sub("^_", "", rest)) else
      1L
    list(k = k, gi = gi, dp = if (is_mix) paste0("mu", k) else "mu")
  }
  out <- list()
  for (nm in omb_par_names(code)) {
    if (grepl("^b(_(mu|disc|hu|theta)[0-9]*)?$", nm)) {
      dp <- if (nm == "b") "mu" else sub("^b_", "", nm)
      b <- coef_of(dp)
      cn <- setdiff(colnames(bX(dp)), if (centered(dp)) "Intercept")
      cn[cn == "Intercept"] <- "(Intercept)"
      out[[nm]] <- array(unname(b[cn]), length(cn))
    } else if (grepl("^Intercept_(disc|hu|theta)[0-9]*$", nm)) {
      dp <- sub("^Intercept_", "", nm)
      b <- coef_of(dp)
      X <- bX(dp)
      cn <- setdiff(colnames(X), "Intercept")
      out[[nm]] <- unname(b[["(Intercept)"]]) +
        sum(colMeans(X)[cn] * unname(b[cn]))
    } else if (grepl("^Intercept(_mu[0-9]+)?(_[0-9]+)?$", nm)) {
      a <- kg(nm, "Intercept")
      out[[nm]] <- as.numeric(slice(comp_tau(a$k), a$gi) - shift_of(a$dp))
    } else if (grepl("^first_Intercept(_mu[0-9]+)?(_[0-9]+)?$", nm)) {
      a <- kg(nm, "first_Intercept")
      out[[nm]] <- slice(comp_tau(a$k), a$gi)[1] - shift_of(a$dp)
    } else if (grepl("^delta(_mu[0-9]+)?(_[0-9]+)?$", nm)) {
      a <- kg(nm, "delta")
      tk <- slice(comp_tau(a$k), a$gi)
      out[[nm]] <- tk[2] - tk[1]
    } else if (grepl("^fixed_Intercept(_[0-9]+)?$", nm)) {
      gi <- if (grepl("_[0-9]+$", nm)) as.integer(sub("^.*_", "", nm)) else
        1L
      kf <- which(!vapply(mx$comps, function(cp) {
        identical(cp$threshold, "sum_to_zero")
      }, NA))[1]
      out[[nm]] <- as.numeric(slice(comp_tau(if (is.na(kf)) 1L else kf),
                                    gi))
    } else if (grepl("^bcs(_mu[0-9]+)?$", nm)) {
      lp <- lpk(if (nm == "bcs") "mu" else sub("^bcs_", "", nm))
      v <- unlist(lapply(lp$cs, function(ct) est[[ct$par]]))
      out[[nm]] <- matrix(v, nrow = length(lp$cs), byrow = TRUE)
    } else if (identical(nm, "theta")) {
      eta <- vapply(seq_len(mx$K), function(k) {
        if (k == fam$mix$ref) return(0)
        unname(coef_of(paste0("theta", k))[["(Intercept)"]])
      }, 0)
      out[[nm]] <- exp(eta - max(eta)) / sum(exp(eta - max(eta)))
    } else if (grepl("^hu[0-9]*$", nm)) {
      out[[nm]] <- stats::plogis(unname(coef_of(nm)[["(Intercept)"]]))
    } else {
      stop("no translation rule for Stan parameter ", nm)
    }
  }
  out
}

omb_check <- function(f_frm, f_brm, fam_frm, fam_brm, data,
                      patch_delta = FALSE, move = 0.1) {
  # a mixture's optimizer can stop on its relative-function test, or
  # with a component at a degenerate boundary (mixture(sratio(), acat())
  # with cs(x) here), neither of which is what this checks: the density
  # is compared at whatever point the fit reached; brms says cs() is
  # experimental
  fit <- allow_warnings(frm(f_frm, family = fam_frm, data = data),
                        c("singular convergence", "degenerate boundary",
                          "Category specific effects for this family"))
  quiet <- function(e) {
    allow_warnings(suppressMessages(e),
                   "Category specific effects for this family")
  }
  pr <- quiet(brms::get_prior(f_brm, data = data, family = fam_brm))
  pr$prior <- ""
  code <- quiet(as.character(
    brms::stancode(f_brm, data = data, family = fam_brm, prior = pr)))
  if (patch_delta) {
    # brms's equidistant mixture declares a bare `delta` per component
    # and reads delta_mu<k> (upstream defect 4): the program its
    # transformed parameters describe, compiled as such
    lines <- strsplit(code, "\n")[[1]]
    di <- grep("^  real<lower=0> delta;|^  real delta;", lines)
    for (j in seq_along(di)) {
      lines[di[j]] <- sub("delta;", paste0("delta_mu", j, ";"), lines[di[j]])
    }
    code <- paste(lines, collapse = "\n")
  }
  sdat <- quiet(brms::standata(f_brm, data = data, family = fam_brm,
                               prior = pr))
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  K <- family(fit)[["mix"]][["K"]] %||% 1L
  lconst <- if ("theta" %in% omb_par_names(code)) lgamma(K) else 0
  for (mv in c(0, move)) {
    p <- fit$opt$par + mv * sin(seq_along(fit$opt$par) * 1.7)
    est <- fit$obj$env$parList(p)
    up <- rstan::unconstrain_pars(sf, omb_translate(fit, est, code, sdat))
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
    ours <- -fit$obj$fn(p) + lconst
    expect_lt(abs(lp - ours), 1e3 * .Machine$double.eps * abs(ours),
              label = paste(deparse1(f_brm$formula), "move", mv))
  }
}

test_that("ordinal mixtures are brms's compiled density", {
  skip_unless_brms_fit()
  d <- omb_data(20261030)
  bb <- function(f) brms::bf(f)
  B <- asNamespace("brms")
  mix <- function(...) suppressMessages(brms::mixture(...))
  omb_check(bf(y ~ x), bb(y ~ x), mixture(cumulative(), cumulative()),
            mix(B$cumulative(), B$cumulative()), d)
  omb_check(bf(y ~ x), bb(y ~ x), mixture(cumulative("probit"), sratio()),
            mix(B$cumulative("probit"), B$sratio()), d)
  omb_check(bf(y ~ x, disc1 ~ 0 + z, theta1 ~ z),
            brms::bf(y ~ x, disc1 ~ 0 + z, theta1 ~ z),
            mixture(cumulative(), cumulative()),
            mix(B$cumulative(), B$cumulative()), d)
  omb_check(bf(y ~ x), bb(y ~ x),
            mixture(acat(), cratio(threshold = "sum_to_zero"), order = "mu"),
            mix(B$acat(), B$cratio(threshold = "sum_to_zero"),
                order = "mu"), d)
  omb_check(bf(y ~ cs(x)), bb(y ~ cs(x)), mixture(sratio(), acat()),
            mix(B$sratio(), B$acat()), d)
  omb_check(bf(y | thres(gr = g) ~ x), bb(y | thres(gr = g) ~ x),
            mixture(cumulative(), sratio(), order = "mu"),
            mix(B$cumulative(), B$sratio(), order = "mu"), d)
  omb_check(bf(y ~ x), bb(y ~ x),
            mixture(cumulative(threshold = "equidistant"),
                    sratio(threshold = "equidistant")),
            mix(B$cumulative(threshold = "equidistant"),
                B$sratio(threshold = "equidistant")), d,
            patch_delta = TRUE)
  omb_check(bf(yh ~ x, hu1 ~ z), brms::bf(yh ~ x, hu1 ~ z),
            mixture(hurdle_cumulative("probit"),
                    hurdle_cumulative("probit")),
            mix(B$hurdle_cumulative("probit"),
                B$hurdle_cumulative("probit")), d)
})

test_that("hurdle_cumulative() with thres(gr = ) and cs() is brms's", {
  skip_unless_brms_fit()
  d <- omb_data(20261031)
  B <- asNamespace("brms")
  omb_check(bf(yh | thres(gr = g) ~ x, hu ~ z),
            brms::bf(yh | thres(gr = g) ~ x, hu ~ z),
            hurdle_cumulative(), B$hurdle_cumulative(), d)
  omb_check(bf(yh | thres(gr = g) ~ x, disc ~ 0 + z),
            brms::bf(yh | thres(gr = g) ~ x, disc ~ 0 + z),
            hurdle_cumulative("probit"), B$hurdle_cumulative("probit"), d)
  # the probit: brms's logit program reads the top category out of
  # range once cs() or a modeled disc puts it on hurdle_cumulative_
  # logit_lpmf (dev/upstream-bugs.md, brms-1)
  # a smaller move: a larger one crosses two thresholds on some row,
  # where both densities are NaN
  omb_check(bf(yh ~ cs(x), disc ~ 0 + z, hu ~ z),
            brms::bf(yh ~ cs(x), disc ~ 0 + z, hu ~ z),
            hurdle_cumulative("probit"), B$hurdle_cumulative("probit"), d,
            move = 0.02)
})

test_that("cs() on cumulative() is brms's program", {
  # lifted at 0.68.0 (dev/rel068-cs-lpcheck.R has six shapes, a
  # mixture among them: on these two-class data a mixture with cs()
  # stops short of an optimum, where the density is NaN); a small
  # move, since a larger one crosses two thresholds on some row, where
  # both densities are NaN
  skip_unless_brms_fit()
  d <- omb_data(20261032)
  B <- asNamespace("brms")
  omb_check(bf(y ~ z + cs(x), disc ~ 0 + z),
            brms::bf(y ~ z + cs(x), disc ~ 0 + z),
            cumulative("probit"), B$cumulative("probit"), d, move = 0.02)
})

