# Reviewer of lane seflags, claim 3: false alarms of the curved-ridge
# correction on IDENTIFIED models, weak but identified. Per fit: code,
# lost set, SE warnings, max |SE / ref - 1| over the fixed coefficients
# (ref: lm()/lmer()/glm() where one exists, else the inverse of the
# exact AD Hessian at the optimum, which is the right answer for an
# identified model without random effects), and, on the lane, every
# se_ridge_curvature() call: eigenvalue, corrected curvature, and
# whether the direction was then removed.
#   Rscript dev/seflags-rev-falarm.R <lane|base> <group> <out.tsv> [seeds]
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
a <- commandArgs(trailingOnly = TRUE)
grp <- a[2]
out <- a[3]
seeds <- if (length(a) > 3) eval(parse(text = a[4])) else 1:20
ridge_log <- list()
if (exists("se_ridge_curvature", ns)) {
  trace("se_ridge_curvature", exit = quote({
    rv <- returnValue()
    ridge_log[[length(ridge_log) + 1L]] <<- c(ev = e$values[k], c0 = rv)
  }), print = FALSE, where = ns)
}
cat("group\tmodel\tseed\tcode\tlost\tse_warn\tother_warn\tmax_rel\tn_ridge\tmin_c0_over_ev\tfit_s\n",
    file = out)
rec <- function(model, seed, expr, ref = NULL, coef_names = NULL) {
  ridge_log <<- list()
  t0 <- proc.time()[["elapsed"]]
  r <- cap(expr)
  t1 <- proc.time()[["elapsed"]] - t0
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(grp, model, seed, "ERROR", gsub("[\t\n]", " ", substr(fit, 1, 80)),
        NA, NA, NA, NA, NA, round(t1, 2), sep = "\t", file = out,
        append = TRUE)
    cat("\n", file = out, append = TRUE)
    return(invisible())
  }
  se <- se_outer(fit)
  on <- names(se)
  if (is.null(ref)) {
    # the exact AD Hessian at the optimum (models without random effects)
    he <- tryCatch(fit$obj$he(fit$opt$par), error = function(e) NULL)
    fit$obj$fn(fit$opt$par)
    if (!is.null(he)) {
      ref <- tryCatch(sqrt(diag(solve(he))), error = function(e) NULL)
      if (!is.null(ref)) names(ref) <- on
    }
  }
  rel <- NA
  if (!is.null(ref)) {
    nmz <- if (is.null(coef_names)) names(ref) else coef_names
    s1 <- se[match(nmz, on)]
    s0 <- if (is.null(coef_names)) ref else ref[seq_along(nmz)]
    rel <- max(abs(s1 / s0 - 1))
  }
  sew <- sum(grepl("^Standard errors are not available", r$w))
  rl <- if (length(ridge_log)) do.call(rbind, ridge_log) else NULL
  mr <- if (is.null(rl)) NA else min(rl[, "c0"] / rl[, "ev"], na.rm = TRUE)
  cat(grp, model, seed, fit$opt$convergence, lost_of(fit), sew,
      length(r$w) - sew, signif(rel, 4), NROW(rl), signif(mr, 4),
      round(t1, 2), sep = "\t", file = out, append = TRUE)
  cat("\n", file = out, append = TRUE)
}
lmer_ref <- function(fo, d) {
  m <- suppressWarnings(suppressMessages(lme4::lmer(fo, data = d,
                                                    REML = FALSE)))
  sqrt(diag(as.matrix(stats::vcov(m))))
}

for (s in seeds) {
  set.seed(s)
  if (grp == "nl") {
    # Emax with doses far below ED50: only emax / ED50 is well determined
    dose <- rep(c(0, 5, 10, 25, 50), each = 4)
    d <- data.frame(dose = dose, y = 2 + 30 * dose / (100 + dose) +
                      rnorm(20, 0, 0.5))
    rec("emax_low_doses", s, frm(bf(y ~ e0 + emax * dose / (exp(led) + dose),
                                    e0 ~ 1, emax ~ 1, led ~ 1, nl = TRUE),
                                 data = d,
                                 start = list(beta = c(2, 30, log(100)))))
    # Michaelis-Menten, 6 points (Puromycin-like concentrations)
    S <- c(0.02, 0.06, 0.11, 0.22, 0.56, 1.1)
    d <- data.frame(S = S, v = 200 * S / (0.06 + S) + rnorm(6, 0, 8))
    rec("mm_6pts", s, frm(bf(v ~ vm * S / (exp(lk) + S), vm ~ 1, lk ~ 1,
                             nl = TRUE), data = d,
                          start = list(beta = c(200, log(0.06)))))
    # logistic growth, only the early part observed
    tt <- seq(0, 3, length.out = 15)
    d <- data.frame(tt = tt, y = 10 / (1 + exp((4 - tt) / 1)) +
                      rnorm(15, 0, 0.1))
    rec("logistic_early", s, frm(bf(y ~ asym / (1 + exp((xmid - tt) /
                                                         exp(ls))),
                                    asym ~ 1, xmid ~ 1, ls ~ 1, nl = TRUE),
                                 data = d,
                                 start = list(beta = c(10, 4, 0))))
    # exponential decay over a short range
    tt <- seq(0, 0.3, length.out = 20)
    d <- data.frame(tt = tt, y = 5 * exp(-1 * tt) + rnorm(20, 0, 0.05))
    rec("decay_short", s, frm(bf(y ~ y0 * exp(-exp(lk) * tt), y0 ~ 1,
                                 lk ~ 1, nl = TRUE), data = d,
                              start = list(beta = c(5, 0))))
    # psychometric, 80 trials over a narrow range
    x <- runif(80, -0.5, 0.5)
    d <- data.frame(x = x, y = rbinom(80, 1, pnorm(x, 0.1, 0.6)))
    rec("psychometric_narrow", s,
        frm(bf(y ~ pnorm(x, m0, exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
            family = bernoulli(link = "identity"), data = d,
            start = list(beta = c(0, log(0.6)))))
    # lapse psychometric, 300 trials
    x <- runif(300, 200, 400)
    d <- data.frame(x = x, y = rbinom(300, 1, 0.04 * 0.5 + 0.96 *
                                        pnorm(x, 300, 25)))
    rec("lapse", s, frm(bf(y ~ lapse * 0.5 + (1 - lapse) *
                             pnorm(x, m0, exp(ls)),
                           lapse ~ 1, m0 ~ 1, ls ~ 1, nl = TRUE),
                        family = bernoulli(link = "identity"), data = d,
                        start = list(beta = c(0.05, 300, log(25)))))
    # calendar-year logistic
    yr <- seq(1900, 2000, length.out = 60)
    d <- data.frame(yr = yr, y = 50 / (1 + exp((1950 - yr) / 12)) +
                      rnorm(60, 0, 3))
    rec("logistic_years", s, frm(bf(y ~ Asym / (1 + exp((xmid - yr) /
                                                         exp(lscal))),
                                    Asym ~ 1, xmid ~ 1, lscal ~ 1,
                                    nl = TRUE), data = d,
                                 start = list(beta = c(50, 1950,
                                                       log(12)))))
    # two-compartment-like biexponential, weak second rate
    tt <- seq(0.1, 8, length.out = 25)
    d <- data.frame(tt = tt, y = 4 * exp(-1.5 * tt) + 1 * exp(-0.2 * tt) +
                      rnorm(25, 0, 0.03))
    rec("biexp", s, frm(bf(y ~ a1 * exp(-exp(l1) * tt) +
                             a2 * exp(-exp(l2) * tt),
                           a1 ~ 1, l1 ~ 1, a2 ~ 1, l2 ~ 1, nl = TRUE),
                        data = d,
                        start = list(beta = c(4, log(1.5), 1, log(0.2)))))
  } else if (grp == "lin") {
    n <- 100
    x1 <- rnorm(n)
    for (r in c(0.999, 0.99999, 0.9999999)) {
      x2 <- r * x1 + sqrt(1 - r^2) * rnorm(n)
      d <- data.frame(x1 = x1, x2 = x2, g = factor(rep(1:10, 10)))
      d$y <- 1 + x1 + x2 + rnorm(10, 0, 0.5)[d$g] + rnorm(n)
      m <- lm(y ~ x1 + x2, data = d)
      rec(paste0("lm_cor", r), s, frm(y ~ x1 + x2, data = d),
          ref = sqrt(diag(vcov(m))), coef_names = c("(Intercept)", "x1", "x2"))
      rec(paste0("lmer_cor", r), s, frm(y ~ x1 + x2 + (1 | g), data = d),
          ref = lmer_ref(y ~ x1 + x2 + (1 | g), d),
          coef_names = c("(Intercept)", "x1", "x2"))
      d$c <- rpois(n, exp(0.5 + 0.3 * x1 + 0.3 * x2))
      m <- glm(c ~ x1 + x2, family = poisson(), data = d)
      rec(paste0("pois_cor", r), s, frm(c ~ x1 + x2, family = poisson(),
                                       data = d),
          ref = sqrt(diag(vcov(m))), coef_names = c("(Intercept)", "x1", "x2"))
    }
    # uncentered covariate, glm and glmm
    x <- rnorm(n, 1000, 1)
    d <- data.frame(x = x, g = factor(rep(1:10, 10)))
    d$y <- rbinom(n, 1, plogis(0.3 * (x - 1000)))
    m <- glm(y ~ x, family = binomial(), data = d)
    rec("logit_uncentered", s, frm(y ~ x, family = bernoulli(), data = d),
        ref = sqrt(diag(vcov(m))), coef_names = c("(Intercept)", "x"))
  } else if (grp == "poly") {
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
    for (deg in 3:6) {
      rhs <- paste(paste0("x", 1:deg), collapse = " + ")
      fo <- stats::as.formula(paste("y ~", rhs))
      cn <- c("(Intercept)", paste0("x", 1:deg))
      rec(paste0("poly", deg, "_lm"), s, frm(fo, data = d),
          ref = sqrt(diag(vcov(lm(fo, data = d)))), coef_names = cn)
      fo2 <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
      rec(paste0("poly", deg, "_lmer"), s, frm(fo2, data = d),
          ref = lmer_ref(fo2, d), coef_names = cn)
    }
  } else if (grp == "mixed") {
    for (sdg in c(0.02, 0.05, 0.1)) {
      d <- data.frame(x = rnorm(300), g = factor(rep(1:30, 10)))
      d$y <- 1 + 0.5 * d$x + rnorm(30, 0, sdg)[d$g] + rnorm(300)
      rec(paste0("ri_sd", sdg), s, frm(y ~ x + (1 | g), data = d),
          ref = lmer_ref(y ~ x + (1 | g), d), coef_names = c("(Intercept)",
                                                             "x"))
      d$y2 <- 1 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g] +
        rnorm(30, 0, sdg)[d$g] * d$x + rnorm(300)
      rec(paste0("rs_sd", sdg), s, frm(y2 ~ x + (1 + x | g), data = d),
          ref = lmer_ref(y2 ~ x + (1 + x | g), d),
          coef_names = c("(Intercept)", "x"))
      d$c <- rpois(300, exp(0.3 + 0.3 * d$x + rnorm(30, 0, sdg)[d$g]))
      m <- suppressWarnings(suppressMessages(lme4::glmer(
        c ~ x + (1 | g), family = poisson(), data = d)))
      rec(paste0("pois_ri_sd", sdg), s, frm(c ~ x + (1 | g),
                                            family = poisson(), data = d),
          ref = sqrt(diag(as.matrix(vcov(m)))),
          coef_names = c("(Intercept)", "x"))
    }
  } else if (grp == "bounded") {
    n <- 150
    x <- rnorm(n)
    yy <- 0.2 + 0.5 * x + rnorm(n)
    d <- data.frame(x = x, y = yy)
    dt <- d[d$y > 0, ]
    rec("trunc_lb0", s, frm(y | trunc(lb = 0) ~ x, data = dt))
    d$cen <- ifelse(d$y > 1, "right", "none")
    d$yc <- pmin(d$y, 1)
    rec("cens_right", s, frm(yc | cens(cen) ~ x, data = d))
    mu <- plogis(-2 + 0.5 * x)
    yb <- rbeta(n, mu * 200, (1 - mu) * 200)
    yb[sample(n, 3)] <- 0
    db <- data.frame(x = x, y = yb)
    rec("zoib_phi200", s, frm(y ~ x, family = zero_one_inflated_beta(),
                              data = db))
    yb2 <- rbeta(n, 0.98 * 20, 0.02 * 20)
    yb2[sample(n, 2)] <- 1
    rec("zoib_near1", s, frm(y ~ 1, family = zero_one_inflated_beta(),
                             data = data.frame(y = yb2)))
    yt <- 1 + 0.5 * x + rt(n, df = 30)
    rec("student_nu30", s, frm(y ~ x, family = student(),
                               data = data.frame(x = x, y = yt)))
    rec("negbin_large_shape", s,
        frm(c ~ x, family = negbinomial(),
            data = data.frame(x = x, c = rpois(n, exp(1 + 0.3 * x)))))
  }
}
cat("done", grp, "\n")
