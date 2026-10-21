# Release reviewer, 0.73.0: frmtmb.sample's loo_*() functions on draws of
# censored and truncated models, negbinomial2 and dirichlet_multinomial.
# Draws are jittered around the ML estimate as in frmtmb.sample's
# test-loo-predict.R (no sampler). With sd = 0 every draw is the
# estimate, so log_lik() must sum to logLik() and posterior_epred() must
# equal fitted().
#   Rscript dev/rel073-rev-loo.R <lib>   (log dev/rel073-rev-log/loo-<tag>.txt)
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib: frmtmb", format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), find.package("frmtmb"), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
mk_draws <- function(fit, nd = 200L, sd = 0.03, seed = 5L) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) {
    v <- fit$estimates[[cp]]
    fx <- fit$frame[["betad_fixed_idx"]]
    if (identical(cp, "betad") && length(fx)) v <- v[-fx]
    v
  }))
  set.seed(seed)
  M <- matrix(rep(est, each = nd) + stats::rnorm(nd * length(est), 0, sd),
              nd, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
tryv <- function(label, expr) {
  r <- tryCatch(q(expr), error = function(e)
    structure(gsub("\\s+", " ", substr(conditionMessage(e), 1, 300)),
              class = "err"))
  if (inherits(r, "err")) cat("  ", label, ": ERROR:", r, "\n")
  invisible(r)
}
set.seed(20261021)
n <- 120
d <- data.frame(x = rnorm(n), N = 15L)
d$k <- rpois(n, exp(1 + 0.4 * d$x))
d$cc <- ifelse(d$k >= 5, "right", "none"); d$kc <- pmin(d$k, 5)
d$nb <- rnbinom(n, mu = exp(1.3 + 0.4 * d$x), size = 3)
d$cnb <- ifelse(d$nb >= 6, "right", ifelse(d$nb <= 1, "left", "none"))
d$yb <- rbeta(n, 3, 4)
P <- t(apply(cbind(0, 0.4 + 0.3 * d$x, -0.3), 1,
             function(e) exp(e) / sum(exp(e))))
d$Y <- t(vapply(seq_len(n), function(i) {
  p <- rgamma(3, P[i, ] * 5); rmultinom(1, 15, p / sum(p))
}, numeric(3)))
dt <- d[d$k >= 1, ]
cases <- list(
  pois_cens = list(bf(kc | cens(cc) ~ x), poisson(), d),
  pois_trunc = list(bf(k | trunc(lb = 1) ~ x), poisson(), dt),
  nb2_cens = list(bf(nb | cens(cnb) ~ x), brmsfamily("negbinomial2"), d),
  nb2_trunc = list(bf(nb | trunc(lb = 1, ub = 40) ~ x),
                   brmsfamily("negbinomial2"), d[d$nb >= 1, ]),
  nb2 = list(bf(nb ~ x), brmsfamily("negbinomial2"), d),
  beta_trunc = list(bf(yb | trunc(lb = 0.05) ~ x), Beta(), d[d$yb > 0.05, ]),
  dm = list(bf(Y | trials(N) ~ x), dirichlet_multinomial(), d))
dr <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- q(frm(cs[[1]], family = cs[[2]], data = cs[[3]]))
  d0 <- mk_draws(fit, nd = 20L, sd = 0)
  ll0 <- tryv("log_lik sd 0", log_lik(d0))
  cat(sprintf("%-10s logLik %.8f | log_lik() at the estimate, row sums: %s\n",
              nm, as.numeric(logLik(fit)),
              if (inherits(ll0, "err")) "ERROR" else
                paste(format(range(rowSums(ll0)), nsmall = 8), collapse = " to ")))
  fe <- tryv("fitted", fitted(fit))
  pe <- tryv("posterior_epred sd 0", posterior_epred(d0))
  if (!inherits(fe, "err") && !inherits(pe, "err")) {
    fv <- if (length(dim(fe)) == 3L) fe[, 1, ] else fe[, 1]
    pv <- if (length(dim(pe)) == 3L) pe[1, , ] else pe[1, ]
    cat(sprintf("   posterior_epred at the estimate vs fitted(): max rel diff %.2e\n",
                max(abs(pv / fv - 1))))
  } else if (inherits(fe, "err") != inherits(pe, "err")) {
    cat("   MISMATCH: fitted()", if (inherits(fe, "err")) "refuses" else "runs",
        "but posterior_epred() of draws",
        if (inherits(pe, "err")) "refuses" else "runs", "\n")
  }
  dd <- mk_draws(fit)
  dr[[nm]] <- dd
  for (fn in c("loo", "loo_epred", "loo_linpred", "loo_predict",
               "loo_predictive_interval", "loo_R2")) {
    r <- tryv(fn, get(fn)(dd))
    if (!inherits(r, "err")) {
      s <- if (fn == "loo") sprintf("elpd_loo %.3f", r$estimates[1, 1]) else
        paste("dim", paste(dim(as.array(r)), collapse = "x"),
              "| finite", all(is.finite(unlist(r))))
      cat("  ", fn, ":", s, "\n")
    }
  }
  if (nm == "pois_trunc" || nm == "nb2_trunc") {
    le <- loo_epred(dd)
    cat(sprintf("   loo_epred / fitted(): range %.4f to %.4f\n",
                min(le / fitted(fit)[, 1]), max(le / fitted(fit)[, 1])))
  }
}
w <- tryv("loo_model_weights nb2 vs nb2_cens", loo_model_weights(dr$nb2,
                                                                  dr$nb2))
if (!inherits(w, "err")) print(w)
cat("END\n")
