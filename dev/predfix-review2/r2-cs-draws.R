source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
suppressPackageStartupMessages(library(frmtmb.sample))
cat("frmtmb.sample from", find.package("frmtmb.sample"),
    " StanHeaders", as.character(packageVersion("StanHeaders")),
    " rstan", as.character(packageVersion("rstan")), "\n")
# cs() in posterior_predict() on draws. Seed 808. Per cell, the share of
# draws in category k against the mean of posterior_epred()'s
# probability, z in binomial Monte Carlo units over draws x rows.
set.seed(808)
n <- 300
ng <- 15
g <- factor(rep(seq_len(ng), each = n / ng))
x <- rnorm(n); w <- runif(n, 1, 3)
sim_ord <- function(cuts, slopes, xx, uu) {
  K <- length(cuts) + 1; out <- rep(K, length(xx)); alive <- rep(TRUE, length(xx))
  for (k in seq_along(cuts)) {
    s <- alive & runif(length(xx)) < plogis(cuts[k] + slopes[k] * xx + uu)
    out[s] <- k; alive <- alive & !s
  }
  out
}
d <- data.frame(x, w, g,
                yo = sim_ord(c(-0.5, 0.3, 0.2), c(1.4, -1.0, 0.3), x,
                             rnorm(ng, 0, 0.3)[g]),
                yo2 = sim_ord(c(0.2, -0.4), c(-0.8, 0.9), log(w), 0),
                y2 = 0.5 * x + rnorm(n))
cellz <- function(lab, pp, ep) {
  # pp: draws x rows (categories 1..K); ep: draws x rows x K
  K <- dim(ep)[3]
  z <- vapply(seq_len(K), function(k) {
    obs <- colMeans(pp == k); ex <- colMeans(ep[, , k, drop = FALSE])
    D <- nrow(pp)
    max(abs(obs - ex) / sqrt(pmax(ex * (1 - ex), 1e-9) / D))
  }, 0)
  pooled <- vapply(seq_len(K), function(k) {
    obs <- mean(pp == k); ex <- mean(ep[, , k])
    abs(obs - ex) / sqrt(ex * (1 - ex) / length(pp))
  }, 0)
  cat(sprintf("%-40s draws=%d rows=%d max cell |z|=%.2f pooled max |z|=%.2f\n",
              lab, nrow(pp), ncol(pp), max(z), max(pooled)))
}
samp <- function(fit) suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 1000, warmup = 300, seed = 8,
             refresh = 0)))
nd <- data.frame(x = rep(c(2.5, -2.5), each = 20), w = 1,
                 g = factor("1", levels = levels(g)))
nd1 <- data.frame(x = 2.5, w = 1, g = factor("1", levels = levels(g)))
for (fam in list(sratio(), acat(), cratio())) {
  lab <- fam$family
  fit <- frm(bf(yo ~ cs(x) + (1 | g)), family = fam, data = d)
  ds <- samp(fit)
  for (cfg in list(list("in sample", NULL, NULL),
                   list("newdata 40 rows", nd, NULL),
                   list("newdata 1 row", nd1, NULL),
                   list("newdata re_formula = NA", nd, NA))) {
    set.seed(5)
    pp <- tryCatch(posterior_predict(ds, newdata = cfg[[2]], re_formula = cfg[[3]]),
                   error = function(e) e)
    ep <- tryCatch(posterior_epred(ds, newdata = cfg[[2]], re_formula = cfg[[3]]),
                   error = function(e) e)
    if (inherits(pp, "error") || inherits(ep, "error")) {
      cat(lab, cfg[[1]], "ERROR:", if (inherits(pp, "error")) conditionMessage(pp),
          if (inherits(ep, "error")) conditionMessage(ep), "\n"); next
    }
    if (is.null(dim(pp))) pp <- matrix(pp, ncol = 1)
    cellz(paste(lab, cfg[[1]]), pp, ep)
  }
  # the absent-term control: what a dropped cs() looks like at x = 2.5
  cat(sprintf("   epred P(Y=1) at x = 2.5: %.3f; share drawn: %.3f\n",
              mean(posterior_epred(ds, newdata = nd1)[, , 1]),
              {set.seed(5); mean(posterior_predict(ds, newdata = nd1) == 1)}))
  up <- tryCatch(posterior_predict(ds, newdata = data.frame(x = 1, w = 1, g = "zz"),
                                   allow_new_levels = TRUE),
                 error = function(e) conditionMessage(e))
  cat("   unseen level with allow_new_levels:", substr(paste(up)[1], 1, 120), "\n")
}

cat("\n-- multivariate draws with two cs() responses --\n")
mf <- frm(mvbf(bf(yo ~ cs(x) + (1 | g), family = sratio()),
               bf(yo2 ~ cs(log(w)), family = acat()),
               bf(y2 ~ x, family = gaussian()), rescor = FALSE), data = d)
dm <- tryCatch(samp(mf), error = function(e) e)
if (inherits(dm, "error")) {
  cat("frm_sample on mv: ERROR", conditionMessage(dm), "\n")
} else {
  ndw <- data.frame(x = rep(c(2.5, -2.5), each = 20), w = rep(c(1, 3), 20),
                    g = factor("1", levels = levels(g)))
  for (r in c("yo", "yo2")) for (cfg in list(list("in sample", NULL),
                                             list("newdata 40 rows", ndw),
                                             list("newdata 1 row", ndw[1, ]))) {
    set.seed(5)
    pp <- tryCatch(posterior_predict(dm, newdata = cfg[[2]], resp = r),
                   error = function(e) e)
    ep <- tryCatch(posterior_epred(dm, newdata = cfg[[2]], resp = r),
                   error = function(e) e)
    if (inherits(pp, "error") || inherits(ep, "error")) {
      cat("mv", r, cfg[[1]], "ERROR:", if (inherits(pp, "error")) conditionMessage(pp),
          if (inherits(ep, "error")) conditionMessage(ep), "\n"); next
    }
    if (is.null(dim(pp))) pp <- matrix(pp, ncol = 1)
    cellz(paste("mv", r, cfg[[1]]), pp, ep)
  }
  set.seed(5)
  ppall <- tryCatch(posterior_predict(dm, newdata = ndw), error = function(e) e)
  cat("mv posterior_predict all responses:",
      if (inherits(ppall, "error")) conditionMessage(ppall) else
        paste(dim(ppall), collapse = "x"), "\n")
  if (!inherits(ppall, "error") && length(dim(ppall)) == 3) {
    epo2 <- posterior_epred(dm, newdata = ndw, resp = "yo2")
    cellz("mv yo2 from the all-responses call", ppall[, , "yo2"], epo2)
  }
}
