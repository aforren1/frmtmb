# Lane latperf: the sparse joint covariance against the dense route on the
# same fit, in one process (dev/rtmb-pitfalls.md item 21). The dense route is
# the one get_joint_cov() took before: jcs_joint_cov() mocked to NULL, with a
# cache that keeps the fit's sdreport and drops the joint covariance.
# Usage: Rscript dev/latperf-agree.R   (lane library first). Seeds below.
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat(sprintf("lib: %s frmtmb %s\n", dirname(dirname(system.file(package =
  "frmtmb"))), format(packageVersion("frmtmb"))))
ns <- asNamespace("frmtmb")

mk <- function(seed, G, K, sd = 0.6) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  d$s <- factor(sample.int(6, n, TRUE))
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.6), K, sd = sd)))
  eta <- 0.6 + 0.4 * d$x + e + c(-0.3, 0.2, 0, 0.1, -0.1, 0.3)[d$s]
  d$cnt <- stats::rpois(n, exp(eta))
  d$tr <- 8L
  d$k <- stats::rbinom(n, d$tr, stats::plogis(eta - 0.6))
  d$zc <- ifelse(stats::runif(n) < 0.2, 0L, d$cnt)
  d$ord <- cut(eta + stats::rlogis(n), c(-Inf, 0, 1, 2, Inf),
               ordered_result = TRUE)
  d
}

dense_copy <- function(fit) {
  f2 <- fit
  keep <- setdiff(ls(fit$cache), c("Vjoint"))
  f2$cache <- list2env(mget(keep, envir = fit$cache))
  f2
}
with_dense <- function(expr) {
  old <- get("jcs_joint_cov", envir = ns)
  unlockBinding("jcs_joint_cov", ns)
  assign("jcs_joint_cov", function(fit) NULL, envir = ns)
  on.exit({
    assign("jcs_joint_cov", old, envir = ns)
    lockBinding("jcs_joint_cov", ns)
  })
  expr
}
rel <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  c(max = max(abs(a[ok] / b[ok] - 1)), n = sum(ok),
    na_mismatch = sum(is.finite(a) != is.finite(b)))
}

cases <- list(
  list("ar1 poisson", "cnt ~ x + ar(t, g, cov = TRUE)", poisson()),
  list("ar1 + (1|s) poisson", "cnt ~ x + (1 | s) + ar(t, g, cov = TRUE)",
       poisson()),
  list("ma1 poisson", "cnt ~ x + ma(t, g, cov = TRUE)", poisson()),
  list("arma1 poisson", "cnt ~ x + arma(t, g, cov = TRUE)", poisson()),
  list("cosy poisson", "cnt ~ x + cosy(t, g)", poisson()),
  list("unstr poisson", "cnt ~ x + unstr(t, g)", poisson()),
  list("condar poisson", "cnt ~ x + ar(t, g, p = 2)", poisson()),
  list("ar1 binomial trials", "k | trials(tr) ~ x + ar(t, g, cov = TRUE)",
       binomial()),
  list("ar1 zip", "zc ~ x + ar(t, g, cov = TRUE)", zero_inflated_poisson()),
  list("ar1 cumulative", "ord ~ x + ar(t, g, cov = TRUE)", cumulative())
)
d <- mk(11, 30, 12)
nd <- data.frame(t = 1:4, g = factor(rep("new", 4)), x = c(-1, 0, 0.5, 1),
                 s = factor(c(1, 2, 3, 4), levels = levels(d$s)), tr = 8L)
for (cs in cases) {
  fit <- suppressWarnings(suppressMessages(
    frm(stats::as.formula(cs[[2]]), family = cs[[3]], data = d)))
  fd <- dense_copy(fit)
  sparse <- frmtmb:::jc_is_sparse(frmtmb:::get_joint_cov(fit))
  a1 <- suppressWarnings(fitted(fit))
  b1 <- with_dense(suppressWarnings(fitted(fd)))
  ee <- function(v) if (is.matrix(v) && "Est.Error" %in% colnames(v)) {
    v[, "Est.Error"]
  } else if (length(dim(v)) == 3L) v[, "Est.Error", ] else v
  r_fit <- rel(ee(a1), ee(b1))
  a2 <- frm_linpred(fit, se.fit = TRUE)$se.fit
  b2 <- with_dense(frm_linpred(fd, se.fit = TRUE)$se.fit)
  r_lp <- rel(a2, b2)
  if (grepl("cumulative", cs[[1]])) {
    r_res <- c(NA, 0, 0)
  } else {
  a3 <- suppressWarnings(residuals(fit))
  b3 <- with_dense(suppressWarnings(residuals(fd)))
  r_res <- rel(ee(a3), ee(b3))
  }
  r_nd <- tryCatch({
    if (grepl("condar", cs[[1]])) c(max = NA, n = 0, na_mismatch = 0) else {
      a4 <- suppressWarnings(fitted(fit, newdata = nd))
      b4 <- with_dense(suppressWarnings(fitted(fd, newdata = nd)))
      rel(ee(a4), ee(b4))
    }
  }, error = function(e) {
    cat("  newdata error:", conditionMessage(e), "\n")
    c(max = NA, n = 0, na_mismatch = NA)
  })
  cat(sprintf(paste("RESULT %-22s sparse %s | fitted %.2e (%d, mis %d) |",
                    "linpred %.2e (%d, mis %d) | resid %.2e (%d, mis %d) |",
                    "newdata %.2e (%d, mis %d)\n"),
              cs[[1]], sparse, r_fit[1], r_fit[2], r_fit[3], r_lp[1],
              r_lp[2], r_lp[3], r_res[1], r_res[2], r_res[3], r_nd[1],
              r_nd[2], r_nd[3]))
}
