# Lane sratio: where sratio()'s thresholds do not cross, making them
# unordered must leave the maximum likelihood fit unchanged. This fits
# the sratio models of the existing tests (their data constructors are
# copied here, with the tests' seeds) and writes, per model, the
# log-likelihood, every fixef() estimate and the smallest gap between
# adjacent thresholds, at 17 significant digits.
#
# Run it once on the build before the change and once after:
#
#   Rscript dev/sratio-invariance.R before /opt/rlib/lane-sratio \
#     > dev/sratio-invariance-before-log.txt 2>&1
#   Rscript dev/sratio-invariance.R after /opt/rlib/lane-sratio \
#     > dev/sratio-invariance-after-log.txt 2>&1
#
# The "after" run reads dev/sratio-invariance-before.csv and prints the
# differences. The before run used the consolidated branch at a897706d
# installed into the lane library, before any edit.
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
.libPaths(c(args[2], "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")

sim_ord_data <- function(seed = 91, n = 600) {       # test-ordinal.R
  set.seed(seed)
  x <- rnorm(n)
  p <- plogis(outer(rep(1, n), c(-1, 0.3, 1.4)) - 1.2 * x)
  data.frame(y = rowSums(runif(n) > cbind(p, 1)) + 1L, x = x)
}
ord3_data <- function(seed, n = 250, z = FALSE) {    # v29, ordinal-fitted
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  if (z) dd$z <- stats::rnorm(n)
  eta <- 0.9 * dd$x
  p <- cbind(plogis(-0.8 - eta), plogis(0.6 - eta) - plogis(-0.8 - eta),
             1 - plogis(0.6 - eta))
  dd$y <- factor(apply(p, 1L, function(pr) sample(3L, 1L, prob = pr)),
                 levels = 1:3, ordered = TRUE)
  dd
}
shapes_data <- function(seed) {                      # brms-shapes(-punch)
  set.seed(seed)
  n <- if (seed == 20260917) 120 else 150
  ng <- n / 10
  dd <- data.frame(x = rnorm(n), z = rnorm(n),
                   g = factor(rep(seq_len(ng), each = 10)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(ng, 0, 0.7)[dd$g], 1)
  dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n),
                       c(-Inf, -0.3, 0.8, Inf), labels = 1:3),
                   ordered = TRUE)
  dd
}
mv_gap_data <- function(seed = 5, n = 150) {         # test-mv-gaps.R
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  data.frame(x = x, z = z, y1 = 1 + 0.5 * x + rnorm(n),
             y2 = -1 + rnorm(n), y3 = 0.3 * x + rnorm(n),
             o = cut(x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE),
             o2 = cut(-x + rlogis(n), c(-Inf, 0, 1, Inf), labels = FALSE))
}
thres_data <- function(seed = 11, n = 240, id = TRUE) {  # test-thres.R
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
  if (id) d$id <- factor(sample(1:20, n, TRUE))
  re <- if (id) rnorm(20, 0, 0.4)[d$id] else 0
  u <- rlogis(n, 0.7 * d$x + re)
  tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6),
              c = c(-1.5, -0.4, 0.5, 1.6))
  d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]),
                1L)
  d
}
lik_data <- function() {                             # brms-likelihood
  set.seed(5)
  n <- 300
  do <- data.frame(x = rnorm(n), z = rnorm(n))
  do$y <- ordered(cut(0.9 * do$x + rlogis(n),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  do
}

d91 <- sim_ord_data()
d300 <- sim_ord_data(n = 300)
dth <- thres_data()
ddm <- mv_gap_data()
dsh <- shapes_data(20260917)
dpu <- shapes_data(20260918)
dre <- d91
set.seed(3)
dre$g <- factor(sample(1:30, nrow(dre), TRUE))
dre$w <- sample(1:3, nrow(dre), TRUE)
data("inhaler", package = "brms")

models <- list(
  logit = quote(frm(bf(y ~ x) + sratio(), data = d91)),
  probit = quote(frm(bf(y ~ x), family = sratio("probit"), data = d300)),
  cloglog = quote(frm(bf(y ~ x), family = sratio("cloglog"), data = d300)),
  cauchit = quote(frm(bf(y ~ x), family = sratio("cauchit"), data = d300)),
  probit_approx = quote(frm(bf(y ~ x), family = sratio("probit_approx"),
                            data = d300)),
  v29_42 = quote(frm(bf(y ~ x) + sratio(), data = ord3_data(42))),
  v29_43_cs = quote(frm(bf(y ~ x + cs(x)) + sratio(), data = ord3_data(43))),
  fitted_102_cs = quote(frm(bf(y ~ x + cs(z)) + sratio(),
                            data = ord3_data(102, z = TRUE))),
  shapes_cs = quote(frm(bf(ord ~ x + cs(z)) + sratio(), data = dsh)),
  punch_cs = quote(frm(bf(ord ~ x + cs(z)) + sratio(), data = dpu)),
  likelihood = quote(frm(bf(y ~ x) + sratio(), data = lik_data())),
  likelihood_cs = quote(frm(bf(y ~ x + cs(z)) + sratio(),
                            data = lik_data())),
  inhaler_cs = quote(frm(bf(rating ~ period + carry + cs(treat)) + sratio(),
                         data = inhaler)),
  ranef = quote(frm(bf(y ~ x + (1 | g)) + sratio(), data = dre)),
  weights = quote(frm(bf(y | weights(w) ~ x) + sratio(), data = dre)),
  mv = quote(frm(bf(o ~ x) + cumulative() + bf(o2 ~ x) + sratio() +
                   bf(y1 ~ x) + gaussian(), data = ddm)),
  thres_gr = quote(frm(y | thres(gr = g) ~ x, data = dth, family = sratio())),
  thres_gr_noid = quote(frm(y | thres(gr = g) ~ x,
                            data = thres_data(id = FALSE),
                            family = sratio())),
  thres_5 = quote(suppressWarnings(frm(y | thres(5) ~ x, data = dth,
                                       family = sratio())))
)

# the smallest gap between adjacent thresholds of each vector: rows are
# Intercept[k], Intercept[g,k] or <resp>_Intercept[k]
min_gap <- function(fe) {
  nm <- rownames(fe)
  i <- grep("Intercept\\[", nm)
  key <- sub("\\[([^,]*,)?[0-9]+\\]$", "[\\1]", nm[i])
  v <- fe[i, "Estimate"]
  g <- tapply(v, key, function(t) if (length(t) > 1) min(diff(t)) else NA)
  paste(sprintf("%s %.3g", names(g), g), collapse = ", ")
}

rows <- list()
for (nm in names(models)) {
  fit <- eval(models[[nm]])
  fe <- fixef(fit)
  rows[[nm]] <- data.frame(model = nm, par = c("logLik", rownames(fe)),
                           value = c(as.numeric(logLik(fit)),
                                     fe[, "Estimate"]),
                           se = c(NA, fe[, "Est.Error"]))
  cat(sprintf("%-14s logLik %.17g conv %d max|gr| %.2e min gap: %s\n",
              nm, as.numeric(logLik(fit)), fit$opt$convergence,
              max(abs(fit$obj$gr(fit$opt$par))), min_gap(fe)))
}
res <- do.call(rbind, rows)
rownames(res) <- NULL

if (identical(arm, "before")) {
  utils::write.csv(res[c("model", "par", "value")],
                   "dev/sratio-invariance-before.csv", row.names = FALSE)
} else {
  b <- utils::read.csv("dev/sratio-invariance-before.csv")
  m <- merge(b, res, by = c("model", "par"), suffixes = c(".before",
                                                          ".after"),
             sort = FALSE)
  stopifnot(nrow(m) == nrow(b), nrow(m) == nrow(res))
  m$diff <- m$value.after - m$value.before
  cat("\nper model: |logLik after - before| / |logLik|, the largest",
      "|estimate after - before| over its fixef() rows, and the largest",
      "of that difference over the row's standard error (after; NaN",
      "where the Hessian is singular)\n")
  for (nm in unique(m$model)) {
    s <- m[m$model == nm, ]
    ll <- s[s$par == "logLik", ]
    e <- s[s$par != "logLik", ]
    cat(sprintf(paste("%-14s rel dlogLik %.2e  max |d estimate| %.2e",
                      " max |d| / se %.2e  (%d rows)\n"),
                nm, abs(ll$diff) / abs(ll$value.before),
                max(abs(e$diff)), max(abs(e$diff) / e$se), nrow(e)))
  }
}
