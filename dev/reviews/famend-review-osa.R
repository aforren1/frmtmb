# Review of lane famend, items 3618 and 5726: censored-count OSA residuals
# on cases the lane's tests do not cover, and fractional weights.
#   Rscript dev/reviews/famend-review-osa.R
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
osa_r <- function(fit, ...) residuals(fit, type = "osa", ...)[, "Estimate"]
dp1 <- function(fit, n) {
  lapply(frmtmb:::eval_dpars(fit)[[1L]], function(v) rep_len(as.numeric(v), n))
}
in_step <- function(r, lo, hi) {
  qlo <- qnorm(lo)
  qhi <- qnorm(pmin(hi, 1))
  wid <- qhi - qlo
  wid[!is.finite(wid)] <- 1
  is.finite(r) & r >= qlo - 1e-8 * wid & r <= qhi + 1e-8 * wid
}
# exact given TMB's uniforms: osa_uniforms() draws one runif per residual
# row after set.seed(123)
exact <- function(F, y, lo, hi, u) {
  U <- local({ set.seed(123); runif(sum(u)) })
  W <- F(hi) - F(lo)
  qnorm((F(y) - F(lo) - U * (F(y) - F(y - 1))) / W)
}
cap <- function(expr) {
  w <- character()
  v <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (length(w)) for (x in unique(w)) cat("   W:", substr(x, 1, 150), "\n")
  v
}
report <- function(label, r, u, ref = NULL, ok = NULL) {
  if (inherits(r, "err")) {
    cat(sprintf("%-46s ERROR %s\n", label, substr(r, 1, 120)))
    return(invisible())
  }
  cat(sprintf("%-46s censored NA %d of %d; uncensored finite %d of %d",
              label, sum(is.na(r[!u])), sum(!u), sum(is.finite(r[u])), sum(u)))
  if (!is.null(ref)) {
    cat(sprintf("; max rel |r - exact| %.3g", max(abs(r[u] - ref)) /
                  max(abs(ref))))
  }
  if (!is.null(ok)) cat(sprintf("; in step %d of %d", sum(ok), length(ok)))
  cat("\n")
}

set.seed(31)
n <- 300
x <- rnorm(n)
mu <- exp(1.2 + 0.4 * x)
yy <- rpois(n, mu)

# 1. trunc(lb = 1) and right censoring at 6: window (0, 6]
d <- data.frame(x = x, yy = yy)
d <- d[d$yy >= 1, ]
d$y <- pmin(d$yy, 6)
d$cc <- ifelse(d$yy > 6, "right", "none")
f <- cap(frm(y | trunc(lb = 1) + cens(cc) ~ x, family = poisson(), data = d))
if (!inherits(f, "err")) {
  r <- cap(osa_r(f))
  u <- d$cc == "none"
  m <- dp1(f, nrow(d))$mu[u]
  report("poisson trunc(lb=1) + right cens at 6", r, u,
         exact(function(q) ppois(q, m), d$y[u], 0, 6, u))
}

# 2. left censoring only at 2, negbinomial: window (2, Inf)
d <- data.frame(x = x, yy = rnbinom(n, size = 2, mu = mu))
d$y <- pmax(d$yy, 2)
d$cc <- ifelse(d$yy <= 2, "left", "none")
f <- cap(frm(y | cens(cc) ~ x, family = negbinomial(), data = d))
r <- cap(osa_r(f))
u <- d$cc == "none"
dp <- dp1(f, n)
Fn <- function(q) pnbinom(q, size = dp$shape[u], mu = dp$mu[u])
report("negbinomial left cens at 2", r, u,
       exact(Fn, d$y[u], 2, Inf, u))

# 3. the explicit oneStepGeneric route on a censored poisson, in step
d <- data.frame(x = x, yy = yy)
d$y <- pmin(pmax(d$yy, 1), 6)
d$cc <- ifelse(d$yy <= 1, "left", ifelse(d$yy > 6, "right", "none"))
f <- cap(frm(y | cens(cc) ~ x, family = poisson(), data = d))
u <- d$cc == "none"
m <- dp1(f, n)$mu[u]
W <- ppois(6, m) - ppois(1, m)
r <- cap(osa_r(f, osa_method = "oneStepGeneric"))
if (!inherits(r, "err")) {
  ok <- in_step(r[u], (ppois(d$y[u] - 1, m) - ppois(1, m)) / W,
                (ppois(d$y[u], m) - ppois(1, m)) / W)
  report("poisson (1, 6] oneStepGeneric", r, u, ok = ok)
} else report("poisson (1, 6] oneStepGeneric", r, u)
r2 <- cap(osa_r(f))
report("poisson (1, 6] default (cdf)", r2, u,
       exact(function(q) ppois(q, m), d$y[u], 1, 6, u))
# the wrong window [lo, hi]: how many rows would an in-step test catch
okw <- in_step(r2[u], (ppois(d$y[u] - 1, m) - ppois(0, m)) /
                 (ppois(6, m) - ppois(0, m)),
               (ppois(d$y[u], m) - ppois(0, m)) / (ppois(6, m) - ppois(0, m)))
cat("   the same residuals in the step of the wrong window [1, 6]:",
    sum(okw), "of", length(okw), "\n")

# 4. an uncensored count at the left censoring point itself
d4 <- d
d4$y[which(u)[1]] <- 1
r <- cap(osa_r(cap(frm(y | cens(cc) ~ x, family = poisson(), data = d4))))
report("uncensored row at lo (refusal expected)", r, d4$cc == "none")

# 5. integer weights 2 on a censored poisson
d5 <- d
d5$w <- rep(c(1, 2), length.out = n)
f5 <- cap(frm(y | cens(cc) + weights(w) ~ x, family = poisson(), data = d5))
r <- cap(osa_r(f5))
report("poisson cens + integer weights 1, 2", r, u)

# 6. fractional weights on a Gamma with (1 | g) and on a lognormal
set.seed(8)
g <- factor(rep(1:30, each = 8))
xg <- rnorm(240)
yg <- rgamma(240, shape = 4, rate = 4 / exp(0.5 + 0.3 * xg + rnorm(30, 0, 0.4)[g]))
dg <- data.frame(y = yg, x = xg, g = g, w = rep(c(0.5, 1), 120))
fg <- cap(frm(y | weights(w) ~ x + (1 | g), family = Gamma("log"), data = dg))
r <- cap(osa_r(fg))
cat(sprintf(paste0("Gamma (1|g) weights 0.5/1: finite %d of %d; w = 0.5 rows ",
                   "mean %.3f sd %.3f KS p %.3g; w = 1 rows mean %.3f sd %.3f ",
                   "KS p %.3g\n"), sum(is.finite(r)), length(r),
            mean(r[dg$w == 0.5]), sd(r[dg$w == 0.5]),
            ks.test(r[dg$w == 0.5], "pnorm")$p.value,
            mean(r[dg$w == 1]), sd(r[dg$w == 1]),
            ks.test(r[dg$w == 1], "pnorm")$p.value))
dl <- data.frame(y = exp(1 + 0.5 * x + rnorm(n, 0, 0.7)), x = x,
                 w = runif(n, 0.1, 0.9))
fl <- cap(frm(y | weights(w) ~ x, family = lognormal(), data = dl))
rl <- cap(osa_r(fl))
dpl <- dp1(fl, n)
refl <- qnorm(plnorm(dl$y, dpl$mu, dpl$sigma))
cat(sprintf("lognormal weights in (0.1, 0.9): max rel |r - qnorm(F(y))| %.3g\n",
            max(abs(rl - refl)) / max(abs(refl))))
# a count with a fractional weight keeps NA, with the warning
dc <- data.frame(y = yy, x = x, w = rep(c(0.5, 1), length.out = n))
rc <- cap(osa_r(cap(frm(y | weights(w) ~ x, family = poisson(), data = dc))))
cat("poisson weight 0.5 rows NA:", sum(is.na(rc[dc$w == 0.5])), "of",
    sum(dc$w == 0.5), "\n")
cat("\ndone\n")
