# Reviewer of lane loglik: accuracy near the parameter bounds against
# closed forms (AR(1) and compound symmetry have exact precision
# matrices), and the edge cases of addition terms. Seed 91.
# Output: dev/review-loglik-log/numerics.txt
LIB <- Sys.getenv("RV_LIB", "C:/Users/adf44/source/r/wt-loglik-lib")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
unfit <- function(f, ...) {
  fit <- frm(f, ..., dry_run = "objective")
  class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
  fit
}

set.seed(91)
G <- 30
k <- 8
d <- data.frame(g = rep(seq_len(G), each = k), t = rep(seq_len(k), G))
d$y <- stats::rnorm(nrow(d))
d <- d[sample(nrow(d)), ]
rownames(d) <- NULL

# AR(1), marginal sd s: interior row has mean phi (e_{t-1} + e_{t+1}) /
# (1 + phi^2) and variance s^2 (1 - phi^2) / (1 + phi^2); an end row
# mean phi e_nb and variance s^2 (1 - phi^2)
ar1_closed <- function(d, mu, s, phi) {
  out <- numeric(nrow(d))
  for (g in unique(d$g)) {
    r <- which(d$g == g)
    r <- r[order(d$t[r])]
    e <- d$y[r] - mu
    kk <- length(r)
    for (j in seq_len(kk)) {
      if (j == 1 || j == kk) {
        nb <- if (j == 1) e[2] else e[kk - 1]
        m <- phi * nb
        v <- s^2 * (1 - phi^2)
      } else {
        m <- phi * (e[j - 1] + e[j + 1]) / (1 + phi^2)
        v <- s^2 * (1 - phi^2) / (1 + phi^2)
      }
      out[r[j]] <- stats::dnorm(e[j], m, sqrt(v), log = TRUE)
    }
  }
  out
}
cat("ar(cov = TRUE), 30 groups x 8, against the AR(1) closed form\n")
for (s in c(1, 1e-6)) {
  for (phi in c(0.9, 0.999, 0.99999, 1 - 1e-8, -0.99999)) {
    fit <- unfit(y ~ 1 + ar(t, g, cov = TRUE), data = d)
    fit$estimates$beta[] <- 0
    fit$estimates$betad[] <- log(s)
    fit$estimates$thetaac[] <- phi / sqrt(1 - phi^2)
    ll <- tryCatch(autocor_row_loglik(fit, eval_dpars(fit), "y"),
                   error = function(e) conditionMessage(e))
    ref <- ar1_closed(d, 0, s, phi)
    if (is.character(ll)) {
      cat(sprintf("  sigma %-6g phi %-12.10g ERROR %s\n", s, phi, ll))
    } else {
      cat(sprintf("  sigma %-6g phi %-12.10g rel = %.3e  (max |ref| %.3g)\n",
                  s, phi, rel(ll, ref), max(abs(ref))))
    }
  }
}

# compound symmetry: Q = (I - c J) / (s^2 (1 - rho)), c = rho / (1 + (k -
# 1) rho)
cs_closed <- function(d, mu, s, rho) {
  out <- numeric(nrow(d))
  for (g in unique(d$g)) {
    r <- which(d$g == g)
    e <- d$y[r] - mu
    kk <- length(r)
    cc <- rho / (1 + (kk - 1) * rho)
    f <- 1 / (s^2 * (1 - rho))
    qd <- f * (1 - cc)
    ge <- f * (e - cc * sum(e))
    out[r] <- 0.5 * log(qd / (2 * pi)) - 0.5 * ge^2 / qd
  }
  out
}
cat("cosy(), 30 groups x 8, against the closed form\n")
fitc <- unfit(y ~ 1 + cosy(t, g), data = d)
dd <- fitc$frame$autocor$y$d
a <- 1 / (dd - 1)
for (rho in c(0.5, 0.999999, -a + 1e-3, -a + 1e-6, -a + 1e-9)) {
  fitc$estimates$beta[] <- 0
  fitc$estimates$betad[] <- 0
  fitc$estimates$thetaac[] <- qlogis((rho + a) / (1 + a))
  ll <- tryCatch(autocor_row_loglik(fitc, eval_dpars(fitc), "y"),
                 error = function(e) conditionMessage(e))
  ref <- cs_closed(d, 0, 1, rho)
  if (is.character(ll)) {
    cat(sprintf("  rho %-14.10g ERROR %s\n", rho, ll))
  } else {
    cat(sprintf("  rho %-14.10g rel = %.3e\n", rho, rel(ll, ref)))
  }
}

# the frame's refusals of addition terms beside a covariance term
cat("\naddition terms beside ar(cov = TRUE) in frm()\n")
d$w <- runif(nrow(d), 0.5, 2)
d$se <- runif(nrow(d), 0.1, 0.3)
d$cc <- 0L
d$lb <- -10
d$ym <- d$y
d$ym[3] <- NA
for (f in list(y | weights(w) ~ ar(t, g, cov = TRUE),
               y | se(se, sigma = TRUE) ~ ar(t, g, cov = TRUE),
               y | cens(cc) ~ ar(t, g, cov = TRUE),
               y | trunc(lb = lb) ~ ar(t, g, cov = TRUE),
               ym | mi() ~ ar(t, g, cov = TRUE),
               y | weights(w) ~ cosy(t, g),
               y | se(se, sigma = TRUE) ~ unstr(t, g))) {
  r <- tryCatch({
    frm(f, data = d, dry_run = "objective")
    "ACCEPTED"
  }, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 90)))
  cat(sprintf("  %-46s %s\n", deparse1(f), r))
}
cat("DONE\n")
