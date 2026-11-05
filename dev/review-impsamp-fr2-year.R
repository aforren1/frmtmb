# Fix round 2, item 1. `year` seed 2 read SE ratio 0.867 to 18.27 against
# accuracy2's reference; null2 seed 3 read 7246. Which parameter, and which is
# wrong? Reference: my exact marginal on the CENTERED column, its Hessian mapped to
# the raw column by the delta rule, and the fit on the raw column.
# Usage: REV_ARM=after|before Rscript dev/review-impsamp-fr2-year.R
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs2.R")
ses <- function(f) tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)), error = function(e) NA)
cat("arm", arm, find.package("frmtmb"), "\n")
gh <- gh_rule(61L)
for (s in 1:3) {
  de <- design_list$year(s)
  d <- de$data
  m <- mean(d$x)
  dl <- split(d, d$g)
  nll_c <- function(p) {
    sd <- exp(p[3])
    tot <- 0
    for (dg in dl) {
      fn <- function(u) {
        vapply(u, function(ui) sum(dbinom(dg$y, 1, plogis(p[1] + p[2] * (dg$x - m) + ui), log = TRUE)), 0) +
          dnorm(u, 0, sd, log = TRUE)
      }
      o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
      q <- o$minimum
      cur <- -(fn(q + 1e-3) - 2 * fn(q) + fn(q - 1e-3)) / 1e-6
      w <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
      nd <- q + sqrt(2) * w * gh$x
      a <- fn(nd) + gh$x^2 + log(gh$w)
      tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * w))
    }
    tot
  }
  ng <- function(f, p, h = 1e-4) {
    vapply(seq_along(p), function(i) {
      e <- replace(numeric(length(p)), i, h)
      (f(p + e) - f(p - e)) / (2 * h)
    }, 0)
  }
  o <- nlminb(c(0, 0.1, 0), nll_c, gradient = function(p) ng(nll_c, p), control = list(rel.tol = 1e-10))
  H <- numDeriv::jacobian(function(p) ng(nll_c, p), o$par, method = "simple")
  V <- solve(0.5 * (H + t(H)))
  J <- rbind(c(1, -m, 0), c(0, 1, 0), c(0, 0, 1))
  Vr <- J %*% V %*% t(J)
  se_r <- sqrt(diag(Vr))
  pr <- c(o$par[1] - m * o$par[2], o$par[2], o$par[3])
  cat(sprintf("\nyear seed %d: exact raw par %s\n  exact SE (delta from centered) %s\n", s,
              paste(format(pr, digits = 6), collapse = " "), paste(format(se_r, digits = 5), collapse = " ")))
  w <- character(0)
  fi <- withCallingHandlers(frm(de$formula, family = de$family, data = d, importance = de$draws),
                            warning = function(x) {
                              w <<- c(w, substr(conditionMessage(x), 1, 80))
                              invokeRestart("muffleWarning")
                            })
  fl <- suppressWarnings(frm(de$formula, family = de$family, data = d))
  cat("  importance SE/exact", format(ses(fi) / se_r, digits = 5), " z", format((unname(fi$opt$par) - pr) / se_r, digits = 3),
      " warns", length(w), "\n")
  cat("  Laplace    SE/exact", format(ses(fl) / se_r, digits = 5), " z", format((unname(fl$opt$par) - pr) / se_r, digits = 3), "\n")
  # accuracy2's reference, reproduced: num_hess with h = 1e-3 on the raw coefficient
  nll_r <- function(p) nll_c(c(p[1] + m * p[2], p[2], p[3]))
  num_hess <- function(f, x, h = 1e-3) {
    k <- length(x)
    H <- matrix(0, k, k)
    for (i in seq_len(k)) {
      ei <- replace(numeric(k), i, h)
      H[i, ] <- ng(f, x + ei, 1e-3) - ng(f, x - ei, 1e-3)
      H[i, ] <- H[i, ] / (2 * h)
    }
    0.5 * (H + t(H))
  }
  V2 <- tryCatch(solve(num_hess(nll_r, pr)), error = function(e) NULL)
  if (!is.null(V2)) {
    cat("  accuracy2's reference SE", format(sqrt(abs(diag(V2))), digits = 5), " fit / it", format(ses(fi) / sqrt(abs(diag(V2))), digits = 5), "\n")
  }
  cat("  summary est error", format(summary(fi)$fixed[, 2], digits = 5), "\n")
}
# null2 seed 3: the correlation's chart
set.seed(3)
g <- rep(1:30, each = 2)
x <- rnorm(60)
y <- rbinom(60, 1, 0.4)
d <- data.frame(y, x, g = factor(g))
w <- character(0)
fi <- withCallingHandlers(frm(bf(y ~ x + (x | g)), family = bernoulli(), data = d, importance = 500L),
                          warning = function(x) {
                            w <<- c(w, substr(conditionMessage(x), 1, 70))
                            invokeRestart("muffleWarning")
                          })
fl <- suppressWarnings(frm(bf(y ~ x + (x | g)), family = bernoulli(), data = d))
cat("\nnull2 seed 3: par", format(unname(fi$opt$par), digits = 5), " SE", format(ses(fi), digits = 5), " warns", paste(w, collapse = ";"), "\n")
cat("Laplace par", format(unname(fl$opt$par), digits = 5), " SE", format(ses(fl), digits = 5), "\n")
vc <- frmtmb:::varcorr_matrices(fi)[[1L]]
cat("cor", vc[1, 2] / sqrt(vc[1, 1] * vc[2, 2]), "\n")
for (ch in c(-1000, -668, -300, -100, -30, -10)) {
  p <- replace(unname(fi$opt$par), 5L, ch)
  o <- nlminb(p[1:4], function(b) fi$obj$fn(c(b, ch)), function(b) fi$obj$gr(c(b, ch))[1:4])
  cat(sprintf("  chart %-6g tape profile %.5f\n", ch, o$objective))
}
print(tryCatch(summary(fi)$random, error = function(e) "no random table"))
