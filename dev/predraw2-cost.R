# Lane predraw2: cost of the draw check and of the reweighting, on lane
# ppdraw's cost fits (dev/ppdraw-fix-cost.R: F1 gaussian with 202
# coefficients and 20000 rows, F2 poisson GLMM with 62 outer parameters
# and 2000 levels, F3 categorical with 49, F4 zinb GLMM with 300 levels)
# and the sweep's ordinal mixture. Per fit: one objective evaluation, the
# side law's build, the check (20 evaluations, draw_logratio()), what a
# fired check costs on 1000 draws (draw_logratio() on 1000 and
# draw_resample()), and the first and second predict() of 100 draws.
# Single runs on a shared machine; each block is a minimum of 3 where it
# is under 1 s. Usage: Rscript dev/predraw2-cost.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
q <- function(expr) suppressMessages(suppressWarnings(expr))
el <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
best <- function(f, k = 3L) {
  t <- el(f())
  if (t < 1) for (i in seq_len(k - 1L)) t <- min(t, el(f()))
  t
}
clear <- function(fit) {
  for (k in c("draw_edge_law", "draw_edge_said", "draw_check",
              "draw_check_said")) {
    if (exists(k, envir = fit$cache, inherits = FALSE)) {
      rm(list = k, envir = fit$cache)
    }
  }
}
cost <- function(label, fit, nd = NULL) {
  ds <- ns$fit_draw_space(fit)
  V <- ds$V
  L <- chol(V)
  map <- ds$map
  v0 <- ns$fit_outer_vector(fit, map)
  P <- length(v0)
  t1fn <- best(function() fit$obj$fn(fit$opt$par))
  clear(fit)
  tlaw <- el(law <- q(ns$draw_edge_law(fit, map, V, colSums(L^2), FALSE,
                                       rownames(V) %||% map$names)))
  dw <- list(L = L, law = law, v0 = v0, held = integer(0),
             off = P - length(fit$opt$par))
  set.seed(1)
  Z0 <- matrix(stats::rnorm(P * 20), P, 20)
  D0 <- v0 + if (is.null(law)) crossprod(L, Z0) else
    ns$draw_edge_apply(law, L, Z0)
  tchk <- best(function() ns$draw_logratio(fit, dw, D0, Z0))
  Z <- matrix(stats::rnorm(P * 1000), P, 1000)
  D <- v0 + if (is.null(law)) crossprod(L, Z) else
    ns$draw_edge_apply(law, L, Z)
  tlw <- el(r <- ns$draw_logratio(fit, dw, D, Z))
  trs <- best(function() ns$draw_resample(r$lw, r$stratum))
  clear(fit)
  set.seed(1)
  tp1 <- el(q(predict(fit, newdata = nd, ndraws = 100L)))
  set.seed(1)
  tp2 <- el(q(predict(fit, newdata = nd, ndraws = 100L)))
  chk <- fit$cache[["draw_check"]]
  cat(sprintf(paste("%-34s np %3d  one fn %.4fs  side law %.2fs  check(20)",
                    "%.2fs  fired: 1000 ratios %.2fs, resample %.3fs |",
                    "first predict %.2fs  second %.2fs  excess %.3f\n"),
              label, length(fit$opt$par), t1fn, tlaw, tchk, tlw, trs, tp1,
              tp2, chk$excess %||% NA_real_))
}

set.seed(1)
n <- 20000
d1 <- data.frame(f = factor(sample(sprintf("L%03d", 1:200), n, TRUE)),
                 x = stats::rnorm(n))
d1$y <- stats::rnorm(n, as.integer(d1$f) / 50 + d1$x)
f1 <- q(frm(bf(y ~ f + x), data = d1))
cost("F1 gaussian y ~ f(200) + x", f1, nd = d1[1:20, ])

set.seed(2)
n <- 20000
d2 <- data.frame(f = factor(sample(sprintf("L%02d", 1:60), n, TRUE)),
                 g = factor(sample(1:2000, n, TRUE)), x = stats::rnorm(n))
u <- stats::rnorm(2000, 0, 0.5)
d2$y <- stats::rpois(n, exp(0.2 + as.integer(d2$f) / 100 + 0.3 * d2$x +
                              u[d2$g]))
f2 <- q(frm(bf(y ~ f + x + (1 | g)), family = poisson(), data = d2))
cost("F2 poisson f(60) + x + (1 | g:2000)", f2, nd = d2[1:20, ])

set.seed(3)
n <- 5000
X <- matrix(stats::rnorm(n * 6), n, 6)
B <- matrix(stats::rnorm(6 * 7, 0, 0.5), 6, 7)
eta <- cbind(0, X %*% B)
pr <- exp(eta) / rowSums(exp(eta))
d3 <- data.frame(X)
d3$y <- factor(apply(pr, 1, function(p) sample.int(8, 1, prob = p)))
f3 <- q(frm(bf(y ~ X1 + X2 + X3 + X4 + X5 + X6), family = categorical(),
            data = d3))
cost("F3 categorical 8 x 7", f3, nd = d3[1:20, ])

set.seed(4)
n <- 6000
d4 <- data.frame(g = factor(sample(1:300, n, TRUE)), x = stats::rnorm(n))
u0 <- stats::rnorm(300, 0, 0.5)
u1 <- stats::rnorm(300, 0, 0.3)
mu <- exp(0.5 + 0.3 * d4$x + u0[d4$g] + u1[d4$g] * d4$x)
d4$y <- ifelse(stats::runif(n) < 0.2, 0L, stats::rnbinom(n, mu = mu, size = 2))
f4 <- q(frm(bf(y ~ x + (1 + x | g), zi ~ x),
            family = zero_inflated_negbinomial(), data = d4))
cost("F4 zinb (1 + x | g:300), zi ~ x", f4, nd = d4[1:20, ])

n <- 300L
source("dev/ppdraw-sweep-gen.R")
d5 <- gen("mixture_cumulative", "interior")
f5 <- q(frm(bf(y ~ x), family = mixture(cumulative(), cumulative()),
            data = d5))
cost("F5 ordinal mixture (fires)", f5)
