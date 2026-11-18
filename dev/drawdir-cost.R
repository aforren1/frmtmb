# Lane drawdir: cost of the per-direction reading (draw_check_dirs()) on
# predraw2's cost fits (dev/predraw2-cost.R: F1 gaussian with 202
# coefficients and 20000 rows, F2 poisson GLMM with 62 outer parameters
# and 2000 levels, F3 categorical with 49, F4 zinb GLMM with 300 levels,
# F5 the sweep's ordinal mixture) and two large-p GLMs: F6 poisson
# y ~ 0 + f with 500 levels of 5 rows, F7 gaussian y ~ 0 + f with 1000
# levels of 5 rows. Per fit: one objective evaluation, the summed check
# (20 evaluations, draw_logratio()), the per-direction reading and its
# evaluation count, and the first and second predict() of 100 draws.
# Single runs on a shared machine; a block under 1 s is the least of 3.
# Usage: Rscript dev/drawdir-cost.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
has_dirs <- exists("draw_check_dirs", envir = ns, inherits = FALSE)
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
real <- ns$draw_reweight
cap <- new.env()
grab <- function(object, dw) {
  dw$off <- length(dw$v0) - length(object$opt$par)
  cap$dw <- dw
  dw$D
}
cost <- function(label, fit, nd = NULL) {
  t1fn <- best(function() fit$obj$fn(fit$opt$par))
  clear(fit)
  assignInNamespace("draw_reweight", grab, ns = "frmtmb")
  set.seed(1)
  tlaw <- el(q(ns$predict_par_drawer(fit, TRUE, 2L)))
  assignInNamespace("draw_reweight", real, ns = "frmtmb")
  dw <- cap$dw
  P <- length(dw$v0)
  Z0 <- ns$draw_check_z(P, 20L)
  D0 <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z0) else
    ns$draw_edge_apply(dw$law, dw$L, Z0)
  tchk <- best(function() ns$draw_logratio(fit, dw, D0, Z0))
  tdir <- NA_real_
  dd <- NULL
  if (has_dirs) tdir <- best(function() dd <<- ns$draw_check_dirs(fit, dw))
  clear(fit)
  set.seed(1)
  tp1 <- el(q(predict(fit, newdata = nd, ndraws = 100L)))
  set.seed(1)
  tp2 <- el(q(predict(fit, newdata = nd, ndraws = 100L)))
  chk <- fit$cache[["draw_check"]]
  cat(sprintf(paste("%-36s np %4d  one fn %.4fs  side law+drawer %.2fs",
                    "check(20) %.2fs  dirs %.2fs (%s evals, stat %s) |",
                    "first predict %.2fs  second %.2fs  fire %s\n"),
              label, length(fit$opt$par), t1fn, tlaw, tchk, tdir,
              if (is.null(dd)) "-" else dd$n,
              if (is.null(dd)) "-" else format(dd$stat, digits = 3), tp1,
              tp2, isTRUE(chk$fire)))
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

set.seed(6)
d6 <- data.frame(f = factor(rep(1:500, each = 5)))
d6$y <- stats::rpois(2500, exp(log(2) + stats::rnorm(500, 0, 0.3))[d6$f])
f6 <- q(frm(bf(y ~ 0 + f), family = poisson(), data = d6))
cost("F6 poisson y ~ 0 + f(500)", f6)

set.seed(7)
d7 <- data.frame(f = factor(rep(1:1000, each = 5)))
d7$y <- stats::rnorm(5000, stats::rnorm(1000)[d7$f])
f7 <- q(frm(bf(y ~ 0 + f), data = d7))
cost("F7 gaussian y ~ 0 + f(1000)", f7)
