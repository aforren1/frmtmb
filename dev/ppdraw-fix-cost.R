# Fix round copy of dev/reviews/ppdraw-review-f.R, timing the side law
# (draw_edge_build()) in place of draw_edge_scan(): cost on fits with
# many outer parameters and on Laplace fits with many group levels, and
# messages per call. Usage: Rscript dev/ppdraw-fix-cost.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
el <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
cost <- function(label, fit, nd = NULL, ndraws = 100L) {
  ds <- frmtmb:::fit_draw_space(fit)
  np <- length(fit$opt$par)
  t1fn <- el(for (i in 1:5) fit$obj$fn(fit$opt$par)) / 5

  ts <- el(h <- frmtmb:::draw_edge_build(fit, ds$map, ds$V, diag(ds$V),
                                         !is.null(ds$root)))
  h <- if (isFALSE(h)) integer(0) else c(h$edge, h$narrowed)
  # a fresh copy of the cache, so the first predict pays the scan
  for (k in c("draw_edge_law", "draw_edge_said")) {
    if (!exists(k, envir = fit$cache, inherits = FALSE)) next
    rm(list = k, envir = fit$cache)
  }
  nm <- 0L
  cnt <- function(m) {
    if (inherits(m, "frmtmb_draw_held")) nm <<- nm + 1L
    invokeRestart("muffleMessage")
  }
  set.seed(1)
  tp1 <- el(withCallingHandlers(
    q(predict(fit, newdata = nd, ndraws = ndraws)), message = cnt))
  set.seed(1)
  tp2 <- el(withCallingHandlers(
    q(predict(fit, newdata = nd, ndraws = ndraws)), message = cnt))
  cat(sprintf(paste("%-44s np %4d  one fn %.4fs  scan %.2fs (%d held)",
                    " first predict %.2fs  second %.2fs  held messages %d\n"),
              label, np, t1fn, ts, length(h), tp1, tp2, nm))
}

# F1: gaussian, 200-level factor, 20000 rows (no random effects)
set.seed(1)
n <- 20000
d1 <- data.frame(f = factor(sample(sprintf("L%03d", 1:200), n, TRUE)),
                 x = stats::rnorm(n))
d1$y <- stats::rnorm(n, as.integer(d1$f) / 50 + d1$x)
t <- el(f1 <- q(frm(bf(y ~ f + x), data = d1)))
cat(sprintf("F1 fit %.1fs\n", t))
cost("F1 gaussian y ~ f(200) + x, n 20000", f1, nd = d1[1:20, ])

# F2: poisson GLMM, 60-level factor fixed, (1 | g) with 2000 levels
set.seed(2)
n <- 20000
d2 <- data.frame(f = factor(sample(sprintf("L%02d", 1:60), n, TRUE)),
                 g = factor(sample(1:2000, n, TRUE)), x = stats::rnorm(n))
u <- stats::rnorm(2000, 0, 0.5)
d2$y <- stats::rpois(n, exp(0.2 + as.integer(d2$f) / 100 + 0.3 * d2$x +
                              u[d2$g]))
t <- el(f2 <- q(frm(bf(y ~ f + x + (1 | g)), family = poisson(), data = d2)))
cat(sprintf("F2 fit %.1fs\n", t))
cost("F2 poisson y ~ f(60) + x + (1 | g:2000)", f2, nd = d2[1:20, ])

# F3: categorical, 8 categories, 6 predictors
set.seed(3)
n <- 5000
X <- matrix(stats::rnorm(n * 6), n, 6)
B <- matrix(stats::rnorm(6 * 7, 0, 0.5), 6, 7)
eta <- cbind(0, X %*% B)
pr <- exp(eta) / rowSums(exp(eta))
d3 <- data.frame(X)
d3$y <- factor(apply(pr, 1, function(p) sample.int(8, 1, prob = p)))
t <- el(f3 <- q(frm(bf(y ~ X1 + X2 + X3 + X4 + X5 + X6),
                    family = categorical(), data = d3)))
cat(sprintf("F3 fit %.1fs\n", t))
cost("F3 categorical 8 cats x 7 coefs", f3, nd = d3[1:20, ])

# F4: zinb GLMM with (1 + x | g), 300 levels (a slow Laplace fit)
set.seed(4)
n <- 6000
d4 <- data.frame(g = factor(sample(1:300, n, TRUE)), x = stats::rnorm(n))
u0 <- stats::rnorm(300, 0, 0.5)
u1 <- stats::rnorm(300, 0, 0.3)
mu <- exp(0.5 + 0.3 * d4$x + u0[d4$g] + u1[d4$g] * d4$x)
d4$y <- ifelse(stats::runif(n) < 0.2, 0L, stats::rnbinom(n, mu = mu, size = 2))
t <- el(f4 <- q(frm(bf(y ~ x + (1 + x | g), zi ~ x), family =
                      zero_inflated_negbinomial(), data = d4)))
cat(sprintf("F4 fit %.1fs\n", t))
cost("F4 zinb (1 + x | g:300), zi ~ x", f4, nd = d4[1:20, ])
