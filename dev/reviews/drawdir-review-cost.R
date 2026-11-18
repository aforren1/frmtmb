# Reviewer of lane drawdir, check 5: the reading's cost on two of the
# lane's cost fits (dev/drawdir-cost.R F2 and F7, same data), and the
# cache across predict() calls that differ in newdata, re_formula,
# ndraws and propagate_error.
# Usage: Rscript dev/reviews/drawdir-review-cost.R
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
q <- function(expr) suppressMessages(suppressWarnings(expr))
el <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
real <- ns$draw_reweight
cap <- new.env()
law_of <- function(fit) {
  for (k in c("draw_edge_law", "draw_edge_said", "draw_check",
              "draw_check_said")) {
    if (exists(k, envir = fit$cache, inherits = FALSE)) {
      rm(list = k, envir = fit$cache)
    }
  }
  assignInNamespace("draw_reweight", function(object, dw) {
    dw$off <- length(dw$v0) - length(object$opt$par)
    cap$dw <- dw
    dw$D
  }, ns = "frmtmb")
  t <- el(q(ns$predict_par_drawer(fit, TRUE, 2L)))
  assignInNamespace("draw_reweight", real, ns = "frmtmb")
  list(dw = cap$dw, t = t)
}
# evaluations counted through the objective
count_dirs <- function(fit, dw) {
  fn0 <- fit$obj$fn
  cnt <- 0L
  fit$obj$fn <- function(...) {
    cnt <<- cnt + 1L
    fn0(...)
  }
  t <- el(dd <- ns$draw_check_dirs(fit, dw))
  fit$obj$fn <- fn0
  list(t = t, n = cnt, dd = dd)
}

set.seed(7)
d7 <- data.frame(f = factor(rep(1:1000, each = 5)))
d7$y <- stats::rnorm(5000, stats::rnorm(1000)[d7$f])
f7 <- q(frm(bf(y ~ 0 + f), data = d7))
l7 <- law_of(f7)
for (i in 1:2) {
  gc(reset = TRUE)
  m0 <- sum(gc()[, 2])
  c7 <- count_dirs(f7, l7$dw)
  m1 <- sum(gc()[, 6])
  cat(sprintf("F7 p %d run %d: dirs %.2f s, %d objective calls, stat %.4g | side law+drawer %.2f s | max Mb used %.0f\n",
              length(f7$opt$par), i, c7$t, c7$n, c7$dd$stat, l7$t, m1))
}
Rprof(tf <- tempfile(), interval = 0.01)
invisible(ns$draw_check_dirs(f7, l7$dw))
Rprof(NULL)
sp <- summaryRprof(tf)$by.total
print(utils::head(sp[, c("total.time", "self.time")], 15))

set.seed(2)
n <- 20000
d2 <- data.frame(f = factor(sample(sprintf("L%02d", 1:60), n, TRUE)),
                 g = factor(sample(1:2000, n, TRUE)), x = stats::rnorm(n))
u <- stats::rnorm(2000, 0, 0.5)
d2$y <- stats::rpois(n, exp(0.2 + as.integer(d2$f) / 100 + 0.3 * d2$x +
                              u[d2$g]))
f2 <- q(frm(bf(y ~ f + x + (1 | g)), family = poisson(), data = d2))
l2 <- law_of(f2)
t1 <- el(f2$obj$fn(f2$opt$par + 1e-6))
c2 <- count_dirs(f2, l2$dw)
cat(sprintf("F2 np %d: dirs %.2f s, %d objective calls (%.3f s each), stat %.4g | one moved fn %.3f s | side law+drawer %.2f s\n",
            length(f2$opt$par), c2$t, c2$n, c2$t / c2$n, c2$dd$stat, t1,
            l2$t))

# the cache: one verdict per fit, whatever the call
set.seed(3)
d <- data.frame(g = factor(rep(1:15, each = 8)), x = stats::rnorm(120))
d$y <- 1 + 0.5 * d$x + stats::rnorm(15, 0, 0.8)[d$g] +
  stats::rnorm(15, 0, 0.1)[d$g] * d$x + stats::rnorm(120)
fit <- q(frm(bf(y ~ x + (1 + x | g)), data = d))
calls <- list(
  pe_false = function() predict(fit, propagate_error = FALSE, ndraws = 10),
  plain = function() predict(fit, ndraws = 50, summary = FALSE),
  newdata = function() predict(fit, newdata = d[1:5, ], ndraws = 500,
                               summary = FALSE),
  newlevel = function() predict(fit, newdata = data.frame(
    x = 2, g = factor("new")), allow_new_levels = TRUE, ndraws = 20,
    summary = FALSE),
  re_na = function() predict(fit, re_formula = NA, ndraws = 100,
                             summary = FALSE))
fn0 <- fit$obj$fn
cnt <- 0L
fit$obj$fn <- function(...) {
  cnt <<- cnt + 1L
  fn0(...)
}
for (nm in names(calls)) {
  before <- cnt
  set.seed(1)
  q(calls[[nm]]())
  chk <- fit$cache[["draw_check"]]
  cat(sprintf("cache %-9s objective calls %4d | verdict %s fire %s dir %s n %s\n",
              nm, cnt - before, !is.null(chk), chk$fire %||% NA,
              format(chk$dir %||% NA, digits = 6), chk$dir_n %||% NA))
}
fit$obj$fn <- fn0
# a fresh fit of the same data gives the same verdict
fit2 <- q(frm(bf(y ~ x + (1 + x | g)), data = d))
set.seed(1)
q(predict(fit2, re_formula = NA, ndraws = 7, summary = FALSE))
c1 <- fit$cache[["draw_check"]]
c2 <- fit2$cache[["draw_check"]]
cat("fresh fit verdict identical:", identical(c1, c2), "\n")
cat("DONE\n")
