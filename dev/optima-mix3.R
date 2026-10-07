# Lane optima, item 2: the three-component ordinal mixture of
# test-nonfinite-gradient.R (data of dev/ordmix-rev-probit-sat.R),
# mixture(cumulative("probit"), sratio("cloglog"), acat()). On 0.68.1
# it ended where the probit's log-odds had underflowed, with an infinite
# gradient. Here: the fit on the arm given, then the optimizer run by
# hand on the taped objective, recording where a gradient first stops
# being finite at a finite objective, and which parameters carry it.
#   Rscript dev/optima-mix3.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
set.seed(20261005 + 31)
n <- 400
x <- rnorm(n)
rnorm(n)
sample(c("a", "b"), n, TRUE)
rbinom(n, 1, 0.4)
rlogis(n)
runif(n)
runif(n, 0.5, 2)
cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
d <- data.frame(y = as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf))), x = x)
famsel <- if (length(commandArgs(TRUE)) > 1) commandArgs(TRUE)[2] else "review"
fam <- function() switch(famsel,
  review = mixture(cumulative("probit"), sratio("cloglog"), acat()),
  cum_acat = mixture(cumulative("probit"), acat()))
r <- tryCatch(withCallingHandlers(
  frm(bf(y ~ x), family = fam(), data = d, verbose = TRUE),
  warning = function(w) {
    cat("WARNING:", substr(conditionMessage(w), 1, 160), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) conditionMessage(e))
if (is.character(r)) cat("ERROR:", substr(r, 1, 120), "\n") else {
  cat("logLik", format(as.numeric(logLik(r)), digits = 10), "code",
      r$opt$convergence, "\n")
}
u <- frm(bf(y ~ x), family = fam(), data = d, dry_run = "objective")
obj <- u$obj
tr <- list()
fn <- function(p) {
  v <- obj$fn(p)
  tr[[length(tr) + 1L]] <<- list(p = p, f = v, k = "fn")
  if (is.nan(v)) Inf else v
}
gr <- function(p) {
  g <- obj$gr(p)
  tr[[length(tr) + 1L]] <<- list(p = p, f = obj$fn(p), g = g, k = "gr")
  g
}
o <- tryCatch(stats::nlminb(obj$par, fn, gr,
                            control = list(eval.max = 1000,
                                           iter.max = 1000)),
              error = function(e) conditionMessage(e))
print(if (is.character(o)) o else o[c("objective", "convergence")])
bad <- which(vapply(tr, function(t) t$k == "gr" && !all(is.finite(t$g)), NA))
if (length(bad)) {
  t <- tr[[bad[1]]]
  cat("first non-finite gradient at evaluation", bad[1], "of", length(tr),
      ": objective", format(t$f, digits = 12), "\n")
  print(stats::setNames(as.numeric(t$g), names(obj$par))[!is.finite(t$g)])
  pl <- obj$env$parList(t$p)
  print(pl[setdiff(names(pl), "beta")])
  print(pl$beta)
}
if (length(bad)) {
  # each component's log-density per row at that point, numerically
  u2 <- u
  u2$estimates <- obj$env$parList(tr[[bad[1]]]$p)
  class(u2) <- setdiff(class(u2), "frmtmb_unfitted")
  resp <- names(u2$spec$responses)[1]
  rspec <- u2$spec$responses[[resp]]
  famo <- rspec$family
  dp <- frmtmb:::eval_dpars(u2)[[resp]]
  av <- u2$frame[["aterm_values"]][[resp]]
  ex <- frmtmb:::fit_extras(u2, resp)
  y <- u2$frame[["y"]][[resp]]
  L <- sapply(1:3, function(k) {
    as.numeric(famo[["mix"]][["comp_lpdf"]](y, dp, av, k, ex))
  })
  lp <- famo[["mix"]][["log_pi"]](dp)
  cat("per component: rows with -Inf", colSums(L == -Inf), "; NaN",
      colSums(is.nan(L)), "; min finite", apply(L, 2, function(v)
        min(v[is.finite(v)])), "\n")
  both <- which(rowSums(L == -Inf) >= 2)
  cat("rows with two or more components at -Inf:", length(both), "\n")
  if (length(both)) print(cbind(y = y[both], L[both, , drop = FALSE])[1:min(5, length(both)), ])
}
if (length(bad)) {
  p0 <- tr[[bad[1]]]$p
  nb <- which(!is.finite(tr[[bad[1]]]$g))
  for (i in nb) {
    h <- 1e-6 * max(1, abs(p0[i]))
    e <- replace(numeric(length(p0)), i, h)
    cat(sprintf("param %d (%s) value %.6g: central difference %.6g\n", i,
                names(obj$par)[i], p0[i],
                (obj$fn(p0 + e) - obj$fn(p0 - e)) / (2 * h)))
  }
  # which rows: the gradient of each row's component-2 density
  # through the tape, by dropping rows from the data
  cat("eta of component 2 at the rows, range:",
      range(dp[["mu2"]]), "; tau2:", format(frmtmb:::ord_tau_from_raw(
        u2$estimates$tau_raw2, TRUE), digits = 17), "
")
}
if (length(bad)) {
  nm <- names(obj$par)
  q <- p0
  q[nm == "tau_raw1"] <- c(0, 0, 0)
  q[which(nm == "beta")[1]] <- 0
  cat("component 1 tamed: gradient finite", all(is.finite(obj$gr(q))), "\n")
  q <- p0
  q[which(nm == "tau_raw2")[3]] <- -30
  cat("tau_raw2[3] at -30: gradient finite", all(is.finite(obj$gr(q))), "\n")
  for (r3 in c(-700, -705, -708, -709, -710)) {
    q <- p0
    q[which(nm == "tau_raw2")[3]] <- r3
    cat("tau_raw2[3] at", r3, ": gradient finite", all(is.finite(obj$gr(q))),
        "\n")
  }
}
if (length(bad)) {
  # component 2 alone, taped as a function of (mu2_x, tau_raw2)
  ex0 <- ex
  f2 <- function(v) {
    dpv <- dp
    dpv[["mu2"]] <- v[1] * d$x
    exv <- ex0
    exv$tau_raw2 <- v[2:4]
    famo[["mix"]][["comp_lpdf"]](y, dpv, av, 2L, exv)
  }
  v0 <- c(u2$estimates$beta[2], u2$estimates$tau_raw2)
  tp <- RTMB::MakeTape(function(v) sum(f2(v)), v0)
  cat("component 2 alone: value", format(tp(v0), digits = 10),
      "; gradient", format(tp$jacobian(v0)), "\n")
  rows_bad <- vapply(seq_along(y), function(i) {
    ti <- RTMB::MakeTape(function(v) f2(v)[i], v0)
    any(!is.finite(ti$jacobian(v0)))
  }, NA)
  cat("rows with a non-finite component-2 gradient:", sum(rows_bad),
      "; their y:", paste(table(y[rows_bad]), collapse = " "), "\n")
}
if (length(bad)) {
  ib <- which(rows_bad)
  tau2 <- frmtmb:::ord_tau_from_raw(u2$estimates$tau_raw2, TRUE)
  for (i in ib) {
    cat("row", i, "y", y[i], "eta", format(dp[["mu2"]][i], digits = 17),
        "M", format(tau2 - dp[["mu2"]][i], digits = 17), "\n")
  }
  lk <- frmtmb:::get_link("cloglog")
  lg <- frmtmb:::ord_log_cdf_pair(lk)
  for (i in ib) for (j in 1:3) {
    m <- tau2[j] - dp[["mu2"]][i]
    a <- RTMB::MakeTape(function(x) lg$lF(x), m)
    b <- RTMB::MakeTape(function(x) lg$l1mF(x), m)
    cat(sprintf("  M_%d = %.17g: lF %.6g (d %.6g), l1mF %.6g (d %.6g)\n", j, m,
                a(m), a$jacobian(m), b(m), b$jacobian(m)))
  }
}
if (length(bad)) {
  for (i in ib) {
    ti <- RTMB::MakeTape(function(v) f2(v)[i], v0)
    cat("row", i, "value", format(ti(v0), digits = 10), "jacobian",
        format(ti$jacobian(v0)), "\n")
  }
  # the same row through the sratio pieces by hand
  i <- ib[1]
  xi <- d$x[i]
  h <- function(v) {
    tau <- frmtmb:::ord_tau_from_raw_ad(v[2:4])
    M <- tau - v[1] * xi
    lg$l1mF(M[1]) + lg$l1mF(M[2]) + lg$lF(M[3])
  }
  th <- RTMB::MakeTape(h, v0)
  cat("by hand: value", format(th(v0), digits = 10), "jacobian",
      format(th$jacobian(v0)), "\n")
}
