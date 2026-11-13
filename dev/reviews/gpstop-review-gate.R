# Review of lane gpstop: the gate of gp_body_relin(). Each case runs arms
# lane and norelin in one process and prints how often the relinearization
# ran, whether a fit keeps the closure, the fit's serialized size, and the
# log-likelihood. Data: dev/review-fitedge-designs.R rv_mk() (E wig 120 2 and
# L wig 120 2, the lane's rows) and set.seed(81900 + case index).
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
source(file.path(RV_WT, "dev/review-fitedge-designs.R"))
which_cases <- if (length(a)) strsplit(a[1], ",")[[1]] else "all"
dE <- rv_mk("E", 2, 120L, 2)
dL <- rv_mk("L", 2, 120L, 2)
sz <- function(f) {
  if (inherits(f, "error")) return(NA)
  length(serialize(f, NULL))
}
both <- function(tag, ...) {
  out <- list()
  for (arm in c("lane", "norelin")) {
    q <- rv_run(tag, arm, ...)
    out[[arm]] <- q
    cat(sprintf("SIZE %s %s %s\n", tag, arm, format(sz(q$fit))))
  }
  invisible(out)
}
fE <- bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE)
cases <- list(
  approx_a = function() both("approx_a k20", bf(y ~ a * exp(b * z),
                                                a ~ gp(x, k = 20), b ~ 1,
                                                nl = TRUE), data = dE),
  approx_b = function() both("approx_b k20", bf(y ~ a * exp(b * z), a ~ 1,
                                                b ~ gp(x, k = 20), nl = TRUE),
                             data = dL),
  by = function() {
    d <- dE
    d$y <- ifelse(d$fby == "p", d$y, 0.6 * d$y)
    both("by", bf(y ~ a * exp(b * z), a ~ gp(x, by = fby), b ~ 1, nl = TRUE),
         data = d)
  },
  sigma = function() {
    set.seed(81903)
    d <- dE
    d$y <- d$y + rnorm(nrow(d), 0, 0.05 * exp(0.3 * sin(d$x)))
    both("sigma_gp", bf(y ~ a * exp(b * z), a ~ 1, b ~ 1, sigma ~ gp(x),
                        nl = TRUE), data = d)
  },
  two = function() {
    set.seed(81904)
    d <- dE
    d$y <- (2 + 0.5 * sin(d$x)) * exp((0.3 + 0.1 * cos(d$x)) * d$z) +
      0.1 * rnorm(nrow(d))
    both("two_gp", bf(y ~ a * exp(b * z), a ~ gp(x), b ~ gp(x), nl = TRUE),
         data = d)
  },
  group = function() both("gp_group", bf(y ~ a * exp(b * z),
                                         a ~ gp(x) + (1 | g), b ~ 1,
                                         nl = TRUE), data = dE),
  start = function() both("start_beta", fE, data = dE,
                          start = list(beta = c(2, 0.3))),
  autoscale = function() {
    d <- dE
    d$xs <- d$x * 1e4
    d$wv <- rnorm(nrow(d)) * 1e5
    both("autoscale", bf(y ~ a * exp(b * z), a ~ gp(xs) + wv, b ~ 1,
                         nl = TRUE), data = d,
         control = frmtmb_control(autoscale = TRUE))
  },
  refit = function() {
    for (arm in c("lane", "norelin")) {
      rv_arm(arm)
      f <- suppressWarnings(frm(fE, data = dE))
      CNT$relin <- 0L
      set.seed(81909)
      ynew <- dE$y + rnorm(nrow(dE), 0, 0.05)
      t0 <- proc.time()[[3]]
      r <- suppressWarnings(refit(f, newresp = ynew))
      t1 <- proc.time()[[3]] - t0
      cat(sprintf("REFIT %s ll %.6f relin %d wall %.2f runs %s\n", arm,
                  as.numeric(logLik(r)), CNT$relin, t1,
                  format(r$opt$gp_second[["runs"]] %||% NA)))
      CNT$relin <- 0L
      t0 <- proc.time()[[3]]
      u <- suppressWarnings(update(f, data = transform(dE, y = ynew)))
      t1 <- proc.time()[[3]] - t0
      cat(sprintf("UPDATE %s ll %.6f relin %d wall %.2f\n", arm,
                  as.numeric(logLik(u)), CNT$relin, t1))
      flush.console()
    }
  })
if (!identical(which_cases, "all")) cases <- cases[which_cases]
for (nm in names(cases)) {
  r <- tryCatch(cases[[nm]](), error = function(e) {
    cat("CASE", nm, "ERROR", conditionMessage(e), "\n")
  })
}
# approximate gp(): is it short where the exact one was? best of theta
# starts as the reference
for (spec in list(list("approx_a k20", bf(y ~ a * exp(b * z), a ~ gp(x, k = 20),
                                          b ~ 1, nl = TRUE), dE),
                  list("approx_b k20", bf(y ~ a * exp(b * z), a ~ 1,
                                          b ~ gp(x, k = 20), nl = TRUE), dL))) {
  if (!identical(which_cases, "all") &&
      !any(c("approx_a", "approx_b") %in% which_cases)) break
  rv_arm("lane")
  f0 <- rv_fit(spec[[2]], data = spec[[3]])$fit
  th <- f0$opt$par[names(f0$opt$par) == "theta"]
  best <- -Inf
  for (rho in c(0.1, 0.3, 1, 3)) for (sd0 in c(0.3, 1.5)) {
    th[] <- c(log(sd0), log(rho))
    r <- rv_fit(spec[[2]], data = spec[[3]], start = list(theta = unname(th)))
    if (!inherits(r$fit, "error")) best <- max(best, as.numeric(logLik(r$fit)))
  }
  cat(sprintf("APPROXREF %s default %.6f best of 8 theta starts %.6f gap %.4f\n",
              spec[[1]], as.numeric(logLik(f0)), best,
              as.numeric(logLik(f0)) - best))
}
cat("DONE\n")
