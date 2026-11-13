# Lane gpstop, item 2 prototype: a student() fit's working marginal
# reweighted at its own optimum. At the default's stop, each row gets the
# EM weight of the t density, w = (nu + 1) / (nu + (r / sigma)^2), and the
# gaussian marginal of y with known noise sigma^2 / w gives a theta. Then a
# descent from the default's point with that theta. Args: arm, rows as
# case:curve:n:seed separated by commas (dev/arstep-ng.R data).
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs(a[1])
rows <- strsplit(a[2], ",")[[1]]
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x,
               wig = function(x) sin(7 * x) + 0.2 * x,
               two = function(x) sin(x) + 0.4 * sin(10 * x),
               smooth = function(x) sin(x), step = function(x) 2 * (x > 5) - 1)
ns <- asNamespace("frmtmb")
em_theta <- function(fit, rounds = 1L) {
  fr <- fit$frame
  lp <- fr$linpreds[[1L]]
  bks <- Filter(function(bk) identical(bk$covstruct, "gp"), fr$re_blocks)
  rb <- bks
  at <- 0L
  for (i in seq_along(rb)) {
    ti <- rb[[i]]$theta_idx
    rb[[i]]$theta_idx <- at + seq_along(ti)
    at <- at + length(ti)
  }
  f0 <- list(re_blocks = rb, par_template = list(theta = numeric(at), beta = 0))
  y <- as.numeric(fr$y[[1L]])
  mu <- as.numeric(fitted(fit)[, "Estimate"])
  sg <- exp(fit$estimates$betad[1L])
  nu <- 1 + exp(fit$estimates$betad[2L])
  w <- (nu + 1) / (nu + ((y - mu) / sg)^2)
  n <- length(y)
  m <- ns$gp_marginal_opt(ns$gp_marginal_data(y, matrix(1, n, 1L), lp$Z, rb,
                                              sg^2 / w, rep(1, n), FALSE, 1L,
                                              f0))
  m$theta
}
for (rw in rows) {
  p <- strsplit(rw, ":")[[1]]
  cs <- p[1]; cn <- p[2]; n <- as.integer(p[3]); seed <- as.integer(p[4])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- curves[[cn]](x)
  d <- data.frame(x = x, y = eta + 0.02 * rt(n, 5))
  t0 <- proc.time()[[3]]
  q <- fe_fit(y ~ gp(x), data = d, family = student())
  el0 <- proc.time()[[3]] - t0
  ll0 <- as.numeric(logLik(q$fit))
  th <- tryCatch(em_theta(q$fit), error = function(e) {
    cat("EMERR", conditionMessage(e), "\n"); NULL })
  if (is.null(th)) next
  e <- q$fit$estimates
  t1 <- proc.time()[[3]]
  r <- fe_fit(y ~ gp(x), data = d, family = student(),
              start = list(beta = e$beta, betad = e$betad, theta = th))
  el1 <- proc.time()[[3]] - t1
  ll1 <- if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit))
  cat(sprintf(paste0("STP %s %s %d %d default %.4f conv %d nw %d (%.1fs) th %s",
                     " | em th %s -> %.4f conv %s nw %d (%.1fs) gain %.4f\n"),
              cs, cn, n, seed, ll0, q$fit$opt$convergence, length(q$warn),
              el0, paste(round(e$theta, 3), collapse = ","),
              paste(round(th, 3), collapse = ","), ll1,
              if (is.na(ll1)) "NA" else r$fit$opt$convergence, length(r$warn),
              el1, ll1 - ll0))
  flush.console()
}
cat("DONE\n")
