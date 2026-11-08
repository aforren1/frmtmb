# Lane bigx, fix round, major 2: raw polynomials of x on [1000, 1010]
# (the review's dev/review-bigx-collin.R construction, data seeds 501 to
# 503), against glm(): which columns the frame keeps, whether the
# optimizer's basis engages, the log-likelihood and the SEs.
#   BIGX_ARM=before|after Rscript dev/bigx-polyfar.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
for (s in 501:503) for (deg in 5:8) {
  set.seed(s)
  n <- 200
  x1 <- rnorm(n); invisible(rnorm(n)); invisible(rnorm(n)); invisible(rnorm(n))
  xx <- runif(n, 1000, 1010)
  d <- data.frame(xx = xx, yy = sin(xx / 2) + rnorm(n, 0, 0.3))
  fo <- as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg))
  msg <- character(0)
  x <- withCallingHandlers(bigx_fit(bf(fo), family = gaussian(), data = d),
                           message = function(m) {
                             msg <<- c(msg, conditionMessage(m))
                             invokeRestart("muffleMessage")
                           })
  g <- glm(fo, data = d)
  if (x$err) {
    cat(sprintf("seed %d deg %d ERROR %s\n", s, deg, conditionMessage(x$f)))
    next
  }
  f <- x$f
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(f)))),
                 error = function(e) NA)
  cat(sprintf(paste0("seed %d deg %d kept %d (glm %d) basis %s code %d ",
                     "ll - glm %+.3e SEs finite %d of %d warn %d %s\n"),
              s, deg, length(f$estimates[["beta"]]), sum(!is.na(coef(g))),
              !is.null(f$obj$env$.frmtmb_opt_basis), f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(g)),
              sum(is.finite(se)), length(se), length(x$w),
              if (length(x$w)) substr(x$w[1], 1, 60) else ""))
  if (length(msg)) cat("   message:", substr(msg[1], 1, 120), "\n")
}
