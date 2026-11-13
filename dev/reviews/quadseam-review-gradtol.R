# Review of lane quadseam: the documented control
# frmtmb_control(grad_tol = 1e-4) on a plain one-block Bernoulli fit
# (500 groups of 4, sd 1.2; seeds 1 and 3 of the review's bern500
# design), uninstrumented, with every warning's class and message, the
# settle record and seconds.
#   Rscript quadseam-review-gradtol.R <base|lane> <seeds>
a <- commandArgs(TRUE)
arm <- a[1]
seeds <- eval(parse(text = a[2]))
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
suppressMessages(library(frmtmb))
cat("lib:", arm, dirname(find.package("frmtmb")), "\n")
for (s in seeds) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:500, each = 4)), x = rnorm(2000))
  u <- rnorm(500, 0, 1.2)
  d$y <- rbinom(2000, 1, plogis(-0.5 + 0.7 * d$x + u[d$g]))
  w <- list()
  t0 <- proc.time()[[3]]
  f <- withCallingHandlers(
    frm(y ~ x + (1 | g), family = bernoulli(), data = d, quadrature = TRUE,
        control = frmtmb_control(grad_tol = 1e-4)),
    warning = function(x) {
      w[[length(w) + 1L]] <<- x; invokeRestart("muffleWarning")
    })
  tt <- proc.time()[[3]] - t0
  qs <- f$opt[["quad_settle"]]
  cat(sprintf("FIT seed %d %s logLik %.10f conv %d settle %s secs %.1f nwarn %d\n",
              s, arm, as.numeric(logLik(f)), f$opt$convergence,
              if (is.null(qs)) "-" else paste(format(qs, digits = 3),
                                              collapse = "/"),
              tt, length(w)))
  for (x in w) {
    cat("  WARN", class(x)[1], ":",
        substr(gsub("\\s+", " ", conditionMessage(x)), 1, 700), "\n")
  }
}
cat("END\n")
