# Review of lane nlstart, short check of the fix round: prior classes the
# unit's precision formula does not cover (nl_prior_precision() gives 0),
# on `a * xs`, xs = x * s. Each row against the same model on rellib-r18.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-fix-prior2.R
source("dev/review-nlstart-common.R")
priors <- c(
  "gamma(2, 1000)", "lognormal(0, 0.001)", "lognormal(0, 1)",
  "exponential(1000)", "cauchy(0, 0.001)", "cauchy(0, 1)",
  "student_t(3, 0, 0.001)", "logistic(0, 0.001)",
  "uniform(-1, 1)", "uniform(0, 1e21)", "normal(0, 0.001)")
cat("ROWHEAD prior s seed ll conv nwarn a units\n")
for (seed in c(6201, 6203)) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  for (s in c(1e-4, 1e-8, 1e-20)) {
    d$xs <- d$x * s
    for (pr in priors) {
      lb <- if (grepl("^(gamma|lognormal|exponential)", pr)) 0 else NA
      sp <- if (grepl("^uniform", pr)) {
        lim <- as.numeric(regmatches(pr, gregexpr("-?[0-9.e+]+", pr))[[1]])
        set_prior(pr, nlpar = "a", lb = lim[1], ub = lim[2])
      } else if (!is.na(lb)) {
        set_prior(pr, nlpar = "a", lb = 0)
      } else set_prior(pr, nlpar = "a")
      # a prior bounded below at 0 cannot start at 0
      st <- if (!is.na(lb)) list(beta = c("a_(Intercept)" = 1))
      f <- tt(frm(bf(y ~ a * xs, a ~ 1, nl = TRUE), data = d, prior = sp,
                  start = st))
      cat(sprintf("ROW %s %g %d %.8f %s %d %s %s\n", gsub(" ", "", pr), s, seed,
                  llv(f), if (iserr(f)) NA else f$opt$convergence,
                  length(attr(f, "w")),
                  if (iserr(f)) gsub("\\s+", "_", substr(f$msg, 1, 60)) else
                    format(f$estimates$beta[[1]], digits = 6),
                  if (iserr(f)) "-" else paste(format(f$par_units %||% 1,
                                                      digits = 3),
                                               collapse = ",")))
    }
  }
  flush.console()
}
cat("DONE\n")
