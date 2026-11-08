# Review of lane loglik, attack 5 and the meaning of log_prior.
# Arg 1: "after" (lane library first) or "before" (rellib-r18 alone).
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "after"
libs <- c("C:/Users/adf44/source/r/wt-loglik-lib",
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "before") libs <- libs[-1]
# the review's mutant build, REV_MUT set by the caller
if (arm == "rev") libs[1] <- "C:/Users/adf44/source/r/rev-loglik-lib"
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")),
    format(packageVersion("frmtmb")), "\n")
has_pen <- exists("logLik_penalized", asNamespace("frmtmb"))
run <- function(lbl, expr) {
  ws <- character(0)
  f <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, paste0("<", class(w)[1], "> ", substr(conditionMessage(w), 1, 90)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  pen <- -f$opt$objective
  cat(sprintf("%-44s pen %.6f  logLik %.6f  warnings %d\n", lbl, pen,
              as.numeric(logLik(f)), length(ws)))
  for (w in ws) cat("    ", w, "\n")
  invisible(f)
}

# A. all-zero counts, a tight prior that puts the penalized value above 0
set.seed(301)
d <- data.frame(y = rep(0L, 50))
run("A poisson all 0, Intercept normal(-5, 0.01)",
    frm(y ~ 1, family = poisson(), data = d,
        prior = set_prior("normal(-5, 0.01)", class = "Intercept")))
# B. bernoulli, near-separation, with (1 | g) and a strong sd prior
set.seed(302)
d <- data.frame(g = factor(rep(1:30, each = 2)), x = rnorm(60))
d$y <- as.integer(d$x > 0)
run("B bernoulli separated x, b normal(0, 0.01)",
    frm(y ~ x, family = bernoulli(), data = d,
        prior = set_prior("normal(20, 0.01)", class = "b")))
run("C bernoulli (1|g) all 0, sd normal(5, 0.01)",
    frm(y ~ 1 + (1 | g), family = bernoulli(), data = transform(d, y = 0L),
        prior = set_prior("normal(5, 0.01)", class = "sd")))
# D. binomial with many trials, a nearly perfect fit
set.seed(303)
d <- data.frame(n = 50L, x = rnorm(40))
d$y <- rbinom(40, d$n, plogis(-6 + 0.1 * d$x))
run("D binomial rare events, b normal(0, 0.001)",
    frm(y | trials(n) ~ x, family = binomial(), data = d,
        prior = set_prior("normal(0, 0.001)", class = "b") +
          set_prior("normal(-6, 0.01)", class = "Intercept")))

if (has_pen) {
  cat("\n## what log_prior holds, against brms's lprior\n")
  set.seed(304)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:12, each = 5)))
  d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.8)[d$g] + rnorm(60)
  f <- frm(y ~ x + (1 | g), data = d,
           prior = set_prior("normal(0, 1)", class = "sd"))
  sdv <- VarCorr(f)$g$sd[1, "Estimate"]
  lp <- attr(logLik_penalized(f), "log_prior")
  cat(sprintf("log_prior %.10f\n", lp))
  cat(sprintf("dnorm(sd)               %.10f\n", dnorm(sdv, 0, 1, log = TRUE)))
  cat(sprintf("dnorm(sd) + log(sd)     %.10f  (Jacobian of log sd)\n",
              dnorm(sdv, 0, 1, log = TRUE) + log(sdv)))
  cat(sprintf("brms half: dnorm + log 2 %.10f\n",
              dnorm(sdv, 0, 1, log = TRUE) + log(2)))
  if (requireNamespace("brms", quietly = TRUE)) {
    sc <- brms::stancode(y ~ x + (1 | g), data = d,
                         prior = brms::prior(normal(0, 1), class = sd))
    cat("brms", format(packageVersion("brms")), "lprior lines:\n")
    cat(grep("lprior", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
    sc <- brms::stancode(y ~ x, data = d,
                         prior = brms::prior(normal(0, 1), class = b, lb = 0))
    cat(grep("lprior", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  }
  fb <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("normal(0, 1)", class = "b", coef = "x", lb = 0))
  bx <- fb$estimates$beta[["x"]]
  cat(sprintf("b lb=0: log_prior %.10f  dnorm(b) %.10f  brms adds -normal_lccdf(0) = %.10f\n",
              attr(logLik_penalized(fb), "log_prior"),
              dnorm(bx, 0, 1, log = TRUE), log(2)))
}
cat("done\n")
