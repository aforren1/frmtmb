# Lane loglik fix round (m8): the review cost script, arms base / prefix /
# after (rellib-r18, rev-loglik-lib with REV_MUT unset, lane library).
# Arg 1: "after" or "before". Minimum and median of per-call times over
# repeated blocks, seed 401.
args <- commandArgs(TRUE)
arm <- args[1]
libs <- c("C:/Users/adf44/source/r/wt-loglik-lib",
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
libs[1] <- switch(arm, base = "C:/Users/adf44/source/r/rellib-r18", prefix = "C:/Users/adf44/source/r/rev-loglik-lib", after = "C:/Users/adf44/source/r/wt-loglik-lib")
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")),
    format(packageVersion("frmtmb")), "\n")
per_call <- function(f, reps, blocks = 7L) {
  f()
  t <- vapply(seq_len(blocks), function(i) {
    t0 <- proc.time()[["elapsed"]]
    for (k in seq_len(reps)) f()
    (proc.time()[["elapsed"]] - t0) / reps
  }, 0)
  sprintf("min %8.2f ms  median %8.2f ms", 1000 * min(t), 1000 * median(t))
}
set.seed(401)
d <- data.frame(x = rnorm(200), g = factor(rep(1:20, each = 10)),
                f = factor(sample(letters[1:20], 200, TRUE)))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(200)
p_small <- set_prior("normal(0, 1)", class = "b") +
  set_prior("student_t(3, 0, 1)", class = "sd") +
  set_prior("normal(1, 1)", class = "Intercept")
p_big <- p_small + set_prior("lkj(2)", class = "cor") +
  set_prior("exponential(1)", class = "sigma")
for (lv in letters[2:20]) {
  p_big <- p_big + set_prior("normal(0, 0.5)", class = "b",
                             coef = paste0("f", lv))
}
cases <- list(
  small_noprior = list(fm = y ~ x + (1 | g), pr = NULL),
  small_prior = list(fm = y ~ x + (1 | g), pr = p_small),
  big_prior = list(fm = y ~ x + f + (1 + x | g), pr = p_big))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit_once <- function() suppressWarnings(suppressMessages(
    frm(cs$fm, data = d, prior = cs$pr)))
  fit <- fit_once()
  cat(sprintf("%-14s frm()          %s\n", nm, per_call(fit_once, 1L, 5L)))
  cat(sprintf("%-14s logLik()       %s\n", nm,
              per_call(function() logLik(fit), 20L)))
  cat(sprintf("%-14s summary()      %s\n", nm,
              per_call(function() summary(fit), 5L)))
  cat(sprintf("%-14s print(fit)     %s\n", nm,
              per_call(function() utils::capture.output(print(fit)), 5L)))
  cat(sprintf("%-14s deviance()     %s\n", nm,
              per_call(function() deviance(fit), 20L)))
}
cat("done\n")
