# The numbers vignette("reinforcement-learning") and
# vignette("bayesian-cognitive-modeling") quote for importance =, on the
# fixed build. Seeds 4201:4206 at 100 draws, 40 subjects, as in the
# vignette; dataset 4206 at 100 and 400 draws; Exams_2 at 512 draws.
# Log dev/impsamp-out/rl.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
sys.source(system.file("rl", "rw-delta.R", package = "frmtmb",
                       mustWork = TRUE), envir = globalenv())
sq <- function(e) suppressWarnings(suppressMessages(e))
form <- bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
           beta ~ 1 + (1 | p | id))
pull <- function(f) {
  b <- fixef(f)[, "Estimate"]
  c(b, log_sd_alpha = log(VarCorr(f)$id$sd[1, "Estimate"]))
}
one <- function(s, nt, N) {
  set.seed(s)
  dd <- rl_simulate(rl_bandit_design(40, nt), seed = s)[[1L]]
  fam <- rw_delta(subject = id, trial = trial)
  lap <- sq(frm(form, family = fam, data = dd))
  t0 <- proc.time()[["elapsed"]]
  imp <- tryCatch(sq(frm(form, family = fam, data = dd, importance = N)),
                  error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(imp, "error")) {
    cat(sprintf("%d trials, seed %d, N=%d: ERROR %s\n", nt, s, N,
                substr(conditionMessage(imp), 1, 300)))
    return(NULL)
  }
  pl <- pull(lap)
  pi <- pull(imp)
  im <- imp$importance
  cat(sprintf(
    "%d trials, seed %d, N=%d: log sd %.4f -> %.4f mcse %.3f essmin %.3f rounds %d capped %s secs %.1f | shifts %s\n",
    nt, s, N, pl[["log_sd_alpha"]], pi[["log_sd_alpha"]], im$mcse, im$ess_min,
    im$rounds, im$capped, el,
    paste(format(pi - pl, digits = 2), collapse = " ")))
  list(shift = pi - pl, mcse = im$mcse, ess = im$ess_min)
}
for (nt in c(100L, 20L)) {
  r <- Filter(Negate(is.null), lapply(4201:4206, one, nt = nt, N = 100L))
  sh <- vapply(r, `[[`, numeric(length(r[[1L]]$shift)), "shift")
  cat(sprintf("%d trials: completed %d of 6\n", nt, length(r)))
  print(round(rowMeans(sh), 3))
  cat("largest abs shift, log sd:",
      round(max(abs(sh["log_sd_alpha", ])), 3), "\n")
  cat("mcse range:", round(range(vapply(r, `[[`, 0, "mcse")), 3),
      "; min ess:", round(min(vapply(r, `[[`, 0, "ess")), 3), "\n")
}
invisible(one(4206L, 20L, 400L))
d <- data.frame(k = c(21, 17, 21, 18, 22, 31, 31, 34, 34, 35, 35, 36,
                      39, 36, 35),
                n = 40L, zero = 0, id = factor(1:15))
fo <- bf(k | trials(n) ~ 1 + (1 | id), mu2 ~ 0 + offset(zero),
         theta1 ~ 0 + offset(zero))
lap <- sq(frm(fo, family = mixture(binomial, binomial), data = d))
imp <- sq(frm(fo, family = mixture(binomial, binomial), data = d,
              importance = 512))
cat("Exams_2: Laplace", format(as.numeric(logLik(lap)), digits = 7),
    "corrected", format(as.numeric(logLik(imp)), digits = 7), "mcse",
    format(imp$importance$mcse, digits = 3), "rounds", imp$importance$rounds,
    "\n")
