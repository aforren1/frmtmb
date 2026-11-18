# Reviewer of lane drawdir: the fits of the reviewer's survey that the
# per-direction reading fires on, measured with more draws and against
# a sampled posterior. Per id (dev/reviews/drawdir-review-survey.R): the
# lane's reweighted draws and the normal law's (verdict's fire set
# FALSE), 20000 draws each, seed 1: in-sample MAE of the predictive mean
# against the true conditional mean, and for designs with a new-level
# oracle its coverage of the central 90% and 98% intervals. Then
# frm_sample() (brms's default priors, 4 chains of 1500 after 1500
# warmup, seed 1): MAE of posterior_epred()'s mean against the truth, and
# for the sleepstudy-like designs the new level's predictive drawn by hand
# from the posterior draws.
# Usage: Rscript dev/reviews/drawdir-review-fires.R <id> [...]
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ex <- parse("dev/reviews/drawdir-review-survey.R")
keep <- vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("fq", "ilogit", "cum_mean", "cum_draw",
                                "design", "seeds_of", "said_of")
}, TRUE)
for (e in ex[keep]) eval(e, globalenv())
cover <- function(pn, o) {
  q90 <- stats::quantile(pn, c(0.05, 0.95), names = FALSE)
  q98 <- stats::quantile(pn, c(0.01, 0.99), names = FALSE)
  sprintf("cov90 %.3f w90 %8.2f cov98 %.3f w98 %8.2f",
          mean(o >= q90[1] & o <= q90[2]), q90[2] - q90[1],
          mean(o >= q98[1] & o <= q98[2]), q98[2] - q98[1])
}
for (id in args) {
  p <- strsplit(id, ":")[[1]]
  r <- design(p[1], as.integer(p[2]))
  fit <- r$fit
  o <- NULL
  if (!is.null(r$onew)) {
    set.seed(99)
    o <- r$onew(2e5)
  }
  cat("==", id, "\n")
  arm <- function(lab) {
    set.seed(1)
    s <- said_of(predict(fit, ndraws = 20000, summary = FALSE))
    pin <- s$val
    chk <- fit$cache[["draw_check"]]
    line <- sprintf("  %-5s mae %.4f | %s", lab,
                    mean(abs(colMeans(pin) - r$mu)), s$said)
    if (!is.null(o)) {
      set.seed(1)
      pn <- fq(predict(fit, newdata = r$nd, allow_new_levels = TRUE,
                       ndraws = 20000, summary = FALSE))[, 1]
      line <- paste(line, "|", cover(pn, o))
    }
    cat(line, sprintf("| fire %s dir %.3f (%s)\n", chk$fire, chk$dir,
                      chk$dir_par))
  }
  arm("lane")
  chk <- fit$cache[["draw_check"]]
  chk$fire <- FALSE
  assign("draw_check", chk, envir = fit$cache)
  arm("norm")
  if (p[1] == "sigid") next
  t0 <- proc.time()[[3]]
  ds <- tryCatch(suppressMessages(suppressWarnings(
    frm_sample(fit, chains = 4, iter = 3000, warmup = 1500, seed = 1,
               refresh = 0, cores = 1))), error = function(e) {
                 cat("  frm_sample error:", conditionMessage(e), "\n")
                 NULL
               })
  if (is.null(ds)) next
  ep <- suppressWarnings(posterior_epred(ds))
  rh <- tryCatch(max(rhat(ds), na.rm = TRUE), error = function(e) NA)
  line <- sprintf("  mcmc  mae %.4f (epred, %d draws, %.0f s, max rhat %.3f)",
                  mean(abs(colMeans(ep) - r$mu)), nrow(ep),
                  proc.time()[[3]] - t0, rh)
  if (!is.null(o) && p[1] %in% c("slp", "slpw", "slp10")) {
    A <- as.matrix(ds)
    col <- function(pat) A[, grep(pat, colnames(A))[1]]
    s0 <- col("^sd_Subject__Intercept$")
    s1 <- col("^sd_Subject__Days$")
    rho <- col("^cor_Subject__Intercept__Days$")
    sg <- col("^sigma$")
    set.seed(1)
    pn <- col("^b_Intercept$") + 9 * col("^b_Days$") +
      stats::rnorm(nrow(A)) * sqrt(s0^2 + 81 * s1^2 + 18 * rho * s0 * s1 +
                                     sg^2)
    line <- paste(line, "|", cover(pn, o))
  }
  cat(line, "\n")
}
cat("DONE\n")
