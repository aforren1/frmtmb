# Review of lane quadseam2, attack 5: FEWBERN seed 3 of
# quadseam2-review-latres.R lost frmtmb_latres_laplace on the lane, with
# latres_quadratic() TRUE on a bernoulli. What step did the gate use, and
# what does it read at a step of 1? Also the step-scaled form: a sweep of
# latres_marginal_sd() through h in 1e-3..1 on the same rows.
#   Rscript quadseam2-review-fewbern.R <base|lib>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
gen <- function(seed, G = 4L, K = 5L, sd = 0.8) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  d$e <- as.vector(replicate(G, arima.sim(list(ar = 0.5), K, sd = sd)))
  d
}
for (seed in c(703, 701, 702, 704)) {
  d <- gen(seed)
  d$yp <- rpois(nrow(d), exp(0.3 * d$x + d$e))
  d$yb <- rbinom(nrow(d), 1, plogis(0.3 * d$x + d$e))
  w <- character()
  f <- withCallingHandlers(suppressMessages(
    frm(bf(yb ~ x + ar(t, g, cov = TRUE)), family = bernoulli(), data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  bk <- ns$frame_latres(f$frame)[[1L]]
  th <- f$estimates$theta[bk$theta_idx]
  h <- ns$latres_marginal_sd(th, bk$aux)
  cat(sprintf("SEED %d theta %s marginal sd %.4g\n", seed,
              paste(format(th, digits = 4), collapse = ","), h))
  for (x in w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 200), "\n")
  if (exists("latres_quadratic", ns)) {
    cat("   gate:", ns$latres_quadratic(f, bk), "\n")
    # the gate's differences at other steps, by replacing the sd
    for (hh in c(1e-3, 1e-2, 0.1, 0.3, 1)) {
      a2 <- bk
      orig <- ns$latres_marginal_sd
      unlockBinding("latres_marginal_sd", ns)
      assign("latres_marginal_sd", function(theta, a) hh, envir = ns)
      g <- ns$latres_quadratic(f, a2)
      assign("latres_marginal_sd", orig, envir = ns)
      lockBinding("latres_marginal_sd", ns)
      cat(sprintf("   gate at step %g: %s\n", hh, g))
    }
  }
}
cat("END\n")
