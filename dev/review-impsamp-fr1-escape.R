# Fix round 1: frm()'s Laplace route has sd_up_escape() ("restarted from a larger
# sd"); importance_fit() warm-starts without it. How often does the importance fit
# stay at a false collapse the Laplace route leaves? Standard columns as well as
# shifted ones, data seeds 1:40. Reports only fits where the two routes disagree on
# whether the sd is collapsed, and what the user was told.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
ord <- function(seed, m) {
  set.seed(seed)
  g <- factor(rep(seq_len(15L), each = 10))
  z <- rnorm(150)
  lat <- 0.7 * z + rnorm(15, 0, 0.7)[g] + rlogis(150)
  list(data = data.frame(y = ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf)))),
                         w = 5 * z + m, g = g),
       formula = bf(y ~ w + (1 | g)), family = cumulative(), draws = 200L)
}
bern_small <- function(seed) {
  set.seed(seed)
  g <- factor(rep(seq_len(25L), each = 3))
  x <- rnorm(75)
  list(data = data.frame(y = rbinom(75, 1, plogis(0.3 + 0.5 * x + rnorm(25, 0, 1)[g])), x, g),
       formula = bf(y ~ x + (1 | g)), family = bernoulli(), draws = 200L)
}
pois_small <- function(seed) {
  set.seed(seed)
  g <- factor(rep(seq_len(20L), each = 2))
  x <- rnorm(40)
  list(data = data.frame(y = rpois(40, exp(-0.5 + 0.3 * x + rnorm(20, 0, 0.8)[g])), x, g),
       formula = bf(y ~ x + (1 | g)), family = poisson(), draws = 200L)
}
gauss_small <- function(seed) {
  set.seed(seed)
  g <- factor(rep(seq_len(8L), each = 3))
  x <- rnorm(24)
  list(data = data.frame(y = 0.2 + 0.3 * x + rnorm(8, 0, 0.5)[g] + rnorm(24), x, g),
       formula = bf(y ~ x + (1 | g)), family = gaussian(), draws = 64L)
}
ds <- list(
  bern_small = bern_small, pois_small = pois_small, gauss_small = gauss_small,
  unbal_bern = design_list$unbal_bern, t_outlier = design_list$t_outlier,
  ordinal_std = function(s) ord(s, 0), ordinal_2000 = function(s) ord(s, 2000),
  smallpois = design_list$smallpois, wtbin = design_list$wtbin)
for (nm in names(ds)) {
  n <- 0L; dis <- 0L
  for (s in 1:40) {
    de <- ds[[nm]](s)
    l <- suppressWarnings(tryCatch(frm(de$formula, family = de$family, data = de$data), error = function(e) NULL))
    w <- character(0)
    i <- withCallingHandlers(
      tryCatch(frm(de$formula, family = de$family, data = de$data, importance = de$draws), error = function(e) e),
      warning = function(x) {
        w <<- c(w, substr(conditionMessage(x), 1, 90))
        invokeRestart("muffleWarning")
      })
    if (is.null(l) || inherits(i, "error")) next
    n <- n + 1L
    tl <- l$opt$par[names(l$opt$par) == "theta"][[1]]
    ti <- i$opt$par[names(i$opt$par) == "theta"][[1]]
    gain <- as.numeric(logLik(l)) - as.numeric(logLik(i))
    if (gain > 0.05 + 3 * i$importance$mcse) {
      dis <- dis + 1L
      cat(sprintf("  %s s %d: Laplace ll is %.3f ABOVE importance (mcse %.4f) log sd lap %.3f imp %.3f rounds %d warns: %s\n",
                  nm, s, gain, i$importance$mcse, tl, ti, i$importance$rounds,
                  if (length(w)) paste(w, collapse = " | ") else "NONE"))
    }
  }
  cat(sprintf("%s: %d of %d fits read below the Laplace fit by more than its noise\n", nm, dis, n))
}
