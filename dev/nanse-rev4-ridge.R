# Reviewer, punch 2b: the risk of inverting a positive-definite
# remainder as it is. On the finite-difference path (models with random
# effects) an exact fixed-effect ridge can come out slightly positive
# from the Hessian's noise. Tier 3 is reached when something else fails
# (an empty theta row, a variance at 0); then the remainder, ridge
# included, may be inverted with no flat test. Designs, 12 seeds each,
# group variance 0 (an empty theta row) and 0.5:
#   D1  y ~ a + c + dd, a ~ 0 + x1, c ~ 1 + (1 | g), dd ~ 1  (c + dd ridge)
#   D2  y ~ a + b, a ~ 0 + f (k = 20), b ~ 1 + (1 | g)        (spread ridge)
#   D3  y ~ a + b, a ~ 0 + f, b ~ 0 + h + (1 | g), k = 10     (crossed ridge)
# Per fit: the ridge parameters with a finite SE (should be 0), whether
# a warning or the deferred report named them, and the tier.
#   Rscript dev/nanse-rev4-ridge.R lane|r2|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  r2 = c("C:/Users/adf44/source/r/nanse-rev-lib",
         "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
phrase <- "not available|are not identified"
mk <- function(design, s, sdg) {
  set.seed(s)
  if (design == "D1") {
    n <- 300
    d <- data.frame(x1 = rnorm(n), g = factor(rep(1:15, 20)))
    d$y <- 1 + 0.5 * d$x1 + rnorm(15, 0, sdg)[d$g] + rnorm(n)
    list(d = d, f = bf(y ~ a + c + dd, a ~ 0 + x1, c ~ 1 + (1 | g), dd ~ 1,
                       nl = TRUE), ridge = c("c_(Intercept)", "dd_(Intercept)"))
  } else if (design == "D2") {
    k <- 20
    f <- factor(rep(seq_len(k), each = 5))
    d <- data.frame(f = f, g = factor(rep(1:10, length.out = 5 * k)))
    d$y <- rnorm(k)[f] + rnorm(10, 0, sdg)[d$g] + rnorm(5 * k, 0, 0.5)
    list(d = d, f = bf(y ~ a + b, a ~ 0 + f, b ~ 1 + (1 | g), nl = TRUE),
         ridge = c(paste0("a_f", seq_len(k)), "b_(Intercept)"))
  } else {
    k <- 10
    d <- expand.grid(f = factor(seq_len(k)), h = factor(seq_len(k)))
    d <- d[sample(nrow(d), 6 * k), ]
    d$g <- factor(rep(1:8, length.out = nrow(d)))
    d$y <- rnorm(k)[d$f] + rnorm(k)[d$h] + rnorm(8, 0, sdg)[d$g] +
      rnorm(nrow(d), 0, 0.5)
    list(d = d, f = bf(y ~ a + b, a ~ 0 + f, b ~ 0 + h + (1 | g), nl = TRUE),
         ridge = c(paste0("a_f", seq_len(k)), paste0("b_h", seq_len(k))))
  }
}
tab <- NULL
for (design in c("D1", "D2", "D3")) {
  for (sdg in c(0, 0.5)) {
    for (s in 1:12) {
      m <- mk(design, s, sdg)
      w <- character()
      fit <- tryCatch(withCallingHandlers(frm(m$f, data = m$d),
        warning = function(x) {
          w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
        }, message = function(x) invokeRestart("muffleMessage")),
        error = function(e) NULL)
      if (is.null(fit)) next
      wv <- character()
      V <- withCallingHandlers(vcov(fit, full = TRUE), warning = function(x) {
        wv <<- c(wv, conditionMessage(x)); invokeRestart("muffleWarning")
      })
      se <- suppressWarnings(sqrt(diag(V)))
      names(se) <- ns$outer_par_names(fit)
      r <- intersect(m$ridge, names(se))
      lost <- ns$sdr_of(fit)$se_lost
      tab <- rbind(tab, data.frame(
        design, sdg, seed = s, code = fit$opt$convergence,
        ridge_finite = sum(is.finite(se[r])), n_ridge = length(r),
        min_finite_se = if (any(is.finite(se[r]))) min(se[r][is.finite(se[r])])
          else NA_real_,
        lost = length(lost),
        warned = any(grepl(phrase, c(w, wv)))))
    }
  }
}
print(tab, row.names = FALSE)
cat("\nfits with a finite SE on a ridge parameter, and of those unwarned:\n")
tab$silent <- tab$ridge_finite > 0 & !tab$warned
print(aggregate(cbind(fits = 1, ridge_finite_fits = ridge_finite > 0, silent)
                ~ design + sdg, tab, sum))
