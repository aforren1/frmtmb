# Lane optrap, punch round 2: the wide design's seed 1, where 0.69.0
# stops 9107.66 short at code 0 in silence with the group sd at 6e-16
# after 3 iterations (dev/optrap-p2-wide1.R). Which start brings it to
# lmer: sigma at the residual spread instead of sd(y), and also the mu
# coefficients at least squares.
#   Rscript dev/optrap-p2-wide1-starts.R <lib> [seed]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
s <- if (length(a) > 1) as.integer(a[2]) else 1L
cat("lib", dirname(find.package("frmtmb")), "seed", s, "\n")
set.seed(s)
n <- 101000
p <- 49
Z <- matrix(rnorm(n * p), n, p)
colnames(Z) <- paste0("z", seq_len(p))
d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                       "+ (1 | g)"))
l4 <- as.numeric(logLik(lmer(fo, data = d, REML = FALSE)))
ls <- lm(stats::reformulate(colnames(Z), "y"), data = d)
rs <- stats::mad(residuals(ls))
arms <- list(default = NULL,
             sigma_resid = list(betad = log(rs)),
             sigma_and_beta_ls = list(betad = log(rs),
                                      beta = unname(coef(ls))))
for (nm in names(arms)) {
  w <- character()
  t0 <- proc.time()[[3]]
  f <- withCallingHandlers(frm(fo, data = d, start = arms[[nm]]),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  cat(sprintf("%-18s code %d | logLik - lmer %.4g | group sd %.3g | iterations %d | warnings %d | secs %.1f\n",
              nm, f$opt$convergence, as.numeric(logLik(f)) - l4,
              exp(f$opt$par[names(f$opt$par) == "theta"]),
              f$opt$iterations %||% NA, length(w), proc.time()[[3]] - t0))
}
