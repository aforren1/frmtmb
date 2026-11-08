# Brief step 10 (readers x latperf): a nonlinear body on `shape` beside a
# linear `mu` with a latent residual (`ar(cov = TRUE)`), the case that
# reaches both lanes' edits of frm_linpred() and predict_mean_se(). Its
# standard errors, in sample and on new rows, against frm_lp_basis()'s
# A V A' (a dense reference independent of the sparse route), and
# fitted() and conditional_effects() run.
#
#   Rscript dev/rel080-latnl.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1L], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
set.seed(9101)
S <- if (length(args) > 1L) as.integer(args[2L]) else 6L; Tn <- 25
d <- expand.grid(t = seq_len(Tn), s = factor(seq_len(S)))
d$x <- rnorm(nrow(d)); d$z <- runif(nrow(d), -1, 1)
e <- unlist(lapply(seq_len(S), function(i) {
  as.numeric(arima.sim(list(ar = 0.6), Tn, sd = 0.3))
}))
mu <- exp(0.8 + 0.4 * d$x + e)
d$y <- rnbinom(nrow(d), mu = mu, size = exp(1 + 0.8 * d$z))
f <- bf(y ~ x + ar(time = t, gr = s, cov = TRUE)) +
  nlf(shape ~ exp(a + b * z)) + lf(a ~ 1, b ~ 1)
ws <- character()
fit <- withCallingHandlers(frm(f, data = d, family = negbinomial()),
  warning = function(w) {
    ws <<- c(ws, class(w)[1L]); invokeRestart("muffleWarning")
  })
cat("fit: logLik", format(as.numeric(logLik(fit)), digits = 9), "conv",
    fit$opt$convergence, "warnings", if (length(ws)) paste(ws, collapse = ",") else "none", "\n")
nd <- data.frame(t = c(Tn + 1, Tn + 2), s = factor(c(1, 2), levels = levels(d$s)),
                 x = c(0, 1), z = c(-0.5, 0.5))
for (where in c("in sample", "new rows")) {
  newd <- if (where == "new rows") nd else NULL
  for (dp in c("mu", "shape")) {
    lp <- frm_linpred(fit, newdata = newd, dpar = dp, se.fit = TRUE)
    B <- frm_lp_basis(fit, newdata = newd, dpar = dp)
    A <- as.matrix(B$A)
    ref <- sqrt(rowSums((A %*% B$V) * A))
    se <- as.numeric(lp$se.fit)
    cat(sprintf("%-9s %-5s rows %3d | max |rel| se.fit vs A V A' %.2e | finite %d of %d\n",
                where, dp, length(se), max(abs(se / ref - 1)), sum(is.finite(se)),
                length(se)))
  }
  fv <- fitted(fit, newdata = newd)
  cat(sprintf("%-9s fitted(): Est.Error finite %d of %d, range %.4g to %.4g\n",
              where, sum(is.finite(fv[, "Est.Error"])), nrow(fv),
              min(fv[, "Est.Error"]), max(fv[, "Est.Error"])))
}
ce <- conditional_effects(fit, effects = "z")
cat("conditional_effects(z): rows", nrow(ce[[1]]), "se__ finite",
    sum(is.finite(ce[[1]]$se__)), "\n")
cat("LATNL DONE\n")
