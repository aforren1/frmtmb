# Item 4, the cov = FALSE form: does mi_rows_joint() keep an unread missing
# row of a latent residual's response read, so the refusal still fires where
# its reason is false? Gamma(log) and poisson-free lognormal, `ar(t, g)`
# (cov = FALSE), 3 missing rows, seed 12301 as test-latres.R. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
set.seed(12301)
dg <- expand.grid(t = 1:8, g = factor(1:10))
dg <- dg[order(dg$g, dg$t), ]
dg$x <- rnorm(nrow(dg))
dg$y <- rgamma(nrow(dg), 5, 5 / exp(0.3 * dg$x))
dg$pos <- exp(0.3 * dg$x + rnorm(nrow(dg), 0, 0.5))
dg$y[c(3L, 20L, 41L)] <- NA
dg$pos[c(3L, 20L, 41L)] <- NA
for (sp in list(list(f = bf(y | mi() ~ x + ar(t, g)), fam = Gamma("log")),
                list(f = bf(pos | mi() ~ x + ar(t, g)), fam = lognormal()),
                list(f = bf(y | mi() ~ x + ma(t, g)), fam = Gamma("log")))) {
  r <- fe_fit(sp$f, family = sp$fam, data = dg)
  fr <- tryCatch(frm(sp$f, family = sp$fam, data = dg, dry_run = "frame"),
                 error = function(e) NULL)
  rd <- if (!is.null(fr)) {
    rn <- names(fr$mi_map)[1]
    c(joint = ns$mi_rows_joint(fr, rn),
      read = length(ns$mi_ml_frame(fr)$mi_map[[rn]]$rows))
  }
  cat(sprintf("MC %s | %s | %s\n", deparse1(sp$f$formula),
              if (inherits(r$fit, "error")) paste("ERROR", substr(conditionMessage(r$fit), 1, 140))
              else sprintf("ll %.6f conv %d", as.numeric(logLik(r$fit)), r$fit$opt$convergence),
              paste(names(rd), rd, collapse = " ")))
}
cat("DONE\n")
