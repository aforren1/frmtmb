# Review of lane illse2: the REML log-likelihood of a fit taped in the
# design's coordinates (beta_whitener()) against lme4, 0.71.0 and the
# Jacobian log|det T|.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-remlll.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
data("sleepstudy", package = "lme4")
set.seed(1)
ssu <- sleepstudy[sort(sample(nrow(sleepstudy), 120)), ]
ldT <- function(f) {
  bw <- f$beta_whiten
  if (is.null(bw)) return(0)
  sum(vapply(bw$blocks, function(b) {
    determinant(b$T, logarithm = TRUE)$modulus[[1]]
  }, 0))
}
for (sh in c(0, 1e3, 1e4, 1e5, 1e6)) {
  d <- ssu
  d$t <- d$Days + sh
  f <- suppressWarnings(frm(Reaction ~ t + (1 | Subject), data = d,
                            REML = TRUE))
  l <- lme4::lmer(Reaction ~ t + (1 | Subject), data = d, REML = TRUE)
  fp <- suppressWarnings(frm(Reaction ~ t + (1 | Subject), data = d,
                             control = frmtmb_control(profile = TRUE)))
  fm <- suppressWarnings(frm(Reaction ~ t + (1 | Subject), data = d))
  cat(sprintf(paste0("shift %-6g | REML logLik frm %.8f lme4 %.8f diff %+.6f",
                     " | log|det T| %+.6f | AIC frm %.4f lme4 %.4f\n",
                     "             | profile logLik %.8f ML %.8f\n"),
              sh, as.numeric(logLik(f)), as.numeric(logLik(l)),
              as.numeric(logLik(f)) - as.numeric(logLik(l)), ldT(f),
              AIC(f), AIC(l), as.numeric(logLik(fp)),
              as.numeric(logLik(fm))))
}
# anova of two REML fits with the same fixed effects, one with a random
# slope: whitened designs share T, so the difference should hold
d <- ssu
d$t <- d$Days + 1e5
f1 <- suppressWarnings(frm(Reaction ~ t + (1 | Subject), data = d,
                           REML = TRUE))
l1 <- lme4::lmer(Reaction ~ t + (1 | Subject), data = d, REML = TRUE)
f0 <- suppressWarnings(frm(Reaction ~ Days + (1 | Subject), data = d,
                           REML = TRUE))
cat(sprintf("same model, Days vs Days + 1e5: REML logLik %.8f vs %.8f\n",
            as.numeric(logLik(f0)), as.numeric(logLik(f1))))
