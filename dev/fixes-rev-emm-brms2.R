# Reviewer of lane fixes, claim 1: the sigma ~ scale(z) case of
# dev/fixes-rev-emm-brms.R with brms's own bf() for the brms side.
#   Rscript dev/fixes-rev-emm-brms2.R <lib>
src <- readLines("dev/fixes-rev-emm-brms.R")
src <- src[seq_len(grep("^run[(]\"ns", src) - 1L)]
eval(parse(text = src))
fit <- frm(bf(yg ~ f, sigma ~ scale(z)), data = d, family = gaussian())
bb <- brms_fixed_fit(brms::bf(yg ~ f, sigma ~ scale(z)), gaussian(), d, fit,
                     ndraws = 4)
a <- summary(emmeans(fit, ~ z, dpar = "sigma", at = list(z = c(1, 3))))$emmean
b <- summary(emmeans(bb, ~ z, dpar = "sigma", at = list(z = c(1, 3))))$emmean
cat("frmtmb", fmt(a), "\nbrms  ", fmt(b), "\n")
cat(sprintf("max rel diff %.3g\n", max(abs(a - b) / abs(b))))
