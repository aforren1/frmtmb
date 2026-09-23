# Does frmtmb accept brms's gr(g, by = f), the one brms spelling in which
# a random-effect SD depends on a group-level (factor) covariate?
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
d <- sim_grouped(101)
d$wf <- factor(d$w > 0)
r <- tryCatch(frm(bf(y ~ x + (1 | gr(g, by = wf))), family = gaussian(),
                  data = d),
              error = function(e) e)
if (inherits(r, "error")) {
  cat("frm gr(by =): REFUSED:", conditionMessage(r), "\n")
} else {
  cat("frm gr(by =): FITTED logLik", format(logLik(r), digits = 12), "\n")
  print(VarCorr(r))
}
f2 <- drmTMB::drmTMB(dbf(y ~ x + (1 | g), sigma ~ 1, sd(g) ~ wf),
                     family = gaussian(), data = d)
cat("drm sd(g) ~ wf: logLik", format(logLik(f2), digits = 12), "\n")
