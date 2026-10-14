# Reviewer, lane dblbar: an independent fit check of (1 + f || g).
# frmtmb (lane library) against lme4 with explicit dummy columns and
# against glmmTMB's diag(). Data: dev/dblbar-data.R, seed 20261013.
# Output: dev/dblbar-rev-log/indep.txt
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(frmtmb); library(lme4); library(glmmTMB)
})
sink(file.path(wt, "dev/dblbar-rev-log/indep.txt"), split = TRUE)
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| lme4", as.character(packageVersion("lme4")), "| glmmTMB",
    as.character(packageVersion("glmmTMB")), "\n")
dd <- dblbar_data()
dd$fb <- as.numeric(dd$f == "b")
dd$fc <- as.numeric(dd$f == "c")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
show <- function(lab, ll_f, ll_o, sd_f, sd_o, b_f, b_o) {
  cat(sprintf("%-26s logLik frm %.10f other %.10f diff %.3e\n", lab, ll_f,
              ll_o, ll_f - ll_o))
  cat(sprintf("%-26s sd frm   %s\n", "", paste(sprintf("%.6f", sd_f),
                                               collapse = " ")))
  cat(sprintf("%-26s sd other %s  max rel %.3e\n", "",
              paste(sprintf("%.6f", sd_o), collapse = " "), rel(sd_f, sd_o)))
  cat(sprintf("%-26s b max rel %.3e\n", "", rel(b_f, b_o)))
}
frm_sd <- function(fit) {
  v <- variables(fit)
  e <- as.numeric(unlist(lapply(grep("^sd_", v, value = TRUE), function(n) {
    summary(fit)$random[[1]][paste0("sd(", sub("^sd_g__", "", n), ")"),
                             "Estimate"]
  })))
  e
}
frm_vc <- function(fit) VarCorr(fit)$g$sd[, 1]

## Gaussian, (1 + f || g)
ff <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
lm4 <- lmer(y ~ f + (1 | g) + (0 + fb | g) + (0 + fc | g), data = dd,
            REML = FALSE,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))
vc4 <- as.data.frame(VarCorr(lm4))
sd4 <- vc4$sdcor[match(c("g", "g.1", "g.2"), vc4$grp)]
show("gauss lme4 dummies", as.numeric(logLik(ff)), as.numeric(logLik(lm4)),
     frm_vc(ff), sd4, fixef(ff), fixef(lm4))
gt <- glmmTMB(y ~ f + diag(1 + f | g), data = dd, REML = FALSE)
show("gauss glmmTMB diag", as.numeric(logLik(ff)), as.numeric(logLik(gt)),
     frm_vc(ff), attr(VarCorr(gt)$cond$g, "stddev"), fixef(ff),
     fixef(gt)$cond)
cat("sigma frm", sigma(ff), "lme4", sigma(lm4), "glmmTMB", sigma(gt), "\n")
cat("frm ranef vs lme4 ranef max rel:",
    rel(as.vector(ranef(ff)$g[, "Intercept"]),
        ranef(lm4)$g[levels(dd$g), "(Intercept)"]), "\n")

## Gaussian, (1 + z + f || g)
fm <- frm(bf(y ~ z + f + (1 + z + f || g)), family = gaussian(), data = dd)
gm <- glmmTMB(y ~ z + f + diag(1 + z + f | g), data = dd, REML = FALSE)
show("mix glmmTMB diag", as.numeric(logLik(fm)), as.numeric(logLik(gm)),
     frm_vc(fm), attr(VarCorr(gm)$cond$g, "stddev"), fixef(fm),
     fixef(gm)$cond)

## Poisson, (1 + f || g), Laplace in all three
fp <- frm(bf(xo ~ f + (1 + f || g)), family = poisson(), data = dd)
gp <- glmmTMB(xo ~ f + diag(1 + f | g), family = poisson(), data = dd)
gl <- glmer(xo ~ f + (1 | g) + (0 + fb | g) + (0 + fc | g),
            family = poisson(), data = dd,
            control = glmerControl(optimizer = "bobyqa"))
vcl <- as.data.frame(VarCorr(gl))
show("pois glmmTMB diag", as.numeric(logLik(fp)), as.numeric(logLik(gp)),
     frm_vc(fp), attr(VarCorr(gp)$cond$g, "stddev"), fixef(fp),
     fixef(gp)$cond)
show("pois lme4 dummies", as.numeric(logLik(fp)), as.numeric(logLik(gl)),
     frm_vc(fp), vcl$sdcor[match(c("g", "g.1", "g.2"), vcl$grp)], fixef(fp),
     fixef(gl))

## (0 + f || g): a column per level, against glmmTMB
f0 <- frm(bf(y ~ f + (0 + f || g)), family = gaussian(), data = dd)
g0 <- glmmTMB(y ~ f + diag(0 + f | g), data = dd, REML = FALSE)
show("f0 glmmTMB diag", as.numeric(logLik(f0)), as.numeric(logLik(g0)),
     frm_vc(f0), attr(VarCorr(g0)$cond$g, "stddev"), fixef(f0),
     fixef(g0)$cond)

## numeric (1 + x || g) against lme4's own ||
fn <- frm(bf(y ~ x + (1 + x || g)), family = gaussian(), data = dd)
ln <- lmer(y ~ x + (1 + x || g), data = dd, REML = FALSE)
vcn <- as.data.frame(VarCorr(ln))
show("num lme4 ||", as.numeric(logLik(fn)), as.numeric(logLik(ln)),
     frm_vc(fn), vcn$sdcor[1:2], fixef(fn), fixef(ln))
sink()
