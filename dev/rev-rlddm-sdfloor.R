.libPaths(c("C:/Users/adf44/source/r/rev-rlddm-lib","C:/Users/adf44/source/r/pinlib","C:/Users/adf44/source/r/rellib-0552","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb);library(frmtmb.eam);library(frmtmb.learn)})
truth <- list(alpha=0.35, drift=2.5, bs=1.5, ndt=0.25, sd_alpha=0.5, sd_drift=1.0, sd_bs=0.2, sd_ndt=0.15)
d <- frm_task_design("bandit2arm", n_subject=100L, n_trial=200L, seed=20260908L)
set.seed(20260908L+2L); i <- as.integer(d$id); ns <- 100L
ua <- rnorm(ns,0,truth$sd_alpha); ud <- rnorm(ns,0,truth$sd_drift); ub <- rnorm(ns,0,truth$sd_bs); un <- rnorm(ns,0,truth$sd_ndt)
s <- frm_task_simulate(rlddm(subject=id, trial=trial), d, pars=list(alpha=plogis(qlogis(truth$alpha)+ua[i]), drift=truth$drift+ud[i], bs=truth$bs*exp(ub[i]), ndt=truth$ndt*exp(un[i]), bias=0.5), seed=20260908L)[[1L]]
lv <- levels(s$id); own <- as.numeric(tapply(s$rt, s$id, min))[match(lv, lv)]; tru <- truth$ndt*exp(un)
k <- sum(own*tru)/sum(own*own)
cat("sd(own floors)        :", format(sd(own), digits=9), "\n")
cat("sd(k * own), k =", format(k, digits=6), ":", format(k*sd(own), digits=9), "\n")
cat("sd(0.828719 * own)    :", format(0.828719*sd(own), digits=9), "\n")
cat("sd(truths)            :", format(sd(tru), digits=9), "\n")
cat("the fit reports        : 0.0333713\n")