# Which of the model features brms's six fixture fits use does frmtmb
# accept? This bounds how much of tests.brmsfit-methods.R can be ported
# at all, so the 2.6b estimate needs it. dry_run = "frame" stops before
# the fit, so the answer is about the grammar and the frame, not about
# convergence.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
load("dev/brms-suite/brms/data/epilepsy.rda")
d <- epilepsy
d$AgeSD <- abs(rnorm(nrow(d))) + 0.1
d$w <- 1
d$vol <- rnorm(nrow(d))
d$ord <- factor(rep(1:4, length.out = nrow(d)), ordered = TRUE)

probe <- function(label, expr) {
  r <- tryCatch({expr; "ACCEPTS"},
                error = function(e) paste("refuses:",
                                          substr(conditionMessage(e), 1, 90)))
  cat(sprintf("%-34s %s\n", label, r))
}
F <- function(form, ...) frm(form, d, ..., dry_run = "frame")

probe("mo(Exp)",            F(count ~ mo(Base), family = poisson()))
probe("s(Age)",             F(count ~ s(Age), family = poisson()))
probe("arma(visit, patient)", F(count ~ Trt + arma(visit, patient),
                                family = poisson()))
probe("arma(cov = TRUE)",   F(count ~ Trt + arma(visit, patient, cov = TRUE)))
probe("me(Age, AgeSD)",     F(count ~ me(Age, AgeSD), family = poisson()))
probe("mm(patient, visit)", F(count ~ (1 | mm(patient, visit)),
                              family = poisson()))
probe("cs(Trt)",            F(ord ~ cs(Trt), family = sratio()))
probe("gp(Age)",            F(count ~ gp(Age), family = poisson()))
probe("| weights(w)",       F(count | weights(w) ~ Trt, family = poisson()))
probe("| se(AgeSD)",        F(Age | se(AgeSD) ~ Trt))
probe("| trunc(lb = 0)",    F(count | trunc(lb = 0) ~ Trt, family = poisson()))
probe("| cens(Trt)",        F(count | cens(Trt) ~ zBase, family = poisson()))
probe("| thres(3)",         F(ord | thres(3) ~ Trt, family = cumulative()))
probe("nl = TRUE",          F(bf(count ~ a * Trt + b, a + b ~ 1, nl = TRUE)))
probe("mvbind(y1, y2)",     F(mvbf(bf(count ~ Trt), bf(Age ~ Trt))))
probe("categorical()",      F(ord ~ Trt, family = categorical()))
probe("mixture(gaussian, 2)",
      F(Age ~ Trt, family = mixture(gaussian, gaussian)))
probe("gr(patient, by = Trt)", F(count ~ (1 | gr(patient, by = Trt)),
                                 family = poisson()))
probe("(1 + Trt | visit)",  F(count ~ Trt + (1 + Trt | visit),
                              family = poisson()))
probe("(1 | ID | patient)", F(bf(count ~ Trt + (1 | q | patient),
                                 sigma ~ (1 | q | patient))))
probe("offset(Age)",        F(count ~ Trt + offset(Age), family = poisson()))
probe("t2(Age)",            F(count ~ t2(Age), family = poisson()))
probe("mmc(Age, vol)",      F(count ~ (1 + mmc(Age, vol) | mm(patient, visit)),
                             family = poisson()))
probe("unstr(visit, patient)", F(count ~ Trt + unstr(visit, patient),
                                 family = poisson()))
