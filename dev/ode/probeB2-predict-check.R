# Probe B follow-up: what predict()/fitted() actually return on the
# custom-family route (predict's response scale looked wrong).
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
library(frmtmb)
library(RTMBode)
source("dev/ode/pk-common.R")
d <- sim_pk(); d$idn <- as.integer(d$id)
source("dev/ode/probeB-family-def.R")
form <- bf(conc | vreal(time, dose) + vint(idn) ~ 1 + (1 | id),
           lke ~ 1 + (1 | id), lV ~ 1, lsigma ~ 1)
fit <- frm(form + pk_family, data = d,
           start = list(beta = 0, betad = c(log(0.25), log(8), 0)))
p <- as.numeric(predict(fit)); f <- as.numeric(fitted(fit))
cat("head predict :", format(head(p, 4), digits = 4), "\n")
cat("head fitted  :", format(head(f, 4), digits = 4), "\n")
cat("head truth   :", format(head(d$mu_true, 4), digits = 4), "\n")
cat("predict == lka linear predictor:",
    isTRUE(all.equal(p, as.numeric(predict(fit, dpar = "lka")))), "\n")
cat("max|fitted - truth| :", format(max(abs(f - d$mu_true)), digits = 3), "\n")
cat("residuals use mean_fn:",
    isTRUE(all.equal(as.numeric(residuals(fit)),
                     as.numeric(d$conc - f), tolerance = 1e-8)), "\n")
# predict()'s default type is "link" (the primary dpar's linear
# predictor): on this family that is lka, not the concentration. The
# response mean needs type = "response".
pr <- as.numeric(predict(fit, type = "response"))
cat("predict(type='response') == fitted:", isTRUE(all.equal(pr, f)), "\n")
nd <- d[1:16, ]
cat("predict(newdata, type='response'):",
    format(head(as.numeric(predict(fit, newdata = nd, type = "response")), 3),
           digits = 4), "\n")
cat("  (truth:", format(head(d$mu_true, 3), digits = 4), ")\n")
# a newdata whose ROW SET regroups the subjects differently: the family's
# solve is driven by the vint() column, so a newdata slice that keeps only
# part of a subject's times still solves that subject alone
nd2 <- d[d$idn <= 2, ]
cat("newdata subset of 2 subjects matches the in-sample rows:",
    isTRUE(all.equal(as.numeric(predict(fit, newdata = nd2,
                                        type = "response")),
                     f[d$idn <= 2], tolerance = 1e-6)), "\n")
