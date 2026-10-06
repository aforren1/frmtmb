# Reviewer of lane fixes, claim 1, follow-ups to dev/fixes-rev-emm.R.
#   Rscript dev/fixes-rev-emm2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans); library(splines)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(7101)
n <- 150
d <- data.frame(x = rnorm(n), z = runif(n, 0.5, 4), w = runif(n, 1, 3),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$yc <- rpois(n, d$w * exp(0.2 + 0.3 * log(d$z) + 0.2 * (d$f == "b")))
P <- poisson()
tr <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) {
    paste("ERROR:", conditionMessage(e), "| call:",
          paste(deparse(conditionCall(e)), collapse = " "))
  })
  cat(sprintf("%-40s %s\n", lab, paste(r, collapse = " ")))
}
fo <- yc ~ ns(z, 3) + f + offset(log(w))
fit <- frm(bf(fo), data = d, family = P)
ref <- glm(fo, data = d, family = P)
tr("glm ns+offset response", format(summary(emmeans(ref, "f",
                                                   type = "response"))$rate))
tr("frm ns+offset response", format(summary(emmeans(fit, "f",
                                                   type = "response"))$rate))
tr("frm ns+offset link", format(summary(emmeans(fit, "f"))$emmean))
tr("glm ns+offset link", format(summary(emmeans(ref, "f"))$emmean))
fo2 <- yc ~ poly(z, 2) + f + offset(log(w))
fit2 <- frm(bf(fo2), data = d, family = P)
tr("frm poly+offset response", format(summary(emmeans(fit2, "f",
                                                     type = "response"))$rate))
fo3 <- yc ~ scale(z) + f + offset(log(w))
fit3 <- frm(bf(fo3), data = d, family = P)
tr("frm scale+offset response", format(summary(emmeans(fit3, "f",
                                                      type = "response"))$rate))
fo4 <- yc ~ ns(z, 3) + f
fit4 <- frm(bf(fo4), data = d, family = P)
tr("frm ns no offset response", format(summary(emmeans(fit4, "f",
                                                      type = "response"))$rate))
fo5 <- yc ~ z + f + offset(log(w))
fit5 <- frm(bf(fo5), data = d, family = P)
tr("frm z+offset response", format(summary(emmeans(fit5, "f",
                                                  type = "response"))$rate))
tr("frm ns+offset regrid", format(summary(regrid(ref_grid(fit)))$rate))
tr("frm ns+offset epred", format(summary(emmeans(fit, "f",
                                                epred = TRUE))$emmean))
tr("glm ns+offset response(regrid)", format(summary(regrid(ref_grid(ref))))$rate)
# NA rows: ns() accepts NA, poly() does not (model.frame evaluates it
# before na.omit, for glm too)
dd <- d; dd$z[c(3, 17)] <- NA; dd$f[40] <- NA
fitn <- frm(bf(yc ~ ns(z, 3) + f), data = dd, family = P)
refn <- glm(yc ~ ns(z, 3) + f, data = dd, family = P)
tr("NA rows ns", {
  a <- summary(emmeans(fitn, "f")); b <- summary(emmeans(refn, "f"))
  sprintf("maxrel %.3g gridz frm %.8f glm %.8f",
          max(abs(a$emmean - b$emmean) / abs(b$emmean)),
          unique(ref_grid(fitn)@grid$z), unique(ref_grid(refn)@grid$z))
})
tr("glm poly NA", {glm(yc ~ poly(z, 2) + f, data = dd, family = P); "ok"})
