# Lane ndaterm, fix round 1: which dpars a 2-component gaussian
# mixture has, and the reporting against the natural scale of theta.
# Seed 7. Usage: Rscript dev/ndaterm-r1-theta.R <arm>
arm <- commandArgs(TRUE)[1]
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(7)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.2)
d$y <- ifelse(k == 1, 6, 0) + d$x + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
rs <- fit$spec$responses[[1]]
cat("dpars:", names(rs$dpars), "\n")
nd <- data.frame(x = 0)
for (dn in grep("theta", names(rs$dpars), value = TRUE)) {
  cat(dn, "response", frm_linpred(fit, newdata = nd, dpar = dn,
                                  type = "response"),
      " link", frm_linpred(fit, newdata = nd, dpar = dn, type = "link"),
      "\n")
}
set.seed(1)
p <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)), summary = FALSE,
             ndraws = 4000)
cat("predict means", colMeans(p), "\nfitted",
    fitted(fit, newdata = data.frame(x = c(-1, 0, 1)))[, 1], "\n")
