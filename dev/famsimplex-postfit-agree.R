# fitted() and predict() of the ML fits against brms 2.23.0's on the
# same data: the category means and their standard errors against the
# posterior means and sds of posterior_epred(), and the predictive sds
# against brms's predict(). Fits: dev/famsimplex-brms-fits2.R (dir,
# dirre, dir2 for fitted; all five for predict). frmtmb's predict()
# draws 4000 replicates (seed 1), brms's uses its 4000 draws.
# Output: dev/famsimplex-postfit-agree.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
loadNamespace("brms")
cat("lib:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
forms <- list(dir = bf(Y ~ x, phi ~ z), dirre = bf(Y ~ x + (1 | ID | g)),
              ln = bf(Y ~ x, sigmab ~ z), lnre = bf(Y ~ x + (1 | ID | g)),
              dir2 = bf(Y ~ x))
fams <- list(dir = dirichlet(), dirre = dirichlet(),
             ln = logistic_normal(), lnre = logistic_normal(),
             dir2 = brmsfamily("dirichlet2"))
q <- function(v) sprintf("%.3f", stats::quantile(v, c(0.05, 0.5, 0.95)))
for (nm in names(forms)) {
  bfit <- readRDS(file.path(wt, paste0("famsimplex-brmsfit-", nm, ".rds")))
  fit <- frm(forms[[nm]], family = fams[[nm]], data = bfit$data)
  cat("\n==", nm, "\n")
  if (!startsWith(nm, "ln")) {
    fb <- stats::fitted(bfit)
    ff <- fitted(fit)
    zm <- (fb[, "Estimate", ] - ff[, "Estimate", ]) / fb[, "Est.Error", ]
    rs <- ff[, "Est.Error", ] / fb[, "Est.Error", ]
    cat("fitted: (brms mean - ML) / brms sd, 5/50/95%:", q(zm), "\n")
    cat("fitted: ML se / brms sd, 5/50/95%:", q(rs), "\n")
  }
  pb <- stats::predict(bfit)
  set.seed(1)
  pf <- predict(fit, ndraws = 4000)
  zp <- (pb[, "Estimate", ] - pf[, "Estimate", ]) / pb[, "Est.Error", ]
  rp <- pf[, "Est.Error", ] / pb[, "Est.Error", ]
  cat("predict: (brms mean - ML mean) / brms sd, 5/50/95%:", q(zp), "\n")
  cat("predict: ML sd / brms sd, 5/50/95%:", q(rp), "\n")
  cat("dims identical:", identical(dim(pb), dim(pf)),
      " layer names identical:", identical(dimnames(pb)[[3]],
                                           dimnames(pf)[[3]]), "\n")
}
