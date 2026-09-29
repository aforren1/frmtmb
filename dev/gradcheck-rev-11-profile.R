## Reviewer, claim 5: confint(method = "profile") and the ordinary Wald
## confint(), on two models, saved for an identical() comparison between
## the builds. R/confint.R was edited, so the whole file's public
## behaviour is in scope, not only diagnose().
## usage: Rscript gradcheck-rev-11-profile.R <core-lib> <out.rds>
a <- commandArgs(TRUE)
LIB <- a[1]; OUT <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n")
out <- list()

set.seed(6001)
n <- 600
d1 <- data.frame(x = rnorm(n), z = rnorm(n))
d1$y <- rnorm(n, 1 + 2 * d1$x - 0.5 * d1$z, 1)
f1 <- frm(bf(y ~ x + z), family = gaussian(), data = d1)
out$glm_wald <- confint(f1)
out$glm_profile <- suppressWarnings(confint(f1, parm = c("x", "z"), method = "profile"))
out$glm_vcov <- vcov(f1)

set.seed(6002)
ng <- 40
d2 <- data.frame(g = factor(rep(seq_len(ng), 25)))
d2$x <- rnorm(nrow(d2))
re <- rnorm(ng, 0, 0.7)
d2$y <- rpois(nrow(d2), exp(0.3 + 0.5 * d2$x + re[d2$g]))
f2 <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = poisson(),
                           data = d2))
out$glmm_wald <- confint(f2)
out$glmm_profile <- suppressWarnings(confint(f2, parm = "x", method = "profile"))
out$glmm_vcov <- vcov(f2)

## a bounded fit, where the new bound handling could plausibly leak into
## the covariance path
set.seed(6003)
n <- 300
d3 <- data.frame(x = rnorm(n))
d3$y <- rnorm(n, 1 + 2 * d3$x, 1)
f3 <- suppressWarnings(frm(bf(y ~ x), family = gaussian(), data = d3,
                           prior = set_prior("", class = "b", ub = 0.1)))
out$bounded_wald <- suppressWarnings(confint(f3))
out$bounded_vcov <- vcov(f3)
out$bounded_sdr_pd <- isTRUE(frmtmb:::sdr_of(f3)$pdHess)

## diagnose()'s pre-existing fields, which must not have moved
dg <- function(f) {
  d <- suppressWarnings(diagnose(f, quiet = TRUE))
  d[intersect(names(d),
              c("convergence", "message", "max_grad", "worst_grad",
                "pdHess", "bad_se", "min_eig", "flat", "singular",
                "separation", "unbounded_dpar", "predictor_scale",
                "nonfinite_trials", "scale"))]
}
out$dg1 <- dg(f1); out$dg2 <- dg(f2); out$dg3 <- dg(f3)
saveRDS(out, OUT)
for (nm in names(out)) {
  cat("--", nm, "\n")
  print(out[[nm]])
}
cat("DONE profile\n")
