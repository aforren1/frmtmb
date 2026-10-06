# Reviewer of lane ordmix: the probit's log-odds form saturates at
# |tau - eta| > 38.2 (R/links.R), giving NaN density and an Inf gradient.
# (1) A plain cumulative probit fit evaluated at a far point, on the arm
#     given (lane or base): is the NaN pre-existing?
# (2) arm lane only: the r_three_mixlink model of dev/ordmix-rev-lpcheck.R
#     (data seed 20261005 + 31) from 12 starts: is its reported optimum,
#     where the gradient is Inf, the maximum?
# Usage: Rscript dev/ordmix-rev-probit-sat.R <lane|base>
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
set.seed(20261005 + 31)
n <- 400
x <- rnorm(n); z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
yh <- ifelse(runif(n) < 0.2, 0L, y)
w <- runif(n, 0.5, 2)
cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
d <- data.frame(y, x, z)

## (1) plain cumulative probit
f1 <- frm(bf(y ~ x), family = cumulative("probit"), data = d)
p <- f1$opt$par
for (b in c(5, 10, 20, 40)) {
  q <- p
  q[names(q) == "beta"] <- b
  cat(sprintf("plain probit slope=%g max|tau-eta|=%.1f fn=%s maxgrad=%s\n",
              b, max(abs(outer(b * d$x, c(-1, 1) * 10, "+"))),
              format(f1$obj$fn(q)), format(max(abs(f1$obj$gr(q))))))
}
if (arm != "lane") quit(save = "no")

## (2) the three-component mixture from several starts
fam <- function() mixture(cumulative("probit"), sratio("cloglog"), acat())
fit0 <- frm(bf(y ~ x), family = fam(), data = d,
            control = frmtmb_control(grad_tol = 1e-8))
cat(sprintf("default start: logLik=%.6f maxgrad=%s conv=%d msg=%s\n",
            as.numeric(logLik(fit0)),
            format(max(abs(fit0$obj$gr(fit0$opt$par)))),
            fit0$opt$convergence, fit0$opt$message))
st0 <- fit0$frame[["par_template"]]
set.seed(7)
res <- NULL
for (s in 1:12) {
  st <- lapply(st0, function(v) v + rnorm(length(v), 0, 0.7))
  fs <- tryCatch(suppressWarnings(frm(bf(y ~ x), family = fam(), data = d,
                                      start = st)),
                 error = function(e) NULL)
  if (is.null(fs)) {
    cat("start", s, "error\n")
    next
  }
  gm <- max(abs(fs$obj$gr(fs$opt$par)))
  cat(sprintf("start %2d logLik=%.6f maxgrad=%s conv=%d\n", s,
              as.numeric(logLik(fs)), format(gm), fs$opt$convergence))
  res <- rbind(res, data.frame(s = s, ll = as.numeric(logLik(fs)), gm = gm))
}
cat(sprintf("best of 12 starts %.6f (finite gradient: %s); default %.6f\n",
            max(res$ll), is.finite(res$gm[which.max(res$ll)]),
            as.numeric(logLik(fit0))))
## the same data with the logit in place of the probit
fl <- frm(bf(y ~ x), family = mixture(cumulative(), sratio("cloglog"),
                                      acat()), data = d)
cat(sprintf("logit in component 1: logLik=%.6f maxgrad=%s\n",
            as.numeric(logLik(fl)), format(max(abs(fl$obj$gr(fl$opt$par))))))
