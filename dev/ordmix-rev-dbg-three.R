# Reviewer of lane ordmix: the r_three_mixlink case of
# dev/ordmix-rev-lpcheck.R (data seed 20261005 + 31), where frmtmb's
# gradient at its own optimum is Inf and its density NaN at a point where
# brms's is finite. Which component, and which rows?
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cases <- 31L
set.seed(20261005 + cases)
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
fit <- frm(bf(y ~ x), family = mixture(cumulative("probit"),
                                       sratio("cloglog"), acat()),
           data = d, control = frmtmb_control(grad_tol = 1e-8))
cat("logLik", format(as.numeric(logLik(fit)), digits = 15), "\n")
cat("opt message:", fit$opt$message, " conv:", fit$opt$convergence, "\n")
gr <- fit$obj$gr(fit$opt$par)
names(gr) <- names(fit$opt$par)
print(gr)
print(fit$estimates[c("tau_raw1", "tau_raw2", "tau_raw3")])
print(fixef(fit))
# per-component log densities at the optimum, rows that are not finite
fam <- fit$spec$responses[[1]]$family
# Which rows: the gradient at the fit's own parameters, over the model
# refitted on each single row's complement and on each row alone
p0 <- fit$opt$par
st <- fit$estimates
bad <- integer(0)
for (i in seq_len(n)) {
  fi <- suppressWarnings(frm(bf(y | thres(3) ~ x),
             family = mixture(cumulative("probit"), sratio("cloglog"),
                              acat()),
             data = d[i, , drop = FALSE], start = st,
             control = frmtmb_control(optCtrl = list(iter.max = 1,
                                                     eval.max = 1))))
  gi <- fi$obj$gr(p0)
  if (any(!is.finite(gi))) bad <- c(bad, i)
}
cat("rows with a non-finite gradient alone:", bad, "\n")
print(d[bad, ])
eta1 <- fixef(fit)["mu1_x", "Estimate"] * d$x[bad]
tau1 <- fixef(fit)[1:3, "Estimate"]
cat("component 1 tau - eta at those rows:\n")
print(outer(-eta1, tau1, "+"))
