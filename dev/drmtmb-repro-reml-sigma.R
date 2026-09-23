# Minimal reproduction for a possible drmTMB issue, drmTMB only. NOT
# filed. Under REML = TRUE, adding a sigma random intercept changes
# which parameters are integrated out (beta_sigma joins beta_mu), so a
# model and its sd -> 0 limit report restricted likelihoods on
# different criteria.
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(drmTMB)
set.seed(1)
ng <- 30
n <- ng * 8
g <- factor(rep(seq_len(ng), each = 8))
x <- rnorm(n)
z <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(ng, 0, 0.6)[g] + rnorm(n, 0, exp(-0.2 + 0.3 * z))
d <- data.frame(y, x, z, g)
show <- function(label, value, digits = 10) {
  cat(formatC(label, width = -44), ":", format(value, digits = digits), "\n")
}
random_blocks <- function(fit) {
  paste(unique(names(fit$obj$env$par)[fit$obj$env$random]), collapse = " ")
}

fA <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z), data = d, REML = TRUE)
fB <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d,
             REML = TRUE)
show("REML logLik, sigma ~ z", as.numeric(logLik(fA)))
show("REML logLik, sigma ~ z + (1 | g)", as.numeric(logLik(fB)))
show("fitted sd of the sigma effect", exp(fB$opt$par[["log_sd_sigma"]]), 4)
cat("integrated (random) blocks, A:", random_blocks(fA), "\n")
cat("integrated (random) blocks, B:", random_blocks(fB), "\n")
show("logLik(B) - logLik(A)",
     as.numeric(logLik(fB)) - as.numeric(logLik(fA)), 8)

# ML nests as expected.
mA <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z), data = d)
mB <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d)
show("ML logLik(B) - logLik(A)",
     as.numeric(logLik(mB)) - as.numeric(logLik(mA)), 4)

# Exact attribution: rebuild the ML template of model A (no sigma
# effect) with beta_sigma marked random as well as beta_mu, and evaluate
# it at B's fitted log_sd_mu. B's criterion, with its sigma SD at 1.6e-5,
# is exactly this joint Laplace integral, so the drop is the change of
# criterion and not a change of model.
e <- mA$obj$env
objJ <- TMB::MakeADFun(data = e$data, parameters = e$parList(),
                       map = e$map,
                       random = c("beta_mu", "beta_sigma", "u_mu"),
                       DLL = e$DLL, silent = TRUE)
cat("rebuilt outer parameters:", paste(names(objJ$par), collapse = " "),
    "\n")
show("joint-Laplace criterion of A at B's log_sd_mu",
     -objJ$fn(fB$opt$par[["log_sd_mu"]]))
show("B's reported REML logLik", as.numeric(logLik(fB)))
