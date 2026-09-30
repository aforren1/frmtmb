# The equated mixture through frmtmb.sample, and against a hand-written
# RTMB objective with one shared sigma. Data seed 11; sampler seed 3.
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
library(frmtmb)
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
fam <- mixture(gaussian(), gaussian())
fit <- frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)

# hand-written: one log sigma shared by both components; theta2 is the
# reference, as in frmtmb (mixture_theta_reference())
X <- cbind(1, d$x)
y <- d$y
nll <- function(p) {
  mu1 <- as.vector(X %*% p$b1)
  mu2 <- as.vector(X %*% p$b2)
  s <- exp(p$ls)
  lw1 <- p$t1 - log(1 + exp(p$t1))
  lw2 <- -log(1 + exp(p$t1))
  a <- lw1 + RTMB::dnorm(y, mu1, s, log = TRUE)
  b <- lw2 + RTMB::dnorm(y, mu2, s, log = TRUE)
  m <- 0.5 * (a + b + abs(a - b))
  -sum(m + log(exp(a - m) + exp(b - m)))
}
est <- fit$estimates
p_frm <- list(b1 = unname(est$beta[1:2]), b2 = unname(est$beta[3:4]),
              ls = unname(est$betad[1]), t1 = unname(est$betad[2]))
obj <- RTMB::MakeADFun(nll, p_frm, silent = TRUE)
cat(sprintf("HAND nll at frm estimates %.12f, -logLik(fit) %.12f\n",
            obj$fn(obj$par), -as.numeric(logLik(fit))))
cat(sprintf("HAND rel diff at same point %.3e\n",
            (obj$fn(obj$par) + as.numeric(logLik(fit))) /
              abs(as.numeric(logLik(fit)))))
st <- obj$par
st[] <- 0
st[c(1, 3)] <- c(-1, 3)
opt <- nlminb(st, obj$fn, obj$gr)
cat(sprintf("HAND optimum %.12f, frm optimum %.12f, rel %.3e\n",
            -opt$objective, as.numeric(logLik(fit)),
            (-opt$objective - as.numeric(logLik(fit))) /
              abs(as.numeric(logLik(fit)))))
se <- sqrt(diag(solve(obj$he(opt$par))))
cat("HAND max |est diff| / se:",
    format(max(abs(opt$par - obj$par) / se)), "\n")

if (length(args) && args[1] == "sample") {
  library(frmtmb.sample)
  dr <- frm_sample(fit, chains = 1, iter = 400, seed = 3, refresh = 0)
  print(variables(dr))
  m <- posterior::as_draws_matrix(dr)
  print(colnames(m))
  if (all(c("sigma1", "sigma2") %in% colnames(m))) {
    cat("DRAWS sigma1 identical sigma2:",
        identical(unname(m[, "sigma1"]), unname(m[, "sigma2"])), "\n")
  }
  print(summary(dr))
}
