# Reviewer of lane optstart: four fits of the lane's fit diff, against
# references that do not use frmtmb.
#   Rscript dev/optstart-rev-verify.R <lib> <label>
# Worse on the lane: test-v14.R's zero_inflated_asym_laplace (exact: rq
# median regression on the nonzero rows, sigma from it, glm for zi) and
# frmtmb.sample test-brms-pins.R's (1 + x + z + w | g) (lme4, ML).
# Newly nonzero on the lane: test-sparsex.R's y ~ f + x + z + (1 | g)
# (lme4, ML) and test-optrap.R's negative binomial at the Poisson limit
# (the Poisson log likelihood is the supremum; MASS::glm.nb).
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
cap <- function(expr) {
  w <- character()
  f <- withCallingHandlers(suppressMessages(expr), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  cat(sprintf("  frmtmb logLik %.10f code %d (%s); %d warnings%s\n",
              as.numeric(logLik(f)), f$opt$convergence, f$opt$message,
              length(w), if (length(w)) paste0(": ", paste(substr(w, 1, 70),
                                                         collapse = " | "))
              else ""))
  f
}
# W1
cat("W1 test-v14.R zero_inflated_asym_laplace\n")
set.seed(14)
n <- 600
x <- rnorm(n)
y_ald <- 0.5 + 0.8 * x + rnorm(n)
zi_true <- stats::rbinom(n, 1, plogis(-1 + 0.5 * x))
dd <- data.frame(x = x, y = ifelse(zi_true == 1, 0, y_ald))
fz <- cap(frm(bf(y ~ x, zi ~ x, quantile = 0.5) +
                zero_inflated_asym_laplace(), data = dd))
nz <- dd$y != 0
rq <- quantreg::rq(y ~ x, tau = 0.5, data = dd[nz, ])
r <- residuals(rq)
s <- 0.5 * mean(abs(r))
ll_ald <- sum(log(0.25) - log(s) - 0.5 * abs(r / s))
gz <- glm(I(y == 0) ~ x, family = binomial(), data = dd)
cat(sprintf("  exact: %.10f (rq + glm)\n", ll_ald + as.numeric(logLik(gz))))
# W2
cat("W2 test-brms-pins.R (1 + x + z + w | g)\n")
set.seed(5)
G <- 10; n <- 300
d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                w = stats::rnorm(n), g = factor(rep(1:G, 30)))
S <- matrix(c(1, .5, -.3, .2, .5, 1, .3, -.2, -.3, .3, 1, .1,
              .2, -.2, .1, 1), 4)
L <- t(chol(S * 0.4))
U <- t(L %*% matrix(stats::rnorm(4 * G), 4))
d$y <- 1 + d$x + U[d$g, 1] + U[d$g, 2] * d$x + U[d$g, 3] * d$z +
  U[d$g, 4] * d$w + stats::rnorm(n)
fp <- cap(frm(bf(y ~ x + (1 + x + z + w | g)), family = gaussian(), data = d))
best <- -Inf
for (opt in c("bobyqa", "Nelder_Mead")) {
  l4 <- suppressWarnings(suppressMessages(lmer(y ~ x + (1 + x + z + w | g),
    data = d, REML = FALSE, control = lmerControl(optimizer = opt,
      optCtrl = list(maxfun = 1e5)))))
  cat(sprintf("  lme4 %s: %.10f (singular %s)\n", opt,
              as.numeric(logLik(l4)), isSingular(l4)))
  best <- max(best, as.numeric(logLik(l4)))
}
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  gt <- suppressWarnings(glmmTMB::glmmTMB(y ~ x + (1 + x + z + w | g),
                                          data = d))
  cat(sprintf("  glmmTMB: %.10f\n", as.numeric(logLik(gt))))
}
# N1
cat("N1 test-sparsex.R y ~ f + x + z + (1 | g)\n")
set.seed(7)
n <- 4000
ds <- data.frame(f = factor(rep(sprintf("s%03d", 1:200), each = 20)),
                 x = rnorm(n), z = rnorm(n),
                 g = factor(sample(1:25, n, TRUE)))
fe <- rnorm(200, 0, 0.8)
ds$y <- rnorm(n, 1 + fe[as.integer(ds$f)] + 0.5 * ds$x - 0.3 * ds$z +
                rnorm(25, 0, 0.4)[ds$g], 1)
form <- bf(y ~ f + x + z + (1 | g)) + gaussian()
ctl <- function(sparse) {
  frmtmb_control(sparse_x = sparse, restarts = 2,
                 optCtrl = list(iter.max = 2000, eval.max = 2000,
                                rel.tol = 1e-12))
}
f_d <- cap(frm(form, data = ds, control = ctl(FALSE)))
f_s <- cap(frm(form, data = ds, control = ctl(TRUE)))
f_def <- cap(frm(form, data = ds))
l4 <- lmer(y ~ f + x + z + (1 | g), data = ds, REML = FALSE)
cat(sprintf("  lme4: %.10f; group sd %.4f (singular %s); frmtmb group sd %s\n",
            as.numeric(logLik(l4)), attr(VarCorr(l4)$g, "stddev"),
            isSingular(l4), format(exp(f_d$opt$par[grep("theta",
                                                       names(f_d$opt$par))]))))
# N2
cat("N2 test-optrap.R negbinomial at the Poisson limit\n")
set.seed(2)
invisible(runif(200)); invisible(rnorm(20)); invisible(rnorm(200))
set.seed(2)
n <- 150
x <- rnorm(n)
invisible(rnorm(n))
mu <- stats::plogis(-2 + 0.5 * x)
invisible(rbeta(n, mu * 200, (1 - mu) * 200)); invisible(sample(n, 3))
invisible(rbeta(n, 0.98 * 20, 0.02 * 20)); invisible(sample(n, 2))
invisible(rt(n, df = 30))
dn <- data.frame(x = x, c = rpois(n, exp(1 + 0.3 * x)))
fn <- cap(frm(c ~ x, family = negbinomial(), data = dn))
cat(sprintf("  shape at the stop: %.4g\n",
            exp(fn$opt$par[grep("shape", names(fn$opt$par))])))
cat(sprintf("  Poisson glm (the supremum): %.10f\n",
            as.numeric(logLik(glm(c ~ x, family = poisson(), data = dn)))))
nb <- tryCatch(suppressWarnings(MASS::glm.nb(c ~ x, data = dn)),
               error = function(e) NULL)
if (!is.null(nb)) {
  cat(sprintf("  MASS::glm.nb: %.10f, theta %.4g\n",
              as.numeric(logLik(nb)), nb$theta))
}
