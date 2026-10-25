# Reviewer: a residual AR(1) term on a fixed ordering, written from the
# help page of frmtmb_register_autocor() alone, then checked against a
# dense mvtnorm density, against frmtmb's own ar(p = 1, cov = TRUE), and
# through fitted/residuals/simulate/predict/autocor_matrix.
# Rscript dev/review-sarfcor-toyar1.R
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(mvtnorm); library(Matrix)})
cat("lib:", find.package("frmtmb"), "\n")
frmtmb_local_registries(environment())  # this script's own frame

ar1_build <- function(spec, rows, n_all, lookup) {
  n <- length(rows)
  S <- sparseMatrix(i = 2:n, j = 1:(n - 1), x = 1, dims = c(n, n))
  first <- c(1, numeric(n - 1))
  list(S = S, first = first, n = n)
}
ar1_whiten <- function(y, mu, sigma, theta, aux) {
  phi <- tanh(theta[1])
  e <- (y - mu) / sigma
  u <- e - phi * as.vector(aux$S %*% e) +
    (sqrt(1 - phi^2) - 1) * e * aux$first
  # core passes sigma per row even when sigma_varies = FALSE (see review)
  lsig <- if (length(sigma) == 1L) length(y) * log(sigma) else sum(log(sigma))
  list(u = u, logjac = 0.5 * log(1 - phi^2) - lsig)
}
ar1_unwhiten <- function(u, mu, sigma, theta, aux) {
  phi <- tanh(theta[1])
  u <- as.numeric(u)
  u[1] <- u[1] / sqrt(1 - phi^2)
  e <- as.numeric(stats::filter(u, phi, method = "recursive"))
  mu + sigma * e
}
frmtmb_register_autocor("toyar", args = "M", build = ar1_build,
                        whiten = ar1_whiten, unwhiten = ar1_unwhiten,
                        npar = 1L,
                        natural = function(theta, aux) c(toyphi = tanh(theta[1])),
                        sigma_varies = FALSE)
set.seed(41)
n <- 80
phi0 <- 0.6
e <- as.numeric(stats::filter(rnorm(n) * c(1 / sqrt(1 - phi0^2), rep(1, n - 1)),
                              phi0, method = "recursive"))
d <- data.frame(x = rnorm(n), t = 1:n)
d$y <- 1 + 0.5 * d$x + 0.8 * e
fit <- frm(y ~ x + toyar(1), data = d)
cf <- fixef(fit)[, "Estimate"]
cat("toy fit: coef", format(cf, digits = 6), " sigma",
    format(sigma(fit), digits = 6), "\n")
print(confint_varcorr(fit))
phi <- tanh(fit$estimates$thetaac)
sg <- sigma(fit)
mu <- as.vector(cbind(1, d$x) %*% cf)
S <- sg^2 / (1 - phi^2) * phi^abs(outer(1:n, 1:n, "-"))
cat(sprintf("logLik toy %.10f dense %.10f diff %.3g\n", logLik(fit),
            dmvnorm(d$y, mu, S, log = TRUE),
            as.numeric(logLik(fit)) - dmvnorm(d$y, mu, S, log = TRUE)))
fb <- frm(y ~ x + ar(time = t, p = 1, cov = TRUE), data = d)
cat(sprintf("logLik built-in ar(cov = TRUE) %.10f diff to toy %.3g\n",
            logLik(fb), as.numeric(logLik(fb) - logLik(fit))))
cat("coef built-in", format(fixef(fb)[, "Estimate"], digits = 6), "\n")
cat(sprintf("fitted = mu: %.3g\n", max(abs(fitted(fit)[, "Estimate"] - mu))))
cat(sprintf("residuals = y - mu: %.3g\n",
            max(abs(residuals(fit)[, "Estimate"] - (d$y - mu)))))
rp <- residuals(fit, type = "pearson")[, "Estimate"]
cat(sprintf("pearson = (y - mu) / marginal sd: %.3g\n",
            max(abs(rp - (d$y - mu) / sqrt(diag(S))))))
am <- autocor_matrix(fit)
cat(sprintf("autocor_matrix = S / sigma^2: %.3g\n",
            max(abs(unname(as.matrix(am)) - S / sg^2))))
set.seed(1)
sims <- simulate(fit, nsim = 4000)
sm <- as.matrix(sims)
if (nrow(sm) != n) sm <- t(sm)
cv <- cov(t(sm))
cat(sprintf("simulate: dim %s; max|mean - mu| %.3f (sd scale %.3f); lag-1 cor sim %.3f vs %.3f; var[1] sim %.3f vs %.3f\n",
            paste(dim(sims), collapse = "x"), max(abs(rowMeans(sm) - mu)),
            sqrt(S[1, 1] / 4000), mean(diag(cov2cor(cv)[-1, -n])), phi,
            cv[1, 1], S[1, 1]))
pp <- tryCatch(predict(fit, summary = FALSE), error = function(e) e)
if (inherits(pp, "error")) {
  cat("predict(summary = FALSE) error:", conditionMessage(pp), "\n")
} else cat("predict draws dim", paste(dim(pp), collapse = "x"), "\n")
pr <- predict(fit)
cat("predict() columns:", colnames(pr), "\n")
cat(sprintf("predict Estimate vs mu: %.3f\n", max(abs(pr[, "Estimate"] - mu))))
cat("variables:", variables(fit), "\n")
ll <- tryCatch(log_lik(fit), error = function(e) conditionMessage(e))
cat("log_lik:", if (is.character(ll)) ll else dim(ll), "\n")
nd <- tryCatch(predict(fit, newdata = d), error = function(e) conditionMessage(e))
cat("predict(newdata = the fitted data):", if (is.character(nd)) nd else "works", "\n")
nd <- tryCatch(fitted(fit, newdata = d), error = function(e) conditionMessage(e))
cat("fitted(newdata = the fitted data):", if (is.character(nd)) nd else "works", "\n")

# a term that returns a wrong shape from whiten: logjac per row
bad_whiten <- function(y, mu, sigma, theta, aux) {
  w <- ar1_whiten(y, mu, sigma, theta, aux)
  w$logjac <- rep(w$logjac / length(y), length(y))
  w
}
frmtmb_register_autocor("badar", args = "M", build = ar1_build,
                        whiten = bad_whiten, unwhiten = ar1_unwhiten,
                        npar = 1L,
                        natural = function(theta, aux) c(badphi = tanh(theta[1])),
                        sigma_varies = FALSE)
r <- tryCatch({f <- frm(y ~ x + badar(1), data = d)
  sprintf("FIT logLik %.6f (correct %.6f)", logLik(f), logLik(fit))},
  error = function(e) paste("ERROR:", conditionMessage(e)))
cat("whiten with a per-row logjac:", r, "\n")
# u of the wrong length
bad2 <- function(y, mu, sigma, theta, aux) {
  w <- ar1_whiten(y, mu, sigma, theta, aux); w$u <- w$u[-1]; w
}
frmtmb_register_autocor("badar2", args = "M", build = ar1_build,
                        whiten = bad2, unwhiten = ar1_unwhiten, npar = 1L,
                        natural = function(theta, aux) c(b2 = tanh(theta[1])),
                        sigma_varies = FALSE)
r <- tryCatch({f <- frm(y ~ x + badar2(1), data = d)
  sprintf("FIT logLik %.6f", logLik(f))},
  error = function(e) paste("ERROR:", conditionMessage(e)))
cat("whiten with u one short:", r, "\n")
# an unwhiten that is NOT the inverse
frmtmb_register_autocor("badar3", args = "M", build = ar1_build,
                        whiten = ar1_whiten,
                        unwhiten = function(u, mu, sigma, theta, aux) mu + sigma * u,
                        npar = 1L,
                        natural = function(theta, aux) c(b3 = tanh(theta[1])),
                        sigma_varies = FALSE)
r <- tryCatch({f <- frm(y ~ x + badar3(1), data = d)
  sprintf("FIT, autocor_matrix[1,2] %.4f (true %.4f)", autocor_matrix(f)[1, 2],
          phi / (1 - phi^2))},
  error = function(e) paste("ERROR:", conditionMessage(e)))
cat("unwhiten not the inverse:", r, "\n")
# registration validation: a correct term written with ... is refused?
r <- tryCatch({frmtmb_register_autocor("dotar", args = "M", build = ar1_build,
                        whiten = function(...) ar1_whiten(...),
                        unwhiten = function(...) ar1_unwhiten(...), npar = 1L,
                        natural = function(theta, aux) c(d = tanh(theta[1])),
                        sigma_varies = FALSE); "accepted"},
  error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("whiten = function(...):", r, "\n")
# student under the toy term vs dense dmvt at nu = 5
fs <- frm(y ~ x + toyar(1), family = student(), data = d)
p <- fs$obj$env$last.par.best
bd <- which(names(p) == "betad")
p[bd[2]] <- log(5 - 1)
b <- p[names(p) == "beta"]
sgs <- exp(p[bd[1]])
ph <- tanh(p[names(p) == "thetaac"])
Ss <- sgs^2 / (1 - ph^2) * ph^abs(outer(1:n, 1:n, "-"))
ref <- dmvt(d$y, as.vector(cbind(1, d$x) %*% b), Ss, df = 5, log = TRUE)
cat(sprintf("toy student nu 5: frmtmb %.10f dense %.10f diff %.3g\n",
            -fs$obj$fn(p), ref, -fs$obj$fn(p) - ref))
