# Reviewer, fix round 1: F1 (map check vs correct terms), F2 (lagsar
# prior on the natural scale), F9 (timings at n = 2500).
# Rscript dev/review-sarfcor-fr1-a.R
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb.spatial); library(mvtnorm); library(Matrix)})
cat("lib:", find.package("frmtmb"), "\n")
frmtmb_local_registries(environment())

## F1: an AR(1) term written from the new help page (scalar sigma,
## matrix-capable unwhiten)
ar1_build <- function(spec, rows, n_all, lookup) {
  n <- length(rows)
  list(S = sparseMatrix(i = 2:n, j = 1:(n - 1), x = 1, dims = c(n, n)),
       first = c(1, numeric(n - 1)))
}
ar1_whiten <- function(y, mu, sigma, theta, aux) {
  phi <- tanh(theta[1]); e <- (y - mu) / sigma
  u <- e - phi * as.vector(aux$S %*% e) + (sqrt(1 - phi^2) - 1) * e * aux$first
  list(u = u, logjac = 0.5 * log(1 - phi^2) - length(y) * log(sigma))
}
ar1_unwhiten <- function(u, mu, sigma, theta, aux) {
  phi <- tanh(theta[1]); U <- as.matrix(u)
  U[1, ] <- U[1, ] / sqrt(1 - phi^2)
  E <- as.matrix(stats::filter(U, phi, method = "recursive"))
  out <- mu + sigma * E
  if (is.matrix(u)) out else as.vector(out)
}
frmtmb_register_autocor("toyar", args = "M", build = ar1_build,
                        whiten = ar1_whiten, unwhiten = ar1_unwhiten,
                        npar = 1L, sigma_varies = FALSE,
                        natural = function(theta, aux) c(toyphi = tanh(theta[1])))
set.seed(41)
n <- 80
e <- as.numeric(stats::filter(rnorm(n) * c(1 / sqrt(1 - 0.36), rep(1, n - 1)),
                              0.6, method = "recursive"))
d <- data.frame(x = rnorm(n), t = 1:n)
d$y <- 1 + 0.5 * d$x + 0.8 * e
fit <- frm(y ~ x + toyar(1), data = d)
phi <- tanh(fit$estimates$thetaac); sg <- sigma(fit)
mu <- as.vector(cbind(1, d$x) %*% fixef(fit)[, "Estimate"])
S <- sg^2 / (1 - phi^2) * phi^abs(outer(1:n, 1:n, "-"))
fb <- frm(y ~ x + ar(time = t, p = 1, cov = TRUE), data = d)
cat(sprintf("F1 toy from the new page: fits; logLik - dense %.3g; - ar(cov=TRUE) %.3g; pearson %.3g; autocor_matrix %.3g\n",
            as.numeric(logLik(fit)) - dmvnorm(d$y, mu, S, log = TRUE),
            as.numeric(logLik(fit) - logLik(fb)),
            max(abs(residuals(fit, type = "pearson")[, "Estimate"] -
                      (d$y - mu) / sqrt(diag(S)))),
            max(abs(unname(as.matrix(autocor_matrix(fit))) - S / sg^2))))
# a correct term whose unwhiten is vector-only (the 0.74.0-era contract)
frmtmb_register_autocor("toyvec", args = "M", build = ar1_build,
                        whiten = ar1_whiten, npar = 1L, sigma_varies = FALSE,
                        unwhiten = function(u, mu, sigma, theta, aux) {
                          phi <- tanh(theta[1]); u <- as.numeric(u)
                          u[1] <- u[1] / sqrt(1 - phi^2)
                          mu + sigma * as.numeric(stats::filter(u, phi, method = "recursive"))
                        },
                        natural = function(theta, aux) c(tv = tanh(theta[1])))
r <- tryCatch({frm(y ~ x + toyvec(1), data = d); "FIT"},
              error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("F1 vector-only unwhiten:", r, "\n")
# a correct term at a start far from 0 (phi = tanh(3) = 0.995)
frmtmb_register_autocor("toyfar", args = "M", build = ar1_build,
                        whiten = ar1_whiten, unwhiten = ar1_unwhiten,
                        npar = 1L, sigma_varies = FALSE, start = 3,
                        natural = function(theta, aux) c(tf = tanh(theta[1])))
r <- tryCatch({f <- frm(y ~ x + toyfar(1), data = d)
  sprintf("FIT logLik diff to toyar %.3g", as.numeric(logLik(f) - logLik(fit)))},
  error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("F1 start at phi = 0.995:", r, "\n")
# fcor with an ill-conditioned but positive definite M (cond ~ 1e10)
Vb <- 0.99999^abs(outer(1:n, 1:n, "-"))
cat("cond(Vb):", format(kappa(Vb, exact = TRUE), digits = 3), "\n")
r <- tryCatch({frm(y ~ x + fcor(Vb), data = d, data2 = list(Vb = Vb)); "FIT"},
              error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("F1 fcor ill-conditioned M:", r, "\n")

## F9 and F1 at n = 2500, sparse W (row-standardized and symmetric)
k <- 50; nn <- k * k
id <- matrix(1:nn, k)
ii <- c(id[-k, ], id[, -k]); jj <- c(id[-1, ], id[, -1])
A <- sparseMatrix(i = c(ii, jj), j = c(jj, ii), x = 1, dims = c(nn, nn))
Wr <- Diagonal(x = 1 / rowSums(A)) %*% A
dg <- 1 / sqrt(rowSums(A)); Ws <- Diagonal(x = dg) %*% A %*% Diagonal(x = dg)
set.seed(8)
dd <- data.frame(x = rnorm(nn))
dd$y <- as.vector(solve(Diagonal(nn) - 0.5 * Ws, 1 + 0.5 * dd$x + rnorm(nn)))
dd$ye <- 1 + 0.5 * dd$x + as.vector(solve(Diagonal(nn) - 0.5 * Wr, rnorm(nn)))
tm <- function(lab, expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(force(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-40s %8.2f s\n", lab, proc.time()[["elapsed"]] - t0))
  invisible(r)
}
f1 <- tm("frm lag, symmetric W, n = 2500", frm(y ~ x + sar(Ws), data = dd, data2 = list(Ws = Ws)))
cat("  ", if (is.character(f1)) f1 else paste("lagsar", format(confint_varcorr(f1)$estimate, digits = 4)), "\n")
f2 <- tm("frm error, row-std W, n = 2500", frm(ye ~ x + sar(Wr, type = "error"), data = dd, data2 = list(Wr = Wr)))
cat("  ", if (is.character(f2)) f2 else paste("errorsar", format(confint_varcorr(f2)$estimate, digits = 4)), "\n")
for (ff in list(f1, f2)) {
  if (is.character(ff)) next
  tm("  residuals(pearson)", rp <- residuals(ff, type = "pearson"))
  tm("  autocor_matrix()", am <- autocor_matrix(ff))
  tm("  simulate(nsim = 20)", simulate(ff, nsim = 20))
  # pearson against the dense marginal of autocor_matrix
  cat(sprintf("  pearson vs dense: %.3g\n",
              max(abs(rp[, "Estimate"] - residuals(ff)[, "Estimate"] /
                        (sigma(ff) * sqrt(diag(as.matrix(am))))))))
}

## F2: a normal(0.3, 0.2) prior on lagsar against a hand-written MAP
set.seed(2)
k <- 7
gr <- expand.grid(r = seq_len(k), c = seq_len(k)); n7 <- nrow(gr)
A7 <- outer(1:n7, 1:n7, function(i, j)
  as.numeric(abs(gr$r[i] - gr$r[j]) + abs(gr$c[i] - gr$c[j]) == 1))
d7 <- 1 / sqrt(rowSums(A7)); W7 <- A7 * outer(d7, d7)
ds <- data.frame(x = rnorm(n7))
ds$y <- as.vector(solve(diag(n7) - 0.5 * W7, 1 + 0.5 * ds$x + rnorm(n7, 0, 0.7)))
ev <- eigen(W7, symmetric = TRUE, only.values = TRUE)$values
lo <- 1 / min(ev); hi <- 1 / max(ev)
X <- cbind(1, ds$x)
ll <- function(rho, b, s) {
  Aw <- diag(n7) - rho * W7
  dmvnorm(ds$y, solve(Aw, X %*% b), s^2 * solve(crossprod(Aw)), log = TRUE)
}
obj <- function(p, jac = TRUE) {
  sp <- plogis(p[4]); rho <- lo + (hi - lo) * sp
  v <- ll(rho, p[1:2], exp(p[3])) + dnorm(rho, 0.3, 0.2, log = TRUE)
  if (jac) v <- v + log(hi - lo) + log(sp) + log(1 - sp)
  -v
}
fp <- frm(y ~ x + sar(W7), data = ds, data2 = list(W7 = W7),
          prior = set_prior("normal(0.3, 0.2)", class = "lagsar"))
f0 <- frm(y ~ x + sar(W7), data = ds, data2 = list(W7 = W7))
pst <- c(fixef(fp)[, "Estimate"], log(sigma(fp)), fp$estimates$thetaac)
o1 <- optim(pst, obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
o0 <- optim(pst, obj, jac = FALSE, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
rho_of <- function(th) lo + (hi - lo) * plogis(th)
cat(sprintf("F2 ML rho %.6f; frm with normal(0.3,0.2) rho %.6f\n",
            confint_varcorr(f0)$estimate, confint_varcorr(fp)$estimate))
cat(sprintf("F2 hand MAP with logistic Jacobian rho %.6f (obj %.8f); without Jacobian rho %.6f\n",
            rho_of(o1$par[4]), -o1$value, rho_of(o0$par[4])))
cat(sprintf("F2 frm objective at its optimum %.8f; hand objective at frm's point %.8f\n",
            -fp$obj$fn(fp$obj$env$last.par.best), -obj(pst)))
cat(sprintf("F2 b: frm %s hand %s\n", paste(format(pst[1:2], digits = 7), collapse = " "),
            paste(format(o1$par[1:2], digits = 7), collapse = " ")))
print(frmtmb_autocor_pars(fp))
