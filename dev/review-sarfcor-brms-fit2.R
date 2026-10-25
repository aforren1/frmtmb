# Reviewer: brms 2.23.0 fits (rstan) for
#  (c) a row-standardized queen W with REAL eigenvalues, data drawn at
#      rho = -0.8 below brms's lower bound (min eigenvalue): frmtmb ML vs
#      brms posterior;
#  (d) brms posterior_epred() for lag SAR against solve(I - rho W, mu)
#      on one draw, and against frmtmb's fitted() formula;
#  (e) defect 2: a row-standardized rook W, complex eigen() values.
# Rscript dev/review-sarfcor-brms-fit2.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
cat("lib:", find.package("frmtmb"), "| brms", format(packageVersion("brms")), "\n")
lattice <- function(k, queen = FALSE) {
  g <- expand.grid(r = 1:k, c = 1:k)
  n <- nrow(g)
  outer(1:n, 1:n, function(i, j) {
    dr <- abs(g$r[i] - g$r[j]); dc <- abs(g$c[i] - g$c[j])
    as.numeric(if (queen) pmax(dr, dc) == 1 else dr + dc == 1)
  })
}
A <- lattice(8, queen = TRUE)
n <- nrow(A)
W <- A / rowSums(A)
ev <- eigen(W)$values
cat("queen 8x8 W: complex", is.complex(ev), "range", format(range(Re(ev)), digits = 5), "\n")
set.seed(7)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) + 0.8 * W, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
f <- frm(y ~ x + sar(W), data = d, data2 = list(W = W))
cv <- confint_varcorr(f)
cat(sprintf("frmtmb ML lagsar %.4f (delta CI %.4f, %.4f); logLik %.4f\n",
            cv$estimate, cv$lwr, cv$upr, logLik(f)))
# the ML restricted to brms's interval: profile at rho = min eigenvalue
lo_b <- min(Re(ev))
pr <- function(r) {
  Aw <- diag(n) - r * W
  ys <- as.vector(Aw %*% d$y)
  X <- cbind(1, d$x)
  b <- qr.coef(qr(X), ys)
  e <- ys - X %*% b
  s2 <- sum(e^2) / n
  -0.5 * n * log(2 * pi * s2) - 0.5 * n + as.numeric(determinant(Aw)$modulus)
}
rs <- seq(lo_b + 1e-6, 0.99, length.out = 400)
pl <- sapply(rs, pr)
cat(sprintf("profile logLik at frmtmb rho %.4f: %.4f; best inside brms's (%.4f, 1): rho %.4f, %.4f; LR = %.3f\n",
            cv$estimate, pr(cv$estimate), lo_b, rs[which.max(pl)], max(pl),
            2 * (pr(cv$estimate) - max(pl))))
b <- brms::brm(y ~ x + sar(W), data = d, data2 = list(W = W), chains = 4,
               iter = 2000, seed = 71, refresh = 0, silent = 2,
               backend = "rstan")
dr <- as.data.frame(b)
q <- quantile(dr$lagsar, c(0.025, 0.5, 0.975))
cat(sprintf("brms lagsar median %.4f (2.5%% %.4f, 97.5%% %.4f); share of draws within 0.02 of the bound %.3f; max Rhat %.4f\n",
            q[2], q[1], q[3], mean(dr$lagsar < lo_b + 0.02),
            max(brms::rhat(b), na.rm = TRUE)))
cat(sprintf("brms b_x median %.4f; frmtmb b_x %.4f (SE %.4f)\n", median(dr$b_x),
            fixef(f)["x", "Estimate"], fixef(f)["x", "Est.Error"]))

## (d) posterior_epred for lag SAR on one draw
pe <- brms::posterior_epred(b)
i <- 17
mu_i <- dr$b_Intercept[i] + dr$b_x[i] * d$x
ref <- as.vector(solve(diag(n) - dr$lagsar[i] * W, mu_i))
cat(sprintf("brms posterior_epred draw %d vs solve(I - rho W, mu): max diff %.3g\n",
            i, max(abs(pe[i, ] - ref))))
ref0 <- mu_i
cat(sprintf("  ... vs mu alone: max diff %.3g\n", max(abs(pe[i, ] - ref0))))
# brms with newdata of the same rows (counterfactual x + 1)
nd <- d; nd$x <- nd$x + 1
pe2 <- tryCatch(brms::posterior_epred(b, newdata = nd),
                error = function(e) conditionMessage(e))
if (is.character(pe2)) cat("brms posterior_epred(newdata = same rows):", pe2, "\n") else {
  ref2 <- as.vector(solve(diag(n) - dr$lagsar[i] * W, mu_i + dr$b_x[i]))
  cat(sprintf("brms posterior_epred(newdata = same 64 rows, x + 1): works, draw %d vs A^-1 mu_new %.3g\n",
              i, max(abs(pe2[i, ] - ref2))))
}
r <- tryCatch({fitted(f, newdata = nd); "works"}, error = function(e) conditionMessage(e))
cat("frmtmb fitted(newdata = same rows):", r, "\n")
ll <- tryCatch(dim(brms::log_lik(b)), error = function(e) conditionMessage(e))
cat("brms log_lik dim:", ll, "\n")

## (e) defect 2: rook row-standardized W
Ar <- lattice(8)
Wr <- Ar / rowSums(Ar)
cat("rook W eigen complex:", is.complex(eigen(Wr)$values), "\n")
d$yr <- as.vector(solve(diag(n) - 0.6 * Wr, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
r <- tryCatch({
  b2 <- update(b, formula. = yr ~ x + sar(Wr), newdata = d,
               data2 = list(Wr = Wr), chains = 1, iter = 400, seed = 3,
               refresh = 0)
  paste("brms fit, draws:", nrow(as.data.frame(b2)))
}, error = function(e) paste("brms ERROR:", conditionMessage(e)))
cat("rook W lag fit:", r, "\n")
fr <- frm(yr ~ x + sar(Wr), data = d, data2 = list(Wr = Wr))
cat(sprintf("frmtmb rook W lagsar %.4f\n", confint_varcorr(fr)$estimate))
# brms after taking Re(): does the same model fit?
Wr2 <- Wr
r <- tryCatch({
  sdat <- brms::standata(yr ~ x + sar(Wr2), data = d, data2 = list(Wr2 = Wr2))
  cat("standata eigenMsar class:", class(sdat$eigenMsar), "\n")
  "ok"
}, error = function(e) conditionMessage(e))
