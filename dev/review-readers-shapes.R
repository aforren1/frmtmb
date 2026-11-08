# Review of lane readers: fitted() at new levels of several block shapes.
# Each reference is computed here, independent of the lane's code: a closed
# form exp(eta + v / 2) on a log link with v from VarCorr(), or a nested
# stats::integrate() for the two-direction hurdle. Data seeds per case below.
source("dev/review-readers-setup.R")
rel <- function(a, b) format(signif(a / b - 1, 3))
eta_na <- function(f, nd, dp = NULL) {
  as.numeric(frm_linpred(f, newdata = nd, re_formula = NA, dpar = dp))
}
show <- function(lbl, f, nd, ref, ...) {
  r <- with_warn(fitted(f, newdata = nd, allow_new_levels = TRUE, ...))
  ff <- r$value
  cat(sprintf("%-34s est %s | ref %s | rel %s | err %s | warn %s\n", lbl,
              paste(format(ff[, "Estimate"], digits = 6), collapse = " "),
              paste(format(ref, digits = 6), collapse = " "),
              paste(rel(ff[, "Estimate"], ref), collapse = " "),
              paste(format(ff[, "Est.Error"], digits = 4), collapse = " "),
              paste(wclass(r$warnings), collapse = ",")))
  invisible(ff)
}
J <- 30; n <- 600

# 1. two new grouping factors at once, seed 5201
set.seed(5201)
d <- data.frame(x = runif(n, -1, 1), g = factor(sample(J, n, TRUE)),
                h = factor(sample(J, n, TRUE)))
ug <- rnorm(J, 0, 0.8); uh <- rnorm(J, 0, 0.6)
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + ug[d$g] + uh[d$h]))
f1 <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = poisson(), data = d)
vc <- VarCorr(f1)
nd <- data.frame(x = c(-0.5, 0.5), g = "new", h = "new")
show("two new factors", f1, nd,
     exp(eta_na(f1, nd) + (vc$g$sd[1, 1]^2 + vc$h$sd[1, 1]^2) / 2))
nd2 <- data.frame(x = c(-0.5, 0.5), g = "new", h = d$h[1])
eh <- as.numeric(frm_linpred(f1, newdata = transform(nd2, g = d$g[1]))) -
  as.numeric(frm_linpred(f1, newdata = transform(nd2, g = d$g[1]),
                         re_formula = ~ (1 | g)))
show("one new, one seen", f1, nd2,
     exp(eta_na(f1, nd2) + eh + vc$g$sd[1, 1]^2 / 2))

# 2. nested (1 | a / b): new b inside a seen a, seed 5202
set.seed(5202)
A <- 10; B <- 6
d <- data.frame(x = runif(n, -1, 1), a = factor(sample(A, n, TRUE)))
d$b <- factor(paste0(d$a, "_", sample(B, n, TRUE)))
ua <- rnorm(A, 0, 0.7); ub <- rnorm(nlevels(d$b), 0, 0.5)
names(ub) <- levels(d$b)
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + ua[d$a] + ub[as.character(d$b)]))
f2 <- frm(bf(y ~ x + (1 | a / b)), family = poisson(), data = d)
vc <- VarCorr(f2)
print(names(vc))
sab <- vc[[grep(":", names(vc))[1]]]$sd[1, 1]
nd <- data.frame(x = c(-0.5, 0.5), a = d$a[1], b = "zz")
ea <- as.numeric(frm_linpred(f2, newdata = nd, allow_new_levels = TRUE))
show("nested: new b in seen a", f2, nd, exp(ea + sab^2 / 2))

# 3. || block, seed 5203
set.seed(5203)
d <- data.frame(x = runif(n, -1, 1), g = factor(sample(J, n, TRUE)))
u0 <- rnorm(J, 0, 0.7); u1 <- rnorm(J, 0, 0.5)
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + u0[d$g] + u1[d$g] * d$x))
f3 <- frm(bf(y ~ x + (1 + x || g)), family = poisson(), data = d)
vc <- VarCorr(f3)$g
nd <- data.frame(x = c(-0.9, 0, 0.9), g = "new")
show("|| block", f3, nd,
     exp(eta_na(f3, nd) + (vc$sd[1, 1]^2 + nd$x^2 * vc$sd[2, 1]^2) / 2))

# 4. correlated slope across two predictors: hurdle (1 | p | g) on mu and hu,
# two directions. seed 5204
set.seed(5204)
d <- data.frame(x = runif(n, -1, 1), g = factor(sample(J, n, TRUE)))
S <- matrix(c(0.64, -0.3, -0.3, 0.49), 2)
U <- MASS::mvrnorm(J, c(0, 0), S)
lam <- exp(0.6 + 0.4 * d$x + U[d$g, 1])
cnt <- rpois(n, lam)
while (any(cnt == 0)) cnt[cnt == 0] <- rpois(sum(cnt == 0), lam[cnt == 0])
d$y <- ifelse(runif(n) < plogis(-0.4 + U[d$g, 2]), 0L, cnt)
f4 <- frm(bf(y ~ x + (1 | p | g), hu ~ (1 | p | g)),
          family = hurdle_poisson(), data = d)
vc <- VarCorr(f4)$g
print(vc)
sd2 <- vc$sd[, 1]
rho <- vc$cor[1, 1, 2]
Sh <- diag(sd2) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(sd2)
Lh <- t(chol(Sh))
nd <- data.frame(x = c(-0.5, 0.5), g = "new")
em <- eta_na(f4, nd, "mu")
ehu <- eta_na(f4, nd, "hu")
hmean <- function(i, p = 1) {
  inner <- function(z1) vapply(z1, function(a) {
    stats::integrate(function(z2) {
      um <- Lh[1, 1] * a
      uh <- Lh[2, 1] * a + Lh[2, 2] * z2
      l <- exp(em[i] + um)
      ((1 - plogis(ehu[i] + uh)) * l / (1 - exp(-l)))^p * dnorm(z2)
    }, -12, 12, rel.tol = 1e-11)$value * dnorm(a)
  }, 0)
  stats::integrate(inner, -12, 12, rel.tol = 1e-10)$value
}
ref <- vapply(1:2, hmean, 0)
ref2 <- vapply(1:2, hmean, 0, p = 2)
ff <- show("hurdle (1|p|g) mu+hu, 2 dirs", f4, nd, ref)
cat("  draw sd ref", format(sqrt(ref2 - ref^2), digits = 6), "\n")
# the same row with the cross-predictor correlation removed (rho = 0), to see
# whether the fit's correlation is read at all
cat("  rho", rho, "\n")

# 5. mm(): two members, seed 5205
set.seed(5205)
d <- data.frame(x = runif(n, -1, 1), g1 = factor(sample(J, n, TRUE)),
                g2 = factor(sample(J, n, TRUE)))
levels(d$g2) <- levels(d$g1) <- as.character(seq_len(J))
um <- rnorm(J, 0, 0.8)
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + 0.5 * um[d$g1] + 0.5 * um[d$g2]))
f5 <- frm(bf(y ~ x + (1 | mm(g1, g2))), family = poisson(), data = d)
vc <- VarCorr(f5)
s5 <- vc[[1]]$sd[1, 1]
cat("mm sd", s5, "\n")
nd <- data.frame(x = c(-0.5, 0.5), g1 = "n1", g2 = "n2")
show("mm two different new", f5, nd, exp(eta_na(f5, nd) + 0.5 * s5^2 / 2))
nd <- data.frame(x = c(-0.5, 0.5), g1 = "n1", g2 = "n1")
show("mm same new twice (one level)", f5, nd, exp(eta_na(f5, nd) + s5^2 / 2))
lb <- frm_lp_basis(f5, newdata = nd, allow_new_levels = TRUE)
cat("  extra_var same-new", format(lb$extra_var, digits = 5),
    " s^2 =", s5^2, " s^2/2 =", s5^2 / 2, "\n")
nd <- data.frame(x = c(-0.5, 0.5), g1 = "1", g2 = "n2")
e1 <- as.numeric(frm_linpred(f5, newdata = nd, allow_new_levels = TRUE))
show("mm one seen, one new", f5, nd, exp(e1 + 0.25 * s5^2 / 2))

# 6. gaussian new level: mean unchanged, Est.Error holds the sd
set.seed(5206)
d <- data.frame(x = runif(n, -1, 1), g = factor(sample(J, n, TRUE)))
d$y <- 1 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g] + rnorm(n, 0, 0.5)
f6 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
nd <- data.frame(x = c(-0.5, 0.5), g = "new")
ff <- show("gaussian new level", f6, nd, eta_na(f6, nd))
pl <- frm_linpred(f6, newdata = nd, allow_new_levels = TRUE, se.fit = TRUE)
cat("  Est.Error / frm_linpred se.fit:",
    format(ff[, "Est.Error"] / pl$se.fit, digits = 10), "\n")

# 7. Q columns at a one-direction log-link new level: quantiles of
# exp(eta + s z), widened by the coefficients only
f7 <- f3
nd <- data.frame(x = c(0), g = "new")
ff <- fitted(f7, newdata = nd, allow_new_levels = TRUE)
vc <- VarCorr(f7)$g
s7 <- vc$sd[1, 1]
cat("|| Q2.5", ff[, "Q2.5"], " exp(eta - 1.96 s)",
    exp(eta_na(f7, nd) - qnorm(0.975) * s7), " Q97.5", ff[, "Q97.5"],
    " exp(eta + 1.96 s)", exp(eta_na(f7, nd) + qnorm(0.975) * s7), "\n")

# 8. mixed rows: seen and new in one call
nd <- data.frame(x = c(-0.5, 0.5), g = c(as.character(d$g[1]), "new"))
ff <- fitted(f6, newdata = nd, allow_new_levels = TRUE)
fs <- fitted(f6, newdata = nd[1, ], allow_new_levels = TRUE)
cat("seen row alone vs in a mixed call:", ff[1, "Estimate"] - fs[1, "Estimate"],
    ff[1, "Est.Error"] - fs[1, "Est.Error"], "\n")
