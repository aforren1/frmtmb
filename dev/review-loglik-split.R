# Review of lane loglik, attack 1: is logLik() the data log-likelihood at
# the MAP estimates for each prior class? The reference is computed here,
# not by the lane's dev/loglik-routes.R:
# - "closed": a closed form from the public accessors (fixef(), VarCorr(),
#   sigma(), fitted()) and the data, with no tape and no prior code;
# - "tape": the no-prior fit's own tape f0$obj$fn() at the MAP fit's
#   optimizer vector, after checking that the two layouts match.
# Each row prints the lane's logLik(), logLik_penalized() and both
# references, and counts conditions logLik() raises.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")

dmvn <- function(y, mu, V) {
  R <- chol(V)
  z <- backsolve(R, y - mu, transpose = TRUE)
  -0.5 * length(y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z^2)
}
# gaussian LMM marginal, one block per group; Zfun(rows) gives Z_g
lmm_marg <- function(y, mu, g, Zfun, Sig, sigma) {
  s <- 0
  for (lv in levels(g)) {
    r <- which(g == lv)
    Z <- Zfun(r)
    V <- Z %*% Sig %*% t(Z) + diag(sigma^2, length(r))
    s <- s + dmvn(y[r], mu[r], V)
  }
  s
}
quiet_ll <- function(f) {
  n <- 0L
  v <- withCallingHandlers(
    tryCatch(as.numeric(logLik(f)),
             error = function(e) paste("ERROR", class(e)[1],
                                       substr(conditionMessage(e), 1, 70))),
    warning = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") },
    message = function(m) { n <<- n + 1L; invokeRestart("muffleMessage") })
  list(v = v, n = n)
}
tape_ref <- function(f0, fp) {
  if (!identical(names(f0$opt$par), names(fp$opt$par))) return(NA_real_)
  v <- -f0$obj$fn(fp$opt$par)
  f0$obj$fn(f0$opt$par)
  v
}
row <- function(case, fp, closed = NA_real_, tape = NA_real_) {
  q <- quiet_ll(fp)
  pen <- logLik_penalized(fp)
  ll <- suppressWarnings(as.numeric(q$v))
  rc <- if (is.na(closed)) NA else abs(ll - closed) / abs(closed)
  rt <- if (is.na(tape)) NA else abs(ll - tape) / abs(tape)
  cat(sprintf("%-34s logLik %s  pen %.10f  log_prior %.6f  closed %.10f rel %.2e  tape %.10f rel %.2e  conds %d\n",
              case, if (is.character(q$v)) q$v else sprintf("%.10f", q$v),
              as.numeric(pen), attr(pen, "log_prior"), closed, rc, tape, rt,
              q$n))
}

# S1. lm: b with a density and a bound, centered Intercept, sigma
set.seed(101)
d <- data.frame(x = rnorm(70), z = rnorm(70))
d$y <- 0.4 + 0.9 * d$x - 0.3 * d$z + rnorm(70, 0, 0.8)
p1 <- set_prior("normal(0, 0.3)", class = "b", coef = "x", lb = 0.1) +
  set_prior("normal(0, 0.2)", class = "b", coef = "z") +
  set_prior("student_t(3, 2, 0.3)", class = "Intercept") +
  set_prior("exponential(3)", class = "sigma")
fp <- frm(y ~ x + z, data = d, prior = p1)
f0 <- frm(y ~ x + z, data = d)
cl <- sum(dnorm(d$y, fitted(fp)[, "Estimate"], sigma(fp), log = TRUE))
row("S1 lm b(lb)+Intercept+sigma", fp, cl, tape_ref(f0, fp))

# S2. (1 + x | g) with LKJ on cor, sd, b, sigma
set.seed(102)
d <- data.frame(x = rnorm(120), g = factor(rep(1:15, each = 8)))
u <- cbind(rnorm(15, 0, 0.8), rnorm(15, 0, 0.4))
d$y <- 1 + 0.5 * d$x + u[d$g, 1] + u[d$g, 2] * d$x + rnorm(120, 0, 0.6)
p2 <- set_prior("lkj(4)", class = "cor") +
  set_prior("normal(0, 0.5)", class = "sd") +
  set_prior("normal(0, 0.2)", class = "b") +
  set_prior("student_t(3, 0, 0.5)", class = "sigma")
fp <- frm(y ~ x + (1 + x | g), data = d, prior = p2)
f0 <- frm(y ~ x + (1 + x | g), data = d)
vc <- VarCorr(fp)$g$cov[, "Estimate", ]
fe <- fixef(fp)[, "Estimate"]
mu <- fe[1] + fe[2] * d$x
cl <- lmm_marg(d$y, mu, d$g, function(r) cbind(1, d$x[r]), vc, sigma(fp))
row("S2 (1+x|g) lkj+sd+b+sigma", fp, cl, tape_ref(f0, fp))

# S3. (1 | g): sd with a density and a lower bound; gamma on sd
set.seed(103)
d <- data.frame(x = rnorm(90), g = factor(rep(1:18, each = 5)))
d$y <- 1 + 0.5 * d$x + rnorm(18, 0, 0.2)[d$g] + rnorm(90)
for (pr in list(set_prior("normal(1, 0.3)", class = "sd", lb = 0.4),
                set_prior("gamma(2, 2)", class = "sd"))) {
  fp <- frm(y ~ x + (1 | g), data = d, prior = pr)
  f0 <- frm(y ~ x + (1 | g), data = d)
  fe <- fixef(fp)[, "Estimate"]
  s2 <- VarCorr(fp)$g$sd[1, "Estimate"]^2
  cl <- lmm_marg(d$y, fe[1] + fe[2] * d$x, d$g,
                 function(r) matrix(1, length(r), 1), matrix(s2), sigma(fp))
  row(paste("S3 (1|g) sd", unclass(pr)[[1]]$prior %||% "?"), fp, cl,
      tape_ref(f0, fp))
}

# S4. REML: the restricted likelihood in closed form, beta at GLS
reml_closed <- function(y, X, g, sd_g, sigma) {
  V <- diag(sigma^2, length(y))
  for (lv in levels(g)) {
    r <- which(g == lv)
    V[r, r] <- V[r, r] + sd_g^2
  }
  Vi <- solve(V)
  XtViX <- t(X) %*% Vi %*% X
  bh <- solve(XtViX, t(X) %*% Vi %*% y)
  r <- y - X %*% bh
  n <- length(y); p <- ncol(X)
  -0.5 * ((n - p) * log(2 * pi) + as.numeric(determinant(V)$modulus) +
            as.numeric(determinant(XtViX)$modulus) +
            as.numeric(t(r) %*% Vi %*% r))
}
X <- cbind(1, d$x)
f0 <- frm(y ~ x + (1 | g), data = d, REML = TRUE)
cl0 <- reml_closed(d$y, X, d$g, VarCorr(f0)$g$sd[1, "Estimate"], sigma(f0))
row("S4 REML no prior (convention)", f0, cl0, NA)
if (requireNamespace("lme4", quietly = TRUE)) {
  lm4 <- lme4::lmer(y ~ x + (1 | g), data = d, REML = TRUE)
  cat(sprintf("   lme4 REML logLik %.10f  frm %.10f\n",
              as.numeric(logLik(lm4)), as.numeric(logLik(f0))))
}
fp <- frm(y ~ x + (1 | g), data = d, REML = TRUE,
          prior = set_prior("gamma(2, 2)", class = "sd") +
            set_prior("exponential(2)", class = "sigma"))
cl <- reml_closed(d$y, X, d$g, VarCorr(fp)$g$sd[1, "Estimate"], sigma(fp))
row("S4 REML sd+sigma prior", fp, cl, tape_ref(f0, fp))

# S5. nonlinear: priors on two nlpars
set.seed(105)
d <- data.frame(x = runif(60, 0, 4))
d$y <- 3 * exp(-0.7 * d$x) + rnorm(60, 0, 0.15)
p5 <- set_prior("normal(2, 0.3)", nlpar = "a") +
  set_prior("normal(1, 0.1)", nlpar = "k", lb = 0)
fp <- frm(bf(y ~ a * exp(-k * x), a + k ~ 1, nl = TRUE), data = d,
          prior = p5)
f0 <- frm(bf(y ~ a * exp(-k * x), a + k ~ 1, nl = TRUE), data = d,
          start = list(beta = unname(fp$estimates$beta)))
cl <- sum(dnorm(d$y, fitted(fp)[, "Estimate"], sigma(fp), log = TRUE))
row("S5 nonlinear nlpar a,k", fp, cl, tape_ref(f0, fp))

# S6. mo() with a Dirichlet on the simplex
set.seed(106)
d <- data.frame(xo = sample(1:5, 100, TRUE))
d$y <- 1 + c(0, 0.1, 0.9, 1.0, 2.0)[d$xo] + rnorm(100, 0, 0.5)
d$xo <- factor(d$xo, ordered = TRUE)
p6 <- set_prior("dirichlet(c(5, 5, 5, 5))", class = "simo",
                coef = "moxo1") +
  set_prior("normal(0, 1)", class = "b")
fp <- tryCatch(frm(y ~ mo(xo), data = d, prior = p6),
               error = function(e) { cat("S6 fit error:", conditionMessage(e), "\n"); NULL })
if (!is.null(fp)) {
  f0 <- frm(y ~ mo(xo), data = d)
  cl <- sum(dnorm(d$y, fitted(fp)[, "Estimate"], sigma(fp), log = TRUE))
  row("S6 mo() simplex dirichlet", fp, cl, tape_ref(f0, fp))
}

# S7. distributional: priors on sigma's coefficients (natural-scale chart)
set.seed(107)
d <- data.frame(x = rnorm(100))
d$y <- 1 + 0.5 * d$x + rnorm(100, 0, exp(-0.3 + 0.4 * d$x))
p7 <- set_prior("normal(0, 0.2)", class = "b", dpar = "sigma") +
  set_prior("normal(0, 0.5)", class = "Intercept", dpar = "sigma")
fp <- frm(bf(y ~ x, sigma ~ x), data = d, prior = p7)
f0 <- frm(bf(y ~ x, sigma ~ x), data = d)
sg <- tryCatch(fitted(fp, dpar = "sigma")[, "Estimate"],
               error = function(e) NULL)
cl <- if (is.null(sg)) NA_real_ else
  sum(dnorm(d$y, fitted(fp, dpar = "mu")[, "Estimate"], sg, log = TRUE))
row("S7 dpar sigma b+Intercept", fp, cl, tape_ref(f0, fp))

# S8. ar(1) gaussian with a prior on ar
set.seed(108)
d <- data.frame(t = rep(1:12, 10), g = factor(rep(1:10, each = 12)),
                x = rnorm(120))
e <- unlist(lapply(1:10, function(i) as.numeric(arima.sim(list(ar = 0.6), 12))))
d$y <- 1 + 0.5 * d$x + 0.7 * e
fp <- frm(y ~ x + ar(time = t, gr = g, p = 1), data = d,
          prior = set_prior("normal(0, 0.2)", class = "ar") +
            set_prior("normal(0, 0.3)", class = "b"))
f0 <- frm(y ~ x + ar(time = t, gr = g, p = 1), data = d)
row("S8 ar(1) ar+b", fp, NA, tape_ref(f0, fp))

# S9. exact gp: lscale inv_gamma, sdgp
set.seed(109)
d <- data.frame(x = sort(runif(40, 0, 5)))
d$y <- sin(d$x) + rnorm(40, 0, 0.3)
p9 <- set_prior("inv_gamma(4, 2)", class = "lscale") +
  set_prior("normal(0, 1)", class = "sdgp")
fp <- tryCatch(frm(y ~ gp(x), data = d, prior = p9),
               error = function(e) { cat("S9 fit error:", conditionMessage(e), "\n"); NULL })
if (!is.null(fp)) {
  f0 <- frm(y ~ gp(x), data = d)
  row("S9 exact gp lscale+sdgp", fp, NA, tape_ref(f0, fp))
}

# S10. ordinal thresholds (ordthres chart)
set.seed(110)
d <- data.frame(x = rnorm(150))
d$y <- cut(0.8 * d$x + rlogis(150), c(-Inf, -1, 0, 1.2, Inf),
           labels = FALSE)
d$y <- factor(d$y, ordered = TRUE)
p10 <- set_prior("normal(0, 0.5)", class = "Intercept") +
  set_prior("normal(0, 0.3)", class = "b")
fp <- frm(y ~ x, data = d, family = cumulative(), prior = p10)
f0 <- frm(y ~ x, data = d, family = cumulative())
pr <- fitted(fp)
cl <- tryCatch({
  P <- if (length(dim(pr)) == 3L) pr[, "Estimate", ] else NULL
  sum(log(P[cbind(seq_len(nrow(d)), as.integer(d$y))]))
}, error = function(e) NA_real_)
row("S10 cumulative thresholds+b", fp, cl, tape_ref(f0, fp))

# S11. poisson (1 | g), Laplace and quadrature
set.seed(111)
d <- data.frame(x = rnorm(100), g = factor(rep(1:20, each = 5)))
d$yc <- rpois(100, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.5)[d$g]))
p11 <- set_prior("normal(0, 0.2)", class = "b") +
  set_prior("normal(0, 0.3)", class = "sd") +
  set_prior("normal(1, 0.3)", class = "Intercept")
for (q in c(FALSE, TRUE)) {
  fp <- frm(yc ~ x + (1 | g), data = d, family = poisson(), prior = p11,
            quadrature = q)
  f0 <- frm(yc ~ x + (1 | g), data = d, family = poisson(), quadrature = q)
  row(paste("S11 poisson (1|g) quadrature", q), fp, NA, tape_ref(f0, fp))
}
cat("done\n")
