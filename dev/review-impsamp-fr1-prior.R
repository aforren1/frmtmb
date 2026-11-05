# Fix round 1, item 2: priors under importance =. Three references:
# (a) gaussian: Laplace with the same prior is exact, so importance must equal it;
# (b) bernoulli (1 | g) and (x | g): my own exact posterior mode, base R densities;
# (c) logLik must read as frm()'s other routes read: compare each with the exact
#     data likelihood at its estimate, with and without the prior term.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
options(width = 170)
sq <- function(e) suppressWarnings(suppressMessages(e))
ng_f <- function(f, p, h = 1e-4) {
  vapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, h)
    (f(p + e) - f(p - e)) / (2 * h)
  }, 0)
}
# (a) gaussian identity
set.seed(5)
G <- 25L
g <- factor(rep(seq_len(G), each = 6))
x <- rnorm(150)
z <- rnorm(150)
a <- rnorm(G, 0, 0.8)
b <- rnorm(G, 0, 0.4)
y <- 0.3 + 0.5 * x + a[g] + b[g] * x + rnorm(150, 0, exp(-0.5 + 0.3 * z))
d <- data.frame(y, x, z, g)
pri <- list(
  b = set_prior("normal(0, 0.1)", class = "b"),
  Intercept = set_prior("normal(2, 0.2)", class = "Intercept"),
  sd = set_prior("normal(0, 0.1)", class = "sd"),
  sd_coef = set_prior("student_t(3, 0, 0.3)", class = "sd", coef = "x", group = "g"),
  cor = set_prior("lkj(4)", class = "cor"),
  cor2 = set_prior("lkj(0.5)", class = "cor"),
  sigma_b = set_prior("normal(1, 0.05)", class = "b", dpar = "sigma"),
  sigma_int = set_prior("normal(-2, 0.1)", class = "Intercept", dpar = "sigma"),
  lb = set_prior("normal(0, 5)", class = "b", lb = 0.7),
  combo = c(set_prior("normal(0, 0.1)", class = "b"), set_prior("lkj(4)", class = "cor"),
            set_prior("normal(0, 0.2)", class = "sd"))
)
f <- bf(y ~ x + (x | g), sigma ~ z)
cat("(a) gaussian: importance at 64 draws against Laplace with the same prior\n")
for (nm in names(pri)) {
  w <- character(0)
  fi <- withCallingHandlers(tryCatch(frm(f, family = gaussian(), data = d, prior = pri[[nm]],
                                         importance = 64L), error = function(e) e),
                            warning = function(x) {
                              w <<- c(w, substr(conditionMessage(x), 1, 60))
                              invokeRestart("muffleWarning")
                            })
  fl <- sq(frm(f, family = gaussian(), data = d, prior = pri[[nm]]))
  f0 <- sq(frm(f, family = gaussian(), data = d))
  if (inherits(fi, "error")) {
    cat(sprintf("  %-10s STOP: %s\n", nm, substr(conditionMessage(fi), 1, 150)))
    next
  }
  cat(sprintf("  %-10s obj-lap %+.3g  maxpar %.3g  moved from no-prior %.3f  logLik-lap %+.3g  obj==-logLik %s  rounds %d warns %d %s\n",
              nm, fi$opt$objective - fl$opt$objective, max(abs(fi$opt$par - fl$opt$par)),
              max(abs(fi$opt$par - f0$opt$par)), as.numeric(logLik(fi)) - as.numeric(logLik(fl)),
              isTRUE(all.equal(fi$opt$objective, -as.numeric(logLik(fi)))),
              fi$importance$rounds, length(w), paste(w, collapse = ";")))
}
# (b) bernoulli (1 | g): exact posterior mode from base R densities
set.seed(7)
G <- 30L
g <- factor(rep(seq_len(G), each = 5))
x <- rnorm(150)
y <- rbinom(150, 1, plogis(0.3 + 0.6 * x + rnorm(G, 0, 1.2)[g]))
d1 <- data.frame(y, x, g)
gh <- gh_rule(61L)
idx <- split(seq_along(y), g)
ml1 <- function(p) {
  sd <- exp(p[3])
  tot <- 0
  for (ii in idx) {
    fn <- function(u) {
      vapply(u, function(ui) sum(dbinom(y[ii], 1, plogis(p[1] + p[2] * x[ii] + ui), log = TRUE)), 0) +
        dnorm(u, 0, sd, log = TRUE)
    }
    o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
    m <- o$minimum
    cur <- -(fn(m + 1e-3) - 2 * fn(m) + fn(m - 1e-3)) / 1e-6
    s <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
    nd <- m + sqrt(2) * s * gh$x
    a <- fn(nd) + gh$x^2 + log(gh$w)
    tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * s))
  }
  tot
}
hand <- list(
  b = list(prior = set_prior("normal(0, 0.1)", class = "b"),
           lp = function(p) dnorm(p[2], 0, 0.1, log = TRUE)),
  sd = list(prior = set_prior("normal(0, 0.3)", class = "sd"),
            lp = function(p) dnorm(exp(p[3]), 0, 0.3, log = TRUE) + log(2)),
  sd_t = list(prior = set_prior("student_t(3, 0, 0.3)", class = "sd"),
              lp = function(p) dt(exp(p[3]) / 0.3, 3, log = TRUE) - log(0.3) + log(2)),
  sd_gamma = list(prior = set_prior("gamma(4, 2)", class = "sd"),
                  lp = function(p) dgamma(exp(p[3]), 4, 2, log = TRUE)))
cat("\n(b) bernoulli (1 | g): importance at 1000 draws, quadrature, my exact mode\n")
for (nm in names(hand)) {
  pr <- hand[[nm]]$prior
  fi <- sq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pr, importance = 1000L))
  fq <- sq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pr, quadrature = TRUE))
  fl <- sq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pr))
  res <- list()
  for (jac in c(FALSE, TRUE)) {
    obj <- function(p) ml1(p) - hand[[nm]]$lp(p) - if (jac) p[3] else 0
    o <- nlminb(unname(fi$opt$par), obj, gradient = function(p) ng_f(obj, p))
    res[[as.character(jac)]] <- o
  }
  o <- res[["FALSE"]]
  oj <- res[["TRUE"]]
  H <- numDeriv::jacobian(function(p) ng_f(function(q) ml1(q) - hand[[nm]]$lp(q), p), o$par, method = "simple")
  se <- sqrt(diag(solve(0.5 * (H + t(H)))))
  p <- unname(fi$opt$par)
  cat(sprintf("  %-9s theta: imp %.4f quad %.4f lap %.4f exact-noJac %.4f exact-Jac %.4f | maxz(imp vs noJac) %.3f | SE imp/exact %s\n",
              nm, p[3], fq$opt$par[[3]], fl$opt$par[[3]], o$par[3], oj$par[3], max(abs(p - o$par) / se),
              paste(format(sqrt(diag(frmtmb:::sdr_of(fi)$cov.fixed)) / se, digits = 3), collapse = " ")))
  # (c) logLik: the data alone, or the data plus the prior?
  dat <- ml1(p)
  cat(sprintf("            -logLik imp %.4f  obj imp %.4f  exact data nll at fit %.4f  exact data+prior %.4f  mcse %.4f\n",
              -as.numeric(logLik(fi)), fi$opt$objective, dat, dat - hand[[nm]]$lp(p), fi$importance$mcse))
  dq <- ml1(unname(fq$opt$par))
  cat(sprintf("            quadrature: -logLik %.4f exact data %.4f data+prior %.4f | Laplace: -logLik %.4f exact data %.4f data+prior %.4f\n",
              -as.numeric(logLik(fq)), dq, dq - hand[[nm]]$lp(unname(fq$opt$par)),
              -as.numeric(logLik(fl)), ml1(unname(fl$opt$par)),
              ml1(unname(fl$opt$par)) - hand[[nm]]$lp(unname(fl$opt$par))))
  cat(sprintf("            AIC imp %.3f, anova: %s\n", AIC(fi), tryCatch({
    fi0 <- sq(frm(bf(y ~ 1 + (1 | g)), family = bernoulli(), data = d1, importance = 200L))
    capture.output(anova(fi0, fi))[1]
  }, error = function(e) substr(conditionMessage(e), 1, 110))))
}
# (d) lkj on a correlated slope, against my own 2-d exact mode
set.seed(9)
G <- 30L
g <- rep(seq_len(G), each = 8)
x <- rnorm(240)
bb <- matrix(rnorm(2 * G), G, 2) %*% chol(matrix(c(1, 0.6, 0.6, 1), 2)) %*% diag(c(1, 0.5))
y <- rbinom(240, 1, plogis(0.2 + 0.5 * x + bb[g, 1] + bb[g, 2] * x))
d2 <- data.frame(y, x, g = factor(g))
idx2 <- split(seq_along(y), g)
g1 <- gh_rule(21L)
gr <- as.matrix(expand.grid(1:21, 1:21))
xn <- matrix(g1$x[gr], nrow(gr), 2L)
lw <- rowSums(matrix(log(g1$w)[gr], nrow(gr), 2L)) + rowSums(xn^2)
modes <- lapply(idx2, function(i) c(0, 0))
ml2 <- function(p) {
  S <- as.matrix(us_sigma(p[3:5], 2L))
  ch <- chol(S)
  Pi <- chol2inv(ch)
  ld <- 2 * sum(log(diag(ch)))
  lpr <- function(U) -0.5 * rowSums((U %*% Pi) * U) - log(2 * pi) - 0.5 * ld
  tot <- 0
  for (k in seq_along(idx2)) {
    ii <- idx2[[k]]
    ll <- function(U) {
      eta <- outer(p[1] + p[2] * x[ii], rep(1, nrow(U))) + cbind(1, x[ii]) %*% t(U)
      colSums(matrix(dbinom(rep(y[ii], nrow(U)), 1, plogis(as.vector(eta)), log = TRUE), length(ii), nrow(U)))
    }
    neg <- function(u) {
      U <- matrix(u, 1L, 2L)
      -(ll(U) + lpr(U))
    }
    o <- optim(modes[[k]], neg, method = "BFGS", control = list(reltol = 1e-10))
    modes[[k]] <<- o$par
    Lc <- backsolve(chol(optimHess(o$par, neg)), diag(2))
    U <- matrix(o$par, nrow(xn), 2L, byrow = TRUE) + sqrt(2) * xn %*% t(Lc)
    a <- ll(U) + lpr(U) + lw
    tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(2) + sum(log(abs(diag(Lc)))))
  }
  tot
}
cat("\n(d) bernoulli (x | g), lkj: importance at 1000 draws against my exact mode\n")
for (eta in c(1, 4, 20)) {
  # lkj(eta) density on a 2 x 2 correlation r: (1 - r^2)^(eta - 1), const dropped
  lpj <- function(p) {
    S <- as.matrix(us_sigma(p[3:5], 2L))
    r <- S[1, 2] / sqrt(S[1, 1] * S[2, 2])
    (eta - 1) * log(1 - r^2)
  }
  pr <- set_prior(sprintf("lkj(%g)", eta), class = "cor")
  fi <- sq(frm(bf(y ~ x + (x | g)), family = bernoulli(), data = d2, prior = pr, importance = 1000L))
  fl <- sq(frm(bf(y ~ x + (x | g)), family = bernoulli(), data = d2, prior = pr))
  # the chart's Jacobian d r / d theta_3: r = t / sqrt(1 + t^2), log |dr/dt| = -1.5 log(1 + t^2)
  for (jac in c(FALSE, TRUE)) {
    ob <- function(p) ml2(p) - lpj(p) - if (jac) -1.5 * log(1 + p[5]^2) else 0
    o <- nlminb(unname(fi$opt$par), ob, gradient = function(p) ng_f(ob, p), control = list(rel.tol = 1e-8, iter.max = 60))
    cat(sprintf("  lkj(%g) %s: exact %s nll %.4f\n", eta, if (jac) "with chart Jacobian" else "no Jacobian        ",
                paste(format(o$par, digits = 5), collapse = " "), o$objective))
  }
  cat(sprintf("  lkj(%g) importance %s rounds %d ess %.3f mcse %.4f\n          Laplace    %s\n", eta,
              paste(format(unname(fi$opt$par), digits = 5), collapse = " "), fi$importance$rounds, fi$importance$ess_min,
              fi$importance$mcse, paste(format(unname(fl$opt$par), digits = 5), collapse = " ")))
}
# (e) the proposal under a strong prior on the sd: curvature and ESS at the fit
cat("\n(e) prior pushing the sd far from the data's: ESS and rounds\n")
for (loc in c(0.05, 0.2, 5, 20)) {
  pr <- set_prior(sprintf("normal(%g, %g)", loc, loc / 20), class = "sd")
  fi <- sq(tryCatch(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pr, importance = 500L), error = function(e) e))
  if (inherits(fi, "error")) {
    cat(sprintf("  sd ~ normal(%g, %g): STOP %s\n", loc, loc / 20, substr(conditionMessage(fi), 1, 160)))
    next
  }
  fq <- sq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pr, quadrature = TRUE))
  cat(sprintf("  sd ~ normal(%g, %g): sd imp %.4f quad %.4f, rounds %d capped %s ess %.3f mcse %.4f obj-quad %+.4f\n", loc, loc / 20,
              exp(fi$opt$par[[3]]), exp(fq$opt$par[[3]]), fi$importance$rounds, fi$importance$capped,
              fi$importance$ess_min, fi$importance$mcse, fi$opt$objective - fq$opt$objective))
}
