# Multi-dimensional blocks against an independent exact marginal: product
# Gauss-Hermite around each group's mode, densities from base R. Checks the
# draws' m and L by regression at a moved theta (q = 2, 3), the estimates, SEs,
# stationarity, and the gradient away from the optimum.
# Usage: Rscript dev/review-impsamp-multid.R <design> <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
args <- commandArgs(TRUE)
nm <- args[[1]]
seed <- as.integer(args[[2]])
out <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/multid"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rmv <- function(n, sds, cor) {
  q <- length(sds)
  C <- matrix(cor, q, q)
  diag(C) <- 1
  matrix(rnorm(n * q), n, q) %*% (chol(C) * sds)
}
us_of <- function(th, q) as.matrix(us_sigma(th, q))

designs <- list(
  bern2 = function() {
    g <- rep(1:60, each = 5)
    x <- rnorm(300)
    b <- rmv(60, c(1, 0.5), 0.3)
    y <- rbinom(300, 1, plogis(-0.5 + 0.8 * x + b[g, 1] + b[g, 2] * x))
    d <- data.frame(y, x, g = factor(g))
    list(data = d, formula = bf(y ~ x + (x | g)), family = bernoulli(),
         draws = 500L, q = 2L, nodes = 25L,
         eta_fix = function(p, d) p[1] + p[2] * d$x,
         zmat = function(d) cbind(1, d$x),
         dens = function(p, d, eta) dbinom(d$y, 1, plogis(eta), log = TRUE),
         Sigma = function(p) us_of(p[3:5], 2L))
  },
  pois2 = function() {
    g <- rep(1:25, each = 6)
    x <- rnorm(150)
    b <- rmv(25, c(0.7, 0.4), -0.4)
    y <- rpois(150, exp(0.5 + 0.4 * x + b[g, 1] + b[g, 2] * x))
    d <- data.frame(y, x, g = factor(g))
    list(data = d, formula = bf(y ~ x + (x | g)), family = poisson(),
         draws = 300L, q = 2L, nodes = 25L,
         eta_fix = function(p, d) p[1] + p[2] * d$x,
         zmat = function(d) cbind(1, d$x),
         dens = function(p, d, eta) dpois(d$y, exp(eta), log = TRUE),
         Sigma = function(p) us_of(p[3:5], 2L))
  },
  # the lane's hard shape at other seeds: null binary, two rows per group
  null2 = function() {
    g <- rep(1:30, each = 2)
    x <- rnorm(60)
    y <- rbinom(60, 1, 0.4)
    d <- data.frame(y, x, g = factor(g))
    list(data = d, formula = bf(y ~ x + (x | g)), family = bernoulli(),
         draws = 500L, q = 2L, nodes = 31L,
         eta_fix = function(p, d) p[1] + p[2] * d$x,
         zmat = function(d) cbind(1, d$x),
         dens = function(p, d, eta) dbinom(d$y, 1, plogis(eta), log = TRUE),
         Sigma = function(p) us_of(p[3:5], 2L))
  },
  # diagonal block on a probit: non-canonical and q = 2
  probit_diag = function() {
    g <- rep(1:30, each = 6)
    x <- rnorm(180)
    b <- rmv(30, c(0.8, 0.5), 0)
    y <- rbinom(180, 1, pnorm(0.2 + 0.5 * x + b[g, 1] + b[g, 2] * x))
    d <- data.frame(y, x, g = factor(g))
    list(data = d, formula = bf(y ~ x + (x || g)),
         family = bernoulli(link = "probit"), draws = 400L, q = 2L,
         nodes = 25L,
         eta_fix = function(p, d) p[1] + p[2] * d$x,
         zmat = function(d) cbind(1, d$x),
         dens = function(p, d, eta) dbinom(d$y, 1, pnorm(eta), log = TRUE),
         Sigma = function(p) diag(exp(2 * p[3:4])))
  },
  # a block in mu and a block in sigma over the same factor: the scalar rate
  # cannot follow two curvatures that move differently
  musig = function() {
    g <- rep(1:20, each = 6)
    x <- rnorm(120)
    u <- rnorm(20, 0, 0.6)
    v <- rnorm(20, 0, 0.5)
    y <- rnorm(120, 0.2 + 0.5 * x + u[g], exp(-0.3 + v[g]))
    d <- data.frame(y, x, g = factor(g))
    list(data = d, formula = bf(y ~ x + (1 | g), sigma ~ 1 + (1 | g)),
         family = gaussian(), draws = 500L, q = 2L, nodes = 31L,
         eta2 = TRUE,
         dens2 = function(p, d, U) {
           # p: b0 b1 s0 th_mu th_sig; U is m x 2
           m <- nrow(U)
           vapply(seq_len(m), function(i) {
             sum(dnorm(d$y, p[1] + p[2] * d$x + U[i, 1], exp(p[3] + U[i, 2]),
                       log = TRUE))
           }, 0)
         },
         Sigma = function(p) diag(exp(2 * p[4:5])))
  },
  # same, on a Poisson mean and a negative binomial shape
  nb_mushape = function() {
    g <- rep(1:20, each = 6)
    x <- rnorm(120)
    u <- rnorm(20, 0, 0.5)
    v <- rnorm(20, 0, 0.4)
    y <- rnbinom(120, size = exp(0.7 + v[g]), mu = exp(1 + 0.4 * x + u[g]))
    d <- data.frame(y, x, g = factor(g))
    list(data = d, formula = bf(y ~ x + (1 | g), shape ~ 1 + (1 | g)),
         family = negbinomial(), draws = 500L, q = 2L, nodes = 31L,
         eta2 = TRUE,
         dens2 = function(p, d, U) {
           m <- nrow(U)
           vapply(seq_len(m), function(i) {
             sum(dnbinom(d$y, size = exp(p[3] + U[i, 2]),
                         mu = exp(p[1] + p[2] * d$x + U[i, 1]), log = TRUE))
           }, 0)
         },
         Sigma = function(p) diag(exp(2 * p[4:5])))
  },
  # q = 3: the back-substitution's inner loop runs only from here up
  pois3 = function() {
    g <- rep(1:30, each = 8)
    x <- rnorm(240)
    z <- rnorm(240)
    b <- rmv(30, c(0.6, 0.3, 0.3), 0.2)
    y <- rpois(240, exp(0.5 + 0.3 * x - 0.2 * z + b[g, 1] + b[g, 2] * x +
                          b[g, 3] * z))
    d <- data.frame(y, x, z, g = factor(g))
    list(data = d, formula = bf(y ~ x + z + (1 + x + z | g)),
         family = poisson(), draws = 400L, q = 3L, nodes = 13L,
         eta_fix = function(p, d) p[1] + p[2] * d$x + p[3] * d$z,
         zmat = function(d) cbind(1, d$x, d$z),
         dens = function(p, d, eta) dpois(d$y, exp(eta), log = TRUE),
         Sigma = function(p) us_of(p[4:9], 3L))
  }
)
set.seed(seed)
de <- designs[[nm]]()
q <- de$q
dl <- split(de$data, de$data$g)
if (is.null(de$dens2)) {
  zl <- lapply(dl, de$zmat)
  ll_nodes <- function(p, d, zm, U) {
    fx <- de$eta_fix(p, d)
    eta <- fx + zm %*% t(U)  # n x m
    colSums(matrix(de$dens(p, d[rep(seq_len(nrow(d)), nrow(U)), , drop = FALSE],
                           as.vector(eta)), nrow(d), nrow(U)))
  }
} else {
  zl <- lapply(dl, function(d) NULL)
  ll_nodes <- function(p, d, zm, U) de$dens2(p, d, U)
}
g1 <- gh_rule(de$nodes)
grid <- as.matrix(expand.grid(rep(list(seq_len(de$nodes)), q)))
xn <- matrix(g1$x[grid], nrow(grid), q)
lw <- rowSums(matrix(log(g1$w)[grid], nrow(grid), q)) + rowSums(xn^2)
modes <- lapply(dl, function(d) numeric(q))
exact_nll <- function(p, each = FALSE) {
  S <- de$Sigma(p)
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch)) return(Inf)
  Pinv <- chol2inv(ch)
  ld <- 2 * sum(log(diag(ch)))
  lpr <- function(U) {
    -0.5 * rowSums((U %*% Pinv) * U) - 0.5 * q * log(2 * pi) - 0.5 * ld
  }
  res <- numeric(length(dl))
  for (g in seq_along(dl)) {
    d <- dl[[g]]
    zm <- zl[[g]]
    neg <- function(u) {
      U <- matrix(u, 1L, q)
      -(ll_nodes(p, d, zm, U) + lpr(U))
    }
    o <- optim(modes[[g]], neg, method = "BFGS",
               control = list(reltol = 1e-10, maxit = 200))
    modes[[g]] <<- o$par
    h <- optimHess(o$par, neg)
    ei <- tryCatch(chol(h), error = function(e) NULL)
    if (is.null(ei)) {
      Lc <- ch
    } else {
      Lc <- backsolve(ei, diag(q))
    }
    U <- matrix(o$par, nrow(xn), q, byrow = TRUE) + sqrt(2) * xn %*% t(Lc)
    a <- ll_nodes(p, d, zm, U) + lpr(U) + lw
    res[g] <- max(a) + log(sum(exp(a - max(a)))) +
      0.5 * q * log(2) + sum(log(abs(diag(Lc))))
  }
  if (each) return(res)
  -sum(res)
}
num_grad <- function(f, x, h = 1e-4) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, 0)
}
warn <- character(0)
t0 <- Sys.time()
fit <- withCallingHandlers(
  tryCatch(frm(de$formula, family = de$family, data = de$data,
               importance = de$draws),
           error = function(e) e),
  warning = function(w) {
    warn <<- c(warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
tfit <- as.numeric(Sys.time() - t0, units = "secs")
if (inherits(fit, "error")) {
  cat("ERROR:", conditionMessage(fit), "\n")
  cat(sprintf("%s\t%d\tERROR\t%s\n", nm, seed,
              substr(conditionMessage(fit), 1, 400)),
      file = file.path(out, sprintf("%s-%d.tsv", nm, seed)))
  quit("no")
}
par <- unname(fit$opt$par)
cat("names:", names(fit$opt$par), "\npar:", format(par, digits = 5), "\n")
lap <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
ref <- nlminb(par, exact_nll, gradient = function(p) num_grad(exact_nll, p),
              control = list(iter.max = 100, rel.tol = 1e-9))
pe <- ref$par
He <- numDeriv::jacobian(function(p) num_grad(exact_nll, p), pe,
                         method = "simple", method.args = list(eps = 1e-3))
He <- 0.5 * (He + t(He))
Ve <- tryCatch(solve(He), error = function(e) matrix(NA, length(pe), length(pe)))
se_exact <- sqrt(abs(diag(Ve)))
sdr <- tryCatch(frmtmb:::sdr_of(fit), error = function(e) NULL)
se_fit <- if (!is.null(sdr)) sqrt(diag(sdr$cov.fixed)) else rep(NA, length(pe))
g_at_fit <- num_grad(exact_nll, par)
newton <- sqrt(abs(sum(g_at_fit * (Ve %*% g_at_fit))))
res <- data.frame(design = nm, seed = seed, draws = de$draws,
                  rounds = fit$importance$rounds, grad = fit$importance$grad,
                  capped = fit$importance$capped,
                  ess_min = fit$importance$ess_min, mcse = fit$importance$mcse,
                  dll = -fit$opt$objective + ref$objective,
                  loss = exact_nll(par) - ref$objective,
                  lap_loss = exact_nll(unname(lap$opt$par)) - ref$objective,
                  max_z_est = max(abs(par - pe) / se_exact),
                  max_z_lap = max(abs(unname(lap$opt$par) - pe) / se_exact),
                  newton = newton,
                  se_ratio_min = min(se_fit / se_exact, na.rm = TRUE),
                  se_ratio_max = max(se_fit / se_exact, na.rm = TRUE),
                  ref_code = ref$convergence,
                  n_warn = length(warn), tfit = tfit,
                  warn = substr(paste(warn, collapse = " | "), 1, 160))
print(t(res))
cat("(est-ref)/se:", format((par - pe) / se_exact, digits = 3), "\n")
cat("(lap-ref)/se:", format((unname(lap$opt$par) - pe) / se_exact, digits = 3),
    "\n")
cat("se_fit/se:", format(se_fit / se_exact, digits = 4), "\n")
write.table(res, file.path(out, sprintf("%s-%d.tsv", nm, seed)), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Regression of the draws on z at a moved theta: m and L as the draws have them,
# against the log-determinant on the tape and against the documented H(theta).
ns <- asNamespace("frmtmb")
prop <- ns$imp_frozen_proposal(fit)
plan <- prop$plan
io <- prop$io
pe2 <- environment(io$fn)
mv <- par + rep_len(c(0.5, -0.5), length(par)) * pmin(se_exact, 0.3)
pl <- ns$imp_par_list(prop$template, setNames(mv, names(fit$opt$par)), prop$map)
pr <- pe2$proposal(pl)
ng <- plan$n_group
nd <- plan$n_draw
err_ld <- 0
for (g in unique(round(seq(1, ng, length.out = 8)))) {
  U <- vapply(pr$u, function(uj) as.numeric(uj)[(seq_len(nd) - 1L) * ng + g],
              numeric(nd))
  Z <- t(plan$z[plan$idx[, g], , drop = FALSE])
  fit_u <- lm(U ~ Z)
  L <- t(coef(fit_u)[-1L, , drop = FALSE])
  err_ld <- max(err_ld, abs(log(abs(det(L))) - as.numeric(pr$logdet)[g]))
  stopifnot(max(abs(residuals(fit_u))) < 1e-8)
}
cat("max |log det L from the draws - tape logdet| at moved theta:",
    format(err_ld, digits = 3), "\n")
# joint density at the moved draws against the plain objective
a <- as.matrix(io$amat(pl))
nll <- ns$build_objective(fit$frame)
err_j <- 0
for (j in c(1L, nd)) {
  joint <- pl
  b <- numeric(length(plan$uhat))
  for (qq in seq_len(q)) {
    b[plan$idx[qq, ]] <- as.numeric(pr$u[[qq]])[(j - 1L) * ng + seq_len(ng)]
  }
  joint[["b"]] <- b
  err_j <- max(err_j, abs(sum(a[, j] - plan$zsq[, j]) + nll(joint)))
}
cat("max |weights' joint - plain joint| at moved theta:",
    format(err_j, digits = 3), "\n")

# away from the optimum: value and gradient
nmv <- names(fit$opt$par)
pts <- list(opt = par)
for (i in seq_along(par)) {
  for (k in c(-3, 3)) {
    pts[[paste0("c", i, "_", k)]] <- replace(par, i, par[i] + k * se_exact[i])
  }
}
rows <- lapply(names(pts), function(pn) {
  p <- pts[[pn]]
  v <- fit$obj$fn(setNames(p, nmv))
  e <- exact_nll(p)
  eg <- ns$imp_ess_at(prop, setNames(p, nmv))
  ad <- as.numeric(fit$obj$gr(setNames(p, nmv)))
  fd <- num_grad(function(z) fit$obj$fn(setNames(z, nmv)), p, h = 1e-5)
  data.frame(design = nm, seed = seed, pt = pn, fn = v, exact = e,
             diff = v - e, mcse_at = eg$mcse, ess_min_at = min(eg$ess),
             z = (v - e) / eg$mcse,
             ad_vs_fd = max(abs(ad - fd)) / max(1e-6, max(abs(fd))))
})
rows <- do.call(rbind, rows)
print(format(rows[, -(1:2)], digits = 4))
write.table(rows, file.path(out, sprintf("%s-%d-away.tsv", nm, seed)),
            sep = "\t", quote = FALSE, row.names = FALSE)
