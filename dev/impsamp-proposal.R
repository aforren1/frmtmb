# Gaussian or t? At a fit's estimates, per group, with the package's own
# mode and Hessian replaced by the independent reference's: the worst and
# the median group's ESS fraction at 512 antithetic draws, and the sum of
# the groups' error in log m_g against adaptive Gauss-Hermite, over draw
# seeds 1:30. Proposals: the Laplace Gaussian, and a multivariate t with
# that mode and scale at 4 and 16 degrees of freedom.
# Log: dev/impsamp-out/proposal.txt
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("dev/impsamp-agh.R")
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)

group_setup <- function(family, shape, d, par) {
  nth <- switch(shape, int = 1L, slope = 3L, twoblock = 2L)
  beta <- par[1:2]
  th <- par[2L + seq_len(nth)]
  ld <- if (agh_has_disp(family)) par[[3L + nth]] else 0
  Sigma <- agh_Sigma(shape, th)
  etaf <- beta[[1L]] + beta[[2L]] * d$x
  q <- nrow(Sigma)
  Pinv <- solve(Sigma)
  ldet <- as.numeric(determinant(Sigma, TRUE)$modulus)
  lapply(levels(d$g), function(lv) {
    w <- d$g == lv
    Zg <- if (shape == "int") matrix(1, sum(w), 1) else cbind(1, d$x[w])
    # vectorized over draws: U is q x K
    fU <- function(U) {
      eta <- etaf[w] + Zg %*% U
      -colSums(matrix(agh_logdens(family, rep(d$y[w], ncol(U)), as.vector(eta),
                                  rep(d$ntr[w], ncol(U)), ld),
                      sum(w), ncol(U))) +
        0.5 * colSums(U * (Pinv %*% U)) + 0.5 * ldet + 0.5 * q * log(2 * pi)
    }
    f1 <- function(u) fU(matrix(u, q, 1L))
    u <- optim(rep(0, q), f1, method = "BFGS",
               control = list(reltol = 1e-12, maxit = 500))$par
    h <- 1e-4
    H <- matrix(0, q, q)
    for (i in seq_len(q)) for (j in seq_len(i)) {
      ei <- replace(numeric(q), i, h)
      ej <- replace(numeric(q), j, h)
      H[i, j] <- H[j, i] <- (f1(u + ei + ej) - f1(u + ei - ej) -
                               f1(u - ei + ej) + f1(u - ei - ej)) / (4 * h^2)
    }
    R <- chol(H)
    list(fU = fU, u = u, L = backsolve(R, diag(q)), q = q)
  })
}

log_m <- function(gs, N, seed, nu, alpha = 1) {
  # alpha: the share of pairs drawn from the t; 1 is a plain t, and with
  # nu = Inf a plain Gaussian. The density of the mixture weighs every draw.
  set.seed(seed)
  q <- gs$q
  K <- N %/% 2L
  z <- matrix(rnorm(q * K), q, K)
  z <- cbind(z, -z)
  if (is.finite(nu)) {
    nt <- ceiling(alpha * K)
    s <- c(rchisq(nt, nu) / nu, rep(1, K - nt))
    s <- c(s, s)
    v <- sweep(z, 2L, sqrt(s), "/")
    logt <- lgamma((nu + q) / 2) - lgamma(nu / 2) - 0.5 * q * log(nu * pi) -
      0.5 * (nu + q) * log1p(colSums(v^2) / nu)
    logn <- -0.5 * q * log(2 * pi) - 0.5 * colSums(v^2)
    logq <- if (alpha >= 1) logt else
      log(alpha * exp(logt) + (1 - alpha) * exp(logn))
  } else {
    v <- z
    logq <- -0.5 * q * log(2 * pi) - 0.5 * colSums(v^2)
  }
  U <- gs$u + gs$L %*% v
  lw <- -gs$fU(U) - logq + sum(log(diag(gs$L)))
  mx <- max(lw)
  w <- exp(lw - mx)
  c(lm = mx + log(mean(w)), ess = sum(w)^2 / sum(w^2) / N)
}

cases <- list()
d1 <- E$eq_dat("binomial", seed = 124)
cases[["p124 binomial"]] <- list(d = d1, family = "binomial", shape = "int",
                                 form = bf(y | trials(ntr) ~ x + (1 | g)),
                                 fam = binomial())
set.seed(3)
ng <- 40L
g <- factor(rep(seq_len(ng), each = 3L))
x <- rnorm(ng * 3L)
u <- rnorm(ng, 0, 2)
d2 <- data.frame(y = rbinom(ng * 3L, 1, plogis(-0.5 + 0.5 * x + u[g])),
                 x = x, g = g, ntr = 1L)
cases[["hard bernoulli int, sd 2"]] <- list(d = d2, family = "bernoulli",
  shape = "int", form = bf(y ~ x + (1 | g)), fam = bernoulli())
set.seed(5)
ng <- 30L
g <- factor(rep(seq_len(ng), each = 2L))
x <- rnorm(2L * ng)
d4 <- data.frame(y = rbinom(2L * ng, 1, 0.4), x = x, g = g, ntr = 1L)
cases[["bernoulli slope, 30 groups of 2"]] <- list(d = d4, family = "bernoulli",
  shape = "slope", form = bf(y ~ x + (x | g)), fam = bernoulli())
d5 <- impsamp_data("poisson", "slope", "30x3", 2)
cases[["poisson slope 30x3"]] <- list(d = d5, family = "poisson", shape = "slope",
  form = bf(y ~ x + (1 + x | g)), fam = poisson())
d6 <- impsamp_data("gaussian", "int", "30x3", 2)
cases[["gaussian int 30x3"]] <- list(d = d6, family = "gaussian", shape = "int",
  form = bf(y ~ x + (1 | g)), fam = gaussian())
N <- 512L
for (nm in names(cases)) {
  cc <- cases[[nm]]
  fit <- suppressWarnings(suppressMessages(frm(cc$form, family = cc$fam,
                                               data = cc$d, importance = 2000L)))
  par <- impsamp_start(cc$family, cc$shape, fit)
  gs <- group_setup(cc$family, cc$shape, cc$d, par)
  exact <- agh_ll(cc$family, cc$shape, cc$d, par,
                  if (cc$shape == "int") 25L else 15L)
  for (cfg in list(c(Inf, 1), c(16, 1), c(4, 1), c(4, 0.1), c(4, 0.25), c(16, 0.25))) {
    nu <- cfg[[1L]]
    alpha <- cfg[[2L]]
    r <- vapply(1:30, function(s) {
      v <- vapply(gs, log_m, c(lm = 0, ess = 0), N = N, seed = s, nu = nu, alpha = alpha)
      c(err = sum(v["lm", ]) - exact, essmin = min(v["ess", ]),
        essmed = median(v["ess", ]))
    }, c(err = 0, essmin = 0, essmed = 0))
    cat(sprintf(
      "%-34s nu=%-4s a=%-4s err mean %+.4f sd %.4f | ESS min(median over seeds) %.3f, worst seed %.3f, group median %.3f\n",
      nm, format(nu), format(alpha), mean(r["err", ]), sd(r["err", ]), median(r["essmin", ]),
      min(r["essmin", ]), median(r["essmed", ])))
  }
}
