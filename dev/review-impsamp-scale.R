# A fit is equivariant under the scale and shift of a covariate. The fixed
# absolute step of imp_anchor_slopes() is not: h = 1e-4 on the coefficient of a
# column of scale s moves the predictor by 1e-4 * s. Reference: an independent
# exact optimum on the standard column, mapped to the column's coefficients.
# Usage: REV_ARM=after|before Rscript dev/review-impsamp-scale.R <family> <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
args <- commandArgs(TRUE)
fam <- args[[1]]
seed <- as.integer(args[[2]])
draws <- if (length(args) >= 3) as.integer(args[[3]]) else 300L
set.seed(seed)
G <- 15L
g <- factor(rep(seq_len(G), each = 10))
x <- rnorm(150)
u <- rnorm(G, 0, 0.8)
if (fam == "poisson") {
  y <- rpois(150, exp(0.5 + 0.5 * x + u[g]))
  family <- poisson()
  dens <- function(y, eta) dpois(y, exp(eta), log = TRUE)
} else {
  y <- rbinom(150, 1, plogis(0.2 + 0.7 * x + u[g]))
  family <- bernoulli()
  dens <- function(y, eta) dbinom(y, 1, plogis(eta), log = TRUE)
}
gh <- gh_rule(61L)
dl <- split(data.frame(y, x, g), g)
exact_nll <- function(p) {
  sd <- exp(p[3])
  tot <- 0
  for (d in dl) {
    fn <- function(v) {
      vapply(v, function(vi) sum(dens(d$y, p[1] + p[2] * d$x + vi)), 0) +
        dnorm(v, 0, sd, log = TRUE)
    }
    o <- optimize(function(v) -fn(v), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
    m <- o$minimum
    h <- 1e-3
    cur <- -(fn(m + h) - 2 * fn(m) + fn(m - h)) / h^2
    s <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
    nodes <- m + sqrt(2) * s * gh$x
    a <- fn(nodes) + gh$x^2 + log(gh$w)
    tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * s))
  }
  tot
}
ng <- function(f, p, h = 1e-4) {
  vapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, h)
    (f(p + e) - f(p - e)) / (2 * h)
  }, 0)
}
d0 <- data.frame(y, x, g)
l0 <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = family, data = d0))
ref <- nlminb(unname(l0$opt$par), exact_nll, gradient = function(p) ng(exact_nll, p),
              control = list(rel.tol = 1e-10))
pe <- ref$par
se <- sqrt(diag(solve(numDeriv::jacobian(function(p) ng(exact_nll, p), pe,
                                         method = "simple"))))
cat("exact:", format(pe, digits = 6), "se", format(se, digits = 3), "nll",
    format(ref$objective, digits = 9), "\n")
cases <- list(
  s1 = c(1, 0), s1e2 = c(1e2, 0), s1e3 = c(1e3, 0), s1e4 = c(1e4, 0),
  s1e5 = c(1e5, 0), s1e6 = c(1e6, 0), s1em3 = c(1e-3, 0),
  shift100 = c(1, 100), shift2000 = c(1, 2000), shift1e5 = c(1, 1e5),
  s100_shift = c(100, 1e4))
rows <- list()
for (cn in names(cases)) {
  a <- cases[[cn]][1]
  b <- cases[[cn]][2]
  d <- data.frame(y, xx = a * x + b, g)
  w <- character(0)
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    tryCatch(frm(bf(y ~ xx + (1 | g)), family = family, data = d,
                 importance = draws),
             error = function(e) e),
    warning = function(x2) {
      w <<- c(w, conditionMessage(x2))
      invokeRestart("muffleWarning")
    })
  tt <- as.numeric(Sys.time() - t0, units = "secs")
  lp <- suppressWarnings(frm(bf(y ~ xx + (1 | g)), family = family, data = d))
  back <- function(p) c(p[1] + p[2] * b, p[2] * a, p[3])
  zl <- (back(unname(lp$opt$par)) - pe) / se
  if (inherits(fit, "error")) {
    rows[[cn]] <- data.frame(arm, fam, seed, case = cn, rounds = NA, moved1 = NA,
                             grad = NA, conv = NA, zmax = NA, ztheta = NA,
                             zlap_theta = zl[3], lossexact = NA, ll_diff = NA,
                             mcse = NA, ess = NA, nw = NA, secs = tt,
                             msg = substr(conditionMessage(fit), 1, 120))
    next
  }
  p <- back(unname(fit$opt$par))
  z <- (p - pe) / se
  rows[[cn]] <- data.frame(arm, fam, seed, case = cn,
                           rounds = fit$importance$rounds,
                           moved1 = fit$importance$moves[[1]],
                           grad = fit$importance$grad,
                           conv = fit$opt$convergence,
                           zmax = max(abs(z)), ztheta = z[3],
                           zlap_theta = zl[3],
                           lossexact = exact_nll(p) - ref$objective,
                           ll_diff = ref$objective - fit$opt$objective,
                           mcse = fit$importance$mcse,
                           ess = fit$importance$ess_min, nw = length(w),
                           secs = tt,
                           msg = substr(paste(w, collapse = " | "), 1, 120))
}
r <- do.call(rbind, rows)
options(width = 200)
print(format(r[, -(1:3)], digits = 3), row.names = FALSE)
write.table(r, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/scale-%s-%s-%d.tsv",
  arm, fam, seed), sep = "\t", quote = FALSE, row.names = FALSE)
