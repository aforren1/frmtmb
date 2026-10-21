# Densities of negbinomial2 and dirichlet_multinomial against brms
# 2.23.0's Stan code, compiled through rstan::expose_stan_functions(),
# and against an independent reference, over a grid with small and
# large sigma / phi. frmtmb is evaluated twice: on doubles (row_lpdf()),
# and through an RTMB tape with the linear predictors riding along as
# the objective hands them (.eta_<dpar>), which is what a fit sees.
#   Rscript dev/nb2dm-density.R > dev/nb2dm-out/density.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(brms)
})
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "brms", format(packageVersion("brms")),
    "rstan", format(packageVersion("rstan")), "\n")

# brms's own Stan function for the dirichlet-multinomial, cut out of the
# program brms writes, so the reference is brms's text and not a copy
set.seed(1)
dd <- data.frame(x = rnorm(6))
Ym <- matrix(1L, 6, 3, dimnames = list(NULL, c("a", "b", "c")))
dd$Y <- Ym
dd$tr <- 3L
sc <- stancode(Y | trials(tr) ~ x, data = dd, family = dirichlet_multinomial())
lines <- strsplit(sc, "\n")[[1]]
i0 <- grep("real dirichlet_multinomial_logit2_lpmf", lines)
i1 <- i0 + which(lines[(i0 + 1):length(lines)] == "  }")[1]
dm_fun <- paste(lines[i0:i1], collapse = "\n")
# brms's likelihood statement for negbinomial2, read from its program
sc2 <- stancode(bf(y ~ x, sigma ~ x), data = data.frame(y = 0:5, x = 1:6),
                family = brmsfamily("negbinomial2"))
cat("brms negbinomial2 statement:",
    grep("target \\+= neg_binomial", strsplit(sc2, "\n")[[1]], value = TRUE),
    "\n")
sc3 <- stancode(y | rate(d) ~ x,
                data = data.frame(y = 0:5, x = 1:6, d = 1:6),
                family = brmsfamily("negbinomial2"))
cat("brms negbinomial2 rate statement:",
    grep("target \\+= neg_binomial", strsplit(sc3, "\n")[[1]], value = TRUE),
    "\n")
code <- paste0("functions {\n", dm_fun, "\n",
  "  real nb2_ref(int y, real eta, real sigma) {\n",
  "    return neg_binomial_2_log_lpmf(y | eta, inv(sigma));\n  }\n",
  "  real nb2_rate_ref(int y, real eta, real sigma, real denom) {\n",
  "    return neg_binomial_2_log_lpmf(y | eta + log(denom), ",
  "inv(sigma) .* denom);\n  }\n",
  "  real dm_ref(array[] int y, vector mu, real phi) {\n",
  "    return dirichlet_multinomial_logit2_lpmf(y | mu, phi);\n  }\n",
  "}\n")
cat(code)
rstan::expose_stan_functions(rstan::stanc(model_code = code))

ulp_rel <- function(a, b) {
  ifelse(a == b, 0, abs(a - b) / pmax(abs(b), .Machine$double.xmin))
}

## negbinomial2 --------------------------------------------------------
fam <- frmtmb::brmsfamily("negbinomial2")
g <- expand.grid(y = c(0, 1, 2, 5, 10, 50, 200, 1000),
                 eta = c(-6, -1, 0, 1, 3, 7),
                 sigma = c(1e-8, 1e-4, 0.01, 0.3, 1, 10, 1e3, 1e6))
stan_nb <- mapply(function(y, e, s) nb2_ref(as.integer(y), e, s),
                  g$y, g$eta, g$sigma)
ref_nb <- dnbinom(g$y, size = 1 / g$sigma, mu = exp(g$eta), log = TRUE)
dbl_nb <- as.numeric(frmtmb::row_lpdf(fam, g$y, g$y,
                              list(mu = exp(g$eta), sigma = g$sigma),
                              list(), NULL))
# the tape, with the linear predictors the objective supplies
tp <- RTMB::MakeTape(function(p) {
  n <- length(g$y)
  eta <- p[1:n]
  ls <- p[n + 1:n]
  dp <- list(mu = exp(eta), sigma = exp(ls), .eta_mu = eta, .eta_sigma = ls)
  frmtmb::row_lpdf(fam, g$y, g$y, dp, list(), NULL)
}, c(g$eta, log(g$sigma)))
tap_nb <- as.numeric(tp(c(g$eta, log(g$sigma))))
cat("\nnegbinomial2:", nrow(g), "grid points\n")
fin <- is.finite(stan_nb)
cat("  finite: stan", sum(fin), "dnbinom", sum(is.finite(ref_nb)),
    "frmtmb doubles", sum(is.finite(dbl_nb)), "tape", sum(is.finite(tap_nb)),
    "\n")
for (nm in c("dbl_nb", "tap_nb")) {
  v <- get(nm)
  r1 <- ulp_rel(v[fin], stan_nb[fin])
  r2 <- ulp_rel(v, ref_nb)
  cat(sprintf("  %s vs stan: max rel %.3g (at %s); vs dnbinom: max rel %.3g\n",
              nm, max(r1), paste(g[fin, ][which.max(r1), ], collapse = "/"),
              max(r2[is.finite(r2)])))
}
r_sd <- ulp_rel(stan_nb[fin], ref_nb[fin])
cat(sprintf("  stan vs dnbinom: max rel %.3g (at %s)\n", max(r_sd),
            paste(g[fin, ][which.max(r_sd), ], collapse = "/")))
# where stan and dnbinom disagree, which one does frmtmb agree with
bad <- which(fin & r_sd > 1e-10)
if (length(bad)) {
  cat("  rows where stan and dnbinom differ by more than 1e-10:\n")
  print(cbind(g[bad, ], stan = stan_nb[bad], dnbinom = ref_nb[bad],
              frm_dbl = dbl_nb[bad], frm_tape = tap_nb[bad]), digits = 15)
}
# an exact reference: the rising factorial as a sum of logs and the
# log(r + mu) split so that no term cancels; it decides between stan
# and dnbinom where they disagree
nb_exact <- function(y, eta, sigma) {
  r <- 1 / sigma
  mu <- exp(eta)
  lrm <- if (r > mu) log(r) + log1p(mu / r) else eta + log1p(r / mu)
  a <- if (r > mu) -r * log1p(mu / r) else r * (log(r) - lrm)
  rf <- if (y > 0) sum(log(r + seq_len(y) - 1)) else 0
  rf - lgamma(y + 1) + a + y * (eta - lrm)
}
ex_nb <- mapply(nb_exact, g$y, g$eta, g$sigma)
cat(sprintf("  vs the exact reference (max rel): frmtmb tape %.3g, stan %.3g, dnbinom %.3g\n",
            max(ulp_rel(tap_nb, ex_nb)), max(ulp_rel(stan_nb, ex_nb)),
            max(ulp_rel(ref_nb, ex_nb))))
for (s in sort(unique(g$sigma))) {
  k <- g$sigma == s
  cat(sprintf("    sigma %-6g max rel to exact: frmtmb %.2g  stan %.2g  dnbinom %.2g\n",
              s, max(ulp_rel(tap_nb[k], ex_nb[k])),
              max(ulp_rel(stan_nb[k], ex_nb[k])),
              max(ulp_rel(ref_nb[k], ex_nb[k]))))
}
# the tape's derivative in log(sigma) against the analytic one,
# -r dL/dr with r = 1 / sigma, where the analytic form does not cancel
J <- tp$jacobian(c(g$eta, log(g$sigma)))
n <- nrow(g)
dls <- diag(J[, n + 1:n])
r <- 1 / g$sigma
mu <- exp(g$eta)
an <- -r * (digamma(g$y + r) - digamma(r) + log(r) - log(r + mu) +
              (mu - g$y) / (r + mu))
# a second reference: a Richardson central difference of the exact
# value, which does not cancel the way the digamma form does at large r
ex_d <- function(h) {
  (mapply(nb_exact, g$y, g$eta, g$sigma * exp(h)) -
     mapply(nb_exact, g$y, g$eta, g$sigma * exp(-h))) / (2 * h)
}
rich <- (4 * ex_d(1e-3) - ex_d(2e-3)) / 3
ok <- g$sigma >= 1e-4 & g$sigma <= 1e3 & abs(an) > 1e-8
rel <- abs(dls[ok] - an[ok]) / abs(an[ok])
rel2 <- abs(dls[ok] - rich[ok]) / abs(rich[ok])
w <- which(ok)[which.max(rel)]
cat(sprintf("  d/dlog(sigma), %d points (sigma in [1e-4, 1e3]): tape vs digamma form max rel %.3g, tape vs Richardson %.3g\n",
            sum(ok), max(rel), max(rel2)))
cat(sprintf("    worst digamma point y %g eta %g sigma %g: tape %.10g digamma %.10g Richardson %.10g\n",
            g$y[w], g$eta[w], g$sigma[w], dls[w], an[w], rich[w]))
cat("  finite tape derivatives over the whole grid:", sum(is.finite(dls)),
    "of", n, "\n")
# rate(): brms multiplies mu and the shape by the exposure
gr <- expand.grid(y = c(0, 3, 40), eta = c(-2, 0, 2), sigma = c(0.01, 0.5, 5),
                  d = c(0.1, 1, 7.5))
stan_r <- mapply(function(y, e, s, d) nb2_rate_ref(as.integer(y), e, s, d),
                 gr$y, gr$eta, gr$sigma, gr$d)
frm_r <- as.numeric(frmtmb::row_lpdf(fam, gr$y, gr$y,
                             list(mu = exp(gr$eta), sigma = gr$sigma),
                             list(rate = gr$d), NULL))
cat(sprintf("  rate(): %d points, max rel vs stan %.3g\n", nrow(gr),
            max(ulp_rel(frm_r, stan_r))))

## dirichlet_multinomial ----------------------------------------------
dm_exact <- function(y, alpha) {
  # sum of logs of rising factorials: no lgamma difference, so exact to
  # rounding at any precision
  phi <- sum(alpha)
  T <- sum(y)
  out <- lgamma(T + 1) - sum(lgamma(y + 1))
  for (k in seq_along(y)) {
    if (y[k] > 0) out <- out + sum(log(alpha[k] + seq_len(y[k]) - 1))
  }
  out - sum(log(phi + seq_len(T) - 1))
}
set.seed(3)
cases <- list()
for (K in c(2L, 3L, 5L)) {
  for (phi in c(1e-4, 0.05, 1, 8, 300, 1e5, 1e8)) {
    for (r in 1:6) {
      T <- sample(c(0, 1, 7, 40, 500), 1)
      mu <- c(0, rnorm(K - 1, 0, 2))
      p <- exp(mu - max(mu)); p <- p / sum(p)
      y <- as.vector(rmultinom(1, T, p))
      cases[[length(cases) + 1L]] <- list(K = K, phi = phi, mu = mu, y = y)
    }
  }
}
dm_res <- do.call(rbind, lapply(cases, function(cs) {
  fm <- frmtmb::dirichlet_multinomial(K = cs$K)
  dpn <- fm$primary_dpars
  dp <- c(stats::setNames(as.list(cs$mu[-1]), dpn), list(phi = cs$phi))
  ym <- matrix(cs$y, 1)
  dbl <- as.numeric(frmtmb::row_lpdf(fm, ym, ym, dp, list(trials = sum(cs$y)), NULL))
  tpe <- RTMB::MakeTape(function(p) {
    dpp <- c(stats::setNames(as.list(p[seq_along(dpn)]), dpn),
             list(phi = exp(p[length(p)]), .eta_phi = p[length(p)]))
    frmtmb::row_lpdf(fm, ym, ym, dpp, list(trials = sum(cs$y)), NULL)
  }, c(cs$mu[-1], log(cs$phi)))
  tap <- as.numeric(tpe(c(cs$mu[-1], log(cs$phi))))
  stan <- dm_ref(as.integer(cs$y), cs$mu, cs$phi)
  p <- exp(cs$mu - max(cs$mu)); p <- p / sum(p)
  ex <- dm_exact(cs$y, cs$phi * p)
  data.frame(K = cs$K, phi = cs$phi, T = sum(cs$y), stan = stan, exact = ex,
             dbl = dbl, tape = tap)
}))
cat("\ndirichlet_multinomial:", nrow(dm_res), "cases\n")
cat(sprintf("  frmtmb doubles vs stan: max rel %.3g; tape vs stan: %.3g\n",
            max(ulp_rel(dm_res$dbl, dm_res$stan)),
            max(ulp_rel(dm_res$tape, dm_res$stan))))
cat(sprintf("  frmtmb doubles vs exact: max abs %.3g; stan vs exact: max abs %.3g\n",
            max(abs(dm_res$dbl - dm_res$exact)),
            max(abs(dm_res$stan - dm_res$exact))))
cat("  by phi, max abs error against the exact reference:\n")
print(aggregate(cbind(frm = abs(dbl - exact), stan = abs(stan - exact)) ~ phi,
                data = dm_res, FUN = max), digits = 3)
cat("  identical(frmtmb doubles, stan):", sum(dm_res$dbl == dm_res$stan),
    "of", nrow(dm_res), "\n")
# K = 2 is the beta-binomial
bb <- frmtmb::brmsfamily("beta_binomial")
k2 <- dm_res$K == 2
cat("  K = 2 equals beta_binomial():\n")
k2c <- cases[k2]
bbv <- vapply(k2c, function(cs) {
  p <- exp(cs$mu[2]) / (1 + exp(cs$mu[2]))
  as.numeric(frmtmb::row_lpdf(bb, cs$y[2], cs$y[2], list(mu = p, phi = cs$phi),
                      list(trials = sum(cs$y)), NULL))
}, 0)
cat(sprintf("    max abs diff %.3g over %d cases\n",
            max(abs(bbv - dm_res$dbl[k2])), sum(k2)))
