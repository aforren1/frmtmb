# Lane boundfam: what se_edge_limit()'s probes read on one fit.
#   Rscript dev/boundfam-probe.R <lib or "base"> <case>
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
fit_case <- function(k) {
  if (k == "wgam3" || k == "wgam1") {
    s <- if (k == "wgam3") 3 else 1
    set.seed(s)
    n <- 60
    d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
    d$y <- 1 + d$x + rnorm(n, 0, 0.5) + rnorm(n, 0, 0.5)
    w <- runif(n, 0.5, 2)
    d$w <- w / mean(w)
    d$yg <- rgamma(n, shape = 5, rate = 5 / exp(0.5 + 0.3 * d$x +
                                                  rnorm(n, 0, 0.3)))
    return(suppressWarnings(frm(yg | weights(w) ~ x + (1 | obs),
                                family = Gamma("log"), data = d)))
  }
  if (k == "eam") {
    suppressPackageStartupMessages(library(frmtmb.eam))
    set.seed(3)
    dat <- ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
    kk <- sample(300, 18)
    dat$rt[kk] <- stats::runif(18, 0.12, 2.5)
    return(suppressWarnings(frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
      family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                       lognormal()), data = dat)))
  }
  if (k == "maling") {
    d <- data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30, 20,
                          6, 44, 44, 27, 25, 17, 14, 27, 35, 30), n = 45L)
    return(frm(bf(k | trials(n) ~ 1),
               family = mixture(beta_binomial, beta_binomial), data = d,
               start = list(beta = qlogis(c(0.55, 0.95)),
                            betad = c(log(c(15, 15)), 0))))
  }
  stop("unknown case")
}
fit <- fit_case(a[2])
p <- fit$opt$par
np <- length(p)
nm <- ns$outer_par_names(fit)
u <- rep_len(fit$par_units %||% 1, np)
h <- ns$fit_outer_hessian(fit)
H <- h$H
cat("code", fit$opt$convergence, " grad", signif(drop(fit$obj$gr(p)), 3),
    "\n")
cat("diag H * u^2:", signif(abs(diag(H)) * u^2, 3), "\n")
cand <- ns$se_limit_candidates(fit)
cat("candidates:", nm[cand], "\n")
f0 <- fit$obj$fn(p)
for (j in cand) {
  at <- function(t) {
    q <- p
    q[j] <- p[j] + t * u[j]
    as.numeric(fit$obj$fn(q)) - f0
  }
  ts <- c(-12, -2, 2, 12)
  cat(sprintf("%-22s est %.4g  dobj at %s: %s\n", nm[j], p[j],
              paste(ts, collapse = ","),
              paste(signif(vapply(ts, at, 0), 4), collapse = " ")))
  cat("   back to zero:", signif(at(-p[j] / u[j]), 4), " dir",
      ns$se_limit_dir(at, p[j] / u[j], 1e-3), "\n")
}
lim <- ns$se_edge_limit(fit, H)
cat("se_edge_limit:", nm[lim], "\n")
print(ns$se_limit_info(fit))
li <- ns$se_limit_info(fit)
if (length(li$dir)) {
  q <- p
  j <- match(names(li$dir), nm)
  q[j] <- p[j] + li$dir * 2 * u[j]
  rest <- setdiff(seq_len(np), j)
  c0 <- ns$se_limit_curv(fit, p, rest, u)
  c1 <- ns$se_limit_curv(fit, q, rest, u)
  print(data.frame(par = nm[rest], c0 = signif(c0, 6), c1 = signif(c1, 6),
                   ratio = signif(c1 / c0, 8)))
}
