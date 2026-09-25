# Lane thres: validation of the thres() addition term.
#
#   FRMTMB_LIB=/opt/rlib/lane-thres Rscript dev/thres-validate.R
#
# 1. brms's density at a shared parameter point. For each of the four
#    ordinal families, grouped thresholds with a different number of
#    categories per group, frmtmb's objective at a chosen point against
#    the sum of brms's own R densities (brms:::dcumulative, dsratio,
#    dcratio, dacat, the functions posterior_epred() uses), each row
#    read through its group's slice of the thresholds, as brms's Stan
#    code does with Jthres. Links: logit and probit (cauchit for
#    cumulative too, which takes the plain-difference branch).
# 2. MASS::polr at the ML optimum: with a group-specific slope the
#    grouped model factorizes into one polr fit per group, so the sum of
#    their log-likelihoods is the grouped fit's. Groups have 3, 4 and 5
#    categories.
# 3. ordinal::clm with nominal = ~ g at the ML optimum: a shared slope
#    and a threshold vector per group, the same model when every group
#    has the same categories.
# 4. thres(x = K) without groups: unobserved top categories, compared
#    with brms's density at a point.
# Seeds are fixed below; every number printed is full precision.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-thres")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

sim_data <- function(seed, n = 240) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
  u <- rlogis(n, 0.7 * d$x)
  tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6), c = c(-1.5, -0.4, 0.5,
                                                          1.6))
  d$y <- vapply(seq_len(n), function(i) {
    1L + sum(u[i] > tau[[d$g[i]]])
  }, 1L)
  d
}

# brms's density for one row, from brms's own R functions
brms_row_lpmf <- function(fam, link, y, eta, thres) {
  f <- get(paste0("d", fam), envir = asNamespace("brms"))
  log(f(y, eta = eta, thres = matrix(thres, 1L), disc = 1, link = link))
}

# a parameter point: thresholds per group, a slope
pt_thres <- list(a = c(-0.8, 0.1, 1.1), b = c(-0.3, 0.9),
                 c = c(-1.2, -0.3, 0.4, 1.9))
pt_b <- 0.55
raw_of <- function(th, ordered) {
  unlist(lapply(th, function(t) if (ordered) c(t[1], log(diff(t))) else t))
}

cat("\n== 1. frmtmb objective vs brms densities at a shared point ==\n")
d <- sim_data(11)
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  links <- switch(fam, cumulative = c("logit", "probit", "cauchit"),
                  acat = "logit", c("logit", "probit"))
  for (lk in links) {
    ff <- get(fam, envir = asNamespace("frmtmb"))
    fit <- frm(y | thres(gr = g) ~ x, data = d, family = ff(link = lk))
    ordered <- fam %in% c("cumulative", "sratio")
    par <- fit$obj$par
    par[names(par) == "beta"] <- pt_b
    par[names(par) == "tau_raw"] <- raw_of(pt_thres, ordered)
    nll <- fit$obj$fn(par)
    ref <- sum(vapply(seq_len(nrow(d)), function(i) {
      brms_row_lpmf(fam, lk, d$y[i], pt_b * d$x[i], pt_thres[[d$g[i]]])
    }, 0))
    cat(sprintf("%-10s %-8s frmtmb %.12f brms %.12f rel %.3e\n", fam, lk,
                -nll, ref, abs(-nll - ref) / abs(ref)))
  }
}

cat("\n== 2. MASS::polr per group, group-specific slopes ==\n")
d <- sim_data(12, n = 600)
fit <- frm(y | thres(gr = g) ~ g:x, data = d, family = cumulative())
ll_polr <- sum(vapply(c("a", "b", "c"), function(gg) {
  dg <- d[d$g == gg, ]
  as.numeric(logLik(MASS::polr(factor(y) ~ x, data = dg)))
}, 0))
cat(sprintf("frmtmb %.10f polr sum %.10f rel %.3e\n", as.numeric(logLik(fit)),
            ll_polr, abs(as.numeric(logLik(fit)) - ll_polr) / abs(ll_polr)))
th_polr <- unlist(lapply(c("a", "b", "c"), function(gg) {
  MASS::polr(factor(y) ~ x, data = d[d$g == gg, ])$zeta
}))
th_frm <- fixef(fit)[grepl("^Intercept", rownames(fixef(fit))), 1]
cat("thresholds, max |frmtmb - polr| / max|polr|:",
    format(max(abs(th_frm - th_polr)) / max(abs(th_polr)), digits = 3),
    "\n")
print(cbind(frmtmb = th_frm, polr = th_polr))

cat("\n== 3. ordinal::clm, shared slope, nominal = ~ g ==\n")
set.seed(13)
n <- 600
d3 <- data.frame(x = rnorm(n), g = sample(c("a", "b"), n, TRUE))
u <- rlogis(n, 0.8 * d3$x)
d3$y <- ifelse(d3$g == "a", 1 + (u > -1) + (u > 0.3) + (u > 1.4),
               1 + (u > -0.2) + (u > 0.5) + (u > 2))
fit3 <- frm(y | thres(gr = g) ~ x, data = d3, family = cumulative())
clm3 <- ordinal::clm(factor(y) ~ x, nominal = ~ g, data = d3)
cat(sprintf("frmtmb %.10f clm %.10f rel %.3e\n", as.numeric(logLik(fit3)),
            as.numeric(logLik(clm3)),
            abs(as.numeric(logLik(fit3)) - as.numeric(logLik(clm3))) /
              abs(as.numeric(logLik(clm3)))))
cat(sprintf("slope frmtmb %.8f clm %.8f\n", fixef(fit3)["x", 1],
            coef(clm3)[["x"]]))
# clm's nominal thresholds: level a is the baseline, level b adds the
# g-contrast per threshold
cf <- coef(clm3)
th_clm <- c(cf[1:3], cf[1:3] + cf[4:6])
th_f3 <- fixef(fit3)[grepl("^Intercept", rownames(fixef(fit3))), 1]
cat("thresholds, max |frmtmb - clm|:", format(max(abs(th_f3 - th_clm)),
                                             digits = 3), "\n")

cat("\n== 4. thres(x = 5) without groups, brms density at a point ==\n")
d4 <- sim_data(14)
d4$y <- pmin(d4$y, 3)
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  ff <- get(fam, envir = asNamespace("frmtmb"))
  fit <- frm(y | thres(5) ~ x, data = d4, family = ff(),
             prior = set_prior("normal(0, 3)", class = "Intercept"))
  ordered <- fam %in% c("cumulative", "sratio")
  th <- c(-1, -0.2, 0.7, 1.5, 2.4)
  par <- fit$obj$par
  par[names(par) == "beta"] <- pt_b
  par[names(par) == "tau_raw"] <- raw_of(list(th), ordered)
  ll <- sum(frmtmb:::row_lpdf(fit$spec$responses$y$family, d4$y, d4$y,
                              list(mu = pt_b * d4$x), list(),
                              list(tau_raw = raw_of(list(th), ordered))))
  ref <- sum(vapply(seq_len(nrow(d4)), function(i) {
    brms_row_lpmf(fam, "logit", d4$y[i], pt_b * d4$x[i], th)
  }, 0))
  cat(sprintf("%-10s nthres %d frmtmb %.12f brms %.12f rel %.3e\n", fam,
              length(fit$estimates$tau_raw), ll, ref,
              abs(ll - ref) / abs(ref)))
}
