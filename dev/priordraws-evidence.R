# Evidence for sample_prior against brms 2.23.0, on the lane's build.
#   Part 1: frm_sample(sample_prior = "yes") draws against the exact
#           (truncated) prior of each column, one-sample KS, frmtmb only.
#   Part 2: the same columns against brms's own prior draws, two-sample
#           KS (models B and D2 of dev/priordraws-brms-names.R).
#   Part 3: sample_prior = "only" against brms's: means and quantiles of
#           the parameters and of the prior predictive, z-scores with
#           Monte Carlo standard errors.
#   Part 4: Evid.Ratio of point hypotheses against brms's, with the
#           between-block Monte Carlo error of each side.
# Seeds: data seed 1; samplers seed 11 (frmtmb) and 11 (brms).
# Usage: Rscript dev/priordraws-evidence.R [part ...]
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
parts <- commandArgs(TRUE)
if (!length(parts)) parts <- c("1", "2", "3", "4")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$yb <- rbinom(60, 1, 0.5)
d$yp <- rbeta(60, 2, 2)
d$yc <- rpois(60, 3)
d$yo <- pmin(d$yc, 4) + 1
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
P <- function(...) set_prior(...)

# an independent reading of a prior string: the CDF of the density brms
# names, truncated to [lb, ub]
ref_cdf <- function(s, lb = -Inf, ub = Inf) {
  k <- sub("[(].*$", "", s)
  a <- as.numeric(strsplit(sub("^[^(]*[(](.*)[)]$", "\\1", s), ",")[[1]])
  F <- switch(k,
    normal = function(x) pnorm(x, a[1], a[2]),
    student_t = function(x) pt((x - a[2]) / a[3], a[1]),
    cauchy = function(x) pcauchy(x, a[1], a[2]),
    exponential = function(x) pexp(x, a[1]),
    gamma = function(x) pgamma(x, a[1], rate = a[2]),
    inv_gamma = function(x) pgamma(1 / x, a[1], rate = a[2],
                                   lower.tail = FALSE),
    beta = function(x) pbeta(x, a[1], a[2]),
    logistic = function(x) plogis(x, a[1], a[2]),
    lognormal = function(x) plnorm(x, a[1], a[2]),
    uniform = function(x) punif(x, a[1], a[2]))
  lo <- if (is.finite(lb)) F(lb) else 0
  hi <- if (is.finite(ub)) F(ub) else 1
  function(x) (pmin(pmax(F(x), lo), hi) - lo) / (hi - lo)
}

# brms's declared bounds, by class
cls_bounds <- function(nm) {
  if (grepl("^(sigma|sd_|sds_|phi|shape)", nm)) return(c(0, Inf))
  if (grepl("^nu", nm)) return(c(1, Inf))
  if (grepl("^(hu|zi)", nm)) return(c(0, 1))
  c(-Inf, Inf)
}

if ("1" %in% parts) {
  cat("\n#### PART 1: prior draws against the exact prior (one-sample KS)\n")
  cases <- list(
    B = list(bf(y ~ x + z + (1 + x | g)), gaussian(),
             P("normal(0, 1)", class = "b") +
               P("normal(0, 2)", class = "b", coef = "z") +
               P("normal(1, 3)", class = "Intercept") +
               P("exponential(2)", class = "sigma") +
               P("lkj(2)", class = "cor") +
               P("cauchy(0, 1)", class = "sd", group = "g", coef = "x")),
    C2 = list(bf(y ~ x + z), gaussian(),
              P("normal(0, 1)", class = "b", lb = 0, ub = 1)),
    D2 = list(bf(y ~ mo(m) + x), gaussian(), P("normal(0, 1)", class = "b")),
    K = list(bf(y ~ x), student()),
    E = list(bf(yp ~ x), Beta()),
    S = list(bf(yc ~ x), hurdle_poisson()),
    M = list(bf(y ~ s(x)), gaussian()),
    O = list(bf(y ~ x) + bf(z ~ x) + set_rescor(TRUE), NULL,
             P("normal(0, 1)", class = "b", resp = "y")),
    R = list(bf(y ~ x), mixture(gaussian(), gaussian()),
             P("normal(0, 3)", class = "Intercept", dpar = "mu1") +
               P("normal(1, 3)", class = "Intercept", dpar = "mu2"))
  )
  # the prior each column is drawn from, read off the effective prior
  # table: the user's rows above, brms's defaults otherwise
  col_prior <- function(ds, nm) {
    pt <- as.data.frame(prior_summary(ds))
    pick <- function(cl, coef = "", group = "", dpar = "", resp = "") {
      r <- pt[pt$class == cl & pt$coef == coef & pt$group == group &
                pt$dpar == dpar & pt$resp == resp, ]
      if (!nrow(r) && nzchar(coef)) return(pick(cl, "", group, dpar, resp))
      if (!nrow(r) && nzchar(group)) return(pick(cl, "", "", dpar, resp))
      r
    }
    r <- switch(nm,
      Intercept = pick("Intercept"), b = pick("b"),
      b_x = pick("b", "x"), b_z = pick("b", "z"), bsp = pick("b"),
      sigma = pick("sigma"), nu = pick("nu"), phi = pick("phi"),
      hu = pick("hu"), sd_g__Intercept = pick("sd", "Intercept", "g"),
      sd_g__x = pick("sd", "x", "g"), sds_sx = pick("sd"),
      b_y = pick("b", resp = "y"), Intercept_y = pick("Intercept",
                                                      resp = "y"),
      Intercept_z = pick("Intercept", resp = "z"),
      sigma_y = pick("sigma", resp = "y"), sigma_z = pick("sigma",
                                                          resp = "z"),
      Intercept_mu1 = pick("Intercept", dpar = "mu1"),
      Intercept_mu2 = pick("Intercept", dpar = "mu2"),
      sigma1 = pick("sigma1"), sigma2 = pick("sigma2"), NULL)
    r
  }
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    ds <- quiet(frm_sample(cs[[1]], data = d, family = cs[[2]],
                           prior = if (length(cs) > 2) cs[[3]],
                           sample_prior = "yes", chains = 1, iter = 4000,
                           seed = 11, refresh = 0))
    pd <- prior_draws(ds)
    for (cn in names(pd)) {
      v <- pd[[cn]]
      if (grepl("^cor_|^rescor$", cn)) {
        eta <- if (identical(nm, "B")) 2 else 1
        a <- eta - 1 + 2 / 2
        ks <- ks.test(v, function(x) pbeta((x + 1) / 2, a, a))
        what <- sprintf("LKJ(%g) marginal, M = 2", eta)
      } else if (grepl("^(simo|theta)", cn)) {
        D <- sum(grepl(paste0("^", sub("[[].*$", "", cn), "[[]"),
                       names(pd)))
        ks <- ks.test(v, function(x) pbeta(x, 1, D - 1))
        what <- sprintf("Dirichlet(1) marginal, D = %d", D)
      } else {
        r <- col_prior(ds, cn)
        if (is.null(r) || !nrow(r)) {
          cat(sprintf("  %-4s %-16s NO REFERENCE\n", nm, cn))
          next
        }
        bd <- cls_bounds(cn)
        lb <- if (!is.na(r$lb[1])) as.numeric(r$lb[1]) else bd[1]
        ub <- if (!is.na(r$ub[1])) as.numeric(r$ub[1]) else bd[2]
        ks <- ks.test(v, ref_cdf(r$prior[1], lb, ub))
        what <- sprintf("%s on [%g, %g]", r$prior[1], lb, ub)
      }
      cat(sprintf("  %-4s %-16s n=%d KS D=%.4f p=%.3f  %s\n", nm, cn,
                  length(v), ks$statistic, ks$p.value, what))
    }
  }
}

if (any(c("2", "3", "4") %in% parts)) {
  loadNamespace("brms")
  brmq <- function(...) quiet(brms::brm(..., refresh = 0, cores = 1))
}

if ("2" %in% parts) {
  cat("\n#### PART 2: prior draws against brms's (two-sample KS)\n")
  pB <- c(brms::prior(normal(0, 1), class = b),
          brms::prior(normal(0, 2), class = b, coef = z),
          brms::prior(normal(1, 3), class = Intercept),
          brms::prior(exponential(2), class = sigma),
          brms::prior(lkj(2), class = cor),
          brms::prior(cauchy(0, 1), class = sd, group = g, coef = x))
  fB <- P("normal(0, 1)", class = "b") +
    P("normal(0, 2)", class = "b", coef = "z") +
    P("normal(1, 3)", class = "Intercept") +
    P("exponential(2)", class = "sigma") + P("lkj(2)", class = "cor") +
    P("cauchy(0, 1)", class = "sd", group = "g", coef = "x")
  two <- list(
    B = list(brms = list(y ~ x + z + (1 + x | g), prior = pB),
             frm = list(bf(y ~ x + z + (1 + x | g)), gaussian(), fB)),
    D2 = list(brms = list(y ~ mo(m) + x,
                          prior = brms::prior(normal(0, 1), class = b)),
              frm = list(bf(y ~ mo(m) + x), gaussian(),
                         P("normal(0, 1)", class = "b"))))
  for (nm in names(two)) {
    tb <- two[[nm]]
    fb <- do.call(brmq, c(tb$brms, list(data = d, sample_prior = "yes",
                                        chains = 2, iter = 4000,
                                        seed = 11)))
    ds <- quiet(frm_sample(tb$frm[[1]], data = d, family = tb$frm[[2]],
                           prior = tb$frm[[3]], sample_prior = "yes",
                           chains = 2, iter = 4000, seed = 11,
                           refresh = 0))
    pb <- prior_draws(fb)
    pf <- prior_draws(ds)
    cat(sprintf("  %s names identical: %s\n", nm,
                identical(names(pb), names(pf))))
    # brms's default Intercept prior reads the response, as frmtmb's
    # does; both print it
    cat("  brms Intercept prior:",
        fb$prior$prior[fb$prior$class == "Intercept"][1], "\n")
    for (cn in intersect(names(pb), names(pf))) {
      ks <- suppressWarnings(ks.test(pb[[cn]], pf[[cn]]))
      cat(sprintf("  %-3s %-16s brms mean %8.4f sd %8.4f | frmtmb mean %8.4f sd %8.4f | KS D=%.4f p=%.3f\n",
                  nm, cn, mean(pb[[cn]]), sd(pb[[cn]]), mean(pf[[cn]]),
                  sd(pf[[cn]]), ks$statistic, ks$p.value))
    }
  }
}

mc_z <- function(a, b) {
  # a, b: iterations x chains matrices; z of the difference of means
  sa <- posterior::mcse_mean(a)
  sb <- posterior::mcse_mean(b)
  (mean(a) - mean(b)) / sqrt(sa^2 + sb^2)
}
mc_zq <- function(a, b, p) {
  sa <- posterior::mcse_quantile(a, probs = p)
  sb <- posterior::mcse_quantile(b, probs = p)
  (quantile(a, p) - quantile(b, p)) / sqrt(sa^2 + sb^2)
}
chains_of <- function(v, nc) matrix(v, ncol = nc)

if ("3" %in% parts) {
  cat("\n#### PART 3: sample_prior = \"only\" against brms\n")
  fb <- brmq(y ~ x + (1 | g), data = d,
             prior = brms::prior(normal(0, 1), class = b),
             sample_prior = "only", chains = 4, iter = 4000, seed = 11)
  ds <- quiet(frm_sample(bf(y ~ x + (1 | g)), data = d, family = gaussian(),
                         prior = P("normal(0, 1)", class = "b"),
                         sample_prior = "only", chains = 4, iter = 4000,
                         seed = 11, refresh = 0))
  cat("  brms Intercept prior:",
      fb$prior$prior[fb$prior$class == "Intercept"][1], "| frmtmb:",
      as.data.frame(prior_summary(ds))$prior[
        as.data.frame(prior_summary(ds))$class == "Intercept"][1], "\n")
  mb <- as.matrix(fb)
  mf <- as.matrix(ds)
  for (v in c("b_Intercept", "b_x", "sigma", "sd_g__Intercept",
              "r_g[1,Intercept]")) {
    a <- chains_of(mb[, v], 4)
    b <- chains_of(mf[, v], 4)
    cat(sprintf("  %-17s brms mean %8.4f q05 %8.4f q50 %8.4f q95 %8.4f | frmtmb mean %8.4f q05 %8.4f q50 %8.4f q95 %8.4f | z mean %6.2f q05 %6.2f q50 %6.2f q95 %6.2f\n",
                v, mean(a), quantile(a, .05), quantile(a, .5),
                quantile(a, .95), mean(b), quantile(b, .05),
                quantile(b, .5), quantile(b, .95), mc_z(a, b),
                mc_zq(a, b, .05), mc_zq(a, b, .5), mc_zq(a, b, .95)))
  }
  set.seed(5)
  yb <- posterior_predict(fb)[, 1]
  set.seed(5)
  yf <- posterior_predict(ds)[, 1]
  a <- chains_of(yb, 4)
  b <- chains_of(yf, 4)
  cat(sprintf("  yrep[,1] brms mean %8.4f q05 %8.4f q50 %8.4f q95 %8.4f | frmtmb mean %8.4f q05 %8.4f q50 %8.4f q95 %8.4f | z mean %6.2f q05 %6.2f q50 %6.2f q95 %6.2f\n",
              mean(a), quantile(a, .05), quantile(a, .5), quantile(a, .95),
              mean(b), quantile(b, .05), quantile(b, .5), quantile(b, .95),
              mc_z(a, b), mc_zq(a, b, .05), mc_zq(a, b, .5),
              mc_zq(a, b, .95)))
  cat("  frmtmb rhat max:", max(rhat(ds), na.rm = TRUE), " divergences:",
      sum(nuts_params(ds)$Value[nuts_params(ds)$Parameter ==
                                  "divergent__"]), "\n")
}

if ("4" %in% parts) {
  cat("\n#### PART 4: Evid.Ratio against brms\n")
  set.seed(4)
  d4 <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  d4$y <- 0.8 + 0.2 * d4$x + rnorm(6, 0, 0.4)[d4$g] + rnorm(60)
  fb <- brmq(y ~ x + (1 | g), data = d4,
             prior = c(brms::prior(normal(0, 1), class = b),
                       brms::prior(student_t(3, 0, 1), class = sd),
                       brms::prior(student_t(3, 0, 1), class = sigma)),
             sample_prior = "yes", chains = 4, iter = 6000, seed = 11)
  ds <- quiet(frm_sample(bf(y ~ x + (1 | g)), data = d4,
                         family = gaussian(),
                         prior = P("normal(0, 1)", class = "b") +
                           P("student_t(3, 0, 1)", class = "sd") +
                           P("student_t(3, 0, 1)", class = "sigma"),
                         sample_prior = "yes", chains = 4, iter = 6000,
                         seed = 11, refresh = 0))
  kde <- frmtmb.sample:::er_kde_at
  blocks <- function(n, k = 8) split(seq_len(n), rep(seq_len(k),
                                                      each = n %/% k))
  # the ratio and its error: the spread of each density estimate across
  # 8 consecutive blocks, combined by the delta method
  er_se <- function(post, pri) {
    num <- kde(post)
    den <- kde(pri)
    sn <- sd(vapply(blocks(length(post)), function(i) kde(post[i]), 1)) /
      sqrt(8)
    sd_ <- sd(vapply(blocks(length(pri)), function(i) kde(pri[i]), 1)) /
      sqrt(8)
    c(er = num / den, se = num / den * sqrt((sn / num)^2 + (sd_ / den)^2))
  }
  hyps <- list(c("x = 0", "b"), c("x = 0.2", "b"),
               c("sigma = 1", NA), c("sd_g__Intercept = 0.5", NA))
  for (h in hyps) {
    cls <- if (is.na(h[2])) NULL else h[2]
    hb <- quiet(brms::hypothesis(fb, h[1], class = cls, seed = 1))
    hf <- quiet(hypothesis(ds, h[1], class = cls, seed = 1))
    sb <- er_se(hb$samples$H1, hb$prior_samples$H1)
    sf <- er_se(hf$samples$H1, hf$prior_samples$H1)
    z <- (sf["er"] - sb["er"]) / sqrt(sf["se"]^2 + sb["se"]^2)
    cat(sprintf("  %-22s brms ER %8.4f (se %.4f, printed %8.4f) | frmtmb ER %8.4f (se %.4f, printed %8.4f) | z %6.2f\n",
                h[1], sb["er"], sb["se"], hb$hypothesis$Evid.Ratio,
                sf["er"], sf["se"], hf$hypothesis$Evid.Ratio, z))
  }
}
