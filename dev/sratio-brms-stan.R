# Lane sratio: brms 2.23.0's COMPILED log density at frmtmb's estimate,
# on fits whose unconstrained optimum has crossing sratio thresholds.
#
#   Rscript dev/sratio-brms-stan.R > dev/sratio-brms-stan-log.txt 2>&1
#
# brms declares sratio's thresholds as a plain `vector`
# (brms:::has_ordered_thres(sratio()) is FALSE), so its mode may have
# them cross. For each case brms writes the Stan program, rstan
# compiles it (about a minute each), and log_prob and its gradient are
# evaluated at frmtmb's estimate with adjust_transform = FALSE. With
# flat priors the log density must equal frmtmb's logLik and its
# gradient must vanish, which says frmtmb's estimate is brms's mode.
# With brms's normal(0, 2) on class Intercept it must equal minus
# frmtmb's penalized objective with NO Jacobian added: an unconstrained
# vector has none, on either side. cratio, whose thresholds were
# already unordered, runs as the control for that convention.
#
# The cases:
#   grouped  y | thres(gr = g) ~ x, data seed 11 without the random
#            effect, n = 240 (dev/thres-sratio-order.R); on the old
#            build level a's second and third thresholds sat 5.6e-8
#            apart, on the ordering boundary
#   inhaler  rating ~ period + carry + cs(treat), brms's own inhaler
#            data; the old build's thresholds 2 and 3 sat 1.8e-8 apart
#   o2       the sratio response of test-mv-gaps.R's data (seed 5),
#            fitted alone; the old build's thresholds sat 6.3e-9 apart
# and, without Stan, the multivariate model of test-mv-gaps.R, whose
# log-likelihood must be the sum of its three univariate fits'.
.libPaths(c("/opt/rlib/stan", "/opt/rlib/lane-sratio", "/opt/rlib/base",
            "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cat("frmtmb from", find.package("frmtmb"), "\n")
cat("rstan", as.character(packageVersion("rstan")), "StanHeaders",
    as.character(packageVersion("StanHeaders")), "\n")

set.seed(11)
n <- 240
dg <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
u <- rlogis(n, 0.7 * dg$x)
tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6),
            c = c(-1.5, -0.4, 0.5, 1.6))
dg$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[dg$g[i]]]), 1L)

data("inhaler", package = "brms")

set.seed(5)
n <- 150
x <- rnorm(n)
z <- rnorm(n)
dm <- data.frame(x = x, z = z, y1 = 1 + 0.5 * x + rnorm(n),
                 y2 = -1 + rnorm(n), y3 = 0.3 * x + rnorm(n),
                 o = cut(x + rlogis(n), c(-Inf, -1, 0, 1, Inf),
                         labels = FALSE),
                 o2 = cut(-x + rlogis(n), c(-Inf, 0, 1, Inf),
                          labels = FALSE))

cases <- list(
  grouped = list(form = y | thres(gr = g) ~ x, data = dg),
  inhaler = list(form = rating ~ period + carry + cs(treat), data = inhaler),
  o2 = list(form = o2 ~ x, data = dm)
)

# brms's parameters at frmtmb's estimate. brms centers X and not Xcs,
# and its ungrouped `Intercept` is the threshold minus means_X' b;
# under thres(gr = ) it does not center at all.
brms_pars <- function(fit, sdat, grouped) {
  fe <- fixef(fit)[, "Estimate"]
  xn <- colnames(sdat$X)
  b <- unname(fe[xn])
  pars <- list(b = array(b, length(b)))
  if (grouped) {
    th <- fit$spec$responses[[1L]]$family[["thres"]]
    for (k in seq_along(th$groups)) {
      nm <- paste0("Intercept[", th$groups[k], ",", seq_len(th$nthres[k]),
                   "]")
      pars[[paste0("Intercept_", k)]] <- as.array(unname(fe[nm]))
    }
    return(pars)
  }
  th <- unname(fe[grep("^Intercept\\[", names(fe))])
  pars$Intercept <- as.array(th - sum(colMeans(sdat$X) * b))
  if (!is.null(sdat$Kcs)) {
    cn <- colnames(sdat$Xcs)
    pars$bcs <- t(vapply(cn, function(v) {
      unname(fe[paste0(v, "[", seq_along(th), "]")])
    }, th))
    dim(pars$bcs) <- c(length(cn), length(th))
  }
  pars
}

stan_at <- function(case, family, prior_int = NULL) {
  bfam <- get(family, envir = asNamespace("brms"))()
  ffam <- get(family, envir = asNamespace("frmtmb"))()
  bf0 <- bf(case$form)
  pr <- get_prior(bf0, data = case$data, family = bfam)
  pr$prior <- ""
  fpr <- NULL
  if (!is.null(prior_int)) {
    pr$prior[pr$class == "Intercept" & pr$coef == "" & pr$group == ""] <-
      prior_int
    fpr <- set_prior(prior_int, class = "Intercept")
  }
  fit <- frm(case$form, data = case$data, family = ffam, prior = fpr)
  code <- make_stancode(bf0, data = case$data, family = bfam, prior = pr)
  sdat <- make_standata(bf0, data = case$data, family = bfam, prior = pr)
  mod <- stan_model(model_code = code)
  sf <- suppressMessages(sampling(mod, data = sdat, chains = 0))
  grouped <- isTRUE(fit$spec$responses[[1L]]$family[["thres"]][["grouped"]])
  up <- unconstrain_pars(sf, brms_pars(fit, sdat, grouped))
  gr <- grad_log_prob(sf, up, adjust_transform = FALSE)
  list(fit = fit, lp = log_prob(sf, up, adjust_transform = FALSE),
       grad = as.numeric(gr))
}

# the smallest gap between adjacent thresholds, per threshold vector
min_gap <- function(fit) {
  fe <- fixef(fit)[, "Estimate"]
  i <- grep("^Intercept\\[", names(fe))
  key <- sub("\\[([^,]*,)?[0-9]+\\]$", "[\\1]", names(fe)[i])
  g <- tapply(fe[i], key, function(t) min(diff(t)))
  paste(sprintf("%s %.4g", names(g), g), collapse = ", ")
}

cat("\n== flat priors: brms log_prob against frmtmb logLik at the MLE ==\n")
for (nm in names(cases)) {
  r <- stan_at(cases[[nm]], "sratio")
  ll <- as.numeric(logLik(r$fit))
  cat(sprintf(paste0("%-8s brms %.17g frmtmb %.17g diff %.3g ",
                     "max|grad| %.2e\n         min gap %s\n"),
              nm, r$lp, ll, r$lp - ll, max(abs(r$grad)), min_gap(r$fit)))
}

cat("\n== normal(0, 2) on class Intercept: brms log_prob against minus",
    "frmtmb's penalized objective, no Jacobian ==\n")
for (fam in c("sratio", "cratio")) {
  for (nm in c("grouped", "inhaler")) {
    r <- stan_at(cases[[nm]], fam, "normal(0, 2)")
    ours <- -r$fit$opt$objective
    cat(sprintf("%-7s %-8s brms %.17g frmtmb %.17g diff %.3g max|grad| %.2e\n",
                fam, nm, r$lp, ours, r$lp - ours, max(abs(r$grad))))
  }
}

cat("\n== multivariate: cumulative o, sratio o2, gaussian y1 ==\n")
# brms is attached after frmtmb, so bf() and the families are brms's
fb <- frmtmb::bf
mv <- frm(fb(o ~ x) + frmtmb::cumulative() + fb(o2 ~ x) +
            frmtmb::sratio() + fb(y1 ~ x) + stats::gaussian(), data = dm)
su <- sum(vapply(list(
  frm(fb(o ~ x) + frmtmb::cumulative(), data = dm),
  frm(fb(o2 ~ x) + frmtmb::sratio(), data = dm),
  frm(fb(y1 ~ x) + stats::gaussian(), data = dm)), function(f) {
    as.numeric(logLik(f))
  }, 0))
cat(sprintf("mv logLik %.17g, sum of univariate %.17g, rel diff %.3g\n",
            as.numeric(logLik(mv)), su,
            abs(as.numeric(logLik(mv)) - su) / abs(su)))
fe <- fixef(mv)[, "Estimate"]
cat("mv o2 thresholds:", format(fe[grep("^o2_Intercept", names(fe))],
                                digits = 10), "\n")
