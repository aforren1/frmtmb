# Reviewer of lane optima, re-check (c) and (d): every reader of a mo()
# model's draws, and the designs beyond one term. For each model: the
# simo_ columns against the chart columns a reader hands back
# (mo_simplex() of them, to rounding), brms's names for the same model
# (brms 2.23.0, a short run), and the readers run on the draws. On the
# gaussian fit1 the linear predictor and log_lik are rebuilt by hand
# from the simo_ draws.
#   Rscript dev/optima-rev2-readers.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
try1 <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("  %-34s ERROR %s\n", label,
                substr(gsub("\n", " ", conditionMessage(r)), 1, 160)))
  }
  invisible(r)
}
q <- function(expr) suppressWarnings(suppressMessages(expr))

set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
strong <- data.frame(income, ls)
fit <- frm(bf(ls ~ mo(income)), data = strong, family = gaussian())
s <- q(frm_sample(fit, chains = 2, iter = 1000, warmup = 500, seed = 11,
                  cores = 1, refresh = 0))
cat("== fit1 readers ==\n")
m <- as.matrix(s)
W <- m[, grep("^simo_", colnames(m))]
im <- frmtmb.sample:::draws_internal_matrix(s)
# the internal matrix keeps the D - 1 chart columns under simo_ names
zc <- im[, grep("^simo_", colnames(im)), drop = FALSE]
cat("internal columns:", colnames(im), "\n")
Wc <- t(apply(zc, 1, frmtmb::mo_simplex))
cat("max |mo_simplex(chart draw) - simo_ draw|:",
    format(max(abs(Wc - W)), digits = 3), "\n")
code <- as.integer(strong$income) - 1L
eta <- t(vapply(seq_len(nrow(m)), function(i) {
  cz <- c(0, cumsum(W[i, ]))
  m[i, "b_Intercept"] + m[i, "bsp_moincome"] * 3 * cz[code + 1L]
}, numeric(100)))
pl <- try1("posterior_linpred", posterior_linpred(s))
if (is.matrix(pl)) {
  cat("posterior_linpred vs hand-built from simo_: max |diff|",
      format(max(abs(pl - eta)), digits = 3), "relative",
      format(max(abs(pl - eta)) / max(abs(eta)), digits = 3), "\n")
}
pe <- try1("posterior_epred", posterior_epred(s))
if (is.matrix(pe)) cat("posterior_epred - linpred max |diff|",
                       format(max(abs(pe - eta)), digits = 3), "\n")
ll <- try1("log_lik", log_lik(s))
if (is.matrix(ll)) {
  llh <- t(vapply(seq_len(nrow(m)), function(i) {
    stats::dnorm(strong$ls, eta[i, ], m[i, "sigma"], log = TRUE)
  }, numeric(100)))
  cat("log_lik vs hand-built: max |diff|", format(max(abs(ll - llh)),
                                                  digits = 3), "\n")
}
nd <- data.frame(income = factor(lev, levels = lev, ordered = TRUE))
pn <- try1("posterior_epred(newdata)", posterior_epred(s, newdata = nd))
if (is.matrix(pn)) {
  hand <- t(vapply(seq_len(nrow(m)), function(i) {
    cz <- c(0, cumsum(W[i, ]))
    m[i, "b_Intercept"] + m[i, "bsp_moincome"] * 3 * cz
  }, numeric(4)))
  cat("posterior_epred(newdata) vs hand: max |diff|",
      format(max(abs(pn - hand)), digits = 3), "\n")
}
f1 <- try1("fitted", fitted(s))
p1 <- try1("predict", q(predict(s)))
ce <- try1("conditional_effects", q(conditional_effects(s, "income")))
if (!inherits(ce, "error")) {
  e <- ce[[1]]$estimate__
  cat("conditional_effects estimates:", format(e, digits = 4), "\n")
}
ad <- try1("as_draws_df", posterior::as_draws_df(s))
cat("variables():", try1("variables", variables(s)), "\n")
rv <- try1("as_draws_rvars", posterior::as_draws_rvars(s))
if (!inherits(rv, "error")) cat("rvars names:", names(rv), "\n")
h <- try1("hypothesis", hypothesis(s, "simo_moincome1[1] > 0.5",
                                   class = NULL))
if (!inherits(h, "error")) print(h$hypothesis)
h2 <- try1("hypothesis (no brackets)",
           hypothesis(s, "bsp_moincome > 15", class = NULL))
if (!inherits(h2, "error")) print(h2$hypothesis)
lo <- try1("loo", q(loo(s)))
if (!inherits(lo, "error")) cat("loo elpd", format(lo$estimates[1, 1]), "\n")
wa <- try1("waic", q(waic(s)))
bs <- try1("bridge_sampler", q(bridge_sampler(s, silent = TRUE)))
if (!inherits(bs, "error")) cat("bridge logml", format(bs$logml), "\n")
ps <- try1("posterior_summary", posterior_summary(s))
fx <- try1("fixef", fixef(s))
su <- try1("summary", summary(s))
lp <- try1("log_posterior", log_posterior(s))
up <- try1("update(newdata)", q(update(s, newdata = strong[1:80, ],
                                         chains = 1, iter = 300,
                                         refresh = 0)))
if (!inherits(up, "error")) {
  cat("update() draws simo_ columns:",
      paste(grep("^simo_", colnames(as.matrix(up)), value = TRUE),
            collapse = " "), "\n")
}

# (d) designs beyond one term
set.seed(5)
n <- 300
d <- data.frame(x1 = sample(0:3, n, TRUE), x2 = sample(0:4, n, TRUE),
                z = rnorm(n), g = factor(sample(1:15, n, TRUE)))
u <- rnorm(15, sd = 0.5)
d$y <- c(0, 1, 1, 2)[d$x1 + 1] + c(0, 0.3, 0.6, 0.6, 1)[d$x2 + 1] +
  0.3 * d$z + u[d$g] + rnorm(n, sd = exp(c(0, 0.2, 0.2, 0.4)[d$x1 + 1]))
d$yn <- 2 * exp(-0.2 * c(0, 1, 1, 3)[d$x1 + 1]) * (1 + 0.3 * d$z) +
  rnorm(n, sd = 0.3)
designs <- list(
  two_terms_interaction = list(f = bf(y ~ mo(x1) * z + mo(x2))),
  sigma = list(f = bf(y ~ z, sigma ~ mo(x1))),
  nlpar = list(f = bf(yn ~ a * (1 + bz * z), a ~ mo(x1), bz ~ 1,
                      nl = TRUE)),
  ranef = list(f = bf(y ~ mo(x1) + (1 | g)))
)
for (nm in names(designs)) {
  cat("\n==", nm, "==\n")
  f <- designs[[nm]]$f
  ft <- try1("frm", q(frm(f, data = d, family = gaussian())))
  if (inherits(ft, "error")) next
  sm <- try1("frm_sample", q(frm_sample(ft, chains = 2, iter = 800,
                                        warmup = 400, seed = 3, cores = 1,
                                        refresh = 0)))
  if (inherits(sm, "error")) next
  mm <- as.matrix(sm)
  sim <- grep("^simo", colnames(mm), value = TRUE)
  cat("frmtmb simo columns:", sim, "\n")
  im <- frmtmb.sample:::draws_internal_matrix(sm)
  for (tm in frmtmb::mo_frame_terms(ft)) {
    zc <- im[, intersect(tm$names, colnames(im)), drop = FALSE]
    Wc <- t(apply(zc, 1, frmtmb::mo_simplex))
    cat(sprintf("  %s -> %s: max |mo_simplex(chart) - simo_| %.3g\n",
                tm$zeta, paste(tm$names, collapse = ","),
                max(abs(Wc - mm[, tm$names, drop = FALSE]))))
  }
  sp <- rstan::get_sampler_params(sm$stanfit, inc_warmup = FALSE)
  rh <- vapply(sim, function(cn) {
    posterior::rhat(posterior::as_draws_array(sm)[, , cn])
  }, 0)
  cat("  divergences", sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
      "| max R-hat on weights", format(max(rh), digits = 4), "\n")
  pe <- try1("posterior_epred", posterior_epred(sm))
  if (is.matrix(pe)) cat("  posterior_epred finite:", all(is.finite(pe)),
                         "\n")
  ll <- try1("log_lik", log_lik(sm))
  lo <- try1("loo", q(loo(sm)))
  ce <- try1("conditional_effects", q(conditional_effects(sm, "x1")))
  b <- try1("brm", q(brms::brm(f, data = d, chains = 1, iter = 200,
                               seed = 1, refresh = 0)))
  if (!inherits(b, "error")) {
    bs <- grep("^simo", posterior::variables(b), value = TRUE)
    cat("brms simo columns:   ", bs, "\n")
    cat("  names identical:", identical(sort(sim), sort(bs)), "\n")
  }
}
