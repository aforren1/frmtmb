# Reviewer: how wide is D1? For every core family whose compat row says
# nl is refused, fit the same predictor written linearly and as an nl
# body, and the intercept-only model. An nl fit whose logLik equals the
# intercept-only fit (and not the linear fit) dropped the body. brms
# 2.23.0's make_stancode() on the nl spelling says whether brms builds it.
# Seed 4901.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); requireNamespace("brms")})
set.seed(4901)
n <- 300
x <- round(rnorm(n), 3)
eta <- 0.3 + 0.6 * x
pr <- plogis(eta)
cat_y <- function(k) {
  p <- cbind(1, exp(0.5 * x), exp(-0.5 * x + 0.2), exp(0.3))[, seq_len(k)]
  p <- p / rowSums(p)
  apply(p, 1, function(q) sample(seq_len(k), 1, prob = q))
}
ord <- ordered(cut(eta + rlogis(n), c(-Inf, -0.2, 0.5, 1.3, Inf),
                   labels = 1:4))
pos <- exp(0.5 * eta + rnorm(n, 0, 0.4))
cnt <- rpois(n, exp(eta))
unit <- rbeta(n, 4 * pr, 4 * (1 - pr))
tr <- rep(10L, n)
bin <- rbinom(n, 10, pr)
simp <- t(apply(cbind(exp(0.4 * x), 1, exp(-0.3 * x)), 1, function(a) {
  g <- rgamma(3, 5 * a); g / sum(g) }))
mult <- t(apply(simp, 1, function(p) rmultinom(1, 10, p)))
d <- data.frame(x = x, ord = ord, pos = pos, cnt = cnt, unit = unit, tr = tr,
                bin = bin, ber = rbinom(n, 1, pr),
                real = eta + rnorm(n),
                circ = ((eta + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi,
                hpos = ifelse(runif(n) < 0.2, 0, pos),
                hcnt = ifelse(runif(n) < 0.25, 0L, 1L + cnt),
                zcnt = ifelse(runif(n) < 0.2, 0L, cnt),
                zunit = ifelse(runif(n) < 0.2, 0, unit),
                zoiunit = ifelse(runif(n) < 0.1, 0, ifelse(runif(n) < 0.1, 1,
                                                          unit)),
                zbin = ifelse(runif(n) < 0.2, 0L, bin),
                hord = ordered(ifelse(runif(n) < 0.2, 0L, as.integer(ord)),
                               levels = 0:4),
                cat = factor(cat_y(3)), xunit = pmin(pmax(unit, 0.01), 0.99),
                tw = ifelse(runif(n) < 0.3, 0, pos))
d$S <- simp
d$M <- mult
# response, family text, addition
F <- list(
  Gamma = c("pos", "Gamma(link = 'log')"),
  lognormal = c("pos", "lognormal()"),
  weibull = c("pos", "weibull()"),
  exponential = c("pos", "exponential()"),
  frechet = c("pos", "frechet()"),
  inverse.gaussian = c("pos", "inverse.gaussian(link = 'log')"),
  shifted_lognormal = c("pos", "shifted_lognormal()"),
  student = c("real", "student()"),
  skew_normal = c("real", "skew_normal()"),
  asym_laplace = c("real", "asym_laplace()"),
  exgaussian = c("real", "exgaussian()"),
  gen_extreme_value = c("real", "gen_extreme_value()"),
  von_mises = c("circ", "von_mises()"),
  poisson = c("cnt", "poisson()"),
  negbinomial = c("cnt", "negbinomial()"),
  negbinomial2 = c("cnt", "negbinomial2()"),
  nbinom1 = c("cnt", "nbinom1()"),
  geometric = c("cnt", "geometric()"),
  compois = c("cnt", "com_poisson()"),
  discrete_weibull = c("cnt", "discrete_weibull()"),
  hurdle_poisson = c("hcnt", "hurdle_poisson()"),
  hurdle_negbinomial = c("hcnt", "hurdle_negbinomial()"),
  zero_inflated_poisson = c("zcnt", "zero_inflated_poisson()"),
  zero_inflated_negbinomial = c("zcnt", "zero_inflated_negbinomial()"),
  binomial = c("bin | trials(tr)", "binomial()"),
  beta_binomial = c("bin | trials(tr)", "beta_binomial()"),
  zero_inflated_binomial = c("zbin | trials(tr)", "zero_inflated_binomial()"),
  zero_inflated_beta_binomial = c("zbin | trials(tr)",
                                  "zero_inflated_beta_binomial()"),
  bernoulli = c("ber", "bernoulli()"),
  beta = c("xunit", "Beta()"),
  xbeta = c("unit", "xbeta()"),
  zero_inflated_beta = c("zunit", "zero_inflated_beta()"),
  zero_one_inflated_beta = c("zoiunit", "zero_one_inflated_beta()"),
  hurdle_gamma = c("hpos", "hurdle_gamma()"),
  hurdle_lognormal = c("hpos", "hurdle_lognormal()"),
  tweedie = c("tw", "tweedie()"),
  cumulative = c("ord", "cumulative()"),
  sratio = c("ord", "sratio()"),
  cratio = c("ord", "cratio()"),
  acat = c("ord", "acat()"),
  hurdle_cumulative = c("hord", "hurdle_cumulative()"),
  categorical = c("cat", "categorical()"),
  multinomial = c("M | trials(tr)", "multinomial()"),
  dirichlet = c("S", "dirichlet()"),
  dirichlet2 = c("S", "dirichlet2()"),
  logistic_normal = c("S", "logistic_normal()"),
  dirichlet_multinomial = c("M | trials(tr)", "dirichlet_multinomial()"))
ev <- function(txt, pkg) eval(parse(text = txt), envir = asNamespace(pkg))
ll <- function(f) if (inherits(f, "error")) NA_real_ else
  as.numeric(logLik(f))
tryf <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
ordinal <- c("cumulative", "sratio", "cratio", "acat", "hurdle_cumulative")
rows <- list()
for (nm in names(F)) {
  resp <- F[[nm]][1]
  fam_f <- tryf(ev(F[[nm]][2], "frmtmb"))
  if (inherits(fam_f, "error")) {
    rows[[nm]] <- data.frame(family = nm, status = NA, lin = NA, nl = NA,
                             int = NA, verdict = "family ctor failed",
                             brms = NA, msg = conditionMessage(fam_f))
    next
  }
  # an ordinal body has no intercept (thresholds), every other has b0
  body <- if (nm %in% ordinal) "b1 * x" else "b0 + b1 * x"
  pars <- if (nm %in% ordinal) "b1 ~ 1" else "b0 + b1 ~ 1"
  f_lin <- tryf(frm(as.formula(paste(resp, "~ x")), family = fam_f, data = d))
  f_nl <- tryf(frm(bf(as.formula(paste(resp, "~", body)),
                      as.formula(pars), nl = TRUE), family = fam_f, data = d))
  f_0 <- tryf(frm(as.formula(paste(resp, "~ 1")), family = fam_f, data = d))
  st <- tryCatch(as.data.frame(frm_compat("nl", nm))$status[1],
                 error = function(e) NA)
  a <- ll(f_lin); b <- ll(f_nl); c0 <- ll(f_0)
  v <- if (inherits(f_nl, "error")) "frm refuses" else
    if (is.na(a)) "linear fit failed" else
    if (abs(a - b) < 1e-4 * max(1, abs(a))) "nl honored" else
    if (!is.na(c0) && abs(c0 - b) < 1e-4 * max(1, abs(c0))) "BODY DROPPED" else
    "nl differs from both"
  bb <- tryf(brms::make_stancode(
    brms::bf(as.formula(paste(resp, "~", body)), as.formula(pars),
             nl = TRUE), family = ev(F[[nm]][2], "brms"), data = d))
  rows[[nm]] <- data.frame(
    family = nm, status = st, lin = round(a, 3), nl = round(b, 3),
    int = round(c0, 3), verdict = v,
    brms = if (inherits(bb, "error")) paste("refuses:",
      substr(conditionMessage(bb), 1, 60)) else "builds",
    msg = if (inherits(f_nl, "error")) substr(conditionMessage(f_nl), 1, 90)
      else if (inherits(f_lin, "error")) substr(conditionMessage(f_lin), 1, 90)
      else "")
}
res <- do.call(rbind, rows)
options(width = 250)
print(res, row.names = FALSE, right = FALSE)
cat("\nverdicts:\n"); print(table(res$verdict, res$status, useNA = "ifany"))
saveRDS(res, "dev/review-equiv-scratch/d1-width.rds")
