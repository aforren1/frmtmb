# Reviewer probe 5: edge cases. Usage: Rscript dev/review-reslot-p5.R
# before|after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
q <- function(expr) suppressWarnings(suppressMessages(expr))
msg <- function(expr) tryCatch({ force(expr); "ok" },
                               error = function(e) conditionMessage(e))
dd <- re_data(1L)
set.seed(1)
dd$yn <- 2 * exp(0.3 * dd$x) + rnorm(nrow(dd), 0, 0.3)
show_e <- function(ri) {
  for (e in ri$entries) {
    cat(sprintf("  entry %-8s %-12s %-8s %-9s flat=%s\n", e$comp,
                paste(e$idx, collapse = ","), e$dist$kind, e$scale,
                isTRUE(e$dist$flat)))
  }
}
cat("== logistic_normal under prior = 'flat'\n")
set.seed(2)
P <- matrix(stats::rgamma(4 * nrow(dd), 2), ncol = 4)
P <- P / rowSums(P)
dd$Y <- P
fl <- tryCatch(q(frm(bf(Y ~ x), data = dd, family = logistic_normal(),
                     dry_run = "objective")), error = function(e) e)
if (inherits(fl, "error")) cat("  FIT ERROR:", conditionMessage(fl), "\n") else {
  r <- q(frmtmb.sample:::sample_resolve_priors(fl, "flat"))
  show_e(r$ri)
}

cat("\n== Intercept with nlpar, edge spellings\n")
fnl <- bf(yn ~ a * exp(k * x), a ~ 0 + x, k ~ 1, nl = TRUE)
cat("  Intercept nlpar a on a ~ 0 + x:",
    msg(validate_prior(set_prior("normal(0, 1)", class = "Intercept",
                                 nlpar = "a"), fnl, data = dd)), "\n")
cat("  printed spec:\n")
print(as.data.frame(set_prior("normal(0, 1)", class = "Intercept",
                              nlpar = "a"))[, c("prior", "class", "coef",
                                                "nlpar")])
fnl2 <- bf(yn ~ a * exp(k * x), a ~ 1, k ~ 1, nl = TRUE)
fit <- q(frm(fnl2, data = dd, prior = set_prior("normal(2, 1)",
                                                class = "Intercept",
                                                nlpar = "a") +
               set_prior("normal(0, 1)", nlpar = "k")))
ps <- as.data.frame(prior_summary(fit))
print(ps[, intersect(c("prior", "class", "coef", "nlpar", "source"),
                     names(ps))], row.names = FALSE)

cat("\n== centered intercept refusal: guard-absent spellings\n")
for (f in list(bf(y ~ x, center = FALSE), bf(y ~ 0 + Intercept + x),
               bf(y ~ x, sigma ~ 0 + Intercept + x),
               bf(y ~ x, sigma ~ x, center = FALSE))) {
  dp <- if (length(f$pforms)) "sigma" else ""
  cat("  ", deparse(f$formula), if (nzchar(dp)) "(dpar sigma)", ":",
      msg(validate_prior(set_prior("normal(0, 1)", class = "b",
                                   coef = "Intercept", dpar = dp),
                         f, data = dd)), "\n")
}

cat("\n== bounds-only class-wide sd vs group-level sd on frm_sample\n")
fit1 <- q(frm(bf(y ~ x + (1 | g)), data = dd, dry_run = "objective"))
for (pl in list(set_prior("", class = "sd", lb = 0.1),
                set_prior("", class = "sd", group = "g", lb = 0.1),
                set_prior("", class = "sd", group = "g",
                          coef = "Intercept", lb = 0.1))) {
  cat("  ", format(pl$class), "group=", pl[[1]]$group, "coef=",
      pl[[1]]$coef, "\n")
  r <- tryCatch(q(frmtmb.sample:::sample_resolve_priors(fit1, pl)),
                error = function(e) e)
  if (inherits(r, "error")) cat("  ERROR", conditionMessage(r), "\n") else
    show_e(r$ri)
}

cat("\n== MAP fit with no flat entries: prior = 'flat' drops its penalty\n")
dd$cnt <- rpois(nrow(dd), exp(0.3 + 0.2 * dd$x))
fm <- q(frm(bf(cnt ~ x), family = poisson(), data = dd,
            prior = set_prior("normal(5, 0.1)", class = "b")))
fml <- q(frm(bf(cnt ~ x), family = poisson(), data = dd))
ds <- q(frm_sample(fm, prior = "flat", chains = 1, iter = 200, warmup = 100,
                   seed = 1, refresh = 0))
p <- fml$obj$par
cat("  objective at the ML mode: ML fit", format(fml$obj$fn(p), digits = 12),
    " flat run's", format(ds$obj$fn(p) %||% NA, digits = 12),
    " MAP fit's", format(fm$obj$fn(p), digits = 12), "\n")
