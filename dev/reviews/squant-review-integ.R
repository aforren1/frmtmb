# Review of lane squant, items 4 and 6: what summary() and confint()
# report, the sandwich through summary(vcov =); brms interop
# (stancode, brmsfamily, default priors); frm_compat(); predict, fitted,
# simulate, log_lik, loo, osa, cens, trunc, mi (response and predictor);
# mixture(quantile_elf(), quantile_elf()); frmtmb.sample's frm_sample();
# the weights check where its condition is absent.
# Run (worktree root):
#   Rscript dev/reviews/squant-review-integ.R > dev/reviews/squant-review-integ.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
try_ <- function(label, expr) {
  ws <- character()
  v <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(gsub("\\s+", " ", conditionMessage(w)), 1, 160))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ws <<- c(ws, paste("MESSAGE:", substr(gsub("\\s+", " ",
                                               conditionMessage(m)), 1, 160)))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  if (inherits(v, "error")) {
    cat(sprintf("[%s] ERROR: %s\n", label,
                substr(gsub("\\s+", " ", conditionMessage(v)), 1, 300)))
  } else {
    cat(sprintf("[%s] ok%s\n", label,
                if (length(ws)) paste0(" | conditions: ",
                                       paste(unique(ws), collapse = " || "))
                else ""))
  }
  invisible(v)
}
set.seed(601)
n <- 300
d <- data.frame(x = runif(n), row = seq_len(n))
d$y <- 1 + 2 * d$x + (0.5 + d$x) * rnorm(n)
fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(), data = d)

cat("\n== summary(fit)\n")
print(summary(fit))
cat("\n== confint(fit) (default)\n")
print(try_("confint", confint(fit)))
cat("\n== summary(fit, vcov = vcov_cluster(fit, ~ row, 'CR0', full = TRUE))\n")
s2 <- try_("summary vcov", summary(fit, vcov = vcov_cluster(fit, ~ row, "CR0",
                                                             full = TRUE)))
if (!inherits(s2, "error")) print(s2)
try_("confint vcov", print(confint(fit, vcov = vcov_cluster(fit, ~ row, "CR0",
                                                            full = TRUE))))
try_("vcov(cluster = ~ row)", print(sqrt(diag(vcov(fit, cluster = ~ row)))))

cat("\n== brms interop\n")
try_("stancode(fit)", stancode(fit))
try_("standata(fit)", standata(fit))
try_("frmtmb::brmsfamily('quantile_elf')",
     print(brmsfamily("quantile_elf")$family))
if (requireNamespace("brms", quietly = TRUE)) {
  try_("brms::brmsfamily('quantile_elf')", brms::brmsfamily("quantile_elf"))
}
try_("default_prior", print(default_prior(bf(y ~ x, quantile = 0.9),
                                          family = quantile_elf(), data = d)))
try_("get_prior", print(get_prior(bf(y ~ x), family = quantile_elf(),
                                  data = d)))

cat("\n== frm_compat()\n")
try_("frm_compat('quantile_elf')", print(head(frm_compat("quantile_elf"))))
ft <- frm_compat_features()
cat("features named quantile_elf:", sum(grepl("quantile_elf", ft$name)),
    "| asym_laplace:", sum(grepl("asym_laplace", ft$name)), "\n")
try_("frm_compat('asym_laplace', 'cens()')",
     print(frm_compat("asym_laplace", "cens()")))

cat("\n== post-fit methods\n")
try_("fitted", print(head(fitted(fit), 2)))
try_("fitted dpar quantile", print(head(fitted(fit, dpar = "quantile"), 2)))
try_("predict", print(head(predict(fit), 2)))
try_("posterior_predict", print(dim(posterior_predict(fit, ndraws = 20))))
try_("posterior_epred", print(dim(posterior_epred(fit, ndraws = 20))))
try_("simulate", print(dim(simulate(fit, nsim = 3, seed = 1))))
try_("log_lik", print(dim(log_lik(fit))))
try_("loo", print(loo(fit)$estimates))
try_("residuals osa", print(summary(residuals(fit, type = "osa")[, 1])))
try_("residuals pearson", print(head(residuals(fit, type = "pearson"), 2)))
try_("residuals deviance", print(head(residuals(fit, type = "deviance"), 2)))
try_("conditional_effects", conditional_effects(fit))
try_("predict newdata", print(predict(fit, newdata = data.frame(x = 0.5))))
try_("hypothesis", print(hypothesis(fit, "x > 0")))

cat("\n== addition terms\n")
d$yt <- d$y
keep <- d$y > 0
dt <- d[keep, ]
try_("trunc(lb = 0)", {
  ft2 <- frm(bf(yt | trunc(lb = 0) ~ x, quantile = 0.5),
             family = quantile_elf(), data = dt)
  cat(" conv", ft2$opt$convergence, "logLik", as.numeric(logLik(ft2)), "\n")
  print(head(predict(ft2), 2))
})
d$ci <- 0
d$ci[d$y > 3.5] <- 1
d$ub <- pmin(d$y, 3.5)
try_("cens right", {
  fc <- frm(bf(ub | cens(ci) ~ x, quantile = 0.5), family = quantile_elf(),
            data = d)
  cat(" conv", fc$opt$convergence, "\n")
  print(summary(residuals(fc, type = "osa")[, 1]))
})
d$lo <- d$y - 0.3
d$hi <- d$y + 0.3
d$cint <- "interval"
try_("cens interval", {
  fi <- frm(bf(lo | cens(cint, hi) ~ x, quantile = 0.5),
            family = quantile_elf(), data = d)
  cat(" conv", fi$opt$convergence, "\n")
})
dm <- d
dm$ym <- dm$y
dm$ym[sample(n, 20)] <- NA
try_("mi() response", {
  fm <- frm(bf(ym | mi() ~ x, quantile = 0.5), family = quantile_elf(),
            data = dm)
  cat(" conv", fm$opt$convergence, "nobs", nobs(fm), "\n")
})
dm$z <- dm$x + rnorm(n, 0, 0.3)
dm$z[sample(n, 25)] <- NA
try_("mi() predictor of quantile_elf", {
  fz <- frm(bf(y ~ mi(z)) + bf(z | mi() ~ x, quantile = 0.5,
                               family = quantile_elf()) + set_rescor(FALSE),
            family = gaussian(), data = dm)
  cat(" conv", fz$opt$convergence, "\n")
})
try_("weights 2, no group term", {
  d$w <- 2
  fw <- frm(bf(y | weights(w) ~ x, quantile = 0.5), family = quantile_elf(),
            data = d)
  cat(" conv", fw$opt$convergence, "\n")
})
dg <- d
dg$g <- factor(rep(seq_len(n / 4), each = 4))
dg$w <- 1.2
try_("weights 1.2, (1 | g) of 4 rows", {
  fg <- frm(bf(y | weights(w) ~ x + (1 | g), quantile = 0.5),
            family = quantile_elf(smoothing = 1), data = dg)
  cat(" conv", fg$opt$convergence, "\n")
})

cat("\n== mixture(quantile_elf(), quantile_elf())\n")
set.seed(602)
dmx <- data.frame(x = runif(400))
z <- rbinom(400, 1, 0.4)
dmx$y <- ifelse(z == 1, 1 + dmx$x + rnorm(400, 0, 0.5),
                5 + 2 * dmx$x + rnorm(400, 0, 0.5))
try_("mixture", {
  fmx <- frm(bf(y ~ x, quantile1 = 0.5, quantile2 = 0.5),
             family = mixture(quantile_elf(), quantile_elf()), data = dmx)
  cat(" conv", fmx$opt$convergence, "\n")
  print(fixef(fmx))
  print(head(fitted(fmx), 2))
})

cat("\n== frmtmb.sample::frm_sample()\n")
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  cat("frmtmb.sample from", find.package("frmtmb.sample"), "\n")
  fs <- try_("frm_sample", frmtmb.sample::frm_sample(fit, chains = 2,
                                                     iter = 800, seed = 1,
                                                     refresh = 0))
  if (!inherits(fs, "error")) {
    print(fixef(fs))
    cat("ML:\n")
    print(fixef(fit))
  }
}
cat("DONE\n")
