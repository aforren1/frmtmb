# Reviewer, item 1 continued: per-threshold rows and densities on a
# multivariate ordinal model and on a mixture with thres(gr = ), against
# brms 2.23.0's default_prior() and a hand density. Seed 20261007.
#   Rscript dev/relrev-mix2.R > dev/relrev-log/mix2.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
set.seed(20261007)
n <- 500
d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
d$y2 <- 1L + (0.5 * d$x + rlogis(n) > 0) + (0.5 * d$x + rlogis(n) > 1)
rows <- function(p) {
  p <- as.data.frame(p)
  keep <- p$class %in% c("Intercept", "delta")
  sort(paste(p$class, p$coef, p$group, p$dpar, p$resp, sep = "|")[keep])
}
q <- function(e) tryCatch(suppressWarnings(suppressMessages(e)),
                          error = function(err) paste("ERR", conditionMessage(err)))
cmp <- function(tag, fr, br) {
  if (is.character(fr) && length(fr) == 1 && startsWith(fr, "ERR")) { cat(tag, "frmtmb", fr, "\n"); return() }
  if (is.character(br) && length(br) == 1 && startsWith(br, "ERR")) { cat(tag, "brms", br, "\n"); return() }
  cat(sprintf("%-30s %s (%d / %d rows)\n", tag, if (identical(fr, br)) "SAME" else "DIFF", length(fr), length(br)))
  if (!identical(fr, br)) cat("  only frmtmb:", setdiff(fr, br), "\n  only brms:", setdiff(br, fr), "\n")
}
cmp("mv two cumulative",
    q(rows(default_prior(mvbf(bf(y ~ x), bf(y2 ~ x)), data = d, family = cumulative()))),
    q(rows(brms::default_prior(brms::mvbf(brms::bf(y ~ x), brms::bf(y2 ~ x)), data = d, family = brms::cumulative()))))
cmp("mixture none, thres(gr = g)",
    q(rows(default_prior(bf(y | thres(gr = g) ~ x), data = d, family = mixture(cumulative(), sratio())))),
    q(rows(brms::default_prior(brms::bf(y | thres(gr = g) ~ x), data = d, family = brms::mixture(brms::cumulative(), brms::sratio())))))
cmp("mixture mu, thres(gr = g)",
    q(rows(default_prior(bf(y | thres(gr = g) ~ x), data = d, family = mixture(cumulative(), sratio(), order = "mu")))),
    q(rows(brms::default_prior(brms::bf(y | thres(gr = g) ~ x), data = d, family = brms::mixture(brms::cumulative(), brms::sratio(), order = "mu")))))
lp_of <- function(fit) {
  ent <- frmtmb:::resolve_prior_input(list(frame = fit$frame, spec = fit$spec), fit$prior)$entries
  -frmtmb:::neg_log_prior_fn(ent)(fit$estimates)
}
# mv: coef 2 on y2's thresholds (centered by y2's own slope)
fit <- q(frm(mvbf(bf(y ~ x), bf(y2 ~ x)), data = d, family = cumulative(),
             prior = set_prior("normal(0.5, 0.3)", class = "Intercept", coef = "2", resp = "y2")))
if (is.character(fit)) cat("mv fit", fit, "\n") else {
  e <- fit$estimates
  cat("mv estimate names:", names(e), "\n")
  nm <- grep("tau_raw", names(e), value = TRUE)
  r2 <- e[[nm[grepl("y2", nm)]]]
  t2 <- c(r2[1], r2[1] + cumsum(exp(r2[-1])))
  b2 <- e$beta[grep("y2", names(e$beta))]
  cat("y2 slope name:", names(b2), "\n")
  ref <- dnorm(t2[2] - mean(d$x) * b2[[1]], 0.5, 0.3, log = TRUE) + sum(r2[-1])
  got <- lp_of(fit)
  cat(sprintf("mv coef 2 resp y2: frmtmb %.15g hand %.15g rel %.3g\n", got, ref, abs(got - ref) / abs(ref)))
}
# mixture none with thres(gr = g): coef 3 of level b on mu1 (cumulative):
# brms does not center a model with grouped thresholds
fit <- q(frm(bf(y | thres(gr = g) ~ x), data = d, family = mixture(cumulative(), sratio()),
             prior = set_prior("normal(1, 0.5)", class = "Intercept", coef = "3", group = "b", dpar = "mu1")))
if (is.character(fit)) cat("gr fit", fit, "\n") else {
  r <- fit$estimates$tau_raw1
  cat("tau_raw1 length", length(r), "\n")
  tb <- c(r[4], r[4] + cumsum(exp(r[5:6])))
  ref <- dnorm(tb[3], 1, 0.5, log = TRUE) + sum(r[5:6])
  got <- lp_of(fit)
  cat(sprintf("mixture gr coef 3 group b mu1: frmtmb %.15g hand %.15g rel %.3g\n", got, ref, abs(got - ref) / abs(ref)))
}
