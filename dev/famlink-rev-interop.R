## Reviewer check for lane wt-famlink, priority 2: does the exact `$`
## method, or the new `link` field, change what emmeans, insight or
## marginaleffects return on a frmtmb fit?
##
## Usage: Rscript dev/famlink-rev-interop.R <lane|base>
## Writes dev/famlink-rev-interop-<arm>.rds (every result, or its error
## message) and, for base, dev/famlink-rev-savedfit-base.rds.
## Seed 20260916. Data: n = 200, 20 groups.
ARM <- commandArgs(trailingOnly = TRUE)[1]
source("dev/famlink-rev-common.R")
suppressMessages({library(emmeans); library(insight); library(marginaleffects)})
set.seed(20260916)
n <- 200
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(rep(1:20, each = 10)))
re <- rnorm(20, 0, 0.4)[d$g]
d$yg <- 1 + 0.5 * d$x + as.numeric(d$f) * 0.3 + rnorm(n)
d$tr <- sample(5:15, n, TRUE)
d$yb <- rbinom(n, d$tr, plogis(-0.2 + 0.6 * d$x + re))
d$y01 <- rbinom(n, 1, plogis(0.3 * d$x))
d$ynb <- rnbinom(n, mu = exp(1 + 0.4 * d$x + re), size = 2)
d$yzi <- ifelse(runif(n) < 0.3, 0, rpois(n, exp(0.8 + 0.3 * d$x)))
fits <- list(
  gaussian = frm(yg ~ x + f, data = d),
  binomial = frm(yb | trials(tr) ~ x + f, data = d, family = binomial()),
  bernoulli = frm(y01 ~ x + f, data = d, family = bernoulli()),
  negbin_re = frm(ynb ~ x + f + (1 | g), data = d, family = negbinomial()),
  zip = frm(bf(yzi ~ x + f, zi ~ 1), data = d,
            family = zero_inflated_poisson())
)
if (identical(ARM, "base")) {
  saveRDS(fits, "dev/famlink-rev-savedfit-base.rds")
}
try1 <- function(expr) {
  tryCatch(withCallingHandlers(expr, warning = function(w)
    invokeRestart("muffleWarning")),
    error = function(e) paste("ERROR:", conditionMessage(e)))
}
out <- list()
for (nm in names(fits)) {
  fit <- fits[[nm]]
  r <- list()
  r$family_link <- try1(family(fit)$link)
  r$model_info <- try1(unclass(insight::model_info(fit)))
  r$link_inverse <- try1(insight::link_inverse(fit)(c(-1, 0, 1)))
  r$find_statistic <- try1(insight::find_statistic(fit))
  r$get_predicted <- try1(as.numeric(insight::get_predicted(fit)))
  r$get_predicted_resp <- try1(as.numeric(insight::get_predicted(
    fit, predict = "expectation")))
  r$emmeans <- try1(as.data.frame(emmeans::emmeans(fit, "f")))
  r$emmeans_resp <- try1(as.data.frame(emmeans::emmeans(fit, "f",
                                                        type = "response")))
  r$ref_grid <- try1(as.data.frame(summary(emmeans::ref_grid(fit))))
  r$predictions <- try1(as.data.frame(marginaleffects::predictions(fit))[
    , c("estimate", "std.error")])
  r$avg_slopes <- try1(as.data.frame(marginaleffects::avg_slopes(fit))[
    , c("term", "contrast", "estimate", "std.error")])
  r$comparisons <- try1(as.data.frame(marginaleffects::avg_comparisons(
    fit, variables = "f"))[, c("contrast", "estimate", "std.error")])
  r$logLik <- as.numeric(logLik(fit))
  r$coef <- try1(fit$estimates)
  r$se <- try1(sqrt(diag(vcov(fit))))
  r$fitted <- try1(as.numeric(as.matrix(fitted(fit))[, 1]))
  out[[nm]] <- r
}
saveRDS(out, sprintf("dev/famlink-rev-interop-%s.rds", ARM))
cat("done", ARM, "\n")
for (nm in names(out)) {
  errs <- vapply(out[[nm]], function(v) is.character(v) && length(v) == 1 &&
                   startsWith(v, "ERROR:"), TRUE)
  cat(nm, ": errors in", paste(names(errs)[errs], collapse = ", "), "\n")
  for (k in names(errs)[errs]) cat("   ", k, out[[nm]][[k]], "\n")
}
