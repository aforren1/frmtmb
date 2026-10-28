# Reviewer probe 7: ndt of shifted_lognormal under the default and under
# prior = "flat"; brms declares real<lower=0, upper=min_Y> ndt. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1)
d <- data.frame(x = rnorm(60))
d$rt <- 0.3 + exp(rnorm(60, -0.5 + 0.2 * d$x, 0.4))
fit <- q(frm(bf(rt ~ x), data = d, family = shifted_lognormal(),
             dry_run = "objective"))
for (pr in list(NULL, "flat")) {
  r <- q(frmtmb.sample:::sample_resolve_priors(fit, pr))
  cat("-- prior =", if (is.null(pr)) "default" else pr, "\n")
  for (e in r$ri$entries) {
    cat(sprintf("  entry %-6s %-3s %-8s %-8s link=%s lower=%s upper=%s lb=%s ub=%s\n",
                e$comp, paste(e$idx, collapse = ","), e$dist$kind, e$scale,
                e$link$name %||% "-", format(e$dist$lower %||% NA),
                format(e$dist$upper %||% NA), format(e$lb %||% NA),
                format(e$ub %||% NA)))
  }
  if (length(r$ri$lower)) cat("  lower:", names(r$ri$lower), r$ri$lower, "\n")
  if (length(r$ri$upper)) cat("  upper:", names(r$ri$upper), r$ri$upper, "\n")
}
cat("min rt:", min(d$rt), "\n")
