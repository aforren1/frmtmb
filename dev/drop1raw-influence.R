# influence(groups = "g") against lme4's influence.merMod() on models
# with a data-dependent basis. lme4 refits with update(subset = ), so
# model.frame() evaluates poly() and scale() on every row first and the
# deletion keeps the fit's basis. Arg: base | lane.
# Log: dev/drop1raw-out/influence-<arm>.txt (seed 5, n = 200, 8 groups)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1,
                g = factor(sample(letters[1:8], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + rnorm(8)[d$g] + rnorm(n)
for (f in list(z ~ poly(x2, 2) + x + (1 | g),
               z ~ scale(x) + log(x2) + (1 | g))) {
  fit <- q(frm(bf(f), data = d))
  ref <- lme4::lmer(f, data = d, REML = FALSE)
  iref <- stats::influence(ref, groups = "g")
  lme <- iref[["fixed.effects[-g]"]]
  for (how in c("default data", "data = d")) {
    inf <- tryCatch(q(if (how == "data = d") influence(fit, groups = "g",
                                                       data = d) else
                        influence(fit, groups = "g")),
                    error = function(e) e)
    if (inherits(inf, "error")) {
      cat(sprintf("%-34s %-13s ERROR %s\n", deparse1(f), how,
                  substr(conditionMessage(inf), 1, 60)))
      next
    }
    fx <- inf$fixed[rownames(lme), seq_len(ncol(lme)), drop = FALSE]
    if (how == "data = d") cat("  columns:", colnames(inf$fixed), "| lme4:", colnames(lme), "
")
    # columns in the same order: intercept then the terms
    full_diff <- abs(fx - lme)
    dfb_ref <- abs(sweep(lme, 2, lme4::fixef(ref)))
    cat(sprintf("%-34s %-13s max |frmtmb - lme4| %.2e, max |dfbeta| %.2e, ratio %.2e\n",
                deparse1(f), how, max(full_diff), max(dfb_ref),
                max(full_diff) / max(dfb_ref)))
  }
}
