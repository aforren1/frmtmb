# Is a one-value mo() refused by frm_sample() on this build? The data of
# dev/mointerop-sample.R (in the mointerop tree), its last case.
#   Rscript dev/rel077-mo1-probe.R [lib]
a <- commandArgs(TRUE)
lib <- if (length(a)) a[1] else "C:/Users/adf44/source/r/rellib-r15"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-rel077/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
d$x1 <- 4
fit <- suppressWarnings(frm(y ~ mo(x1) + z, data = d))
for (tm in frmtmb:::mo_frame_terms(fit)) cat("term", tm$label, "D", tm$D, "\n")
for (route in c("formula", "fit")) {
  r1 <- tryCatch(suppressWarnings(suppressMessages(
    if (route == "formula") {
      frm_sample(y ~ mo(x1) + z, data = d, chains = 1, iter = 200, seed = 3,
                 refresh = 0)
    } else {
      frm_sample(fit, chains = 1, iter = 200, seed = 3, refresh = 0)
    })), error = function(e) e)
  cat(route, "route:", if (inherits(r1, "error")) {
    substr(conditionMessage(r1), 1, 120)
  } else paste("sampled; variables", paste(variables(r1), collapse = " ")),
  "\n")
}
