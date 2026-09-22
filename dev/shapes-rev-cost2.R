# Reviewer, priority 5, part 2: the paths that refit MANY times. Each
# probe starts from a FRESH fit so that no earlier call can have warmed
# the cache, which is what hid the cost in part 1.
#
#   Rscript dev/shapes-rev-cost2.R lane|base

arg <- commandArgs(trailingOnly = TRUE)
which_lib <- if (identical(arg[1], "base")) "base" else "lane"
lib <- if (which_lib == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", lib, "\n")

CNT <- new.env(parent = emptyenv()); CNT$n <- 0L
ns <- asNamespace("frmtmb")
orig <- get("autoscale_sdreport", envir = ns)
unlockBinding("autoscale_sdreport", ns)
assign("autoscale_sdreport", function(...) { CNT$n <- CNT$n + 1L
  orig(...) }, envir = ns)
lockBinding("autoscale_sdreport", ns)

probe <- function(lab, e) {
  CNT$n <- 0L
  t0 <- proc.time()
  v <- tryCatch(suppressWarnings(e), error = function(c)
    paste0("ERROR: ", substr(conditionMessage(c), 1, 70)))
  cat(sprintf("  %-44s sdreports=%4d  %7.2f s  %s\n", lab, CNT$n,
              (proc.time() - t0)[3],
              if (is.character(v)) v else class(v)[1]))
}

set.seed(99)
n <- 600L
d <- data.frame(x = rnorm(n), g = factor(rep(1:30, length.out = n)))
d$y <- rnorm(n, 1 + 0.5 * d$x + rnorm(30, 0, 0.7)[d$g], 1)
mk <- function() frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)

cat("\n== one accessor, on a FRESH fit each time ==\n")
probe("fixef(fresh)", fixef(mk()))
probe("fixef(fresh, flatten = TRUE)",
      if (which_lib == "lane") fixef(mk(), flatten = TRUE) else fixef(mk()))
probe("ranef(fresh)", ranef(mk()))
probe("VarCorr(fresh)", VarCorr(mk()))
probe("coef(fresh)", coef(mk()))
probe("simulate(fresh, nsim = 2)", simulate(mk(), nsim = 2))
probe("logLik(fresh)", logLik(mk()))
probe("print(fresh)", { f <- mk()
  invisible(utils::capture.output(print(f))) })

cat("\n== the refitting paths ==\n")
probe("frm_bootstrap(fit, B = 8)",
      frm_bootstrap(mk(), R = 8, seed = 1))
probe("frm_allfit(fit)", frm_allfit(mk()))
probe("influence(fit)", influence(mk()))
probe("diagnose(fit)", diagnose(mk()))
