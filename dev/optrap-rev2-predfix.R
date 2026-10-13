# Reviewer of lane optrap, re-check: test-predfix.R's poisson y ~ 0 + xs
# at x scale 1e-6 with autoscale = FALSE, against glm(), per build.
#   Rscript dev/optrap-rev2-predfix.R <lib or none>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(23)
n <- 250
xt <- stats::rnorm(n)
d <- data.frame(y = stats::rpois(n, exp(0.5 + 0.4 * xt)),
                g = 1 + 0.4 * xt + stats::rnorm(n), xs = xt * 1e-6)
gl <- as.numeric(logLik(glm(y ~ 0 + xs, family = poisson(), data = d)))
w <- character()
f <- withCallingHandlers(frm(bf(y ~ 0 + xs) + poisson(), data = d,
                             control = frmtmb_control(autoscale = FALSE)),
  warning = function(x) {w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")})
cat("lib", dirname(find.package("frmtmb")), "| code", f$opt$convergence,
    "| logLik - glm", signif(as.numeric(logLik(f)) - gl, 4), "| polish",
    paste(f$opt$polish %||% "none", collapse = " "), "| warnings", length(w), "\n")
