## Reviewer: the warning TEXT. Its first clause is kept verbatim from
## the old check ("Large maximum absolute gradient at the optimum (N)"),
## but N is no longer the maximum absolute gradient: grad_warning_msg()
## formats v$proj. On a fit with a bound holding a component the two
## differ, so the sentence names one quantity and prints another.
## usage: Rscript gradcheck-rev-16-msg.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w)
}
## a bound that binds on one coefficient, and a loosened optimizer so
## the components no bound holds are genuinely short
set.seed(9301)
n <- 3000
dd <- data.frame(x = rnorm(n), z = rnorm(n))
dd$y <- rpois(n, exp(0.4 + 0.9 * dd$x - 0.6 * dd$z))
r <- gw(frm(bf(y ~ x + z), family = poisson(), data = dd,
            prior = set_prior("", class = "b", coef = "x", ub = 0.1),
            control = frmtmb_control(
              restarts = 0,
              optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                             iter.max = 2000, eval.max = 2000))))
d <- diagnose(r$fit, quiet = TRUE)
cat("diagnose()$max_grad :", format(d$max_grad, digits = 6), "at",
    d$worst_grad, "\n")
cat("diagnose()$grad_proj:", format(d$grad_proj %||% NA, digits = 6), "\n")
cat("held                :", paste(d$grad_bound_held %||% "", collapse = ","),
    "\n")
cat("headroom            :", format(d$grad_headroom %||% NA, digits = 6),
    "\n\nwarnings:\n")
for (w in r$w) cat("  ", w, "\n")
cat("\nprinted diagnose():\n")
writeLines(paste0("   ", capture.output(diagnose(r$fit))))
cat("DONE msg\n")
