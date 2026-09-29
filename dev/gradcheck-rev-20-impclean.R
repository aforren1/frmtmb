## Reviewer re-check (c), the ABSENT case of the new importance guard.
##
## R/confint.R:1251 sets gv <- NULL on any importance fit, and line 1482
## reads `(degenerate || (!is.null(gv) && !isTRUE(gv$warn)))`. With gv
## NULL and degenerate FALSE that conjunct is FALSE for EVERY importance
## fit, so the clean line is withheld unconditionally rather than only
## when the verdict would have complained. Base gated the same line on
## `out$max_grad < 1e-3`, which a well-converged importance fit passes.
##
## Construct the case the guard is supposed to leave alone: an importance
## fit that converges cleanly and whose max|grad| is under grad_tol.
## usage: Rscript gradcheck-rev-20-impclean.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w)
}
row <- function(tag, r) {
  f <- r$fit
  if (inherits(f, "condition")) {
    cat(sprintf("%-34s ERROR %s\n", tag, conditionMessage(f))); return(NULL)
  }
  d <- diagnose(f, quiet = TRUE)
  out <- capture.output(diagnose(f))
  cat(sprintf(paste0("%-34s conv=%d max_grad=%-11s pdHess=%-5s nbadse=%d",
                     " under_tol=%-5s CLEAN LINE=%-5s nwarn=%d\n"),
              tag, f$opt$convergence, format(d$max_grad, digits = 5),
              d$pdHess, length(d$bad_se),
              d$max_grad < (f$control$grad_tol %||% 1e-3),
              any(grepl("No convergence problems", out)), length(r$w)))
  cat("     other clean conjuncts: singular=", is.null(d$singular),
      " separation=", is.null(d$separation),
      " unbounded=", is.null(d$unbounded_dpar),
      " predictor_scale=", is.null(d$predictor_scale), "\n", sep = "")
  invisible(d)
}

## a small, well-specified GLMM with plenty of draws, so the fit lands on
## its optimum and the Monte Carlo gradient is small
for (spec in list(list(ng = 25, per = 40, nd = 2000L, seed = 9501),
                  list(ng = 30, per = 30, nd = 4000L, seed = 9502),
                  list(ng = 20, per = 50, nd = 2000L, seed = 9503))) {
  set.seed(spec$seed)
  dd <- data.frame(g = factor(rep(seq_len(spec$ng), spec$per)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(spec$ng, 0, 0.5)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  r <- gw(frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd,
              importance = spec$nd, se = TRUE))
  row(sprintf("importance nd=%d seed=%d", spec$nd, spec$seed), r)
  ## the SAME data and control with no importance, as the reference for
  ## what the clean line does on an otherwise identical fit
  r0 <- gw(frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd,
               se = TRUE))
  row(sprintf("  same data, NO importance seed=%d", spec$seed), r0)
}
cat("DONE impclean\n")
