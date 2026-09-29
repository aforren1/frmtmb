## Reviewer, claim 6: the block that is pinned to FLIP when the
## bound-aware standard errors are fixed.
##
## The block asserts expect_false(diagnose(fit)$pdHess) on ONE bounded
## fit. The guarded condition is "the covariance machinery is not
## bound-aware, so the UNCONSTRAINED Hessian at a point on a bound is
## indefinite". Constructing the case where that condition is ABSENT
## needs a bounded fit whose unconstrained Hessian at the constrained
## optimum happens to be positive definite. If one exists, the block
## does not pin bound-awareness: it pins a property of that one design.
## usage: Rscript gradcheck-rev-10-flip.R <core-lib>
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
  list(fit = val, w = w,
       grad = any(grepl("Large maximum absolute gradient", w, fixed = TRUE)))
}
set.seed(101)
n <- 300
dd <- data.frame(x = rnorm(n))
dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
for (cap in c(0.1, 0.5, 1, 1.5, 1.8, 1.9, 1.95, 1.99)) {
  r <- gw(frm(bf(y ~ x), family = gaussian(), data = dd,
              prior = set_prior("", class = "b", ub = cap)))
  d <- diagnose(r$fit, quiet = TRUE)
  H <- r$fit$obj$he(r$fit$opt$par)
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf(paste0("ub=%-5s x=%-9s pdHess=%-5s held=[%s] gradwarn=%-5s",
                     "  eigenvalues %s\n"),
              format(cap), format(r$fit$opt$par[[2]], digits = 6),
              d$pdHess, paste(d$grad_bound_held %||% "", collapse = ","),
              r$grad,
              paste(format(ev, digits = 4), collapse = " ")))
  cat(sprintf("        bad_se=[%s]  expect_false(pdHess) would %s\n",
              paste(d$bad_se, collapse = ","),
              if (isTRUE(d$pdHess)) "FAIL" else "pass"))
}
cat("DONE flip\n")
