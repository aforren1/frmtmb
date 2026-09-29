## Reviewer re-check 2: the ABSENT case of the new fallback gate.
##
## With no verdict the clean line now falls back to
## `out$max_grad < (fit$control$grad_tol %||% 1e-3)`. That claims to
## honour the setting. Constructing the case where the guarded thing is
## absent means MOVING grad_tol across each fit's own max|grad| and
## seeing the line appear and disappear. A gate that ignored the setting
## would give the same answer in every column.
##
## Seed 9501 max|grad| is 0.0034419 and seed 9503's is 0.00094754
## (dev/gradcheck-rev-20-impclean.R), so 1e-4 is below both, 1e-3 sits
## between them and 1e-2 is above both.
## usage: Rscript gradcheck-rev-21-impgate.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x

mk <- function(ng, per, seed) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(seq_len(ng), per)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.5)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  dd
}
specs <- list(list(ng = 25, per = 40, nd = 2000L, seed = 9501),
              list(ng = 30, per = 30, nd = 4000L, seed = 9502),
              list(ng = 20, per = 50, nd = 2000L, seed = 9503))
tols <- c(1e-4, 1e-3, 1e-2, 1e-1)

cat(sprintf("%-6s %-12s %s\n", "seed", "max|grad|",
            paste(sprintf("gt=%-7s", format(tols)), collapse = " ")))
for (sp in specs) {
  dd <- mk(sp$ng, sp$per, sp$seed)
  line <- character(0)
  gm <- NA_real_
  for (gt in tols) {
    f <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = poisson(),
                              data = dd, importance = sp$nd, se = TRUE,
                              control = frmtmb_control(grad_tol = gt)))
    d <- diagnose(f, quiet = TRUE)
    out <- capture.output(diagnose(f))
    if (is.na(gm)) gm <- d$max_grad
    ## every other conjunct must stay clean, or the column below is not
    ## measuring the gradient gate
    others <- d$convergence == 0 && isTRUE(d$pdHess) &&
      !length(d$bad_se) && is.null(d$singular) &&
      is.null(d$separation) && is.null(d$unbounded_dpar) &&
      is.null(d$predictor_scale)
    cl <- any(grepl("No convergence problems", out))
    line <- c(line, sprintf("%-5s%s", cl, if (others) " " else "*"))
    ## the verdict must stay withheld at every tolerance
    stopifnot(is.na(d$grad_proj %||% NA_real_),
              is.na(d$grad_headroom %||% NA_real_),
              !length(d$grad_bound_held %||% character(0)))
  }
  cat(sprintf("%-6d %-12s %s\n", sp$seed, format(gm, digits = 6),
              paste(sprintf("%-10s", line), collapse = " ")))
  cat(sprintf("%-6s %-12s expected: %s\n", "", "",
              paste(sprintf("%-10s", sprintf("%-5s ", gm < tols)),
                    collapse = " ")))
}
cat("\n* = some other conjunct was not clean, so that cell does not",
    "measure the gradient gate\n")
cat("The line must follow max|grad| < grad_tol in every cell, or the",
    "fallback\nis not honouring frmtmb_control(grad_tol =).\n")
cat("DONE impgate\n")
