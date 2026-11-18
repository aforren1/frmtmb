# rel084-kt copy of dev/reviews/cdfname-review-rds-b.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed); the .rds is read from the cdfname lane tree by absolute path.
# Review probe, arm B: load arm A's 0.83.0 fits on the lane's core and
# call each post-fit method. A method must either refuse with the
# rename class or give the same answer as on 0.83.0; a silent change
# (a censored row read as if uncensored) is the defect looked for.
#   Rscript dev/reviews/cdfname-review-rds-b.R [base]
args <- commandArgs(trailingOnly = TRUE)
base <- length(args) && args[[1L]] == "base"
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
o <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-cdfname/dev/reviews/cdfname-review-old-fits.rds")
# where the fit keeps family objects carrying the old element
find_old <- function(x, path = "fit", depth = 0L) {
  if (depth > 6L || !is.list(x)) return(character(0))
  out <- if ("lcdf" %in% names(x)) path else character(0)
  nms <- names(x) %||% rep("", length(x))
  for (i in seq_along(x)) {
    nm <- if (nzchar(nms[i])) paste0("$", nms[i]) else paste0("[[", i, "]]")
    out <- c(out, find_old(.subset2(x, i), paste0(path, nm), depth + 1L))
  }
  out
}
show <- function(nm, what, expr) {
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(force(expr)))
    if (is.numeric(v)) sprintf("OK %s", paste(format(head(as.numeric(v), 3),
                                                     digits = 10),
                                              collapse = " "))
    else paste("OK", class(v)[1L])
  }, error = function(e) {
    sprintf("ERROR [%s] %s", paste(class(e)[1:2], collapse = "/"),
            substr(gsub("\n", " ", conditionMessage(e)), 1, 140))
  })
  cat(sprintf("B %-14s %-22s %s\n", nm, what, r))
}
for (nm in names(o$fits)) {
  fit <- o$fits[[nm]]
  dd <- if (grepl("weibull", nm)) o$d else o$dp
  cat("PATHS", nm, paste(head(find_old(unclass(fit)), 6), collapse = " "),
      "\n")
  show(nm, "logLik", logLik(fit))
  show(nm, "predict()", predict(fit))
  show(nm, "fitted()", fitted(fit))
  show(nm, "predict(newdata)", predict(fit, newdata = dd[1:5, ]))
  show(nm, "residuals(osa)", residuals(fit, type = "osa"))
  show(nm, "simulate()", simulate(fit, nsim = 1, seed = 1)[[1]])
  show(nm, "update()", logLik(update(fit, data = dd)))
}
