# Review of lane ordstart, two probes on the lane build.
#   Rscript ordstart-review-paths.R autoscale out.txt
#     the 40-fit data with x rescaled by 1e-4 (sd below autoscale's
#     1e-3 cut, so frm() runs its autoscale pre-fit): does the step-start
#     warning, the mixture_starts record and the answer survive?
#   Rscript ordstart-review-paths.R permute out.txt
#     the same fits with the data rows permuted (10 permutations): does
#     the chosen maximum depend on row order?
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("mode", mode, "frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
ns <- asNamespace("frmtmb")
rec <- new.env()
rec$picks <- 0L
suppressMessages(trace("mixture_starts_pick", where = ns, print = FALSE,
  exit = quote(assign("picks", get("picks", envir = .GlobalEnv$rec) + 1L,
                      envir = .GlobalEnv$rec))))
assign("rec", rec, envir = .GlobalEnv)
fam <- function(nm) switch(nm, cum_cum = mixture(cumulative(), cumulative()),
                           cum_sratio = mixture(cumulative(), sratio()))
fitw <- function(f, d, nm, ...) {
  w <- character()
  rec$picks <- 0L
  t0 <- proc.time()[["elapsed"]]
  ft <- withCallingHandlers(frm(f, family = fam(nm), data = d, ...),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(fit = ft, w = w, picks = rec$picks,
       el = proc.time()[["elapsed"]] - t0)
}
desc <- function(r) {
  sprintf("ll %.4f picks %d ms %s stepw %s degen %s hint %s el %.2f",
          as.numeric(logLik(r$fit)), r$picks,
          if (is.null(r$fit$opt$mixture_starts)) "NULL" else "set",
          any(grepl("jittered starts reached a higher", r$w)),
          any(grepl("degenerate boundary", r$w)),
          if (any(grepl("Compare several starts", r$w))) "compare" else
            if (any(grepl("None of", r$w))) "none-of" else "-",
          r$el)
}
cases <- list(cum_cum = c(3, 4, 10, 11, 14, 15, 16, 18),
              cum_sratio = c(3, 4, 6, 9, 14, 15, 18))
if (mode == "autoscale") {
  for (nm in names(cases)) for (s in cases[[nm]]) {
    d <- mk(s, mix = TRUE)
    d$xs <- d$x * 1e-4
    a <- fitw(bf(y ~ x + z), d, nm)
    b <- fitw(bf(y ~ xs + z), d, nm)
    cat(sprintf("AUTOSCALE %s %2d | plain %s | scaled %s | autoscaled %s\n",
                nm, s, desc(a), desc(b), isTRUE(b$fit$autoscaled)))
  }
} else if (mode == "permute") {
  for (nm in names(cases)) for (s in cases[[nm]]) {
    d <- mk(s, mix = TRUE)
    res <- character()
    lls <- numeric()
    for (p in 0:10) {
      dp <- if (p == 0) d else {
        set.seed(500 + p)
        d[sample(nrow(d)), ]
      }
      r <- fitw(bf(y ~ x + z), dp, nm)
      lls <- c(lls, as.numeric(logLik(r$fit)))
      res <- c(res, paste0(sprintf("%.4f", lls[length(lls)]),
                           if (any(grepl("jittered starts reached", r$w)))
                             "W" else "",
                           if (any(grepl("degenerate boundary", r$w)))
                             "D" else ""))
    }
    cat(sprintf("PERMUTE %s %2d distinct %d spread %.4f | %s\n", nm, s,
                length(unique(round(lls, 4))), diff(range(lls)),
                paste(res, collapse = " ")))
  }
}
cat("DONE\n")
