# Release review 0.80.0: the nlstart reviewer's fresh row binomial mm 22,
# which warned on the nlstart lane build and is silent on the merged build
# while 4.61 short. Prints every condition and the fit's state, by arm.
#
#   Rscript dev/review-rel-mm22.R <lib> [<row spec fam,body,seed> ...]
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
src <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-nlstart/dev/review-nlstart-fresh.R")
a <- grep("^F <- list", src)
b <- grep("^cat\\(\"ROWHEAD", src)
eval(parse(text = src[a:(b - 1L)]))
rows <- if (length(args) > 1L) strsplit(args[-1L], ",") else
  list(c("binomial", "mm", "22"))
one <- function(fm, fam, d, opts, start = NULL) {
  cs <- character()
  r <- withr::with_options(opts, tryCatch(withCallingHandlers(
    if (is.null(start)) frm(fm, family = F[[fam]], data = d) else
      frm(fm, family = F[[fam]], data = d, start = start),
    warning = function(w) {
      cs <<- c(cs, paste0("  W[", paste(class(w), collapse = ","), "] ",
                          substr(conditionMessage(w), 1, 400)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cs <<- c(cs, paste0("  M[", class(m)[1L], "] ",
                          substr(conditionMessage(m), 1, 400)))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      cs <<- c(cs, paste("  ERROR", conditionMessage(e)))
      NULL
    }))
  if (is.null(r)) {
    cat(cs, sep = "\n")
    return(invisible(NULL))
  }
  bb <- r$estimates$beta
  cat(sprintf("  ll %.6f conv %d msg %s\n  beta %s\n",
              as.numeric(logLik(r)), r$opt$convergence, r$opt$message,
              paste(sprintf("%s=%.6g", names(bb), bb), collapse = "; ")))
  ob <- r$obj$env$.frmtmb_opt_basis
  cat("  basis on tape:", if (is.null(ob)) "none" else
        paste("nl", isTRUE(ob$nl)), "\n")
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(r)))), error = function(e)
    paste("vcov error:", conditionMessage(e)))
  cat("  se:", format(se, digits = 4), "\n")
  cat(if (length(cs)) cs else "  (no conditions)", sep = "\n")
  invisible(r)
}
for (rw in rows) {
  fam <- rw[1]; body <- rw[2]; seed <- as.numeric(rw[3])
  d <- gen(fam, body, seed * 1000 + match(body, names(B)))
  fm <- form(fam, body)
  cat("ROW", fam, body, seed, "\n")
  cat(" default\n"); r <- one(fm, fam, d, list())
  cat(" probe off\n"); one(fm, fam, d, list(frmtmb.nl_sign_probe = FALSE))
  cat(" opt_basis off\n"); one(fm, fam, d, list(frmtmb.opt_basis = FALSE))
}
cat("MM22 DONE\n")
