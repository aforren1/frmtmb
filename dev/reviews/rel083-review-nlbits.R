# Release review of 0.83.0: nl_starts. Fits lane nlstop2's mm cases
# (dev/nlstop-fresh.R generator, the starts grid of nlstop2-pole.R) and
# prints the hex of objective and par and the condition classes, so that
# rellib-r22 (nl_starts TRUE and FALSE) and wt-nlstop2-lib compare bitwise.
#   Rscript dev/reviews/rel083-review-nlbits.R <lib> <TRUE|FALSE> <out>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(frmtmb))
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel083")
on <- as.logical(a[2])
out <- file(a[3], "w")
say <- function(...) cat(..., "\n", sep = "", file = out)
say("lib: ", dirname(find.package("frmtmb")), " ", as.character(packageVersion("frmtmb")),
    " nl_starts ", on)
src <- parse("dev/nlstop-fresh.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
ctl <- if (on) NULL else frmtmb_control(nl_starts = FALSE)
for (cs in list(c("bernoulli", 55), c("binomial", 11), c("binomial", 44))) {
  fam <- cs[1]; seed <- as.numeric(cs[2])
  d <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
  fm <- form(fam, "mm")
  starts <- list(default = NULL)
  for (K in c(-6, -5.05, -4.95, -2.5, -0.15, -0.05, 0.5, 50)) {
    for (Vm in c(-2, 0.5, 5)) starts[[sprintf("Vm%g_K%g", Vm, K)]] <- list(beta = c(Vm, K))
  }
  for (sn in names(starts)) {
    cl <- character()
    args <- list(fm, family = F[[fam]], data = d, start = starts[[sn]])
    if (!is.null(ctl)) args$control <- ctl
    r <- tryCatch(withCallingHandlers(do.call(frm, args),
      warning = function(w) { cl <<- c(cl, paste0("W:", class(w)[1])); invokeRestart("muffleWarning") },
      message = function(m) { cl <<- c(cl, paste0("M:", class(m)[1])); invokeRestart("muffleMessage") }),
      error = function(e) e)
    if (inherits(r, "error")) { say("ROW ", fam, " ", seed, " ", sn, " ERR ", substr(conditionMessage(r), 1, 80)); next }
    say("ROW ", fam, " ", seed, " ", sn, " ", sprintf("%a", r$opt$objective), " ",
        paste(sprintf("%a", unname(r$opt$par)), collapse = ","), " code ",
        r$opt$convergence, " ", paste(cl, collapse = "+"))
  }
}
close(out)
