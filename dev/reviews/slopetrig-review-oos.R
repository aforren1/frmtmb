# Review of lane slopetrig: out-of-sample check of the ratio trigger on
# the designs of dev/reviews/slopetrig-review-designs.R. Per design and
# seed, four arms in one process: ref (the same model on the affine
# twin of the column, default control), F (autoscale = FALSE), T
# (autoscale = TRUE), D (the library's default). Writes one tsv row per
# case; short is read later against the best arm.
# Usage: Rscript slopetrig-review-oos.R <base|lane> <out.tsv> [seeds]
a <- commandArgs(TRUE)
arm <- a[1]
out <- a[2]
seeds <- if (length(a) > 2) eval(parse(text = a[3])) else 1:10
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")), "\n")
source("dev/reviews/slopetrig-review-designs.R")
run1 <- function(expr) {
  nw <- 0L
  wt <- character(0)
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      nw <<- nw + 1L
      wt <<- c(wt, substr(conditionMessage(w), 1, 60))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    return(list(ll = NA_real_, nw = nw, conv = NA_integer_, t = el,
                as = NA, wt = paste("ERROR", substr(conditionMessage(f), 1, 80))))
  }
  list(ll = as.numeric(logLik(f)), nw = nw,
       conv = as.integer(f$opt$convergence), t = el,
       as = isTRUE(f$autoscaled), wt = paste(gsub("[\t\n]", " ", wt),
                                              collapse = " || "))
}
hdr <- TRUE
for (nm in names(designs)) for (seed in seeds) {
  ds <- designs[[nm]]
  d <- ds$data(seed)
  fam <- ds$family
  rr <- if (!is.null(ds$ref)) run1(frm(ds$form, data = ds$ref(d), family = fam))
  rF <- run1(frm(ds$form, data = d, family = fam,
                 control = frmtmb_control(autoscale = FALSE)))
  rT <- run1(frm(ds$form, data = d, family = fam,
                 control = frmtmb_control(autoscale = TRUE)))
  rD <- run1(frm(ds$form, data = d, family = fam))
  row <- data.frame(design = nm, seed = seed,
    ll_ref = rr$ll %||% NA, ll_F = rF$ll, ll_T = rT$ll, ll_D = rD$ll,
    w_ref = rr$nw %||% NA, w_F = rF$nw, w_T = rT$nw, w_D = rD$nw,
    c_F = rF$conv, c_T = rT$conv, c_D = rD$conv, as_D = rD$as,
    t_F = rF$t, t_T = rT$t, t_D = rD$t,
    wt_D = rD$wt, wt_F = rF$wt)
  utils::write.table(row, out, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = hdr, append = !hdr)
  hdr <- FALSE
  cat(nm, seed, "done\n")
}
cat("DONE\n")
