# Short check of predraw2's fix round, item 1(c): the three mixture
# designs of dev/predraw2-mixcov.R, seeds 1 to 40, verdict only (it does
# not depend on ndraws), against the lane's final-build tsv
# (dev/predraw2-mc42-<g|p|c>-lane.tsv).
# Usage (from the lane tree): Rscript dev/reviews/predraw2-short-c.R
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
n <- 300L
source("dev/ppdraw-sweep-gen.R", local = TRUE)
des <- c(g = "mixture_gaussian:boundary", p = "mixture_poisson:boundary",
         c = "mixture_cumulative:interior")
for (k in names(des)) {
  ref <- utils::read.delim(sprintf("dev/predraw2-mc42-%s-lane.tsv", k))
  ds <- strsplit(des[[k]], ":")[[1]]
  agree <- 0L; differ <- character(0); nf <- 0L; maxd <- 0
  for (sd in ref$seed) {
    d <- gen(ds[1], ds[2], sd)
    s <- spec(ds[1])
    fit <- tryCatch(q(frm(s[[1]], family = s[[2]], data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) { differ <- c(differ, paste(sd, "fit failed")); next }
    set.seed(1)
    invisible(tryCatch(q(ns$predict_par_drawer(fit, TRUE, 2L)),
                       error = function(e) NULL))
    chk <- fit$cache[["draw_check"]]
    fire <- isTRUE(chk$fire)
    r <- ref[ref$seed == sd, ]
    nf <- nf + fire
    ex <- chk$excess %||% NA_real_
    if (fire == isTRUE(as.logical(r$fire))) agree <- agree + 1L else
      differ <- c(differ, sprintf("%d lane %s now %s", sd, r$fire, fire))
    if (is.finite(ex) && is.finite(r$excess)) maxd <- max(maxd, abs(ex - r$excess))
  }
  cat(sprintf("%s %-28s seeds %d fired now %d (tsv %d) agree %d | max |excess - tsv| %.3g | differ: %s\n",
              k, des[[k]], nrow(ref), nf, sum(as.logical(ref$fire)), agree,
              maxd, paste(differ, collapse = "; ")))
}
})
