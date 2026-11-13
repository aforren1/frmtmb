# Lane ordstart, item 4: what the default starts cost. In one process,
# interleaved, starts = 0 (the 0.81.0 path) against the default (10
# starts), on the two-class data of dev/optima-csmix.R: n = 300 and
# n = 3000 with y ~ x + z, and n = 1000 with a random intercept over 50
# groups. Per arm: the minimum elapsed seconds over 3 rounds, and the
# objective and gradient evaluations the optimizer runs asked for
# (fit$opt$evals), which do not depend on the load.
#   Rscript dev/ordstart-cost.R out.txt [seeds]
args <- commandArgs(trailingOnly = TRUE)
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:3
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
designs <- list(
  n300 = list(n = 300, f = bf(y ~ x + z)),
  n3000 = list(n = 3000, f = bf(y ~ x + z)),
  re1000 = list(n = 1000, f = bf(y ~ x + z + (1 | g))))
quiet <- function(e) suppressWarnings(suppressMessages(e))
for (s in seeds) {
  for (nm in names(designs)) {
    ds <- designs[[nm]]
    d <- mk(s, n = ds$n, mix = TRUE)
    d$g <- factor(rep_len(1:50, ds$n))
    res <- list()
    for (r in 1:3) {
      for (arm in c("off", "on")) {
        ctl <- if (arm == "off") frmtmb_control(starts = 0) else
          frmtmb_control()
        t0 <- proc.time()[["elapsed"]]
        f <- quiet(frm(ds$f, family = mixture(cumulative(), cumulative()),
                       data = d, control = ctl))
        el <- proc.time()[["elapsed"]] - t0
        res[[arm]] <- rbind(res[[arm]], data.frame(
          el = el, evals = f$opt$evals %||% NA, ll = as.numeric(logLik(f)),
          moved = isTRUE(f$opt$mixture_starts[["gain"]] != 0)))
      }
    }
    cat(sprintf(paste("COST seed %d %-7s off %.2fs %5d evals | on %.2fs %5d",
                      "evals | ratio time %.1f evals %.1f | ll off %.4f on",
                      "%.4f\n"),
                s, nm, min(res$off$el), res$off$evals[1], min(res$on$el),
                res$on$evals[1], min(res$on$el) / min(res$off$el),
                res$on$evals[1] / res$off$evals[1], res$off$ll[1],
                res$on$ll[1]))
  }
}
