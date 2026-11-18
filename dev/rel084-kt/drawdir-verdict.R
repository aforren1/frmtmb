# rel084-kt copy of dev/drawdir-verdict.R: rellib-r22 removed from .libPaths() (lib = argument 1, then user lib). No other change.
# Lane drawdir: the draw check's verdict on a build, per fit. Fits each
# case (dev/drawdir-cases.R), calls predict(ndraws = 200, seed 5), and
# records the cached verdict (summed: fire, excess, rel, n; per
# direction: dir, dir_par, dir_n), what was said (classes), and the time
# of the first predict() and of a second one.
# Usage: Rscript dev/drawdir-verdict.R <lib or "base"> <out.tsv> <group>
#        [id filter regex]
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1]
.libPaths(c(lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
out <- args[2]
source("dev/drawdir-cases.R")
ids <- dc_ids(args[3])
if (length(args) > 3) ids <- grep(args[4], ids, value = TRUE)
rows <- list()
for (id in ids) {
  r <- tryCatch(dc_fit(id), error = function(e) NULL)
  if (is.null(r) || is.null(r$fit)) {
    cat(id, "fit failed\n")
    next
  }
  fit <- r$fit
  said <- character(0)
  one <- function() {
    set.seed(5)
    withCallingHandlers(
      predict(fit, ndraws = 200, summary = FALSE),
      warning = function(w) {
        said <<- c(said, class(w)[1])
        invokeRestart("muffleWarning")
      }, message = function(m) {
        said <<- c(said, class(m)[1])
        invokeRestart("muffleMessage")
      })
  }
  t0 <- proc.time()[[3]]
  ok <- !is.null(tryCatch(one(), error = function(e) NULL))
  t1 <- proc.time()[[3]]
  tryCatch(one(), error = function(e) NULL)
  t2 <- proc.time()[[3]]
  chk <- fit$cache[["draw_check"]]
  g <- function(x, d) if (is.null(x)) d else x
  row <- data.frame(id = id, np = length(fit$opt$par), ok = ok,
                    verdict = !is.null(chk), fire = isTRUE(chk$fire),
                    excess = g(chk$excess, NA), rel = g(chk$rel, NA),
                    n = g(chk$n, NA), dir = g(chk$dir, NA),
                    dir_par = g(chk$dir_par, NA), dir_n = g(chk$dir_n, NA),
                    t_first = t1 - t0, t_second = t2 - t1,
                    said = paste(unique(said), collapse = ","))
  rows[[length(rows) + 1L]] <- row
  cat(sprintf("%-40s np %3d fire %-5s ex %7.3f rel %6.3f n %2s | dir %8.3f (%s, %s evals) | t %.2f %.2f | %s\n",
              id, row$np, row$fire, row$excess, row$rel, row$n, row$dir,
              row$dir_par, row$dir_n, row$t_first, row$t_second, row$said))
}
utils::write.table(do.call(rbind, rows), out, sep = "\t", row.names = FALSE,
                   quote = FALSE)
cat("DONE", length(rows), "\n")
