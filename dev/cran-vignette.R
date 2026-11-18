# Lane cran: knit one vignette as R CMD check rebuilds it (NOT_CRAN
# unset), with every chunk timed, to show which chunks to precompute or
# set eval = FALSE.
#
#   Rscript dev/cran-vignette.R <vignette.Rmd> <scratch dir>
#
# The vignette directory is copied into the scratch directory first, so
# the knit writes nothing into the worktree. Prints one CHUNK line per
# chunk (label, CPU s, elapsed s) and one VIGNETTE line with the totals.
local({
t_proc0 <- proc.time()
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(!nzchar(Sys.getenv("NOT_CRAN")))
a <- commandArgs(trailingOnly = TRUE)
src <- normalizePath(a[1], winslash = "/")
work <- file.path(a[2], sub("[.]Rmd$", "", basename(src)))
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
file.copy(list.files(dirname(src), full.names = TRUE), work,
          recursive = TRUE)
ctl <- system.time({
  s <- 0
  for (i in seq_len(4e6)) s <- s + sqrt(i)
})
rows <- list()
t_start <- NULL
knitr::knit_hooks$set(cran_timer = function(before, options) {
  if (before) {
    t_start <<- proc.time()
  } else {
    d <- proc.time() - t_start
    rows[[length(rows) + 1L]] <<- data.frame(
      label = options$label, cpu = sum(d[c(1, 2)]), real = d[[3]])
  }
  NULL
})
knitr::opts_chunk$set(cran_timer = TRUE)
t0 <- proc.time()
ok <- tryCatch({
  rmarkdown::render(file.path(work, basename(src)), quiet = TRUE,
                    envir = new.env(parent = globalenv()))
  TRUE
}, error = function(e) {
  cat("RENDERERROR ", conditionMessage(e), "\n", sep = "")
  FALSE
})
t1 <- proc.time() - t0
d <- do.call(rbind, rows)
if (!is.null(d)) {
  d <- d[order(-d$real), ]
  for (j in seq_len(nrow(d))) {
    cat("CHUNK\t", basename(src), "\t", d$label[j], "\t",
        sprintf("%.2f", d$cpu[j]), "\t", sprintf("%.2f", d$real[j]), "\n",
        sep = "")
  }
}
tp <- proc.time() - t_proc0
cat("VIGNETTE\t", basename(src), "\tok=", ok,
    "\tchunks=", if (is.null(d)) 0 else nrow(d),
    "\trender_cpu=", sprintf("%.2f", sum(t1[c(1, 2)])),
    "\trender_real=", sprintf("%.2f", t1[[3]]),
    "\tproc_real=", sprintf("%.2f", tp[[3]]),
    "\tctl_real=", sprintf("%.2f", ctl[[3]]), "\n", sep = "")
})
