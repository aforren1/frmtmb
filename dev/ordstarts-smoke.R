# Lane ordstarts, build round: a quick look at the build on a few known
# cases. arm = lane (wt-ordstarts-lib first) or base (rellib-r22).
#   Rscript dev/ordstarts-smoke.R <arm>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
nobs <- 1000
G <- 50
re <- "int"
run <- readLines("dev/ordstarts-run.R")
eval(parse(text = run[grep("^mkg <- ", run):(grep("^form <- ", run) - 1L)]))
one <- function(lab, expr) {
  w <- character()
  m <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, class(x)[1])
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, class(x)[1])
    invokeRestart("muffleMessage")
  })
  el <- proc.time()[["elapsed"]] - t0
  ms <- f$opt$mixture_starts
  cat(sprintf("%-28s ll %.4f code %d el %.1fs | ms %s | warn %s | msg %s\n",
              lab, as.numeric(logLik(f)), f$opt$convergence, el,
              if (is.null(ms)) "NULL" else
                paste(names(ms), signif(ms, 4), sep = "=", collapse = " "),
              paste(unique(w), collapse = ","),
              paste(unique(m), collapse = ",")))
  invisible(f)
}
cs <- mixture(cumulative(), sratio())
cc <- mixture(cumulative(), cumulative())
for (cfg in list(c(1000, 13), c(1000, 34), c(300, 28))) {
  d <- mk(cfg[2], n = cfg[1], mix = TRUE)
  lab <- paste0("noRE cs n", cfg[1], " s", cfg[2])
  one(paste(lab, "default"), frm(bf(y ~ x + z), family = cs, data = d))
  one(paste(lab, "starts0"), frm(bf(y ~ x + z), family = cs, data = d,
                                 control = frmtmb_control(starts = 0)))
  one(paste(lab, "starts10"), frm(bf(y ~ x + z), family = cs, data = d,
                                  control = frmtmb_control(starts = 10)))
}
d <- mkg(13)
one("RE cs n1000 s13 default", frm(bf(y ~ x + z + (1 | g)), family = cs,
                                   data = d))
d <- mk(4, n = 300, mix = TRUE)
one("noRE cc n300 s4 default", frm(bf(y ~ x + z), family = cc, data = d))
