# Reviewer of lane ordstarts, build round. Probe 4 (questions 1 and 6):
# the opt-in route. The build's frmtmb_control(starts = 10), which the
# degenerate warning and the docs name, against 0.83.0's starts = 10 (its
# default without random effects), no random effects. Also the compat cell
# 2 of test-compat.R under both.
#   Rscript ordstarts-build-review-p4.R <base|lane> <fam> <n> <seeds> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
which_lib <- args[1]
famname <- args[2]
nobs <- as.integer(args[3])
seeds <- eval(parse(text = args[4]))
out <- args[5]
.libPaths(c(if (which_lib == "lane") "C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
compat <- function() {
  set.seed(20261005)
  n <- 400
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d
}
fam <- function() {
  switch(famname,
    cc = mixture(cumulative(), cumulative()),
    cs = mixture(cumulative(), sratio()),
    compat2 = mixture(cumulative(threshold = "sum_to_zero"), sratio()))
}
form <- if (famname == "compat2") bf(y ~ x) else bf(y ~ x + z)
grab <- function(expr) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage")),
    error = function(e) e)
  list(v = v, w = w, el = proc.time()[["elapsed"]] - t0)
}
first <- !file.exists(out)
for (seed in seeds) {
  d <- if (famname == "compat2") compat() else mk(seed, nobs)
  for (arm in c("s0", "s10")) {
    ctl <- frmtmb_control(starts = if (arm == "s0") 0 else 10)
    r <- grab(frm(form, family = fam(), data = d, control = ctl))
    f <- r$v
    row <- data.frame(lib = which_lib, fam = famname, n = nobs, seed = seed,
                      arm = arm, el = r$el, err = inherits(f, "error"),
                      ll = NA, code = NA, nw = NA, snw = NA, degen = NA,
                      ran = NA)
    if (!inherits(f, "error")) {
      s <- grab(summary(f))
      ms <- f$opt$mixture_starts
      row$ll <- as.numeric(logLik(f))
      row$code <- f$opt$convergence
      row$nw <- length(r$w)
      row$snw <- length(s$w)
      row$degen <- any(grepl("degenerate boundary", c(r$w, s$w),
                             fixed = TRUE))
      row$ran <- if (is.null(ms)) 0 else ms[["starts"]]
      if (famname == "compat2") {
        saveRDS(f$estimates, sub("[.]tsv$", paste0("-", which_lib, "-", arm,
                                                   "-est.rds"), out))
      }
    }
    utils::write.table(row, out, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = first, append = !first)
    first <- FALSE
  }
  cat(which_lib, famname, nobs, "seed", seed, "done\n")
}
cat("DONE p4\n")
