# Lane ordstarts, build round (backlog 6136, the user's decision C0 + D):
# the build against the measurement. Per seed, in one process, the lane
# build's fit with starts = 0 (A0: no jittered and no swapped start, the
# 0.83.0 default with random effects) and its default (D: the swapped
# start where the families differ). For D also the warnings at its first
# summary(). Data: dev/ordstarts-run.R's mkg() for the random-effect
# designs, mk(seed, n, mix = TRUE) for re = none, and the K = 3 data of
# dev/ordstarts-swap.R (k3, sep 2) for family k3.
#   Rscript dev/ordstarts-after.R <family> <re> <n> <G> <seeds> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
re <- args[2]
nobs <- as.integer(args[3])
G <- as.integer(args[4])
seeds <- eval(parse(text = args[5]))
out <- args[6]
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
stopifnot(grepl("wt-ordstarts-lib", find.package("frmtmb")))
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
run <- readLines("dev/ordstarts-run.R")
eval(parse(text = run[grep("^mkg <- ", run):(grep("^form <- ", run) - 1L)]))
sep <- 2
sw <- readLines("dev/ordstarts-swap.R")
eval(parse(text = sw[grep("^mk3 <- ", sw):(grep("^# re = \"none\"", sw) -
                                              1L)]))
form <- switch(re,
  none = bf(y ~ x + z),
  int = bf(y ~ x + z + (1 | g)),
  slope = bf(y ~ x + z + (1 + x | g)))
fam <- function() {
  switch(famname,
    cum_cum = mixture(cumulative(), cumulative()),
    cum_sratio = mixture(cumulative(), sratio()),
    k3 = mixture(cumulative(), sratio(), acat()))
}
clean <- function(s) paste(substr(gsub("[\t\n\r]", " ", s), 1, 120),
                           collapse = " | ")
grab <- function(expr) {
  w <- character()
  m <- character()
  mc <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      mc <<- c(mc, class(x)[1])
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  list(v = v, w = w, m = m, mc = mc, el = proc.time()[["elapsed"]] - t0)
}
q <- function(e) suppressWarnings(suppressMessages(e))
first <- !file.exists(out)
for (seed in seeds) {
  d <- if (famname == "k3") mk3(seed) else if (re == "none") {
    mk(seed, n = nobs, mix = TRUE)
  } else mkg(seed)
  for (arm in c("A0", "D")) {
    ctl <- if (arm == "A0") frmtmb_control(starts = 0) else frmtmb_control()
    r <- grab(frm(form, family = fam(), data = d, control = ctl))
    row <- data.frame(family = famname, re = re, n = nobs, G = G, seed = seed,
                      arm = arm, elapsed = r$el)
    f <- r$v
    if (inherits(f, "error")) {
      row$err <- clean(conditionMessage(f))
    } else {
      resp <- names(f$spec$responses)[1]
      ms <- f$opt$mixture_starts
      cl <- tryCatch(ns$mixture_ord_start_class(f, resp),
                     error = function(e) NULL)
      sm <- grab(summary(f))
      row <- cbind(row, data.frame(
        err = "", ll = as.numeric(logLik(f)), code = f$opt$convergence,
        evals = f$opt$evals %||% NA,
        swaps = if (is.null(ms)) 0 else ms[["swaps"]],
        swapped = if (is.null(ms)) 0 else ms[["swapped"]],
        cls = if (is.null(cl)) "error" else if (length(cl$step)) "step" else
          if (length(cl$boundary)) "flat" else "proper",
        nwarn = length(r$w), s_nwarn = length(sm$w),
        degW = any(grepl("ended at a degenerate boundary", r$w)),
        lostW = any(grepl("Standard errors are not available",
                          c(r$w, sm$w))),
        swapMsg = any(r$mc == "frmtmb_mixture_swap"),
        stepMsg = any(r$mc == "frmtmb_mixture_step_start"),
        warns = clean(r$w), s_warns = clean(sm$w)))
    }
    cols <- c("family", "re", "n", "G", "seed", "arm", "elapsed", "err",
              "ll", "code", "evals", "swaps", "swapped", "cls", "nwarn",
              "s_nwarn", "degW", "lostW", "swapMsg", "stepMsg", "warns",
              "s_warns")
    for (cn in setdiff(cols, names(row))) row[[cn]] <- NA
    utils::write.table(row[cols], out, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = first, append = !first)
    first <- FALSE
  }
  cat(famname, re, nobs, "seed", seed, "done\n")
}
cat("DONE after\n")
