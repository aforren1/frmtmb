# Lane ordstart, item 1: the 40-fit table of dev/optima-ordmix-starts.R,
# ported to the round's libraries and recorded start by start. Per seed
# (two-class latent data of dev/optima-csmix.R, n = 300): the default fit
# (run 0), then 10 starts jittered about the default's start (every outer
# coordinate by N(0, 0.5), seed 1000 + s, the generator of the original
# script). Each run keeps its logLik, optimizer code, the warnings it
# raised and mixture_ord_degeneracy()'s row per component, so the
# better points can be classified afterwards.
#   Rscript dev/ordstart-starts.R base|lane seeds out.tsv [family] [nstart]
#     [n] [starts]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seeds <- eval(parse(text = args[2]))
out <- args[3]
famname <- if (length(args) > 3) args[4] else "cum_cum"
nstart <- if (length(args) > 4) as.integer(args[5]) else 10L
nobs <- if (length(args) > 5) as.integer(args[6]) else 300L
# a count of jittered starts for the default fit (lane arm only)
starts_n <- if (length(args) > 6) as.integer(args[7])
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r20",
  lane = c("C:/Users/adf44/source/r/wt-ordstart-lib",
           "C:/Users/adf44/source/r/rellib-r20"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
fam <- function() {
  switch(famname,
    cum_cum = mixture(cumulative(), cumulative()),
    cum_sratio = mixture(cumulative(), sratio()))
}
fit1 <- function(d, start = NULL, control = NULL) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  args <- list(bf(y ~ x + z), family = fam(), data = d, start = start)
  if (!is.null(control)) args$control <- control
  f <- tryCatch(withCallingHandlers(do.call(frm, args),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    row <- data.frame(ll = NA, code = NA, degen_warn = NA, step_warn = NA,
                      starts_gain = NA, conv_warn = NA, maxabs = NA,
                      elapsed = el, err = conditionMessage(f))
    for (k in 1:2) {
      for (v in c("sharpen", "near", "collapsed", "weight", "reach",
                  "cover")) {
        row[[paste0(v, k)]] <- NA
      }
    }
    row$warn <- ""
    return(list(row = row, fit = NULL))
  }
  dg <- tryCatch(frmtmb:::mixture_ord_degeneracy(
    f, names(f$spec$responses)[1]), error = function(e) NULL)
  row <- data.frame(ll = as.numeric(logLik(f)), code = f$opt$convergence,
    degen_warn = any(grepl("ended at a degenerate boundary", w)),
    step_warn = any(grepl("jittered starts reached a higher", w)),
    starts_gain = if (is.null(f$opt$mixture_starts)) NA else
      f$opt$mixture_starts[["gain"]],
    conv_warn = any(grepl("converge|gradient", w)),
    maxabs = max(abs(f$opt$par)), elapsed = el, err = "")
  for (k in 1:2) {
    for (v in c("sharpen", "near", "collapsed", "weight", "reach",
                "cover")) {
      row[[paste0(v, k)]] <- if (is.null(dg)) NA else dg[[v]][k]
    }
  }
  row$warn <- paste(substr(gsub("[\t\n]", " ", w), 1, 60), collapse = " | ")
  list(row = row, fit = f)
}
rows <- list()
for (s in seeds) {
  d <- mk(s, n = nobs, mix = TRUE)
  r0 <- if (is.null(starts_n)) fit1(d) else
    fit1(d, control = frmtmb_control(starts = starts_n))
  rows[[length(rows) + 1L]] <- cbind(seed = s, run = 0L, r0$row)
  tpl <- frm(bf(y ~ x + z), family = fam(), data = d,
             dry_run = "objective")$estimates
  set.seed(1000 + s)
  for (i in seq_len(nstart)) {
    st <- lapply(tpl, function(v) v + stats::rnorm(length(v), 0, 0.5))
    # disc is held at 1 (mapped), so its slot keeps its value
    dk <- grepl("^disc", names(tpl$betad))
    st$betad[dk] <- tpl$betad[dk]
    r <- fit1(d, st)
    rows[[length(rows) + 1L]] <- cbind(seed = s, run = i, r$row)
  }
  cat("seed", s, "done\n")
}
X <- do.call(rbind, rows)
X$family <- famname
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(X), "rows to", out, "\n")
