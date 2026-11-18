# Lane ordstarts (backlog 6136): how often 10 jittered starts move an
# ordinal mixture with random effects, and what they cost. Per seed, in
# one process and interleaved, the default fit (no starts with random
# effects) and frmtmb_control(starts = 10). Data: the two-class latent
# data of dev/optima-csmix.R (mk(), mix = TRUE) plus a group effect on
# the latent scale, drawn after the rest so the fixed part is the same
# data as dev/silent2-recost.R: u0 ~ N(0, 0.5) per group, and for the
# slope design also u1 ~ N(0, 0.3) per group times x. Groups are
# rep_len(1:G, n), as in dev/silent2-recost.R.
# Per fit: logLik, code, every warning and message with its class,
# elapsed seconds, optimizer evaluations, the starts record, and (after
# the clock stops) the degeneracy rows, the threshold standard errors
# the fit-end check saw and those of the forced covariance, the
# smallest eigenvalue of the outer Hessian and the largest gradient.
#   Rscript dev/ordstarts-run.R <family> <re> <n> <G> <seeds> <out.tsv>
#     family: cum_cum | cum_sratio ; re: int | slope | slopeu (no cor)
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
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
mkg <- function(seed) {
  d <- mk(seed, n = nobs, mix = TRUE)
  # the latent draw of mk() again, so the group effect enters before the
  # cut points; same RNG order as mk()
  set.seed(seed)
  x <- rnorm(nobs)
  z <- rnorm(nobs)
  cls <- rbinom(nobs, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(nobs)
  stopifnot(identical(x, d$x), identical(z, d$z))
  g <- rep_len(seq_len(G), nobs)
  u0 <- rnorm(G, 0, 0.5)
  lat <- lat + u0[g]
  if (re %in% c("slope", "slopeu")) {
    u1 <- rnorm(G, 0, 0.3)
    lat <- lat + u1[g] * x
  }
  d$y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  d$g <- factor(g)
  d
}
form <- switch(re,
  int = bf(y ~ x + z + (1 | g)),
  slope = bf(y ~ x + z + (1 + x | g)),
  slopeu = bf(y ~ x + z + (1 + x || g)))
fam <- function() {
  switch(famname,
    cum_cum = mixture(cumulative(), cumulative()),
    cum_sratio = mixture(cumulative(), sratio()))
}
clean <- function(s) substr(gsub("[\t\n\r]", " ", s), 1, 160)
fit1 <- function(d, ctl) {
  w <- character()
  wc <- character()
  m <- character()
  mc <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(form, family = fam(), data = d, control = ctl),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      wc <<- c(wc, class(x)[1])
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      mc <<- c(mc, class(x)[1])
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  list(f = f, el = el, w = w, wc = wc, m = m, mc = mc)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
row_of <- function(r) {
  f <- r$f
  base <- data.frame(elapsed = r$el,
    nwarn = length(r$w), nmsg = length(r$m),
    degen_warn = any(grepl("ended at a degenerate boundary", r$w)),
    re_hint = any(grepl("runs no jittered starts", r$w)),
    step_msg = any(r$mc == "frmtmb_mixture_step_start") ||
      any(grepl("jittered starts reached a higher", c(r$w, r$m))),
    conv_warn = any(grepl("converge|gradient", r$w, ignore.case = TRUE)),
    warn_classes = paste(unique(r$wc), collapse = ","),
    msg_classes = paste(unique(r$mc), collapse = ","),
    warns = paste(clean(r$w), collapse = " | "),
    msgs = paste(clean(r$m), collapse = " | "))
  if (inherits(f, "error")) {
    base$err <- clean(conditionMessage(f))
    return(base)
  }
  base$err <- ""
  resp <- names(f$spec$responses)[1]
  base$ll <- as.numeric(logLik(f))
  base$code <- f$opt$convergence
  base$evals <- f$opt$evals %||% NA
  base$np <- length(f$opt$par)
  base$maxabs <- max(abs(f$opt$par))
  ms <- f$opt$mixture_starts
  for (v in c("starts", "reached", "gain", "degenerate_above", "tier")) {
    base[[paste0("ms_", v)]] <- if (is.null(ms)) NA else ms[[v]]
  }
  base$maxgrad <- tryCatch(max(abs(quiet(f$obj$gr(f$opt$par)))),
                           error = function(e) NA)
  cl <- tryCatch(ns$mixture_ord_start_class(f, resp), error = function(e) e)
  base$cls <- if (inherits(cl, "error")) "error" else
    if (length(cl$step)) "step" else if (length(cl$boundary)) "flat" else
      "proper"
  dg <- tryCatch(quiet(ns$mixture_ord_degeneracy(f, resp)),
                 error = function(e) NULL)
  for (k in 1:2) {
    for (v in c("sharpen", "near", "cover", "collapsed", "weight")) {
      base[[paste0(v, k)]] <- if (is.null(dg)) NA else dg[[v]][k]
    }
  }
  base$se_deferred <- is.null(f$cache$se_fit_an) && is.null(f$cache$sdr)
  se0 <- tryCatch(quiet(ns$mixture_ord_thres_se(f, resp)),
                  error = function(e) c(NA, NA))
  base$tse_fit1 <- se0[1]
  base$tse_fit2 <- se0[2]
  an <- tryCatch(quiet(ns$se_fit_analysis(f, force = TRUE)),
                 error = function(e) NULL)
  se1 <- tryCatch(quiet(ns$mixture_ord_thres_se(f, resp)),
                  error = function(e) c(NA, NA))
  base$tse1 <- se1[1]
  base$tse2 <- se1[2]
  V <- an$V
  base$maxse <- if (is.matrix(V)) {
    s <- sqrt(pmax(diag(V), 0))
    if (any(is.finite(s))) max(s[is.finite(s)]) else NA
  } else NA
  base$nlost <- length(an$lost)
  H <- f$cache$hessian_fixed$H
  base$mineig <- if (is.matrix(H) && all(is.finite(H))) {
    min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  } else NA
  base
}
first <- !file.exists(out)
for (s in seeds) {
  d <- mkg(s)
  for (arm in c("default", "starts10")) {
    ctl <- if (arm == "default") frmtmb_control() else
      frmtmb_control(starts = 10)
    r <- fit1(d, ctl)
    row <- cbind(family = famname, re = re, n = nobs, G = G, seed = s,
                 arm = arm, row_of(r))
    cols <- c("family", "re", "n", "G", "seed", "arm", "elapsed", "ll",
      "code", "evals", "np", "maxabs", "ms_starts", "ms_reached", "ms_gain",
      "ms_degenerate_above", "ms_tier", "maxgrad", "cls", "degen_warn",
      "re_hint", "step_msg", "conv_warn", "nwarn", "nmsg", "se_deferred",
      "tse_fit1", "tse_fit2", "tse1", "tse2", "maxse", "nlost", "mineig",
      paste0(rep(c("sharpen", "near", "cover", "collapsed", "weight"), 2),
             rep(1:2, each = 5)),
      "warn_classes", "msg_classes", "err", "warns", "msgs")
    for (cn in setdiff(cols, names(row))) row[[cn]] <- NA
    row <- row[cols]
    utils::write.table(row, out, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = first, append = !first)
    first <- FALSE
    cat(sprintf("%s %s n%d seed %d %-8s ll %.4f code %s el %.1fs cls %s\n",
                famname, re, nobs, s, arm, row$ll, row$code, row$elapsed,
                row$cls))
  }
}
cat("DONE", famname, re, nobs, "\n")
