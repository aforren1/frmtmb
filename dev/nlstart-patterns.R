# Per pattern: where does each sign pattern of nl_sign_probe() end, on the
# 450 rows of dev/nlstart-scan.R, and what does it cost? The default fit
# is made with the probe off; each pattern then runs on its tape as
# nl_sign_probe() runs it. Feeds the choice of how many patterns to run.
# Usage: NLSTART_ARM=after [NLSTART_BASE=opt] Rscript dev/nlstart-patterns.R
#   <family>; NLSTART_BASE=opt starts the rest of the model at the fit
fam <- commandArgs(trailingOnly = TRUE)[1]
source("dev/nlstart-common.R")
src <- readLines("dev/nlstart-scan.R")
eval(parse(text = src[grep("^llv <-", src):(grep("^cat[(]\"ROWHEAD", src) - 1L)]))
ns <- asNamespace("frmtmb")
options(frmtmb.nl_sign_probe = FALSE)
for (seed in 1:5) for (slope in c(0.6, -0.6)) for (bn in names(bodies)) {
  b <- bodies[[bn]]
  d <- gen(seed * 100 + 1, slope)
  fm <- do.call(bf, c(list(as.formula(paste("y ~", deparse1(b$f)))),
                      lapply(b$p, function(p) as.formula(paste(p, "~ 1"))),
                      list(nl = TRUE)))
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(suppressMessages(frm(fm, family = F[[fam]],
                                               data = d)))
  t_fit <- proc.time()[["elapsed"]] - t0
  frame <- fit$frame
  obj <- fit$obj
  pn <- names(obj$par)
  at_obj <- integer(0)
  for (lp in frame[["linpreds"]]) {
    if (is.null(lp[["nl_body"]])) next
    bc <- ns$nl_body_coefs(frame, lp)
    added <- ns$nl_body_additive(frame, lp, fit$estimates, bc$coefs)
    k <- 0L
    for (par in names(bc$coefs)) for (at in bc$coefs[[par]]) {
      k <- k + 1L
      if (!added[k]) at_obj <- c(at_obj, which(pn == par)[at])
    }
  }
  pats <- ns$nl_sign_patterns(length(at_obj))
  bounds <- ns$fit_outer_box(fit)
  out <- character(0)
  for (s in pats) {
    p <- if (identical(Sys.getenv("NLSTART_BASE"), "opt")) fit$opt$par else obj$par
    p[at_obj] <- s
    t1 <- proc.time()[["elapsed"]]
    o <- tryCatch(suppressWarnings(ns$optimize_obj(obj, fit$control, bounds,
                                                   fit$par_units,
                                                   start_par = p)),
                  error = function(e) NULL)
    tp <- proc.time()[["elapsed"]] - t1
    out <- c(out, sprintf("%s:%.6f:%.3f",
                          paste(ifelse(s > 0, "+", "-"), collapse = ""),
                          if (is.null(o)) NA else -o$objective, tp))
  }
  cat(sprintf("PAT %s %d %.1f %s %.6f %.3f %s %s\n", fam, seed, slope, bn,
              -fit$opt$objective, t_fit,
              paste(ifelse(fit$opt$par[at_obj] > 0, "+", "-"), collapse = ""),
              paste(out, collapse = " ")))
  flush.console()
}
