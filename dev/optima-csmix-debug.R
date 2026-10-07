# Lane optima, item 3: where does the NaN gradient of a cs() ordinal
# mixture come from? One seed of dev/optima-csmix.R's mix_cum_sratio
# (seed 3 died with "NA/NaN gradient evaluation" on 0.68.1); the
# objective is taped unfitted, the optimizer run by hand, and the
# components' per-row densities read at the point where the gradient
# first stops being finite.
#   Rscript dev/optima-csmix-debug.R base|lane [seed] [design]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seed <- if (length(args) > 1) as.integer(args[2]) else 3L
design <- if (length(args) > 2) args[3] else "cum_sratio"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
d <- mk(seed, mix = TRUE)
fam <- switch(design,
  cum_sratio = mixture(cumulative(), sratio()),
  cum_cum = mixture(cumulative(), cumulative()))
u <- suppressWarnings(frm(bf(y ~ x + cs(z)), family = fam, data = d,
                          dry_run = "objective"))
obj <- u$obj
cat("outer parameters:", length(obj$par), "\n")
print(table(names(obj$par)))
trace <- list()
fn <- function(p) {
  v <- obj$fn(p)
  trace[[length(trace) + 1L]] <<- list(p = p, f = v, kind = "fn")
  if (is.nan(v)) Inf else v
}
gr <- function(p) {
  g <- obj$gr(p)
  trace[[length(trace) + 1L]] <<- list(p = p, g = g, kind = "gr",
                                       f = obj$fn(p))
  g
}
o <- tryCatch(stats::nlminb(obj$par, fn, gr,
                            control = list(eval.max = 1000,
                                           iter.max = 1000)),
              error = function(e) conditionMessage(e))
print(if (is.character(o)) o else o[c("objective", "convergence",
                                       "message")])
# frm()'s restart: one more run from where the first stopped
if (!is.character(o)) {
  cat("objective at the returned par:", format(obj$fn(o$par), digits = 10),
      "; reported:", format(o$objective, digits = 10), "\n")
  pb <- obj$env$last.par.best
  cat("objective at last.par.best:", format(obj$fn(pb), digits = 10), "\n")
  o2 <- tryCatch(stats::nlminb(o$par, fn, gr,
                               control = list(eval.max = 1000,
                                              iter.max = 1000)),
                 error = function(e) conditionMessage(e))
  print(if (is.character(o2)) o2 else o2[c("objective", "convergence",
                                           "message")])
}
kinds <- vapply(trace, `[[`, "", "kind")
gfin <- vapply(trace, function(t) {
  if (t$kind == "gr") all(is.finite(t$g)) else NA
}, NA)
ffin <- vapply(trace, function(t) is.finite(t$f), NA)
cat("evaluations: fn", sum(kinds == "fn"), " gr", sum(kinds == "gr"),
    "; non-finite fn values", sum(!ffin), "; non-finite gradients",
    sum(!gfin, na.rm = TRUE), "\n")
bad <- which(kinds == "gr" & !gfin)
if (length(bad)) {
  t <- trace[[bad[1]]]
  cat("first non-finite gradient at evaluation", bad[1], ": objective",
      format(t$f, digits = 10), "\n")
  print(stats::setNames(t$g, names(obj$par))[!is.finite(t$g)])
  # the cumulative component's per-row thresholds, tau_k - cs_k z - eta
  pl0 <- obj$env$parList(t$p)
  tr <- pl0$tau_raw1
  tau <- cumsum(c(tr[1], exp(tr[-1])))
  eta <- pl0$beta[1] * d$x
  A <- outer(-eta, tau, "+") - outer(d$z, pl0$bcs3)
  lo <- cbind(-Inf, A)
  up <- cbind(A, Inf)
  gap <- up[cbind(seq_len(nrow(d)), d$y)] - lo[cbind(seq_len(nrow(d)), d$y)]
  cat("component 1 (cumulative): rows whose own category's thresholds",
      "cross:", sum(gap < 0), "; touch (gap < 1e-12):",
      sum(gap >= 0 & gap < 1e-12), "; min gap", format(min(gap), digits = 4),
      "\n")
  cross_any <- rowSums(t(apply(A, 1, diff)) < 0)
  cat("component 1: rows with any crossing pair (any category):",
      sum(cross_any > 0), "\n")
  pl <- obj$env$parList(t$p)
  print(pl[setdiff(names(pl), "beta")])
  print(pl$beta)
  # the cs() thresholds per row of component 1 and 2
  saveRDS(list(p = t$p, pl = pl), sprintf("dev/optima-log/csmix-bad-%s-%d.rds",
                                          design, seed))
}
