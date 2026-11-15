# Review of lane predraw2: (1) the log Jacobian anchor on the fired
# ordinal-mixture fits of the coverage corpus that also have a side law;
# (2) loo::psis() on log ratios with -Inf (a draw where the objective is
# not finite).
# Usage: Rscript dev/reviews/predraw2-review-e.R (from the lane root)
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")

cat("\n== loo::psis with -Inf among the log ratios\n")
set.seed(1)
lw <- c(stats::rnorm(999), -Inf)
r <- tryCatch(loo::psis(lw, r_eff = 1), error = function(e) e)
cat("psis(c(rnorm(999), -Inf)):", if (inherits(r, "error")) {
  paste("ERROR", conditionMessage(r))
} else paste("k", r$diagnostics$pareto_k), "\n")
r <- tryCatch(ns$draw_resample(lw, character(1000)), error = function(e) e)
cat("draw_resample on it:", if (inherits(r, "error")) {
  paste("ERROR", conditionMessage(r))
} else paste("k", r$k, "picked the -Inf draw", 1000 %in% r$idx), "\n")

cat("\n== anchor on fired mixture_cumulative interior fits with a side law\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
orig <- ns$draw_logratio
one <- function(v) length(v) == 1L && v == 1
anchor <- function(law, X) {
  add <- numeric(ncol(X))
  if (is.null(law)) return(add)
  for (j in setdiff(which(!vapply(law$lo, one, TRUE) |
                          !vapply(law$hi, one, TRUE)), law$edge)) {
    neg <- X[j, ] < 0
    add[neg] <- add[neg] + log(ns$draw_edge_slope(0, law$lo[[j]]))
    add[!neg] <- add[!neg] + log(ns$draw_edge_slope(0, law$hi[[j]]))
  }
  add
}
patched <- function(object, dw, D, Z) {
  r <- orig(object, dw, D, Z)
  if (is.null(r)) return(r)
  X <- crossprod(dw$law$Lc %||% dw$L, Z)
  r$lw <- r$lw + anchor(dw$law, X)
  r
}
cuts <- c(-Inf, -1, 0.5, 2, Inf)
truthP <- function(x) {
  0.6 * t(vapply(x, function(xi) diff(stats::plogis(cuts - xi)), numeric(4))) +
    0.4 * t(vapply(x, function(xi) diff(stats::plogis(cuts - xi - 3)), numeric(4)))
}
seeds <- as.integer(strsplit(Sys.getenv("REVE_SEEDS",
  "3,4,8,9,10,11,13,15,18,25,26,28,29,31,32,35,36,39,40"), ",")[[1]])
B <- 1000L
for (seed in seeds) {
  d <- gen("mixture_cumulative", "interior", seed)
  s <- spec("mixture_cumulative")
  fit <- tryCatch(suppressMessages(suppressWarnings(
    frm(s[[1]], family = s[[2]], data = d))), error = function(e) NULL)
  if (is.null(fit)) { cat("seed", seed, "fit failed\n"); next }
  sh <- list()
  for (arm in c("lane", "absolute")) {
    assignInNamespace("draw_logratio", if (arm == "lane") orig else patched,
                      ns = "frmtmb")
    set.seed(1)
    p <- suppressWarnings(suppressMessages(
      predict(fit, ndraws = B, summary = FALSE)))
    sh[[arm]] <- t(apply(p, 2L, function(col) tabulate(col, 4L) / B))
  }
  assignInNamespace("draw_logratio", orig, ns = "frmtmb")
  law <- fit$cache[["draw_edge_law"]]
  nms <- ns$outer_par_names(fit)
  nar <- if (is.list(law)) setdiff(which(!vapply(law$lo, one, TRUE) |
                                          !vapply(law$hi, one, TRUE)), law$edge)
  s0 <- if (length(nar)) paste(sprintf("%s lo %.3f hi %.3f", nms[nar],
    vapply(nar, function(j) ns$draw_edge_slope(0, law$lo[[j]]), 0),
    vapply(nar, function(j) ns$draw_edge_slope(0, law$hi[[j]]), 0)),
    collapse = "; ") else "none"
  P <- truthP(d$x)
  cat(sprintf("seed %2d fire %s edge %s | narrowed (slope at 0): %s | mae lane %.4f absolute %.4f\n",
              seed, format(fit$cache[["draw_check"]]$fire),
              paste(nms[if (is.list(law)) law$edge], collapse = ","), s0,
              mean(abs(sh$lane - P)), mean(abs(sh$absolute - P))))
}
