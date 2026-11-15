# Review of lane famend: the nested-quadrature GEV fit of
# test-families-extreme.R ("importance and nested quadrature fit a GEV off
# its support"), which fails under OpenBLAS 0.3.26 LAPACK on the lane.
#   Rscript dev/reviews/famend-review-gevquad.R <base|lane> [seed]
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "| R", R.home(), "| LAPACK",
    La_library(), "\n")
rg <- function(n, mu, s, k) mu + s * ((-log(stats::runif(n)))^(-k) - 1) / k
sim <- function(seed, xi) {
  set.seed(seed)
  g <- rep(seq_len(15), each = 8)
  u <- stats::rnorm(15, 0, 0.7)
  d <- data.frame(g = factor(g), y = rg(120, 1 + u[g], 1, xi))
  d$h <- factor(paste(d$g, rep(1:2, length.out = 120)))
  d
}
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 3
# a third argument forces the cap of exp_neg_capped() for every family in
# this process, to isolate the cap from the rest of the build
if (length(a) > 2) {
  capv <- as.numeric(a[3])
  ns <- asNamespace("frmtmb")
  f <- get("exp_neg_capped", ns)
  body(f) <- substitute({
    if (!inherits(w, "advector") && exact_support()) return(exp(-w))
    exp(-ad_max(w, -CAP))
  }, list(CAP = capv))
  assignInNamespace("exp_neg_capped", f, ns)
  cat("exp_neg_capped forced to cap", capv, "\n")
}
for (s in seeds) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  fq <- tryCatch(withCallingHandlers(
    frm(y ~ 1 + (1 | g) + (1 | h), family = gen_extreme_value(),
        data = sim(s, 0.5), quadrature = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") }),
    error = function(e) structure(conditionMessage(e), class = "err"))
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(fq, "err")) {
    cat(sprintf("seed %d ERROR %s\n", s, substr(fq, 1, 200)))
    next
  }
  cat(sprintf("seed %d: %.1f s code %d objective %.10g logLik %.10g\n", s, el,
              fq$opt$convergence, fq$opt$objective, as.numeric(logLik(fq))))
  cat("   par", paste(sprintf("%.8g", fq$opt$par), collapse = " "), "\n")
  for (x in unique(w)) cat("   W:", substr(gsub("\\s+", " ", x), 1, 300), "\n")
}
cat("done\n")
