# Reviewer short check of lane ordstarts (build fix round 1), items 3 and
# 5: does the inner start (mixture_swap_inner(), swapped random-effect
# modes) change any returned fit beyond speed? Each fit twice in one
# process: as built, and with mixture_swap_inner() a no-op. Also records
# whether the settled point is ever below the starts = 0 fit.
#   Rscript ordstarts-short-p8.R > ordstarts-short-p8.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
orig <- ns$mixture_swap_inner
gen <- function(seed, n, G) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  g <- rep_len(seq_len(G), n)
  u0 <- rnorm(G, 0, 0.5)
  lat <- lat + u0[g]
  y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  data.frame(y = y, x = x, z = z, g = factor(g))
}
q <- function(e) suppressWarnings(suppressMessages(e))
rec <- function(f) {
  ms <- f$opt$mixture_starts
  if (is.null(ms)) "NULL" else
    paste0("swaps=", ms[["swaps"]], " swapped=", ms[["swapped"]])
}
fit <- function(d, fam, ctl = frmtmb_control()) {
  t0 <- proc.time()[["elapsed"]]
  f <- q(frm(bf(y ~ x + z + (1 | g)), family = fam, data = d, control = ctl))
  list(f = f, el = proc.time()[["elapsed"]] - t0)
}
for (des in list(list("cs", 1000, 50), list("cc", 300, 20),
                 list("cc", 1000, 50))) {
  fam <- if (des[[1]] == "cs") mixture(cumulative(), sratio()) else
    mixture(cumulative(), cumulative())
  for (seed in 1:5) {
    d <- gen(seed, des[[2]], des[[3]])
    a0 <- fit(d, fam, frmtmb_control(starts = 0))
    utils::assignInNamespace("mixture_swap_inner", orig, "frmtmb")
    with_inner <- fit(d, fam)
    utils::assignInNamespace("mixture_swap_inner",
                             function(obj, b) invisible(NULL), "frmtmb")
    without <- fit(d, fam)
    utils::assignInNamespace("mixture_swap_inner", orig, "frmtmb")
    l0 <- as.numeric(logLik(a0$f))
    l1 <- as.numeric(logLik(with_inner$f))
    l2 <- as.numeric(logLik(without$f))
    cat(sprintf(paste0("%s n%d seed %d | A0 %.6f | D %.6f (%s, %.1f s) |",
                       " D without inner start %.6f (%s, %.1f s) | dll %.2g",
                       " max|dpar| %.2g | D below A0: %s\n"),
                des[[1]], des[[2]], seed, l0, l1, rec(with_inner$f),
                with_inner$el, l2, rec(without$f), without$el, l1 - l2,
                max(abs(with_inner$f$opt$par - without$f$opt$par)),
                l1 < l0 - 1e-6))
  }
}
cat("DONE p8\n")
