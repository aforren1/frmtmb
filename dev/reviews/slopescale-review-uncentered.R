# Review of lane slopescale, follow-up of slopescale-review-rms.R (a): an
# uncentered random-slope column of small spread (x = c0 + s z, a pH near
# 7 or a log mass near 5). The plan measures the column by its root mean
# square (near c0), so the default does not engage; 0.82.0 engaged on
# y ~ x + (1 + x | g) by sd(x). Seeds 1..10, c0 in {5, 7}, s in {0.03,
# 0.01}, both forms; base default, lane default, lane autoscale = TRUE.
# "short" is against the best of the three. One process.
# Usage: Rscript slopescale-review-uncentered.R
LANE <- "C:/Users/adf44/source/r/wt-slopescale-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(BASE, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
use <- function(lib) {
  if ("frmtmb" %in% loadedNamespaces()) unloadNamespace("frmtmb")
  suppressPackageStartupMessages(library(frmtmb, lib.loc = lib))
}
mk <- function(seed, s, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  data.frame(g = g, x = c0 + s * z, y = eta + rnorm(300))
}
cases <- list()
for (c0 in c(5, 7)) for (s in c(0.03, 0.01)) for (seed in 1:10)
  for (fo in c("y ~ x + (1 + x | g)", "y ~ 1 + (1 + x | g)")) {
    cases[[length(cases) + 1L]] <- list(
      nm = sprintf("c%g s%g %d %s", c0, s, seed, fo),
      fo = stats::as.formula(fo), d = mk(seed, s, c0))
  }
run1 <- function(cs, ctl) {
  nw <- 0L
  f <- withCallingHandlers(
    tryCatch(frm(cs$fo, data = cs$d, control = ctl), error = function(e) e),
    warning = function(w) {
      nw <<- nw + 1L
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(f, "error")) return(c(ll = NA, nw = nw, as = NA))
  c(ll = as.numeric(logLik(f)), nw = nw, as = isTRUE(f$autoscaled))
}
use(BASE)
rb <- lapply(cases, run1, ctl = frmtmb_control())
use(LANE)
rl <- lapply(cases, run1, ctl = frmtmb_control())
rt <- lapply(cases, run1, ctl = frmtmb_control(autoscale = TRUE))
tab <- do.call(rbind, lapply(seq_along(cases), function(i) {
  best <- max(rb[[i]]["ll"], rl[[i]]["ll"], rt[[i]]["ll"], na.rm = TRUE)
  data.frame(case = cases[[i]]$nm,
             short_base = signif(best - rb[[i]][["ll"]], 3),
             short_lane = signif(best - rl[[i]][["ll"]], 3),
             short_T = signif(best - rt[[i]][["ll"]], 3),
             as_base = as.logical(rb[[i]][["as"]]),
             nw_base = rb[[i]][["nw"]], nw_lane = rl[[i]][["nw"]],
             nw_T = rt[[i]][["nw"]])
}))
options(width = 200)
print(tab, row.names = FALSE)
k <- 1e-3
cat("\nfits:", nrow(tab), "; short >", k, ": base", sum(tab$short_base > k),
    "(silent", sum(tab$short_base > k & tab$nw_base == 0), ")",
    "; lane", sum(tab$short_lane > k),
    "(silent", sum(tab$short_lane > k & tab$nw_lane == 0), ")",
    "; lane TRUE", sum(tab$short_T > k), "\n")
cat("lane short where base was not:", sum(tab$short_lane > k & tab$short_base <= k),
    "; base short where lane is not:", sum(tab$short_base > k & tab$short_lane <= k),
    "; max lane short", max(tab$short_lane), "; max TRUE short", max(tab$short_T),
    "\n")
