# Reviewer of lane optstart, round 2: do the start rule and
# newton_polish()'s gate decide the same side of 1e5 on the design
# bisected to the bound (dev/optstart-rev-kappa.R), under this R's BLAS?
#   Rscript dev/optstart-rev2-kappa.R <lib> <label> [deltas]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), " R", R.home(), "\n")
ns <- asNamespace("frmtmb")
set.seed(42)
n <- 300
base <- data.frame(x1 = rnorm(n), e = rnorm(n), g = factor(rep(1:15, 20)))
base$y <- 2 + base$x1 + rnorm(15, 0, 1)[base$g] + rnorm(n)
mk <- function(delta) {
  d <- base
  d$x2 <- d$x1 + delta * d$e
  d
}
fo <- y ~ x1 + x2 + (1 | g)
side <- function(delta) {
  d <- mk(delta)
  fr <- frm(fo, data = d, dry_run = "frame")
  k_rule <- ns$lp_design_kappa(fr$linpreds[[1]])
  k_pol <- ns$predictor_design_kappa(fr)
  st <- par_template(fo, data = d, family = gaussian())
  applies <- !identical(unname(st$betad), log(stats::sd(d$y)))
  c(k_rule = k_rule, k_pol = k_pol, applies = applies,
    polish_ill = k_pol > 1e5)
}
kap <- function(delta) side(delta)[["k_rule"]]
lo <- 1e-6; hi <- 1e-3
for (i in 1:200) {
  mid <- (lo + hi) / 2
  if (kap(mid) > 1e5) lo <- mid else hi <- mid
  if (hi - lo <= 2 * .Machine$double.eps * hi) break
}
ds <- c(lo, hi)
if (length(a) > 2) ds <- c(ds, as.numeric(strsplit(a[3], ",")[[1]]))
for (dl in ds) {
  v <- side(dl)
  cat(sprintf(paste("delta %.17g: rule kappa %.17g, gate kappa %.17g,",
                    "rule applies %s, gate ill %s, consistent %s\n"),
              dl, v[["k_rule"]], v[["k_pol"]], as.logical(v[["applies"]]),
              as.logical(v[["polish_ill"]]),
              as.logical(v[["applies"]]) != as.logical(v[["polish_ill"]])))
}
