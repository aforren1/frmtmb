# Short check of lane ordstart: every step-to-proper move of the fix
# build that this check did not profile (sweep: acat_cumcll_K7 seeds 4
# and 6, theta_formula seed 4; 40-fit table: cum_cum 14 and 15,
# cum_sratio 4, 9, 14, 15, 18). At each returned point: the largest
# |Est.Error| of fixef(), the smallest Hessian eigenvalue, the smallest
# category cover of mixture_ord_degeneracy(), the warnings other than
# the step-start one, and the gain of a refit from the point.
#   Rscript ordstart-short-moves.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/reviews/ordstart-review-sweep.R")
eval(parse(text = src[grep("^gen <- ", src):(grep("^cases <- list", src) - 1L)]))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
cc <- quote(mixture(cumulative(), cumulative()))
cs <- quote(mixture(cumulative(), sratio()))
cases <- list(
  acat_cumcll_K7_seed4 = list(d = quote(gen(4, K = 7)), f = bf(y ~ x + z),
                              fam = quote(mixture(acat(), cumulative("cloglog")))),
  acat_cumcll_K7_seed6 = list(d = quote(gen(6, K = 7)), f = bf(y ~ x + z),
                              fam = quote(mixture(acat(), cumulative("cloglog")))),
  theta_formula_seed4 = list(d = quote(gen(4, K = 4)),
                             f = bf(y ~ x + z, theta2 ~ z), fam = cc),
  cum_cum_14 = list(d = quote(mk(14, mix = TRUE)), f = bf(y ~ x + z), fam = cc),
  cum_cum_15 = list(d = quote(mk(15, mix = TRUE)), f = bf(y ~ x + z), fam = cc),
  cum_sratio_4 = list(d = quote(mk(4, mix = TRUE)), f = bf(y ~ x + z), fam = cs),
  cum_sratio_9 = list(d = quote(mk(9, mix = TRUE)), f = bf(y ~ x + z), fam = cs),
  cum_sratio_14 = list(d = quote(mk(14, mix = TRUE)), f = bf(y ~ x + z), fam = cs),
  cum_sratio_15 = list(d = quote(mk(15, mix = TRUE)), f = bf(y ~ x + z), fam = cs),
  cum_sratio_18 = list(d = quote(mk(18, mix = TRUE)), f = bf(y ~ x + z), fam = cs))
for (nm in names(cases)) {
  ca <- cases[[nm]]
  d <- eval(ca$d)
  fam <- eval(ca$fam)
  w <- character()
  f <- withCallingHandlers(frm(ca$f, family = fam, data = d),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  fe <- fixef(f)
  H <- f$obj$he(f$opt$par)
  ev <- min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  dg <- frmtmb:::mixture_ord_degeneracy(f, "y")
  g <- q(frm(ca$f, family = fam, data = d, start = f$estimates,
             control = frmtmb_control(starts = 0)))
  other <- w[!grepl("jittered starts reached", w)]
  cat(sprintf(paste("MOVE %-22s ll %.4f code %d max|est| %.3g max SE %.3g",
                    "min eig %.3g min cover %.3g sharpen %s refit gain %.4f",
                    "| other warnings: %s\n"),
              nm, as.numeric(logLik(f)), f$opt$convergence,
              max(abs(fe[, "Estimate"])), max(fe[, "Est.Error"]), ev,
              min(dg$cover, na.rm = TRUE),
              paste(signif(dg$sharpen, 3), collapse = ","),
              as.numeric(logLik(g)) - as.numeric(logLik(f)),
              if (length(other)) paste(substr(other, 1, 70), collapse = " | ")
              else "none"))
}
cat("DONE\n")
