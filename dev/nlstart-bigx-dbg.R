# The review's bigx body (b1 * exp(b2 * xb), xb in 0 to 3000, cumulative):
# what start size nl_probe_scale() reads, and where starts of each size end.
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fresh.R")
eval(parse(text = src[grep("^F <- list", src):(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
ns <- asNamespace("frmtmb")
for (seed in c(22, 44, 55)) {
  d <- gen("cumulative", "bigx", seed * 1000 + match("bigx", names(B)))
  fm <- form("cumulative", "bigx")
  f0 <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                            suppressWarnings(suppressMessages(
                              frm(fm, family = cumulative(), data = d))))
  lp <- Filter(function(l) !is.null(l[["nl_body"]]), f0$frame$linpreds)[[1]]
  pr <- list(list(par = "beta", at = 1L), list(par = "beta", at = 2L))
  sc <- ns$nl_probe_scale(f0$frame, lp, f0$estimates, pr)
  cat("seed", seed, "scale mag", sc["mag", ], "key", sc["key", ], "\n")
  for (m2 in c(2^-11, 2^-10, 1e-3, 2e-3)) {
    f <- tt(frm(fm, family = cumulative(), data = d,
                start = list(beta = c("b1_(Intercept)" = -1,
                                      "b2_(Intercept)" = -m2))))
    cat(sprintf("  start (-1, -%.3g): %.4f conv %s\n", m2, llv(f),
                f$opt$convergence))
  }
  # the probe's own route: optimize on the fit's tape from that start
  obj <- f0$obj
  for (m2 in c(2^-11, 1e-3)) {
    p <- obj$par
    p[1:2] <- c(-1, -m2)
    o <- ns$optimize_obj(obj, f0$control, ns$fit_outer_box(f0), f0$par_units,
                         start_par = p)
    cat(sprintf("  tape from (-1, -%.3g): %.4f conv %d\n", m2, -o$objective,
                o$convergence))
  }
}
