# Review of lane quadseam2, attack 4: the lane runs the rule's
# stationary escape only when quad_fit()'s Laplace-side escape did not
# (item 5408) and removes stationary_escape_keep(). Is there a
# quadrature fit where the Laplace escape ran AND the rule's escape on
# base then gained more than rounding? escape_stationary() is wrapped to
# log each call with starts (rule calls carry a tally), as
# dev/quadseam2-escdup.R does. Skew normals with (1 | g) started on
# alpha = 0, small groups, a large group sd, weak true skew.
#   Rscript quadseam2-review-esc2.R <base|lib>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
calls <- new.env()
orig <- get("escape_stationary", ns)
unlockBinding("escape_stationary", ns)
assign("escape_stationary", function(obj, opt, frame, control, bounds,
                                     par_units = NULL, verbose = 0L,
                                     tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  n_st <- length(ns$stationary_escapes(obj, opt, frame))
  o <- orig(obj, opt, frame, control, bounds, par_units, verbose, tally)
  if (n_st) {
    calls$log <- rbind(calls$log, data.frame(
      where = where, gain = o$stationary_escape[["gain"]]))
  }
  o
}, envir = ns)
sn_draw <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * (abs(rnorm(n)) - sqrt(2 / pi)) + sqrt(1 - dl^2) * rnorm(n)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
run <- function(lbl, s, fo, d, start) {
  calls$log <- NULL
  f <- tryCatch(suppressWarnings(suppressMessages(
    frm(fo, family = skew_normal(), data = d, quadrature = TRUE,
        start = start))), error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%s %2d ERROR %s\n", lbl, s, substr(conditionMessage(f), 1, 80)))
    return(invisible())
  }
  lg <- calls$log %||% data.frame(where = character(), gain = numeric())
  la <- lg$gain[lg$where == "laplace"]
  ru <- lg$gain[lg$where == "rule"]
  al <- f$estimates$betad[grep("^alpha", names(f$estimates$betad))]
  cat(sprintf("ROW %s %2d laplace %s rule %s ll %.10f alpha %s\n", lbl, s,
              if (length(la)) format(la[1], digits = 4) else "-",
              if (length(ru)) format(ru[1], digits = 4) else "-",
              as.numeric(logLik(f)),
              paste(format(al, digits = 4), collapse = ",")))
}
for (s in 1:15) {
  set.seed(1000 + s)
  # symmetric, 15 groups of 4
  d <- data.frame(x = rnorm(60), g = factor(rep(1:15, each = 4)))
  d$y <- 0.5 * d$x + rnorm(15, 0, 0.8)[d$g] + rnorm(60)
  run("SYMSMALL", s, bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ 1), d,
      list(betad = c(0, 0)))
  # weak skew, large group sd, 10 groups of 5
  d2 <- data.frame(x = rnorm(50), g = factor(rep(1:10, each = 5)))
  d2$y <- 0.5 * d2$x + rnorm(10, 0, 1.5)[d2$g] + sn_draw(50, 1.5)
  run("WEAK", s, bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ 1), d2,
      list(betad = c(0, 0)))
  # alpha ~ x, both coefficients at 0
  d3 <- data.frame(x = rnorm(80), g = factor(rep(1:16, each = 5)))
  d3$y <- 0.5 * d3$x + rnorm(16, 0, 1)[d3$g] + sn_draw(80, 2)
  run("ALPHAX", s, bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x), d3,
      list(betad = c(0, 0, 0)))
}
cat("END\n")
