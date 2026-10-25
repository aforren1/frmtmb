# Reviewer, fix round 1: how brms 2.23.0 and the fixed frmtmb read each
# bhaz() argument spelling.
#   Rscript dev/review-bhaz-names.R > dev/review-bhaz-log/names.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb.survival)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.survival"),
    "brms", format(packageVersion("brms")), "\n")
set.seed(11); n <- 200; x <- rnorm(n)
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = factor(rep(c("a", "b"), 100)))
strip <- function(m) matrix(as.numeric(m), nrow(m))
cases <- c(
  "bhaz(d = 6)", "bhaz(g = g)", "bhaz(kn = c(0.5, 1.1))",
  "bhaz(k = c(0.5, 1.1))", "bhaz(kno = c(0.5, 1.1))", "bhaz(deg = 2)",
  "bhaz(degr = 2, intercept = FALSE)", "bhaz(de = 2)", "bhaz(NA)",
  "bhaz(NA, 6)", "bhaz(6, gr = NA)", "bhaz(df = 6, NA)",
  "bhaz(gr = NA, d = 6)", "bhaz(warn.outside = FALSE)",
  "bhaz(warn = FALSE)", "bhaz(w = FALSE)", "bhaz(intercept = FALSE)",
  "bhaz(int = FALSE)", "bhaz(int = TRUE)", "bhaz(inter = FALSE)",
  "bhaz(i = FALSE)", "bhaz(df = 6, i = TRUE)",
  "bhaz(Boundary.knots = c(0, 30))", "bhaz(Boundary = c(0, 30))",
  "bhaz(Bound = c(0, 30))", "bhaz(B = c(0, 30))", "bhaz(periodic = TRUE)",
  "bhaz(per = TRUE)", "bhaz(p = TRUE)", "bhaz(derivs = 1)",
  "bhaz(integral = TRUE)", "bhaz(foo = 1)", "bhaz(Knots = c(0.5, 1.1))",
  "bhaz(x = 1)", "bhaz(NA, 6, c(0.5, 1.1))", "bhaz(df = 6, d = 7)",
  "bhaz(degree = 2, deg = 3)", "bhaz(knots = c(0.5, 1.1), kn = c(1, 2))",
  "bhaz(d = 6, degree = 2)", "bhaz(D = 6)", "bhaz(df = 7, deg = 2)")
desc <- function(Z, K) {
  a <- attributes(Z)
  sprintf("K %d deg %s int %s knots %s bnd %s%s", K, a$degree, a$intercept,
          paste(format(a$knots, digits = 4), collapse = ","),
          paste(format(a$Boundary.knots, digits = 4), collapse = ","),
          paste0("", if (isTRUE(a$periodic)) " PERIODIC",
                 if (isTRUE(a$integral)) " INTEGRAL",
                 if (!is.null(a$derivs) && a$derivs > 0) " DERIVS"))
}
for (cs in cases) {
  fo <- stats::as.formula(paste("time | cens(cens) +", cs, "~ x"))
  sd <- tryCatch(suppressWarnings(standata(bf(fo), data = dd,
                                           family = brms::cox())),
                 error = function(e) conditionMessage(e))
  fr <- tryCatch(frm(frmtmb::bf(fo), family = frmtmb.survival::cox(),
                     data = dd, dry_run = "frame"),
                 error = function(e) conditionMessage(e))
  b <- if (is.character(sd)) paste("ERR", substr(sd, 1, 70)) else {
    paste0(desc(sd$Zbhaz, sd$Kbhaz),
           if (!is.null(sd$ngrbhaz)) paste(" ngr", sd$ngrbhaz))
  }
  f <- if (is.character(fr)) paste("REFUSED", substr(fr, 1, 150)) else {
    sp <- fr$aterm_values$time$bhaz_spec
    sprintf("K %d deg %d int %s knots %s bnd %s",
            ncol(fr$aterm_values$time$Zbhaz), sp$degree, sp$intercept,
            paste(format(sp$internal, digits = 4), collapse = ","),
            paste(format(sp$boundary, digits = 4), collapse = ","))
  }
  dz <- if (!is.character(sd) && !is.character(fr) &&
            ncol(fr$aterm_values$time$Zbhaz) == sd$Kbhaz) {
    sprintf("%.1e", max(abs(fr$aterm_values$time$Zbhaz -
                              strip(sd$Zbhaz))) / max(abs(sd$Zbhaz)))
  } else "-"
  cat(sprintf("%-42s\n   brms: %s\n   frm : %s\n   rel max|dZ| %s\n", cs, b,
              f, dz))
}
