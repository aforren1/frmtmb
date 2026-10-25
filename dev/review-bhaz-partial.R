# Reviewer: partial argument names in bhaz(), brms 2.23.0 against frmtmb.
#   Rscript dev/review-bhaz-partial.R > dev/review-bhaz-log/partial.txt 2>&1
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
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
strip <- function(m) matrix(as.numeric(m), nrow(m))
cases <- c("bhaz(Boundary = c(0, 30))", "bhaz(Bound = c(0, 30))",
           "bhaz(inter = FALSE)", "bhaz(int = FALSE)", "bhaz(kn = c(0.5, 1.1))",
           "bhaz(deg = 2)", "bhaz(d = 6)", "bhaz(warn.outside = FALSE)",
           "bhaz(6)", "bhaz(Knots = c(0.5, 1.1))")
for (cs in cases) {
  fo <- stats::as.formula(paste("time | cens(cens) +", cs, "~ x"))
  sd <- tryCatch(suppressWarnings(standata(bf(fo), data = dd,
                                           family = brms::cox())),
                 error = function(e) conditionMessage(e))
  fr <- tryCatch(frm(frmtmb::bf(fo), family = frmtmb.survival::cox(),
                     data = dd, dry_run = "frame"),
                 error = function(e) conditionMessage(e))
  b <- if (is.character(sd)) paste("ERR", sd) else
    sprintf("K %d bnd %s int %s", sd$Kbhaz,
            paste(format(attr(sd$Zbhaz, "Boundary.knots"), digits = 5),
                  collapse = ","), attr(sd$Zbhaz, "intercept"))
  f <- if (is.character(fr)) paste("ERR", substr(fr, 1, 110)) else {
    sp <- fr$aterm_values$time$bhaz_spec
    sprintf("K %d bnd %s int %s", ncol(fr$aterm_values$time$Zbhaz),
            paste(format(sp$boundary, digits = 5), collapse = ","),
            sp$intercept)
  }
  same <- if (!is.character(sd) && !is.character(fr) &&
              ncol(fr$aterm_values$time$Zbhaz) == sd$Kbhaz) {
    sprintf("%.2e", max(abs(fr$aterm_values$time$Zbhaz - strip(sd$Zbhaz))))
  } else "-"
  cat(sprintf("%-28s brms: %s\n%-28s frm : %s | max|dZ| %s\n", cs, b, "",
              f, same))
}
