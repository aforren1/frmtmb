# Reviewer, claim 7: row standata:620 reads sdata$idxl_y_x_1. Under the
# lane's helper the read is a real vector; under the base helper it is
# NULL, and all(NULL %in% 9:5) is TRUE, which the harness's hollow-read
# check does not see (dev/aterms2-rev-log-08-portcmp.txt).
# Log: dev/aterms2-rev-log-08b-row620.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
for (h in c("tests/testthat/helper-brms-suite.R",
            "dev/aterms2-rev-tt-basehelper/helper-brms-suite.R")) {
  e <- new.env()
  sys.source(file.path("C:/Users/adf44/source/r/frmtmb-wt-aterms2", h), e)
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = c(rnorm(9), NA), z = rnorm(10),
                    g1 = sample(1:5, 10, TRUE), g2 = 10:1, g3 = 1:10,
                    s = c(FALSE, rep(TRUE, 9)))
  fr <- suppressWarnings(frm(bf(y ~ mi(x, idx = g1)) +
                               bf(x | mi() + index(g2) + subset(s) ~ 1),
                             data = dat, family = gaussian(),
                             dry_run = "frame"))
  sd <- e$brms_standata_view(fr)
  v <- sd$idxl_y_x_1
  cat(h, ": idxl_y_x_1 =", if (is.null(v)) "NULL" else v,
      " all(v %in% 9:5) =", all(v %in% 9:5), "\n")
  bs <- brms::standata(brms::bf(y ~ mi(x, idx = g1)) +
                         brms::bf(x | mi() + index(g2) + subset(s) ~ 1) +
                         brms::set_rescor(FALSE), dat)
  cat("   brms idxl_y_x_1 =", bs$idxl_y_x_1, " identical:",
      identical(as.integer(v), as.integer(bs$idxl_y_x_1)), "\n")
}
