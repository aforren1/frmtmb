# Round 3: a bare addition-term name that is a full-length vector outside
# the data, and a scalar constant, on both builds. Seed 3.
# FORMROBUST_LIB="" for the base arm.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(3)
d <- data.frame(x = rnorm(40), n = 10L)
d$y <- rnorm(40)
d$yb <- rbinom(40, 10, 0.4)
w_out <- runif(40, 0.5, 2)
k <- 10L
tr <- function(lab, e) {
  r <- tryCatch({ e; "fits" }, error = function(err) conditionMessage(err))
  cat(sprintf("%-26s %s\n", lab, r))
}
tr("weights(w_out) env vector", frm(bf(y | weights(w_out) ~ x), data = d))
tr("trials(k) env scalar", frm(bf(yb | trials(k) ~ x), data = d,
                               family = binomial()))
tr("trials(k + 0L)", frm(bf(yb | trials(k + 0L) ~ x), data = d,
                         family = binomial()))
tr("weights(x^2 + k)", frm(bf(y | weights(x^2 + k) ~ x), data = d))
tr("trials(n) data column", frm(bf(yb | trials(n) ~ x), data = d,
                                family = binomial()))
