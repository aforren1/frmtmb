# Lane wt-resmooth. Run the `fosr-fig` chunk of vignette("case-studies")
# as edited, because that chunk called
# frm_linpred(newdata = , re_formula = NA) with NO subject column, which
# this change turns into a named refusal. Checks that the edited call
# answers and that it IS the population coefficient function: the fs
# term contributes nothing at an unseen level, so the curve must equal
# mgcv's own exclusion of that term.
#   Rscript dev/resmooth-vignette-chunk.R > dev/resmooth-vignette-chunk.txt
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
set.seed(101)
N <- 40L
tt <- seq(0, 1, length.out = 21)
x <- rbinom(N, 1, 0.5)
b0 <- function(t) 1 + 2 * sin(2 * pi * t)
b1 <- function(t) 1.5 - 12 * (t - 0.5)^2
bi <- matrix(rnorm(N * 2, 0, 0.6), N, 2)
Y <- outer(rep(1, N), b0(tt)) + outer(x, b1(tt)) +
  bi[, 1] + outer(bi[, 2], tt - 0.5) * 2 +
  matrix(rnorm(N * length(tt), 0, 0.35), N, length(tt))
fos <- data.frame(subject = factor(rep(seq_len(N), each = length(tt))),
                  t = rep(tt, N), x = rep(x, each = length(tt)),
                  y = as.vector(t(Y)))
ffs <- suppressWarnings(frm(bf(y ~ s(t, k = 10) + s(t, by = x, k = 10) +
                                 s(t, subject, bs = "fs", k = 5)),
                            family = gaussian(), data = fos))
xg <- data.frame(t = seq(0, 1, length.out = 100), x = 0,
                 subject = factor("population",
                                  levels = c(levels(fos$subject),
                                             "population")))
xg1 <- transform(xg, x = 1)
f0 <- as.numeric(frm_linpred(ffs, newdata = xg, allow_new_levels = TRUE))
f1 <- as.numeric(frm_linpred(ffs, newdata = xg1,
                             allow_new_levels = TRUE)) - f0
cat("f0 range", sprintf("%.4f", range(f0)), "\n")
cat("f1 range", sprintf("%.4f", range(f1)), "\n")
cat("max |f0 - b0(t)|", sprintf("%.4f", max(abs(f0 - b0(xg$t)))), "\n")
cat("max |f1 - b1(t)|", sprintf("%.4f", max(abs(f1 - b1(xg$t)))), "\n")
gfs <- suppressWarnings(mgcv::gam(y ~ s(t, k = 10) + s(t, by = x, k = 10) +
                                    s(t, subject, bs = "fs", k = 5),
                                  data = fos, method = "ML"))
ref <- as.numeric(predict(gfs,
                          newdata = transform(xg,
                                              subject = fos$subject[1]),
                          exclude = "s(t,subject)"))
cat("max relative gap to mgcv exclude = s(t,subject):",
    sprintf("%.3e", max(abs(f0 - ref)) / max(abs(ref))), "\n")
# and the old call is now the named refusal the vignette documents
old <- tryCatch(frm_linpred(ffs,
                            newdata = xg[, c("t", "x")],
                            re_formula = NA),
                error = function(e) conditionMessage(e))
cat("old call:", substr(old, 1, 110), "\n")
