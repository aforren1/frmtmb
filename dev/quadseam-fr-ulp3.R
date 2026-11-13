# Fix round: the mi() runaway on ulp draw k (dev/quadseam-ulp.R
# "runaway"), where the lane ended without frmtmb_quad_unsettled on k = 3,
# 4, 7 and 14. The settle record, the warnings and the verbose stages.
#   Rscript dev/quadseam-fr-ulp3.R <lib> <k>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
k <- as.integer(a[2])
e <- 1 + 2e-15 * k
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
d$y <- (1 + 0.5 * d$x + rnorm(n)) * e
d$y[c(5, 17, 33)] <- NA
d$w <- 1.5
w <- character()
f <- withCallingHandlers(
  frm(bf(y | weights(w) + mi() ~ x + (1 | obs)), family = gaussian(),
      data = d, quadrature = TRUE, verbose = TRUE),
  warning = function(x) {
    w <<- c(w, paste(class(x)[1], substr(conditionMessage(x), 1, 200)))
    invokeRestart("muffleWarning")
  })
print(f$opt$quad_settle)
cat(f$opt$message, "\n")
for (x in w) cat("WARN", x, "\n")
