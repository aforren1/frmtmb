# Punch round 1, minor 1: the newdata exposure check, against a mutant
# of aterms_for_newdata() without it (the lane's build before the fix).
# Seed 24. Output: dev/aterms2-log-mutant-newrate.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-aterms2/R/predict.R")
a <- grep("^aterms_for_newdata <- function", src)
b <- a + which(src[(a + 1):length(src)] == "}")[1]
fun <- src[a:b]
k <- grep('if (nm == "rate" && any(v <= 0, na.rm = TRUE)) {', fun,
          fixed = TRUE)
stopifnot(length(k) == 1L)
mut <- fun[-(k:(k + 3L))]
report <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-44s %s\n", label, paste(format(r), collapse = " ")))
}
set.seed(24)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
f <- frm(y | rate(time) ~ x, data = d, family = poisson())
nd <- d[1:3, ]
for (tm in c(-2, 0)) {
  nd$time[2] <- tm
  report(sprintf("lane, time = %g", tm),
         fitted(f, newdata = nd)[2, "Estimate"])
}
old <- get("aterms_for_newdata", envir = ns)
unlockBinding("aterms_for_newdata", ns)
eval(parse(text = mut), envir = ns)   # the definition assigns itself
for (tm in c(-2, 0)) {
  nd$time[2] <- tm
  report(sprintf("mutant without the check, time = %g", tm),
         fitted(f, newdata = nd)[2, "Estimate"])
}
assign("aterms_for_newdata", old, envir = ns)
lockBinding("aterms_for_newdata", ns)
