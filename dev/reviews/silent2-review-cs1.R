# Review of lane silent2: full warnings of cs_sr_acat seed 1 (sweep set of
# dev/silent2-ordmeasure.R), base and lane.
#   Rscript dev/reviews/silent2-review-cs1.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
gen <- function(seed, n = 300, K = 4, b1 = 1.5, b2 = -0.8, c1 = 1, c2 = -1, p = 0.4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, p)
  lat <- ifelse(cls == 1, b1 * d$x + c1, b2 * d$x + c2) + rlogis(n)
  cuts <- stats::quantile(lat, seq_len(K - 1) / K)
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}
w <- character()
f <- withCallingHandlers(frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()),
                             data = gen(1)),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
  message = function(m) invokeRestart("muffleMessage"))
cat("arm", arm, "logLik", as.numeric(logLik(f)), "code", f$opt$convergence, "\n")
for (x in w) cat("WARNING:", x, "\n\n")
s <- tryCatch(utils::capture.output(summary(f)), warning = function(x) paste("summary warning:", conditionMessage(x)))
cat(head(s, 60), sep = "\n")
