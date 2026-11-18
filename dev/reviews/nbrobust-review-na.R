# Review of lane nbrobust: can a missing response reach nb_lpmf()? Fits
# with an NA count through mi(), subset() and plain NA rows, on one build.
#   Rscript dev/reviews/nbrobust-review-na.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(3)
d <- data.frame(x = stats::rnorm(80))
d$y <- stats::rnbinom(80, size = 2, mu = exp(1 + 0.3 * d$x))
d$z <- stats::rnorm(80, d$x)
d$s <- rep(c(TRUE, FALSE), 40)
d$y[c(3, 10)] <- NA
try1 <- function(label, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) {
    cat(sprintf("%-28s %s\n", label, substr(r, 1, 200)))
  } else {
    cat(sprintf("%-28s nobs %d logLik %.8f\n", label, stats::nobs(r),
                as.numeric(logLik(r))))
  }
}
fams <- list(nb = negbinomial(), zinb = zero_inflated_negbinomial(),
             hnb = hurdle_negbinomial(), geo = geometric(),
             nb1 = nbinom1())
for (nm in names(fams)) {
  try1(paste(nm, "plain NA rows"),
       frm(y ~ x, family = fams[[nm]], data = d))
  try1(paste(nm, "mi()"),
       frm(bf(y | mi() ~ x), family = fams[[nm]], data = d))
  try1(paste(nm, "mv subset()"),
       frm(bf(y | subset(s) ~ x, family = fams[[nm]]) +
             bf(z ~ x, family = gaussian()), data = d))
}
