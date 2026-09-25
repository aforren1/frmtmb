# Punch round 2, item 10: the minimal reproductions for the two upstream
# reports drafted in dev/phase3b-rtmb-report-pnorm.md (RTMB) and
# dev/phase3b-rtmb-report-tanh.md (TMB). Plain RTMB, no frmtmb. The
# truth is the Mills-ratio expansion, not a second evaluation of
# dnorm/pnorm, since exp(dnorm(log) - pnorm(log.p)) shares the defect.
# Usage: Rscript dev/phase3b-rtmb-repro.R
# Output: dev/phase3b-log/rtmb-repro.txt
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
out <- sprintf("RTMB %s, TMB %s, %s", packageVersion("RTMB"),
               packageVersion("TMB"), R.version.string)

# (1) pnorm(x, log.p = TRUE): first and second derivatives, x << 0
tp <- MakeTape(function(x) pnorm(x, log.p = TRUE), 0)
hp <- tp$jacfun()
d1 <- function(x) -x - 1 / x + 2 / x^3
d2 <- function(x) -1 + 1 / x^2 - 6 / x^4
out <- c(out, "", "(1) d/dx and d2/dx2 of pnorm(x, log.p = TRUE)")
for (x in c(-1e3, -1e5, -1e6, -1e7, -2e8)) {
  g <- tp$jacobian(x)
  out <- c(out, sprintf(
    "x = %6.0e  d tape %.9g  truth %.9g  rel %.2e | d2 tape %.6g truth %.6g",
    x, g, d1(x), abs(g - d1(x)) / abs(d1(x)), hp$jacobian(x), d2(x)))
}
out <- c(out, sprintf("x^2 eps / 2 at -1e5, -1e7: %.2e, %.2e",
                      1e10 * .Machine$double.eps / 2,
                      1e14 * .Machine$double.eps / 2))

# (2) tanh: second derivative past |x| = 710 (TMB's TanhOp::reverse)
t4 <- MakeTape(function(x) tanh(x), 0)
h4 <- t4$jacfun()
out <- c(out, "", "(2) derivatives of tanh(x); truth 0 to double for both")
for (x in c(700, 710, 711, 1e4, -711)) {
  out <- c(out, sprintf("x = %6g  tanh d %g d2 %g", x, t4$jacobian(x),
                        h4$jacobian(x)))
}
dir.create("dev/phase3b-log", showWarnings = FALSE)
writeLines(out, "dev/phase3b-log/rtmb-repro.txt")
cat(out, sep = "\n")
