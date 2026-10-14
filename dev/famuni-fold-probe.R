# Which NaN-guard spellings survive RTMB's taping. Each guard g(x) should
# be exactly 0 with derivatives 0 for x > 0 and NaN for x <= 0. RTMB folds
# an operation with a constant-zero operand when it records the tape, so
# a guard can be NaN on doubles and 0 on the tape.
# Output: dev/famuni-out/fold-probe.txt
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
cat("RTMB", format(packageVersion("RTMB")), "\n")
guards <- list(
  "0 * log(x)" = function(x) 0 * log(x),
  "log(x) * 0" = function(x) log(x) * 0,
  "0 / (x + |x|)" = function(x) 0 / (x + abs(x)),
  "x * 0 + 0 * sqrt(x)" = function(x) 0 * sqrt(x),
  "(x + |x|) / (x + |x|) - 1" = function(x) (x + abs(x)) / (x + abs(x)) - 1,
  "log(x) - log(x)" = function(x) log(x) - log(x),
  "sqrt(x) - sqrt(x)" = function(x) sqrt(x) - sqrt(x)
)
cat(sprintf("%-28s %10s %10s | %10s %10s %10s %10s\n", "guard", "dbl x=-1",
            "dbl x=2", "tape x=-1", "tape x=0", "tape x=2", "grad x=2"))
for (nm in names(guards)) {
  g <- guards[[nm]]
  f <- function(p) g(p[1]) + p[1]^2
  tp <- MakeTape(f, 2)
  v <- function(x) tp(x) - x^2
  gr <- tp$jacobian(2)[1, 1] - 4
  h <- tp$jacfun()$jacobian(2)[1, 1] - 2
  cat(sprintf("%-28s %10s %10s | %10s %10s %10s %10s  hess %s\n", nm,
              format(suppressWarnings(g(-1))), format(g(2)),
              format(suppressWarnings(v(-1))), format(suppressWarnings(v(0))),
              format(v(2)), format(gr), format(h)))
}
# and through MakeADFun, the objective a fit uses
cat("\nthrough MakeADFun, objective = -(guard(p) - p^2):\n")
for (nm in c("0 * log(x)", "log(x) - log(x)")) {
  g <- guards[[nm]]
  obj <- MakeADFun(function(pr) -(g(pr$p) - pr$p^2), list(p = 2),
                   silent = TRUE)
  cat(sprintf("%-28s fn(-1) = %s, fn(2) = %s, gr(2) = %s\n", nm,
              format(suppressWarnings(obj$fn(-1))), format(obj$fn(2)),
              format(obj$gr(2))))
}
