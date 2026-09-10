# Which of the guard's operations RTMB can tape, and whether their
# derivatives are the ones the guard needs.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-07-adops.R
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
library(RTMB)

try_op <- function(nm, f) {
  z <- try(MakeTape(f, c(0.3, 2.0)), silent = TRUE)
  if (inherits(z, "try-error")) {
    cat(sprintf("  %-28s TAPE FAILED: %s\n", nm,
                sub("\n.*", "", conditionMessage(attr(z, "condition")))))
    return(invisible(NULL))
  }
  v <- try(z(c(0.3, 2.0)), silent = TRUE)
  g <- try(z$jacobian(c(0.3, 2.0)), silent = TRUE)
  cat(sprintf("  %-28s value %s  jac %s\n", nm,
              if (inherits(v, "try-error")) "ERR" else
                paste(format(v, digits = 6), collapse = " "),
              if (inherits(g, "try-error")) "ERR" else
                paste(format(as.numeric(g), digits = 6), collapse = " ")))
}

cat("scalar ops on an advector:\n")
try_op("abs(x[1]-x[2])",   function(x) abs(x[1] - x[2]))
try_op("pmin(x[1], x[2])", function(x) pmin(x[1], x[2]))
try_op("pmax(x[1], x[2])", function(x) pmax(x[1], x[2]))
try_op("min(x)",           function(x) min(x))
try_op("max(x)",           function(x) max(x))
try_op("pmin(x, 0.5)",     function(x) sum(pmin(x, 0.5)))
try_op("pmax(x, 0.5)",     function(x) sum(pmax(x, 0.5)))
try_op("abs-built min",    function(x) (x[1] + x[2] - abs(x[1] - x[2])) / 2)
try_op("abs-built max",    function(x) (x[1] + x[2] + abs(x[1] - x[2])) / 2)
try_op("abs-built clamp",  function(x) {
  lo <- 0; hi <- 0.99
  a <- (x[1] + hi - abs(x[1] - hi)) / 2
  (a + lo + abs(a - lo)) / 2
})

cat("\nvector clamp built from abs, on a length-3 advector:\n")
f <- try(MakeTape(function(x) {
  lo <- 0; hi <- 0.99
  a <- (x + hi - abs(x - hi)) / 2
  b <- (a + lo + abs(a - lo)) / 2
  sum(b * b)
}, c(-0.5, 0.4, 5)), silent = TRUE)
if (inherits(f, "try-error")) {
  cat("  FAILED:", conditionMessage(attr(f, "condition")), "\n")
} else {
  cat("  value at (-0.5, 0.4, 5):", f(c(-0.5, 0.4, 5)),
      " expect", sum(c(0, 0.4, 0.99)^2), "\n")
  cat("  jacobian:", format(as.numeric(f$jacobian(c(-0.5, 0.4, 5)))),
      " expect 0 0.8 0\n")
}
