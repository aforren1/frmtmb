# Probe E6: a one-state ode() lost its gradient in E5 but not in E2. The
# two differ only in how the solution matrix is subset (sol[, 2] drops to
# a vector; sol[, -1] stays a matrix) and in the time grid.
library(RTMB)
library(RTMBode)
dyn <- function(t, y, q) list(-q * y)

case <- function(label, times, pick) {
  f <- function(p) {
    sol <- RTMBode::ode(y = 1, times = times, func = dyn,
                        parms = exp(p$lr), method = "lsoda")
    sum(pick(sol)^2)
  }
  o <- MakeADFun(f, list(lr = log(0.3)), silent = TRUE)
  g <- suppressWarnings(tryCatch(format(o$gr(0)), error = function(e)
    paste("ERROR:", conditionMessage(e))))
  cat(sprintf("%-44s fn = %-10s gr = %s\n", label,
              format(o$fn(0), digits = 6), g))
}
case("times=c(0,1,2),  sol[, 2]  (drops to vector)", c(0, 1, 2),
     function(s) s[, 2])
case("times=c(0,1,2),  sol[, -1] (stays matrix)", c(0, 1, 2),
     function(s) s[, -1])
case("times=c(0,1,2),  sol[, 2, drop=FALSE]", c(0, 1, 2),
     function(s) s[, 2, drop = FALSE])
case("9 times,         sol[, 2]", seq(0, 5, length.out = 9),
     function(s) s[, 2])
case("9 times,         sol[, -1]", seq(0, 5, length.out = 9),
     function(s) s[, -1])

cat("\n-- 2-state system, same question --\n")
dyn2 <- function(t, y, q) { "c" <- RTMB::ADoverload("c")
                            list(c(-q[1] * y[1], q[1] * y[1] - q[2] * y[2])) }
case2 <- function(label, pick) {
  f <- function(p) {
    sol <- RTMBode::ode(y = c(1, 0), times = c(0, 1, 2), func = dyn2,
                        parms = exp(p$lr), method = "lsoda")
    sum(pick(sol)^2)
  }
  o <- MakeADFun(f, list(lr = c(log(0.3), log(0.1))), silent = TRUE)
  g <- suppressWarnings(tryCatch(paste(format(o$gr(o$par), digits = 5),
                                       collapse = " "),
                                 error = function(e)
                                   paste("ERROR:", conditionMessage(e))))
  cat(sprintf("%-44s fn = %-10s gr = %s\n", label,
              format(o$fn(o$par), digits = 6), g))
}
case2("sol[, 3] (one column, drops)", function(s) s[, 3])
case2("sol[, -1] (both columns)", function(s) s[, -1])
case2("sol[, 2:3]", function(s) s[, 2:3])
case2("sol[cbind(2:3, c(2,3))] (matrix index)",
      function(s) s[cbind(2:3, c(2, 3))])
