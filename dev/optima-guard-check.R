# Lane optima, item 3: does a guard that is NaN for a non-positive gap
# and exactly 0 otherwise survive RTMB's tape (CppAD folds 0 * x)?
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
g1 <- function(g) 0 * log(g + abs(g))
g2 <- function(g) {
  h <- log(g + abs(g))
  h - h
}
x <- c(-1, 0, 1e-300, 1e-10, 2)
for (nm in c("g1", "g2")) {
  f <- get(nm)
  tp <- MakeTape(function(g) f(g), x)
  cat(nm, "value:", format(tp(x)), "\n")
  J <- MakeTape(function(g) f(g), x)$jacobian(x)
  cat(nm, "diag jacobian:", format(diag(J)), "\n")
}
g3 <- function(g) log(g + abs(g)) - log(2 * g)
g4 <- function(g) log(g) - log(g)
for (nm in c("g3", "g4")) {
  f <- get(nm)
  tp <- MakeTape(function(g) f(g), x)
  cat(nm, "value:", format(tp(x)), "\n")
  J <- MakeTape(function(g) f(g), x)$jacobian(x)
  cat(nm, "diag jacobian:", format(diag(J)), "\n")
}
g7 <- function(g) log(g + abs(g)) - log(2 * abs(g))
for (nm in c("g7")) {
  f <- get(nm)
  tp <- MakeTape(function(g) f(g), x)
  cat(nm, "value:", format(tp(x)), "\n")
  J <- MakeTape(function(g) f(g), x)$jacobian(x)
  cat(nm, "diag jacobian:", format(diag(J)), "\n")
  cat(nm, "numeric, warnings:\n")
  print(withCallingHandlers(f(x), warning = function(w) {
    cat("  WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }))
}
