# Probe E3 cleared every structural axis, so the stacked-PK failure is in
# one of the two remaining differences: deriving ns from length(y) inside
# the dynamics, or building y with the AD-overloaded c().
library(RTMB)
library(RTMBode)
times <- seq(0, 5, length.out = 9)
ns <- 12

run <- function(label, yexpr, dyn) {
  f <- function(par) {
    "c" <- RTMB::ADoverload("c")
    y0 <- yexpr()
    sol <- RTMBode::ode(y = y0, times = times, func = dyn,
                        parms = exp(par$lp), method = "lsoda",
                        atol = 1e-8, rtol = 1e-8)
    sum(sol[, -1]^2)
  }
  obj <- try(MakeADFun(f, list(lp = rep(log(0.3), 3 * ns)), silent = TRUE),
             silent = TRUE)
  v <- if (inherits(obj, "try-error")) NA else
    suppressWarnings(try(obj$fn(obj$par), silent = TRUE))
  cat(sprintf("%-46s fn = %s\n", label,
              if (inherits(v, "try-error") || is.na(v)) "FAIL/NaN" else
                format(v, digits = 8)))
}

dyn_fixed <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  h <- 1:ns
  A <- y[h]; C <- y[ns + h]
  list(c(-p[h] * A, p[h] * A / p[2 * ns + h] - p[ns + h] * C))
}
dyn_derived <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  m <- length(y) / 2
  h <- 1:m
  A <- y[h]; C <- y[m + h]
  list(c(-p[h] * A, p[h] * A / p[2 * m + h] - p[m + h] * C))
}
plain_y <- function() c(rep(100, ns), rep(0, ns))
ad_y <- function() { "c" <- RTMB::ADoverload("c")
                     c(rep(100, ns), rep(0, ns)) }

run("ns from closure, plain y", plain_y, dyn_fixed)
run("ns from length(y)/2, plain y", plain_y, dyn_derived)
run("ns from closure, AD-overloaded c() y", ad_y, dyn_fixed)
run("ns from length(y)/2, AD-overloaded c() y", ad_y, dyn_derived)

cat("\n--- what length does the dynamics actually see? ---\n")
seen <- c()
dyn_spy <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  seen <<- c(seen, length(y))
  h <- 1:ns
  A <- y[h]; C <- y[ns + h]
  list(c(-p[h] * A, p[h] * A / p[2 * ns + h] - p[ns + h] * C))
}
f <- function(par) {
  sol <- RTMBode::ode(y = c(rep(100, ns), rep(0, ns)), times = times,
                      func = dyn_spy, parms = exp(par$lp), method = "lsoda",
                      atol = 1e-8, rtol = 1e-8)
  sum(sol[, -1]^2)
}
obj <- MakeADFun(f, list(lp = rep(log(0.3), 3 * ns)), silent = TRUE)
invisible(obj$fn(obj$par))
cat("length(y) values passed to func:", paste(unique(seen), collapse = ", "),
    "  (calls:", length(seen), ")\n")
