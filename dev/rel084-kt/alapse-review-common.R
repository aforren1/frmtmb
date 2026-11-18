# rel084-kt copy of dev/reviews/alapse-review-common.R: lane_lib = rellib-r23 and the lane arm drops rellib-r22 (rellib-r23 then user lib).
# Review of lane alapse: shared setup.
#   lib_arg "base" = rellib-r22 (0.83.0), "lane" = wt-alapse-lib first.
base_lib <- "C:/Users/adf44/source/r/rellib-r22"
user_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lane_lib <- "C:/Users/adf44/source/r/rellib-r23"
if (!exists("lib_arg")) lib_arg <- commandArgs(TRUE)[1]
.libPaths(if (identical(lib_arg, "base")) c(base_lib, user_lib) else
  c(lane_lib, user_lib))
suppressMessages(library(frmtmb))
suppressMessages(library(quantreg))
options(frmtmb.notices = FALSE, width = 120)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
`%||%` <- function(a, b) if (is.null(a)) b else a

cap <- function(expr) {
  w <- character(); wc <- character(); m <- character(); mc <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(z) {
      w <<- c(w, gsub("\\s+", " ", conditionMessage(z)))
      wc <<- c(wc, class(z)[1])
      invokeRestart("muffleWarning")
    },
    message = function(z) {
      m <<- c(m, gsub("\\s+", " ", conditionMessage(z)))
      mc <<- c(mc, class(z)[1])
      invokeRestart("muffleMessage")
    })
  list(v = v, w = w, wc = wc, m = m, mc = mc)
}

say <- function(lab, r) {
  if (inherits(r$v, "error")) cat("  ", lab, " ERROR [",
                                  class(r$v)[1], "]: ",
                                  conditionMessage(r$v), "\n", sep = "")
  for (i in seq_along(r$w)) {
    cat("  ", lab, " WARN [", r$wc[i], "]: ", substr(r$w[i], 1, 400),
        "\n", sep = "")
  }
  for (i in seq_along(r$m)) {
    cat("  ", lab, " MSG [", r$mc[i], "]: ", substr(r$m[i], 1, 300),
        "\n", sep = "")
  }
}

# weighted ALD log-likelihood
ald_ll <- function(y, mu, s, p, w = 1) {
  r <- y - mu
  sum(w * (log(p) + log1p(-p) - log(s) - r * (p - (r < 0)) / s))
}
chk <- function(r, p) sum(r * (p - (r < 0)))

# Largest decrease of the objective (negative log-lik) found by random
# perturbations of the outer parameters around the fit, per step size:
# a positive value means a nearby point has a higher likelihood.
local_probe <- function(f, nd = 400, seed = 1,
                        steps = c(1e-8, 1e-6, 1e-4, 1e-3, 1e-2)) {
  set.seed(seed)
  p0 <- f$opt$par
  f0 <- f$obj$fn(p0)
  out <- numeric(length(steps))
  for (k in seq_along(steps)) {
    best <- -Inf
    for (i in seq_len(nd)) {
      u <- stats::rnorm(length(p0))
      u <- u / sqrt(sum(u^2))
      v <- f0 - f$obj$fn(p0 + steps[k] * u)
      if (is.finite(v) && v > best) best <- v
    }
    out[k] <- best
  }
  f$obj$fn(p0)
  stats::setNames(out, format(steps))
}
