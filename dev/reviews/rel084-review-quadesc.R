# rel084 release review, quadesc fix round 1: in escape_stationary(), an
# abandoned restart (escape_track()) is vetoed only by a restart that ended
# on its own, gained AND left the point's window; the abandoned result
# carries evals and nonfinite_trials. New data (seed 7071, n 150, alpha 4)
# and a second mock where the finished restart leaves the window:
#   A: +2 runs to the maximum and raises; -2 returns the stall 1e-12 lower
#      (back on the point): the abandoned point must be kept.
#   B: +2 as A; -2 returns a point at alpha 0.8 (off the window), lower than
#      the stall and higher than the abandoned point: the veto must hold.
#   C: no mock (the real fit from alpha = 0).
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-quadesc.R
source("dev/reviews/rel084-review-common.R")
suppressMessages(library(testthat))
ns <- asNamespace("frmtmb")
set.seed(7071)
dl <- 4 / sqrt(17)
dd <- data.frame(y = 1 + dl * abs(rnorm(150)) + sqrt(1 - dl^2) * rnorm(150))
real <- ns$optimize_obj
alpha_of <- function(p) unname(p[names(p) == "betad"][2])
run_case <- function(mode) {
  st <- new.env()
  mock <- function(obj, control, bounds = list(lower = -Inf, upper = Inf),
                   par_units = NULL, verbose = 0L, start_par = obj$par,
                   tally = NULL) {
    a <- alpha_of(start_par)
    r <- real(obj, control, bounds, par_units, verbose, start_par, tally)
    if (isTRUE(all.equal(a, 2))) {
      st$max_obj <- r$objective
      stop("NA/NaN gradient evaluation")
    }
    if (isTRUE(all.equal(a, -2))) {
      r <- st$stall
      if (mode == "A") {
        r$objective <- r$objective - 1e-12
      } else {
        p <- r$par
        p[which(names(p) == "betad")[2]] <- 0.8
        r$par <- p
        r$objective <- as.numeric(obj$fn(p))
        st$off_obj <- r$objective
      }
    }
    if (is.null(st$stall)) st$stall <- r
    r
  }
  fit <- conds_of(with_mocked_bindings(
    frm(bf(y ~ 1), family = skew_normal(), data = dd,
        start = list(betad = c(0, 0))),
    optimize_obj = mock, .package = "frmtmb"))
  f <- fit$v
  cat(sprintf(paste0("CASE %s | stall alpha %.3g obj %.6f | max obj %.6f | off obj %s\n",
                     "   fit alpha %.4f code %s gain %s evals %s nonfinite %s | %s\n"),
              mode, alpha_of(st$stall$par), st$stall$objective,
              st$max_obj %||% NA, format(st$off_obj %||% NA, digits = 10),
              alpha_of(f$opt$par), f$opt$convergence,
              format(f$opt$stationary_escape[["gain"]], digits = 6),
              format(f$opt$evals %||% "NULL"),
              format(f$opt$nonfinite_trials %||% "NULL"), cls_line(fit)))
}
run_case("A")
run_case("B")
r <- conds_of(frm(bf(y ~ 1), family = skew_normal(), data = dd,
                  start = list(betad = c(0, 0))))
cat(sprintf("CASE C (no mock) | alpha %.4f code %d ll %.6f | %s\n",
            alpha_of(r$v$opt$par), r$v$opt$convergence, as.numeric(logLik(r$v)),
            cls_line(r)))
cat("DONE\n")
