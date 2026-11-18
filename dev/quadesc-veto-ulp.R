# Lane quadesc, fix round 1: the mock of test-quad-escape-record.R's
# "a restart back on the point does not veto an abandoned one" over 21
# roundings of its response (y times 1 + 2e-15 k, k = 0..20), one BLAS per
# process. The restart from alpha = +2 raises after it ends; the one from
# -2 returns the stall with an objective 1e-12 lower.
#   Rscript dev/quadesc-veto-ulp.R <lib|base>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] == "lib") "C:/Users/adf44/source/r/wt-quadesc-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "LAPACK:", La_library(),
    "\n")
ns <- asNamespace("frmtmb")
real <- ns$optimize_obj
st <- new.env()
alpha_of <- function(p) unname(p[names(p) == "betad"][2])
unlockBinding("optimize_obj", ns)
assign("optimize_obj", function(obj, control,
                                bounds = list(lower = -Inf, upper = Inf),
                                par_units = NULL, verbose = 0L,
                                start_par = obj$par, tally = NULL) {
  a <- alpha_of(start_par)
  r <- real(obj, control, bounds, par_units, verbose, start_par, tally)
  if (isTRUE(all.equal(a, 2))) stop("NA/NaN gradient evaluation")
  if (isTRUE(all.equal(a, -2))) {
    r <- st$stall
    r$objective <- r$objective - 1e-12
  }
  if (is.null(st$stall)) st$stall <- r
  r
}, envir = ns)
set.seed(4041)
dl <- 4 / sqrt(17)
y0 <- 1 + dl * abs(stats::rnorm(200)) + sqrt(1 - dl^2) * stats::rnorm(200)
ok <- 0L
for (k in 0:20) {
  st$stall <- NULL
  dd <- data.frame(y = y0 * (1 + 2e-15 * k))
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ 1), family = skew_normal(), data = dd,
        start = list(betad = c(0, 0))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  pass <- abs(alpha_of(st$stall$par)) < 0.05 &&
    abs(alpha_of(fit$opt$par)) > 1 && identical(fit$opt$convergence, 1L) &&
    grepl("restart off a stationary point", fit$opt$message, fixed = TRUE) &&
    fit$opt$stationary_escape[["gain"]] > 0 && !is.null(fit$opt$evals) &&
    !is.null(fit$opt$nonfinite_trials) && length(w) >= 1L &&
    all(grepl("Optimizer did not report convergence", w, fixed = TRUE))
  ok <- ok + pass
  cat(sprintf("K %2d pass %s alpha %.6g code %d gain %.6g evals %s nw %d\n",
              k, pass, alpha_of(fit$opt$par), fit$opt$convergence,
              fit$opt$stationary_escape[["gain"]],
              format(fit$opt$evals %||% NA), length(w)))
}
cat("PASS", ok, "of 21\n")
