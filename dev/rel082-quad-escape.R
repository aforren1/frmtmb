# The double stationary escape of a quadrature fit (release review of
# 0.82.0, minor 6): quad_fit() escapes on the Laplace objective, then
# fit_assembled() escapes again on the rule. On skew_normal() fits of
# symmetric data with (1 | g), where alpha = 0 is reached legitimately,
# count each escape's restarts, its time and its gain, and read the two
# records the fit keeps.
#
#   Rscript dev/rel082-quad-escape.R <lib> [seeds]
args <- commandArgs(TRUE)
lib <- args[1]
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:10
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
calls <- new.env()
orig <- get("escape_stationary", ns)
unlockBinding("escape_stationary", ns)
assign("escape_stationary", function(obj, opt, frame, control, bounds,
                                     par_units = NULL, verbose = 0L,
                                     tally = NULL) {
  # quad_fit() passes no tally, fit_assembled() does
  where <- if (is.null(tally)) "laplace" else "rule"
  # the escape runs only where the stop is at the stationary point; the
  # record an earlier escape left on `opt` is not this call's
  n_st <- length(ns$stationary_escapes(obj, opt, frame))
  t0 <- proc.time()[["elapsed"]]
  o <- orig(obj, opt, frame, control, bounds, par_units, verbose, tally)
  el <- proc.time()[["elapsed"]] - t0
  r <- if (n_st) o$stationary_escape
  calls$log <- rbind(calls$log, data.frame(
    where = where, starts = n_st,
    gain = if (is.null(r)) NA else r[["gain"]], sec = el))
  o
}, envir = ns)
tot <- list(rule_runs = 0, rule_sec = 0, fit_sec = 0, lap_runs = 0)
for (s in seeds) {
  set.seed(s)
  n <- 200
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(n)
  calls$log <- NULL
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(suppressMessages(
    frm(bf(y ~ x + (1 | g)), family = skew_normal(), data = d,
        quadrature = TRUE)))
  el <- proc.time()[["elapsed"]] - t0
  lg <- calls$log
  ru <- lg[lg$where == "rule" & lg$starts > 0, , drop = FALSE]
  la <- lg[lg$where == "laplace" & lg$starts > 0, , drop = FALSE]
  tot$rule_runs <- tot$rule_runs + nrow(ru)
  tot$lap_runs <- tot$lap_runs + nrow(la)
  tot$rule_sec <- tot$rule_sec + sum(ru$sec)
  tot$fit_sec <- tot$fit_sec + el
  se <- f$opt$stationary_escape
  sl <- f$opt$stationary_escape_laplace
  cat(sprintf(paste0("SEED %2d fit %.2f s | laplace escape %s (gain %s) | ",
                     "rule escape %s restarts %s %.2f s (gain %s) | ",
                     "record gain %s, laplace record %s\n"),
              s, el, if (nrow(la)) "ran" else "no",
              if (nrow(la)) format(la$gain[1], digits = 3) else "-",
              if (nrow(ru)) "ran" else "no",
              if (nrow(ru)) ru$starts[1] else 0,
              if (nrow(ru)) ru$sec[1] else 0,
              if (nrow(ru)) format(ru$gain[1], digits = 3) else "-",
              if (is.null(se)) "none" else format(se[["gain"]], digits = 3),
              if (is.null(sl)) "none" else format(sl[["gain"]], digits = 3)))
}
cat(sprintf(paste0("SUMMARY %d fits: laplace escapes %d, rule escapes %d, ",
                   "rule escape time %.2f s of %.2f s fit time (%.1f%%)\n"),
            length(seeds), tot$lap_runs, tot$rule_runs, tot$rule_sec,
            tot$fit_sec, 100 * tot$rule_sec / tot$fit_sec))
# The stall case of test-quadrature-defects.R ("a quadrature fit started
# at the stall escapes it"): the start puts alpha at its stationary point
for (s in seeds) {
  set.seed(s)
  n <- 200
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
  m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
  bad <- list(betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3))
  calls$log <- NULL
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(suppressMessages(
    frm(bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1), family = skew_normal(),
        data = dd, quadrature = TRUE, start = bad)))
  el <- proc.time()[["elapsed"]] - t0
  lg <- calls$log
  ru <- lg[lg$where == "rule" & lg$starts > 0, , drop = FALSE]
  la <- lg[lg$where == "laplace" & lg$starts > 0, , drop = FALSE]
  se <- f$opt$stationary_escape
  sl <- f$opt$stationary_escape_laplace
  cat(sprintf(paste0("STALL %2d fit %.2f s | laplace escape %s (gain %s) | ",
                     "rule escape %s %.2f s (gain %s) | record gain %s, ",
                     "laplace record %s\n"),
              s, el, if (nrow(la)) "ran" else "no",
              if (nrow(la)) format(la$gain[1], digits = 3) else "-",
              if (nrow(ru)) "ran" else "no", if (nrow(ru)) ru$sec[1] else 0,
              if (nrow(ru)) format(ru$gain[1], digits = 3) else "-",
              if (is.null(se)) "none" else format(se[["gain"]], digits = 3),
              if (is.null(sl)) "none" else format(sl[["gain"]], digits = 3)))
}
