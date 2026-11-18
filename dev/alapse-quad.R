# Lane alapse, backlog 5449: why quadrature = TRUE breaks down on
# asym_laplace() with (1 | g). The review's g4 design
# (dev/reviews/wcheck2-review-al.R: 80 rows, groups of 4, seeds
# 1 to 10). Per seed: the quadrature fit's outcome and, when it errs,
# which shape quad_breakdown_message() names; and where TMB's Laplace
# objective (the warm start quad_fit() anchors on) ends when nlminb
# runs it from the default start: the log sd, and the objective at log
# sd + 0, 2 and 6 with the rest re-optimized.
#   Rscript dev/alapse-quad.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(80), g = factor(rep(1:20, each = 4)))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[as.integer(d$g)] + rnorm(80)
  d
}
for (s in 1:10) {
  d <- mk(s)
  r <- cap(frm(y ~ x + (1 | g), family = asym_laplace(), data = d,
               quadrature = TRUE))
  out <- if (inherits(r$v, "error")) {
    m <- conditionMessage(r$v)
    paste("ERROR:", if (grepl("singular variance", m)) {
      sub(".*(on a singular variance component [(]sd [^)]*[)]).*", "\\1", m)
    } else if (grepl("calibration expired", m)) {
      "the calibration expired at every point tried"
    } else substr(m, 1, 120))
  } else {
    sprintf("fit, logLik %.4f, log sd %.3f, warnings %s",
            as.numeric(logLik(r$v)),
            r$v$opt$par[names(r$v$opt$par) == "theta"][1],
            paste(unique(r$wc), collapse = ","))
  }
  ob <- suppressWarnings(suppressMessages(
    frm(y ~ x + (1 | g), family = asym_laplace(), data = d,
        dry_run = "objective")))
  obj <- if (!is.null(ob$obj)) ob$obj else ob
  f0 <- function(v) {
    v <- tryCatch(obj$fn(v), error = function(e) NA_real_)
    if (!is.finite(v)) 1e10 else v
  }
  op <- nlminb(obj$par, f0, control = list(iter.max = 500, eval.max = 800))
  it <- which(names(op$par) == "theta")[1]
  prof <- vapply(c(0, 2, 6), function(t) {
    p <- op$par
    p[it] <- p[it] + t
    fr <- setdiff(seq_along(p), it)
    o <- nlminb(p[fr], function(r) { q <- p; q[fr] <- r; f0(q) },
                control = list(iter.max = 500, eval.max = 800))
    -o$objective
  }, 0)
  cat(sprintf("seed %2d | quad: %s | Laplace stop: log sd %.3f, -fn %.3f, code %d | -fn at log sd +0, +2, +6: %s\n",
              s, out, op$par[it], -op$objective, op$convergence,
              paste(sprintf("%.3f", prof), collapse = " ")))
}
cat("DONE\n")
