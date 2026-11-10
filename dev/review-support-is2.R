# Reviewer probe (lane support): the GEV evaluated on doubles inside the
# fitting paths that read the objective off the tape: importance
# sampling's freeze/verify (R/importance.R) and nested quadrature's
# verify (R/quadrature-nested.R). FRMTMB_PORT_LIB picks the build.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
rgev <- function(n, mu, sigma, xi) mu + sigma * ((-log(runif(n)))^(-xi) - 1) / xi
run <- function(label, expr) {
  t0 <- proc.time()[[3]]
  r <- tryCatch(suppressWarnings(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-34s ERROR: %s\n", label,
                substr(conditionMessage(r), 1, 220)))
  } else {
    cat(sprintf("%-34s logLik %.10g  beta %s  (%.1fs)\n", label,
                as.numeric(logLik(r)),
                paste(format(r$opt$par[names(r$opt$par) == "beta"],
                             digits = 8), collapse = " "),
                proc.time()[[3]] - t0))
  }
  invisible(r)
}
for (seed in 1:4) for (xi in c(1, 1.5)) {
  set.seed(seed)
  ng <- 15
  m <- 8
  g <- rep(seq_len(ng), each = m)
  u <- rnorm(ng, 0, 0.7)
  d <- data.frame(g = factor(g),
                  y = rgev(ng * m, 1 + u[g], 1, xi))
  cat("seed", seed, "xi", xi, "\n")
  run("  laplace", frm(y ~ 1 + (1 | g), family = gen_extreme_value(), data = d))
  run("  importance = 500", frm(y ~ 1 + (1 | g), family = gen_extreme_value(),
                                data = d, importance = 500L))
}
