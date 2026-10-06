# Reviewer of lane fixes, re-check: the cost of check_nl_identified()
# on large nonlinear models, against the fit it guards. The guard is
# timed alone on the built frame (min of 3), and the whole frm() with
# the guard and with it replaced by a no-op (the control), interleaved.
#   Rscript dev/fixes-rev2-nlcost.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
guard <- get("check_nl_identified", ns)
noop <- function(spec, frame, prior) invisible(NULL)
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
mk <- function(n, nlev) {
  set.seed(31)
  d <- data.frame(x = runif(n), z = rnorm(n),
                  f = factor(sample(seq_len(nlev), n, TRUE)))
  af <- rnorm(nlev, 2, 0.3)
  d$y <- af[d$f] * exp((0.5 + 0.2 * d$z) * d$x) + 1 + rnorm(n, 0, 0.3)
  d
}
cases <- list(
  list("n=200000, a ~ 1 + z, b ~ 1 + z", 2e5, 2,
       bf(y ~ a * exp(b * x) + c0, a ~ 1 + z, b ~ 1 + z, c0 ~ 1, nl = TRUE)),
  list("n=50000, a ~ 0 + f (300 levels), b ~ 1 + z", 5e4, 300,
       bf(y ~ a * exp(b * x) + c0, a ~ 0 + f, b ~ 1 + z, c0 ~ 1,
          nl = TRUE)))
for (cs in cases) {
  d <- mk(cs[[2]], cs[[3]])
  fr <- suppressMessages(frm(cs[[4]], data = d, dry_run = "frame"))
  g <- sapply(1:3, function(i) {
    gc()
    tm(guard(fr$spec, fr, NULL))
  })
  st <- if (cs[[3]] > 2) list(beta = c(rep(2, cs[[3]]), 0.5, 0, 1)) else
    list(beta = c(2, 0, 0.5, 0, 1))
  fits <- matrix(NA_real_, 2, 2, dimnames = list(c("guard", "noop"), NULL))
  for (r in 1:2) {
    for (arm in c("guard", "noop")) {
      assignInNamespace("check_nl_identified",
                        if (arm == "guard") guard else noop, ns = "frmtmb")
      gc()
      fits[arm, r] <- tm(suppressMessages(suppressWarnings(
        frm(cs[[4]], data = d, start = st))))
    }
  }
  assignInNamespace("check_nl_identified", guard, ns = "frmtmb")
  cat(sprintf("%-44s guard alone %s s | frm with guard %s s, with no-op %s s\n",
              cs[[1]], paste(format(g, digits = 3), collapse = " "),
              paste(format(fits["guard", ], digits = 3), collapse = " "),
              paste(format(fits["noop", ], digits = 3), collapse = " ")))
}
cat("peak memory (Mb, gc max used):", sum(gc()[, 6]), "\n")
