# Lane quadesc: which calibration candidate quad_fit() keeps on the three
# skew normals of item 5015 (cases as in dev/quadesc-probe.R), and the
# escape calls. Logs every quad_keep_best() candidate (stationary flag,
# objective, alpha) and every escape_stationary() call.
#   Rscript dev/quadesc-cands.R <lib|base> <case> [k]
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-quadesc-lib"
.libPaths(c(if (a[1] == "lib") lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
case <- a[2]
k <- if (length(a) > 2) as.numeric(a[3]) else 0
cat("lib:", dirname(find.package("frmtmb")), "CORETYPE:",
    Sys.getenv("OPENBLAS_CORETYPE"), "LAPACK:", La_library(), "\n")
ns <- asNamespace("frmtmb")
alpha_of <- function(par) {
  b <- par[names(par) == "betad"]
  paste(format(b[-1], digits = 4), collapse = ",")
}
wrap <- function(name, f) {
  unlockBinding(name, ns)
  assign(name, f, envir = ns)
}
o_keep <- ns$quad_keep_best
wrap("quad_keep_best", function(best, a) {
  if (!is.null(a) && !is.null(a$obj)) {
    g <- max(abs(a$obj$gr(a$opt$par)))
    cat(sprintf("  CAND stationary %s obj %.8f code %d grad %.3g alpha %s\n",
                a$stationary, a$opt$objective, a$opt$convergence, g,
                alpha_of(a$opt$par)))
  } else if (!is.null(a$retry)) {
    cat("  CAND retry at alpha", alpha_of(a$retry), "\n")
  } else cat("  CAND NULL\n")
  r <- o_keep(best, a)
  cat(sprintf("  KEEP obj %.8f\n", r$opt$objective))
  r
})
o_settle <- ns$quad_settle
wrap("quad_settle", function(obj, opt, tape_at, control, bounds, par_units,
                             vb = 0L, rounds = ns$quad_settle_rounds,
                             warn = TRUE) {
  r <- o_settle(obj, opt, tape_at, control, bounds, par_units, vb, rounds,
                warn)
  cat(sprintf("  SETTLE in %.8f alpha %s -> out %.8f alpha %s rounds %s\n",
              opt$objective, alpha_of(opt$par), r$opt$objective,
              alpha_of(r$opt$par),
              paste(format(r$opt$quad_settle, digits = 3), collapse = "/")))
  r
})
o_esc <- ns$escape_stationary
wrap("escape_stationary", function(obj, opt, frame, control, bounds,
                                   par_units = NULL, verbose = 0L,
                                   tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  st <- ns$stationary_escapes(obj, opt, frame)
  o <- o_esc(obj, opt, frame, control, bounds, par_units, verbose, tally)
  cat(sprintf("  ESC %s starts %d in %.8f out %.8f alpha_in %s alpha_out %s\n",
              where, length(st), opt$objective, o$objective,
              alpha_of(opt$par), alpha_of(o$par)))
  o
})
src <- readLines("dev/quadesc-probe.R")
# the data builders of the probe, without its run loop
i0 <- grep("^sn <- function", src)
i1 <- grep("^fit1 <- function", src) - 1L
eval(parse(text = src[i0:i1]))
m <- make(case)
d <- m$d
d$y <- d$y * (1 + 2e-15 * k)
w <- character()
f <- withCallingHandlers(
  frm(m$fo, family = skew_normal(), data = d, quadrature = TRUE,
      start = list(betad = m$start)),
  warning = function(x) {
    w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
  })
cat(sprintf("FIT %s k %g ll %.10f alpha %s esc %s lap %s warn %s\n", case, k,
            as.numeric(logLik(f)), alpha_of(f$opt$par),
            paste(format(f$opt$stationary_escape, digits = 4), collapse = "/"),
            paste(format(f$opt$stationary_escape_laplace, digits = 4),
                  collapse = "/"),
            paste(w, collapse = ",")))
cat("END\n")
