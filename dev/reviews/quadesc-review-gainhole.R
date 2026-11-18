# Review of lane quadesc: the second build keeps an abandoned restart's
# point only when no restart that ended on its own "gained", and a gain is
# any strict decrease of the objective. A restart that returns to the
# stationary point (SYMSMALL under OpenBLAS 0.3.32: alpha -0.0017) with an
# objective lower by rounding noise would count as a gain and drop the
# abandoned point. This prints, for each restart of the rule's escape, its
# outcome and its objective minus the stall's, to full precision.
#   Rscript dev/reviews/quadesc-review-gainhole.R <lib|base> <case> <k-expr>
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-quadesc-lib"
.libPaths(c(if (a[1] == "lib") lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
case <- a[2]
ks <- eval(parse(text = a[3]))
cat("lib:", dirname(find.package("frmtmb")), "BLAS:",
    extSoftVersion()[["BLAS"]], "LAPACK:", La_library(), "\n")
ns <- asNamespace("frmtmb")
src <- readLines("dev/quadesc-probe.R")
i0 <- grep("^sn <- function", src)
i1 <- grep("^fit1 <- function", src) - 1L
eval(parse(text = src[i0:i1]))
st <- new.env()
st$on <- FALSE
put <- function(name, f) {
  unlockBinding(name, ns)
  assign(name, f, envir = ns)
}
al_of <- function(p) {
  paste(format(p[names(p) == "betad"][-1], digits = 4), collapse = ",")
}
o_opt <- ns$optimize_obj
put("optimize_obj", function(obj, control, bounds = list(lower = -Inf,
                                                         upper = Inf),
                             par_units = NULL, verbose = 0L,
                             start_par = obj$par, tally = NULL) {
  if (!isTRUE(st$on)) {
    return(o_opt(obj, control, bounds, par_units, verbose, start_par,
                 tally))
  }
  r <- tryCatch(o_opt(obj, control, bounds, par_units, verbose, start_par,
                      tally), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("    restart from %s RAISED\n", al_of(start_par)))
    stop(r)
  }
  # QUADESC_REVIEW_NUDGE=1: a restart that ends back within the point's
  # window (|alpha| < 0.05) reads 1e-12 below the stall instead of above
  # it, the other side of the rounding that decides "gained"
  if (identical(Sys.getenv("QUADESC_REVIEW_NUDGE"), "1")) {
    ap <- r$par[names(r$par) == "betad"][-1]
    if (all(abs(ap) < 0.05)) r$objective <- st$f0 - 1e-12
  }
  cat(sprintf("    restart from %s ended code %d alpha %s obj-stall %.3e\n",
              al_of(start_par), r$convergence, al_of(r$par),
              r$objective - st$f0))
  r
})
o_esc <- ns$escape_stationary
put("escape_stationary", function(obj, opt, frame, control, bounds,
                                  par_units = NULL, verbose = 0L,
                                  tally = NULL) {
  rule <- !is.null(tally) && length(obj$env$random) == 0L
  if (rule && length(ns$stationary_escapes(obj, opt, frame))) {
    st$on <- TRUE
    st$f0 <- opt$objective
    cat(sprintf("  rule escape from alpha %s obj %.10f\n", al_of(opt$par),
                opt$objective))
  }
  o <- tryCatch(o_esc(obj, opt, frame, control, bounds, par_units, verbose,
                      tally), finally = st$on <- FALSE)
  o
})
m <- make(case)
for (k in ks) {
  d <- m$d
  d$y <- d$y * (1 + 2e-15 * k)
  cat(sprintf("== %s k %d\n", case, k))
  w <- character()
  f <- withCallingHandlers(
    frm(m$fo, family = skew_normal(), data = d, quadrature = TRUE,
        start = list(betad = m$start)),
    warning = function(x) {
      w <<- c(w, class(x)[1])
      invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage"))
  al <- f$opt$par[names(f$opt$par) == "betad"][-1]
  cat(sprintf("FIT %s k %d ll %.10f alpha %s code %d on_point %s warn %s\n",
              case, k, as.numeric(logLik(f)), al_of(f$opt$par),
              f$opt$convergence, all(abs(al) < 1),
              paste(sort(unique(w)), collapse = ",")))
}
cat("END\n")
