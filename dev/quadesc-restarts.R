# Lane quadesc: item 5015's cause. Inside the rule's stationary escape
# (fit_assembled()), each escape start is optimized two ways:
#   frozen     on the fit's rule, frozen where quad_fit() anchored it
#              (what escape_stationary() does)
#   reanchored on a rule taped at the start, optimized and settled
#              (quad_fit()'s `restart`, as quad_sd_up_escape() uses)
# plus the Laplace escape's own endpoint as a third start. Each endpoint
# is also scored by the exact marginal log-likelihood (one-dimensional
# integrals by stats::integrate(), split at 400 points over +-12 sd), which
# depends on no rule.
#   Rscript dev/quadesc-restarts.R <lib|base> <case> [k]
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
src <- readLines("dev/quadesc-probe.R")
i0 <- grep("^sn <- function", src)
i1 <- grep("^fit1 <- function", src) - 1L
eval(parse(text = src[i0:i1]))
m <- make(case)
d <- m$d
d$y <- d$y * (1 + 2e-15 * k)
ax <- case == "alphax"
# the exact marginal log-likelihood at an outer vector (beta, betad, theta)
exact_ll <- function(par) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  sdg <- exp(par[names(par) == "theta"][1])
  mu <- b[1] + b[2] * d$x
  sig <- exp(bd[1])
  al <- if (ax) bd[2] + bd[3] * d$x else rep(bd[2], nrow(d))
  tot <- 0
  for (g in levels(d$g)) {
    i <- which(d$g == g)
    lf <- function(u) {
      vapply(u, function(v) {
        sum(RTMBdist::dskewnorm2(d$y[i], mu[i] + v, sig, al[i], log = TRUE))
      }, 1) + stats::dnorm(u, 0, sdg, log = TRUE)
    }
    grid <- seq(-12 * sdg, 12 * sdg, length.out = 4001)
    lv <- lf(grid)
    mx <- max(lv[is.finite(lv)])
    br <- seq(-12 * sdg, 12 * sdg, length.out = 401)
    s <- 0
    for (j in seq_len(400)) {
      s <- s + stats::integrate(function(u) exp(lf(u) - mx), br[j],
                                br[j + 1], rel.tol = 1e-10,
                                abs.tol = 0)$value
    }
    tot <- tot + mx + log(s)
  }
  tot
}
alpha_of <- function(par) {
  b <- par[names(par) == "betad"]
  paste(format(b[-1], digits = 4), collapse = ",")
}
wrap <- function(name, f) {
  unlockBinding(name, ns)
  assign(name, f, envir = ns)
}
stash <- new.env()
o_qf <- ns$quad_fit
wrap("quad_fit", function(...) {
  r <- o_qf(...)
  stash$restart <- r$restart
  r
})
o_esc <- ns$escape_stationary
wrap("escape_stationary", function(obj, opt, frame, control, bounds,
                                   par_units = NULL, verbose = 0L,
                                   tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  st <- ns$stationary_escapes(obj, opt, frame)
  o <- o_esc(obj, opt, frame, control, bounds, par_units, verbose, tally)
  cat(sprintf("ESC %s starts %d in %.8f out %.8f alpha_in %s alpha_out %s\n",
              where, length(st), opt$objective, o$objective,
              alpha_of(opt$par), alpha_of(o$par)))
  if (where == "laplace") {
    stash$lap <- o$par
    cat(sprintf("  laplace endpoint exact %.6f\n", exact_ll(o$par)))
  }
  if (where == "rule" && length(st)) {
    saved <- ns$obj_state_save(obj)
    starts <- c(st, list(stash$lap))
    lbl <- c(paste0("from", seq_along(st)), "lapend")
    for (j in seq_along(starts)) {
      p <- starts[[j]]
      vb0 <- obj$env$value.best
      fr <- tryCatch(ns$optimize_obj(obj, control, bounds, par_units,
                                     start_par = p),
                     error = function(e) e)
      if (inherits(fr, "error")) {
        pb <- obj$env$last.par.best[obj$env$lfixed()]
        fb <- tryCatch(obj$fn(pb), error = function(e) NA)
        gb <- tryCatch(obj$gr(pb), error = function(e) NA)
        cat(sprintf(paste0("  %s frozen ERROR %s | value.best before %.6f",
                           " after %.6f | best alpha %s fn %.6f gr finite %s",
                           " exact %.6f\n"),
                    lbl[j], conditionMessage(fr), vb0, obj$env$value.best,
                    alpha_of(pb), fb, all(is.finite(gb)), exact_ll(pb)))
        fr <- NULL
      } else {
        cat(sprintf(paste0("  %s frozen obj %.6f code %d alpha %s",
                           " exact %.6f\n"),
                    lbl[j], fr$objective, fr$convergence, alpha_of(fr$par),
                    exact_ll(fr$par)))
      }
      ra <- tryCatch(stash$restart(p), error = function(e) NULL)
      if (is.null(ra)) {
        cat(sprintf("  %s reanchored NULL\n", lbl[j]))
      } else {
        cat(sprintf(paste0("  %s reanchored obj %.6f code %d alpha %s",
                           " exact %.6f settle %s\n"),
                    lbl[j], ra$opt$objective, ra$opt$convergence,
                    alpha_of(ra$opt$par), exact_ll(ra$opt$par),
                    paste(format(ra$info$settled), ra$info$rounds)))
      }
    }
    ns$obj_state_restore(obj, saved)
  }
  o
})
w <- character()
f <- withCallingHandlers(
  frm(m$fo, family = skew_normal(), data = d, quadrature = TRUE,
      start = list(betad = m$start)),
  warning = function(x) {
    w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
  })
cat(sprintf("FIT %s k %g ll %.10f exact %.6f alpha %s warn %s\n", case, k,
            as.numeric(logLik(f)), exact_ll(f$opt$par), alpha_of(f$opt$par),
            paste(w, collapse = ",")))
cat("END\n")
