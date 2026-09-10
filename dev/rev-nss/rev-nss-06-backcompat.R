# REVIEW of lane nss, attack 7: is `ss_extrapolate = FALSE` really the
# base commit, bit for bit?
#
# The lane checked five compartment schedules. This widens it to 26,
# adding the shapes the correction's guards were written for and the
# ones most likely to route differently: an infusion, `addl`, an IV
# bolus with no depot, three compartments, an AUC compartment with NO
# steady state, an already-converged system whose differences are at
# the solver's noise floor, an oral depot whose differences are exactly
# zero, n_ss below three where `keep` is not full, and n_ss = 0 where
# there is no run-in at all.
#
# The comparison is `identical()` on the double vectors and, where that
# fails, a ULP count. A printed zero is not a measured zero.
#
# Run twice. REV_ARM=ref writes the base commit's answers, then the
# default arm reads them back and compares. `frm_ode()` absorbs unknown
# arguments into `...`, so the reference build accepts
# `ss_extrapolate` and ignores it; that is what makes one script serve
# both arms, and it is also the trap the lane recorded in section 9.
#
# Script path: dev/rev-nss/rev-nss-06-backcompat.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss"

one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
one_iv <- function(t, y, p) list(c(-p[1L] * y[1L]))
one_auc <- function(t, y, p) list(c(-p[2L] * y[1L],
                                    p[2L] * y[1L] - p[1L] * y[2L],
                                    y[2L]))
two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
three_iv <- function(t, y, p)
  list(c(-(p[1L] + p[2L] + p[4L]) * y[1L] + p[3L] * y[2L] +
           p[5L] * y[3L],
         p[2L] * y[1L] - p[3L] * y[2L],
         p[4L] * y[1L] - p[5L] * y[3L]))

cases <- list()
add <- function(tag, dyn, ns, pv, ev, out, n_ss = 20L, tt = NULL,
                atol = 1e-8, rtol = 1e-8) {
  cases[[tag]] <<- list(dyn = dyn, ns = ns, pv = pv, ev = ev,
                        out = out, n_ss = n_ss,
                        tt = tt %||% seq(0, ev$ii[[1L]] %||% 24,
                                         length.out = 17),
                        atol = atol, rtol = rtol)
}
evb <- function(ii, state = 1L, dur = NULL, addl = NULL) {
  e <- data.frame(time = 0, state = state, value = 100, ii = ii,
                  ss = TRUE)
  if (!is.null(dur)) e$duration <- dur
  if (!is.null(addl)) e$addl <- addl
  e
}
add("1 oral fast q8",      one_oral, 2L, list(0.2, 1.1),  evb(8),  2L)
add("1 oral slow q12",     one_oral, 2L, list(0.005, 1), evb(12), 2L)
add("1 oral flipflop q24", one_oral, 2L, list(0.03, 0.024), evb(24), 2L)
add("1 iv q8",             one_iv,   1L, list(0.2),      evb(8),  1L)
add("1 iv slow q24",       one_iv,   1L, list(0.01),     evb(24), 1L)
add("1 iv converged q24",  one_iv,   1L, list(2.0),      evb(24), 1L)
add("1 oral AUC q12",      one_auc,  3L, list(0.05, 1),  evb(12), 3L)
add("2 oral 107h q24",     two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L)
add("2 oral 107h q12",     two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(12), 2L)
add("2 oral 265h q24",     two_oral, 3L, list(0.1, 0.2, 0.008, 1),
    evb(24), 2L)
add("2 oral 670h q24",     two_oral, 3L, list(0.08, 0.15, 0.003, 1),
    evb(24), 2L)
add("2 oral 23h q8",       two_oral, 3L, list(0.2, 0.4, 0.1, 1.1),
    evb(8), 2L)
add("3 iv q24",            three_iv, 3L,
    list(0.15, 0.3, 0.02, 0.05, 0.004), evb(24), 1L)
add("2 oral infusion 2h",  two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24, dur = 2), 2L)
add("1 iv infusion 4h",    one_iv,   1L, list(0.02),
    evb(24, dur = 4), 1L)
add("1 oral addl 3",       one_oral, 2L, list(0.05, 1),
    evb(12, addl = 3L), 2L, tt = seq(0, 48, length.out = 17))
add("2 oral addl 2 q24",   two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24, addl = 2L), 2L, tt = seq(0, 72, length.out = 17))
add("n_ss 0",              two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, n_ss = 0L)
add("n_ss 1",              two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, n_ss = 1L)
add("n_ss 2",              two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, n_ss = 2L)
add("n_ss 3",              two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, n_ss = 3L)
add("n_ss 4",              two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, n_ss = 4L)
add("n_ss 137",            two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, n_ss = 137L)
add("tight tol 1e-12",     two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, atol = 1e-12, rtol = 1e-12)
add("loose tol 1e-6",      two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24), 2L, atol = 1e-6, rtol = 1e-6)
add("2 oral dose central",  two_oral, 3L, list(0.15, 0.3, 0.02, 1),
    evb(24, state = 2L), 2L)

runcase <- function(cs, ext) {
  args <- list(cs$dyn, init = rep(list(0), cs$ns), times = cs$tt,
               parms = cs$pv, events = cs$ev, output = cs$out,
               n_ss = cs$n_ss, ss_tol = 1e-6, atol = cs$atol,
               rtol = cs$rtol)
  if (!is.null(ext)) args$ss_extrapolate <- ext
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(do.call(frm_ode, args), error = function(e)
      paste("ERROR:", conditionMessage(e))),
    warning = function(e) { w <<- c(w, conditionMessage(e))
                            invokeRestart("muffleWarning") })
  list(v = if (is.character(v)) v else as.numeric(v), w = w)
}

if (identical(arm, "ref")) {
  res <- lapply(cases, runcase, ext = NULL)
  saveRDS(res, file.path(OUT, "rev-nss-06-ref.rds"))
  cat("\nwrote", length(res), "reference cases\n")
} else {
  ref <- readRDS(file.path(OUT, "rev-nss-06-ref.rds"))
  ulp <- function(a, b) {
    if (length(a) != length(b)) return(NA_integer_)
    max(vapply(seq_along(a), function(j) {
      x <- a[[j]]; y <- b[[j]]
      if (identical(x, y)) return(0L)
      if (!is.finite(x) || !is.finite(y)) return(NA_integer_)
      n <- 0L
      z <- x
      while (n < 2000L && z != y) {
        z <- nextafter <- if (z < y) z + max(.Machine$double.xmin,
                                             abs(z) * 2^-52) else
          z - max(.Machine$double.xmin, abs(z) * 2^-52)
        n <- n + 1L
      }
      if (z == y) n else NA_integer_
    }, 0L))
  }
  cat(sprintf("\n%-22s %6s %8s %12s %12s\n", "case", "ident", "ulp",
              "maxabsdiff", "default-arm gain"))
  nid <- 0L
  for (tag in names(cases)) {
    a <- ref[[tag]]$v
    b <- runcase(cases[[tag]], FALSE)$v
    e <- runcase(cases[[tag]], TRUE)$v
    if (is.character(a) || is.character(b)) {
      cat(sprintf("%-22s %6s (error path)\n", tag,
                  identical(a, b)))
      if (identical(a, b)) nid <- nid + 1L
      next
    }
    id <- identical(a, b)
    nid <- nid + id
    md <- max(abs(a - b))
    gn <- if (is.character(e)) NA_real_ else
      max(abs(as.numeric(e) - a)) / max(1e-300, md)
    cat(sprintf("%-22s %6s %8s %12.5e %12s\n", tag, id,
                if (id) "0" else format(ulp(a, b)), md,
                if (is.na(gn)) "-" else format(signif(gn, 3))))
  }
  cat(sprintf("\nidentical on %d of %d cases\n", nid, length(cases)))
  cat("\n-- warnings, base against lane FALSE arm --\n")
  for (tag in names(cases)) {
    wa <- ref[[tag]]$w
    wb <- runcase(cases[[tag]], FALSE)$w
    if (length(wa) != length(wb))
      cat(sprintf("%-22s base %d warning(s), lane FALSE %d\n", tag,
                  length(wa), length(wb)))
  }
}
cat("\ndone\n")
