## Aggregate dev/rtmbp-results/*.rds into the tables for benchmarks.md.
## Per-call costs pool every pass (3 interleaved calls-passes plus the
## full-mode run) because this laptop's thermal drift is large compared
## with the effect being measured; min is the least contaminated
## estimator and is reported next to the median.
out <- file("dev/rtmbp-tables.txt", open = "wt")
p <- function(...) cat(..., "\n", file = out, sep = "")

fs <- list.files("dev/rtmbp-results", full.names = TRUE)
R <- lapply(fs, readRDS)
names(R) <- basename(fs)

key <- function(x) sprintf("%s t%d ap%s", x$backend, x$threads, x$autopar)
shp <- function(x) x$shape

ord <- c("RTMB t1 apFALSE", "RTMBp t1 apFALSE", "RTMBp t1 apTRUE",
         "RTMBp t2 apTRUE", "RTMBp t4 apTRUE", "RTMBp t8 apTRUE",
         "RTMBp t16 apTRUE")

pool <- function(shape, k, field) {
  v <- unlist(lapply(R, function(x)
    if (shp(x) == shape && key(x) == k) x[[field]] else NULL))
  if (is.null(v)) numeric(0) else v
}
## fit / sdreport exist only in full-mode files
full1 <- function(shape, k, field) {
  v <- unlist(lapply(R, function(x)
    if (shp(x) == shape && key(x) == k && !is.null(x[[field]])) x[[field]] else NULL))
  if (is.null(v)) numeric(0) else v
}

fmt <- function(v, d = 3) if (!length(v)) "-" else
  sprintf(paste0("%.", d, "f / %.", d, "f / %.", d, "f"), min(v), median(v), max(v))

for (shape in c("A", "B")) {
  p("\n================ SHAPE ", shape, " ================")
  p("\n--- per-call, min / median / max (seconds), n = ",
    length(pool(shape, ord[1], "fn")), " pooled calls ---")
  p(sprintf("%-18s %-22s %-22s %-22s", "config", "fn", "gr (cached)", "gr (fresh)"))
  for (k in ord) {
    p(sprintf("%-18s %-22s %-22s %-22s", k,
              fmt(pool(shape, k, "fn")), fmt(pool(shape, k, "gr_cached")),
              fmt(pool(shape, k, "gr_fresh"))))
  }
  p("\n--- tape build, fit, sdreport: min / median / max (s) ---")
  p(sprintf("%-18s %-22s %-22s %-22s", "config", "tape", "fit", "sdreport"))
  for (k in ord) {
    p(sprintf("%-18s %-22s %-22s %-22s", k,
              fmt(pool(shape, k, "tape")), fmt(full1(shape, k, "fit"), 2),
              fmt(full1(shape, k, "sdreport"), 2)))
  }
  ## speedups vs the RTMB baseline, on mins
  b_fn <- min(pool(shape, ord[1], "fn"))
  b_gr <- min(pool(shape, ord[1], "gr_fresh"))
  b_ft <- if (length(full1(shape, ord[1], "fit"))) min(full1(shape, ord[1], "fit")) else NA
  b_sd <- if (length(full1(shape, ord[1], "sdreport"))) min(full1(shape, ord[1], "sdreport")) else NA
  p("\n--- speedup vs RTMB baseline (on minima) ---")
  p(sprintf("%-18s %8s %8s %8s %8s", "config", "fn", "gr", "fit", "sdreport"))
  for (k in ord) {
    f <- pool(shape, k, "fn"); g <- pool(shape, k, "gr_fresh")
    ft <- full1(shape, k, "fit"); sd <- full1(shape, k, "sdreport")
    p(sprintf("%-18s %8s %8s %8s %8s", k,
              if (length(f)) sprintf("%.2fx", b_fn / min(f)) else "-",
              if (length(g)) sprintf("%.2fx", b_gr / min(g)) else "-",
              if (length(ft)) sprintf("%.2fx", b_ft / min(ft)) else "-",
              if (length(sd)) sprintf("%.2fx", b_sd / min(sd)) else "-"))
  }
  ## whole-fit total: tape + fit + sdreport. Charges autopar's extra tape
  ## cost against its evaluation gains, which is what a user actually pays.
  p("\n--- whole fit = tape + fit + sdreport (s), and speedup ---")
  p(sprintf("%-18s %10s %10s %10s %10s", "config", "min", "median",
            "x (min)", "x (med)"))
  tot <- function(k, f) f(pool(shape, k, "tape")) + f(full1(shape, k, "fit")) +
    f(full1(shape, k, "sdreport"))
  bmin <- tot(ord[1], min); bmed <- tot(ord[1], median)
  for (k in ord) {
    if (!length(full1(shape, k, "fit"))) next
    p(sprintf("%-18s %10.2f %10.2f %9.2fx %9.2fx", k, tot(k, min), tot(k, median),
              bmin / tot(k, min), bmed / tot(k, median)))
  }

  ## optimizer path identity
  p("\n--- optimizer path / objective ---")
  for (k in ord) {
    x <- Filter(function(z) shp(z) == shape && key(z) == k && !is.null(z$fit_nll), R)
    if (!length(x)) next
    x <- x[[1]]
    p(sprintf("%-18s %d it, %d fn, %d gr, nll %.6f", k,
              x$fit_iter[length(x$fit_iter)], x$fit_fev[length(x$fit_fev)],
              x$fit_gev[length(x$fit_gev)], x$fit_nll[length(x$fit_nll)]))
  }
}

## ---- identity: RTMB vs RTMBp at identical parameter vectors ----------
p("\n================ IDENTITY ================")
for (shape in c("A", "B")) {
  base <- Filter(function(z) shp(z) == shape && key(z) == ord[1], R)[[1]]
  p("\nshape ", shape, ": max abs difference vs RTMB at 4 fixed parameter vectors")
  p(sprintf("%-18s %14s %14s %14s %14s", "config", "max|dfn|", "max rel dfn",
            "max|dgr|", "max rel dgr"))
  for (k in ord[-1]) {
    x <- Filter(function(z) shp(z) == shape && key(z) == k, R)
    if (!length(x)) next
    x <- x[[1]]
    dfn <- abs(x$id_fn - base$id_fn)
    rfn <- dfn / abs(base$id_fn)
    dgr <- max(abs(unlist(x$id_gr) - unlist(base$id_gr)))
    rgr <- max(abs((unlist(x$id_gr) - unlist(base$id_gr)) /
                     pmax(abs(unlist(base$id_gr)), 1e-8)))
    p(sprintf("%-18s %14.3e %14.3e %14.3e %14.3e", k, max(dfn), max(rfn), dgr, rgr))
  }
  ## fitted estimates and standard errors
  p("\nshape ", shape, ": fitted estimates / SEs vs RTMB (full-mode runs)")
  bf <- Filter(function(z) shp(z) == shape && key(z) == ord[1] && !is.null(z$fit_par), R)[[1]]
  p(sprintf("%-18s %14s %14s %14s", "config", "max|dpar|", "max|dest|", "max|dSE|"))
  for (k in ord[-1]) {
    x <- Filter(function(z) shp(z) == shape && key(z) == k && !is.null(z$fit_par), R)
    if (!length(x)) next
    x <- x[[1]]
    p(sprintf("%-18s %14.3e %14.3e %14.3e", k,
              max(abs(x$fit_par - bf$fit_par)),
              max(abs(x$sd_est - bf$sd_est)),
              max(abs(x$sd_se - bf$sd_se))))
  }
}

## reference fit
ref <- readRDS("dev/rtmbp-ref.rds")
p("\n================ REFERENCE ================")
p(sprintf("frm() InstEval logLik        : %.6f", ref$frm_logLik))
p(sprintf("hand-rolled RTMB shape A     : %.6f", -ref$nll_A))
p(sprintf("abs difference               : %.3e", abs(ref$frm_logLik + ref$nll_A)))
p(sprintf("frm() wall (se = FALSE)      : %.1f s", ref$frm_wall))
p(sprintf("shape B hand-rolled nll      : %.6f", ref$nll_B))
close(out)
cat(readLines("dev/rtmbp-tables.txt"), sep = "\n")
