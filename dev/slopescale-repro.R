# Lane slopescale, item 1018: the triage's 33/36-fit table, before and
# after. A random slope whose column is not a fixed effect,
# y ~ 1 + (0 + x | g) or y ~ 1 + (1 + x | g), at column spreads 1e-3,
# 1e-4 and 1e-5, against the same data with the column at spread 1.
# Copied from round-briefs/triage-082/scratch-A/slope-unscaled.R, with
# the library arm as an argument and no time cap.
# Usage: Rscript slopescale-repro.R gau|pois before|after
a <- commandArgs(TRUE)
fam_arg <- a[1]
arm <- a[2]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "family", fam_arg, "arm", arm, "\n")
t0 <- proc.time()[[3]]
ll <- function(f) if (inherits(f, "try-error")) NA else as.numeric(logLik(f))
conv <- function(f) if (inherits(f, "try-error")) NA else f$opt$convergence
cnt <- new.env()
cnt$w <- 0L
cnt$m <- character(0)
fit <- function(fo, fm, d) {
  try(withCallingHandlers(
    frm(fo, family = fm, data = d),
    warning = function(w) {
      cnt$w <- cnt$w + 1L
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cnt$m <- c(cnt$m, class(m)[1])
      invokeRestart("muffleMessage")
    }),
    silent = TRUE)
}
fm <- if (fam_arg == "gau") gaussian() else poisson()
rows <- list()
for (form in c("0", "1")) for (s in c(1e-3, 1e-4, 1e-5)) for (seed in 1:6) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * x0
  y <- if (fam_arg == "gau") eta + rnorm(300) else rpois(300, exp(eta))
  d <- data.frame(g = g, x = x0 * s, xs = x0, y = y)
  fo_s <- if (form == "0") y ~ 1 + (0 + x | g) else y ~ 1 + (1 + x | g)
  fo_r <- if (form == "0") y ~ 1 + (0 + xs | g) else y ~ 1 + (1 + xs | g)
  w0 <- cnt$w
  cnt$m <- character(0)
  f_s <- fit(fo_s, fm, d)
  ws <- cnt$w - w0
  bnd <- "frmtmb_boundary_fit" %in% cnt$m
  f_r <- fit(fo_r, fm, d)
  rows[[length(rows) + 1L]] <- data.frame(
    form = form, s = s, seed = seed, ll_small = ll(f_s), ll_ref = ll(f_r),
    short = ll(f_r) - ll(f_s), conv_small = conv(f_s), warn_small = ws,
    boundary_msg = bnd)
}
r <- do.call(rbind, rows)
options(width = 160)
print(r, digits = 10)
cat("\nshort > 1e-6:", sum(r$short > 1e-6, na.rm = TRUE), "of", nrow(r),
    "; short > 1e-6 with no warning:",
    sum(r$short > 1e-6 & r$warn_small == 0, na.rm = TRUE),
    "; boundary message:", sum(r$boundary_msg),
    "; NA:", sum(is.na(r$short)), "; max short", max(r$short, na.rm = TRUE),
    "; max |short|", max(abs(r$short), na.rm = TRUE),
    "\nelapsed", proc.time()[[3]] - t0, "\n")
