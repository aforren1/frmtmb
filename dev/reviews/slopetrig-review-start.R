# Review of lane slopetrig, check 7: the alternative the lane did not
# build. Start each random slope's log sd at the predictor unit over the
# column's root mean square (the default start shifted by -log(rms)),
# under autoscale = FALSE, with no pre-fit; the start is passed by the
# caller. Arms: F (autoscale = FALSE), C (F's own initial vector passed
# as start =, the control: must equal F), R (C with the slope's log sd
# moved by -log(rms)), D (the lane's default), Z (the centered,
# standardized column; the reference).
# Cases: the 5033 designs that 0.83.0 left short (dev/slopetrig-log/
# repro-base.log), then the out-of-sample designs where autoscale = TRUE
# itself stays short (two_slopes, nested_year, cs_year, grby_year).
# Usage: Rscript slopetrig-review-start.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
source("dev/reviews/slopetrig-review-designs.R")
mk <- function(seed, s, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  data.frame(g = g, x = c0 + s * z, y = eta + rnorm(300))
}
q <- function(expr) {
  nw <- 0L
  f <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      nw <<- nw + 1L; invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (!inherits(f, "error")) attr(f, "nw") <- nw
  f
}
ll <- function(f) if (inherits(f, "error")) NA_real_ else as.numeric(logLik(f))
nw <- function(f) if (inherits(f, "error")) NA_integer_ else attr(f, "nw")
shifted_start <- function(f) {
  p <- f$obj$par
  th <- which(names(p) == "theta")
  for (bk in f$frame$re_blocks) {
    cols <- ns$autoscale_block_cols(f$frame, bk)
    sdi <- as.integer(ns$block_sd_idx(bk))
    if (length(sdi) != bk$dim) next
    for (pos in seq_len(bk$dim)) {
      r <- cols[[pos]]$rms
      if (is.na(r)) next
      k <- th[bk$theta_idx[sdi[pos]]]
      p[k] <- p[k] - log(r)
    }
  }
  p
}
as_start <- function(p) split(unname(p), factor(names(p), unique(names(p))))
cases <- list()
for (seed in 1:5) cases[[length(cases) + 1]] <-
  list(paste("year int", seed), mk(seed, 3, 2000), y ~ 1 + (1 + x | g), gaussian(), NULL)
for (seed in c(1, 3, 4, 5)) cases[[length(cases) + 1]] <-
  list(paste("5+0.1z int", seed), mk(seed, 0.1, 5), y ~ 1 + (1 + x | g), gaussian(), NULL)
cases[[length(cases) + 1]] <-
  list("5+0.1z fix 2", mk(2, 0.1, 5), y ~ x + (1 + x | g), gaussian(), NULL)
for (dn in c("two_slopes", "nested_year", "cs_year", "grby_year")) {
  for (seed in 1:3) {
    ds <- designs[[dn]]
    cases[[length(cases) + 1]] <- list(paste(dn, seed), ds$data(seed), ds$form,
                                       ds$family, ds$ref)
  }
}
for (cs in cases) {
  d <- cs[[2]]; fo <- cs[[3]]; fam <- cs[[4]]
  dz <- if (is.null(cs[[5]])) { dd <- d; dd$x <- std(d$x); dd } else cs[[5]](d)
  fF <- q(frm(fo, data = d, family = fam, control = frmtmb_control(autoscale = FALSE)))
  p0 <- fF$obj$par
  fC <- q(frm(fo, data = d, family = fam, start = as_start(p0),
              control = frmtmb_control(autoscale = FALSE)))
  pr <- shifted_start(fF)
  fR <- q(frm(fo, data = d, family = fam, start = as_start(pr),
              control = frmtmb_control(autoscale = FALSE)))
  fD <- q(frm(fo, data = d, family = fam))
  fZ <- q(frm(fo, data = dz, family = fam))
  lz <- ll(fZ)
  best <- max(ll(fF), ll(fC), ll(fR), ll(fD), lz, na.rm = TRUE)
  cat(sprintf("%-16s short: F %8.3g C %8.3g R %8.3g D %8.3g Z %8.3g | warn F %s R %s D %s | start shift %s\n",
              cs[[1]], best - ll(fF), best - ll(fC), best - ll(fR), best - ll(fD),
              best - lz, nw(fF), nw(fR), nw(fD),
              paste(signif((pr - p0)[pr != p0], 3), collapse = ",")))
}
cat("DONE\n")
