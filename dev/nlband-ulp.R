# Lane nlband: the hole start over roundings of the response a few ulps
# apart (y * (1 + 2e-15 k), k = 0 to 23), in one process. Two sets:
#   band: the band data (seeds 1 to 3); each default fit against the fit
#     from (1, -1.9) on the same draw;
#   worse: the truth beyond the data, b = 0.5, with the same hole (seeds 1
#     to 3); each default fit against the same fit with nl_pole_holes()
#     returning NULL, which is 0.83.0's escape, compared bitwise (logLik,
#     estimates, standard errors, conditions).
# NLBAND_ARM=after Rscript dev/nlband-ulp.R
source("dev/nlband-common.R")
holes_fn <- ns$nl_pole_holes
set_holes <- function(on) {
  f <- if (on) holes_fn else function(w) NULL
  unlockBinding("nl_pole_holes", ns)
  assign("nl_pole_holes", f, envir = ns)
  lockBinding("nl_pole_holes", ns)
}
sig <- function(r) {
  if (iserr(r)) return(paste("ERR", r$msg))
  se <- tryCatch(sqrt(diag(suppressWarnings(vcov(r)))), error = function(e) NA)
  list(ll = llv(r), beta = unname(r$estimates$beta), se = unname(se),
       w = attr(r, "wc"), m = attr(r, "mc"))
}
for (s in 1:3) {
  d0 <- band_data(s)
  for (k in 0:23) {
    d <- d0
    d$y <- d$y * (1 + 2e-15 * k)
    ref <- frm(band_form, data = d, start = list(beta = c(1, -1.9)))
    r <- tt(frm(band_form, data = d))
    cat(sprintf("ROW band seed %d k %d | ll %.10f ref %.10f diff %.3g b %.6f code %s msg %s warn %s\n",
                s, k, llv(r), llv(ref), llv(ref) - llv(r),
                r$estimates$beta[[2L]], r$opt$convergence, msum(r), wsum(r)))
  }
}
for (s in 1:3) {
  d0 <- band_data(s)
  d0$y <- 1 / (d0$x + 0.5) + rnorm(nrow(d0), 0, 0.2)
  for (k in 0:23) {
    d <- d0
    d$y <- d$y * (1 + 2e-15 * k)
    set_holes(TRUE)
    on <- tt(frm(band_form, data = d))
    set_holes(FALSE)
    off <- tt(frm(band_form, data = d))
    set_holes(TRUE)
    cat(sprintf("ROW worse seed %d k %d | identical %s ll %.12f b %.8f msg %s warn %s\n",
                s, k, identical(sig(on), sig(off)), llv(on),
                on$estimates$beta[[2L]], msum(on), wsum(on)))
  }
}
cat("DONE\n")
