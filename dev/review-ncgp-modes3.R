# Item 1 and 2: a proposed scan, patched in this process, against the shipped one.
# The trial sd is y's residual scale, the trial sigma a ladder below it, and rho the
# lane's four values. All in y's own units, so the scan is equivariant. Read on the set
# of review-ncgp-modes2.R, on y * 1000, and for equivariance.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
orig <- ns$gp_start_scan
scan2 <- function(obj, template, frame, start) {
  bks <- Filter(function(bk) identical(bk[["covstruct"]], "gp"), frame[["re_blocks"]] %||% list())
  if (!length(bks) || !is.null(start[["theta"]])) return(obj)
  val <- function(tpl) {
    p <- ns$outer_from_template(tpl, obj, frame)
    if (is.null(p)) return(NA_real_)
    v <- tryCatch(obj$fn(p), error = function(e) NA_real_)
    if (isTRUE(is.finite(v))) v else NA_real_
  }
  y <- as.numeric(unlist(frame$y))
  sy <- stats::sd(y[is.finite(y)])
  has_sigma <- length(template[["betad"]]) >= 1L && is.finite(sy) && sy > 0
  best <- template
  fb <- val(best)
  if (is.na(fb)) fb <- Inf
  for (bk in bks) {
    pos <- bk[["theta_idx"]][-1L]
    d <- sqrt(max(Reduce(`+`, bk[["aux_D2"]])))
    if (!isTRUE(d > 0)) next
    for (r in d / 4^(0:3)) for (nr in if (has_sigma) c(1, 0.1, 1e-2, 1e-3) else 1) {
      tpl <- best
      tpl[["theta"]][pos] <- log(r)
      if (has_sigma) {
        tpl[["theta"]][bk[["theta_idx"]][1L]] <- log(sy)
        tpl[["betad"]][1L] <- log(sy * nr)
        tpl[["beta"]][1L] <- mean(y)
      }
      v <- val(tpl)
      if (!is.na(v) && v < fb) { fb <- v; best <- tpl }
    }
  }
  p <- ns$outer_from_template(best, obj, frame)
  if (!is.null(p)) obj$par <- p
  obj
}
set_scan <- function(f) { unlockBinding("gp_start_scan", ns); assign("gp_start_scan", f, envir = ns); lockBinding("gp_start_scan", ns) }
fns <- list(mix = function(x) 2 * sin(3 * x) + x, mix2 = function(x) sin(4 * x) + 0.5 * x,
            mix3 = function(x) 3 * cos(2 * x) - 2 * x, mix4 = function(x) sin(5 * x) + x^2 / 10,
            two = function(x) sin(x) + 0.3 * sin(10 * x), sin3 = function(x) sin(3 * x), sin10 = function(x) sin(10 * x))
tot <- 0; sh <- c(0, 0); sh2 <- c(0, 0)
for (seed in 1:2) for (fn in names(fns)) for (sg in c(0.3, 0.05, 1e-2)) for (n in c(40L, 120L)) {
  set.seed(31000 + seed * 7 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  d <- data.frame(y = y, x = x)
  set_scan(orig); r1 <- rev_fit(y ~ gp(x), d)
  set_scan(scan2); r2 <- rev_fit(y ~ gp(x), d)
  d3 <- d; d3$y <- d$y * 1000
  r3 <- rev_fit(y ~ gp(x), d3)
  set_scan(orig)
  g1 <- rev_gp_read(r1$fit); g2 <- rev_gp_read(r2$fit); g3 <- rev_gp_read(r3$fit)
  s <- cf_setup(y, x)
  m <- cf_max(s, extra = list(c(log(g1$rho), 2 * (log(g1$sigma) - log(g1$sd))), c(log(g2$rho), 2 * (log(g2$sigma) - log(g2$sd)))), ngrid = 26L)
  tot <- tot + 1
  a1 <- g1$ll - m$ll; a2 <- g2$ll - m$ll; a3 <- g3$ll + n * log(1000) - m$ll
  sh <- sh + c(a1 < -1e-3, a2 < -1e-3); sh2 <- sh2 + c(0, a3 < -1e-3)
  cat(sprintf("RY %s s=%g n %d seed %d shipped %.2e | scan2 %.2e conv %d w %d | scan2 y*1000 %.2e | equi ll %.1e sd %.1e rho %.1e sigma %.1e\n", fn, sg, n, seed, a1, a2, g2$conv, length(r2$warn), a3,
              abs(a3 - a2), abs(g3$sd / 1000 / g2$sd - 1), abs(g3$rho / g2$rho - 1), abs(g3$sigma / 1000 / g2$sigma - 1)))
  flush.console()
}
cat(sprintf("RYSUM total %d short shipped %d scan2 %d scan2 at y*1000 %d\n", tot, sh[1], sh[2], sh2[2]))
cat("DONE\n")
