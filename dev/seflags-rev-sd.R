# Reviewer of lane seflags, claim 2: an sd at zero whose term the
# likelihood barely reads on se_sd_gain_up()'s grid (0.01 to 1 times the
# residual sd). Per case: the boundary message, the SE warning, the
# lost set, se_sd_gain_up()'s gain and `moved`, and the objective's
# change at a far larger sd (is the sd read at all?).
#   Rscript dev/seflags-rev-sd.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
show <- function(lab, r) {
  f <- r$v
  if (inherits(f, "rev_err")) { cat(sprintf("%-40s ERROR %s\n", lab, f)); return() }
  bm <- sum(grepl("[Bb]oundary", r$m))
  sw <- sum(grepl("^Standard errors are not available", r$w))
  on <- ns$outer_par_names(f)
  th <- which(grepl("^theta", on))
  info <- vapply(seq_along(th), function(k) {
    g <- tryCatch(ns$se_sd_gain_up0(f, k), error = function(e) NA)
    if (is.null(g)) g <- tryCatch(ns$se_sd_gain_up(f, k), error = function(e) NA)
    p <- f$opt$par
    f0 <- f$obj$fn(p)
    far <- vapply(c(5, 10), function(h) {
      q <- p; q[th[k]] <- p[th[k]] + h; f0 - f$obj$fn(q)
    }, 0)
    f$obj$fn(p)
    sprintf("%s=%.3g gain %.2g moved %s far(+5,+10) %.3g %.3g", on[th[k]],
            exp(p[th[k]]), as.numeric(g),
            format(attr(g, "moved") %||% NA, digits = 3), far[1], far[2])
  }, "")
  cat(sprintf("%-40s code %d | boundary msg %d | SE warn %d | lost %s\n      %s\n",
              lab, f$opt$convergence, bm, sw, lost_of(f),
              paste(info, collapse = " ; ")))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
for (seed in 1:3) {
  for (G in c(3, 5, 10)) {
    set.seed(seed)
    n <- 20 * G
    d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = 20)))
    d$y <- 1 + 0.5 * d$x + rnorm(n, 0, 1)
    for (k in c(1, 1e-2, 1e-3, 1e-4)) {
      d$kk <- k
      show(sprintf("seed %d G %d y ~ a + %g * u, u ~ 0 + (1|g)", seed, G, k),
           cap(frm(bf(y ~ a + kk * u, a ~ 1 + x, u ~ 0 + (1 | g),
                      nl = TRUE), data = d)))
    }
    show(sprintf("seed %d G %d y ~ x + (1|g)", seed, G),
         cap(frm(y ~ x + (1 | g), data = d)))
    # a group effect on sigma, the response in units of 1e-3 and 1e3
    for (u in c(1e-3, 1e3)) {
      d$yu <- d$y * u
      show(sprintf("seed %d G %d sigma ~ 1 + (1|g), units %g", seed, G, u),
           cap(frm(bf(yu ~ x, sigma ~ 1 + (1 | g)), data = d)))
    }
  }
}
