# Fix round 1, item 9 and 11: "se_row_curved() and the rho^2 floor do not touch a gp fit at ordinary theta",
# and M2. In process, one build: the two pieces reverted against as shipped, digit for digit, on 40 fits.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
patch <- function(nm, f) { environment(f) <- ns; unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
o_curved <- ns$se_row_curved; o_corr <- ns$gp_corr
old_corr <- function(theta, blk) {
  Q <- 0
  for (j in seq_along(blk[["aux_D2"]])) {
    rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
    Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
  }
  exp(-Q) + diag(gp_nugget, nrow(blk[["aux_D2"]][[1]]))
}
calls <- 0L; hits <- 0L
cnt <- function(j, fit, p, u) { calls <<- calls + 1L; v <- o_curved(j, fit, p, u); if (v) hits <<- hits + 1L; v }
never <- function(j, fit, p, u) FALSE
read <- function(fml, d, ...) {
  w <- character()
  f <- withCallingHandlers(frm(fml, d, ...), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  nd <- d[seq_len(min(5, nrow(d))), , drop = FALSE]; nd$x <- nd$x + 0.0137
  list(par = f$opt$par, ll = f$opt$objective, se = suppressWarnings(sqrt(diag(as.matrix(vcov(f, full = TRUE))))), fit = suppressWarnings(fitted(f)), new = tryCatch(suppressWarnings(fitted(f, newdata = nd)), error = function(e) NA),
       ran = tryCatch(unlist(ranef(f)), error = function(e) NA), w = w)
}
same <- function(a, b) identical(lapply(a, function(v) if (is.character(v)) v else sprintf("%a", as.numeric(unlist(v)))), lapply(b, function(v) if (is.character(v)) v else sprintf("%a", as.numeric(unlist(v)))))
nd <- 0L; rows <- 0L
for (seed in 1:4) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (sc in c(1, 3, 1000)) {
  set.seed(60000 + seed)
  n <- 50L
  x <- sort(round(runif(n, 0, 8), 4)); z <- rnorm(n); g <- factor(sample(letters[1:5], n, TRUE))
  y <- sc * (sin(2 * x) + 0.3 * z + sg * rnorm(n))
  d <- data.frame(x, z, g, y)
  fml <- switch(seed, y ~ gp(x), y ~ z + gp(x), y ~ gp(x) + (1 | g), y ~ gp(x, z))
  patch("se_row_curved", cnt); patch("gp_corr", o_corr)
  a <- read(fml, d)
  patch("se_row_curved", never); patch("gp_corr", old_corr)
  b <- read(fml, d)
  patch("se_row_curved", o_curved); patch("gp_corr", o_corr)
  s <- same(a, b); rows <- rows + 1L; if (!s) nd <- nd + 1L
  cat(sprintf("RCL %s s=%g c=%g identical %s; curved calls so far %d hits %d; SE NaN %d/%d\n", deparse(fml), sg, sc, s, calls, hits, sum(is.nan(a$se)), sum(is.nan(b$se))))
  flush.console()
}
cat(sprintf("RCLSUM rows %d differ %d se_row_curved calls %d true %d\n", rows, nd, calls, hits))
cat("== floor at extreme theta\n")
set.seed(1); x <- sort(runif(30, 0, 6)); d <- data.frame(x, y = sin(x) + rnorm(30, sd = .2))
f <- suppressWarnings(frm(y ~ gp(x), d)); bk <- f$frame$re_blocks[[1]]
for (lr in c(0, -100, -300, -354, -356, -400, -800)) {
  a <- tryCatch(ns$gp_corr(c(0, lr), bk), error = function(e) NULL); b <- old_corr; environment(b) <- ns; bb <- tryCatch(b(c(0, lr), bk), error = function(e) NULL)
  cat(sprintf("RCF log rho %g: identical to no floor %s; NaN new %d old %d; off-diagonal max %.2e\n", lr, identical(a, bb), sum(is.nan(a)), sum(is.nan(bb)), max(a[upper.tri(a)])))
}
cat("== M2: the rows behind the flat verdict at y * 1000, seed 14, noise 1e-4\n")
set.seed(14); x <- sort(runif(40, 0, 8)); y0 <- sin(x) + 1e-4 * rnorm(40)
for (c0 in c(1, 1000, 1e6, 1e9)) {
  f <- suppressWarnings(frm(y ~ gp(x), data.frame(y = y0 * c0, x = x)))
  fe <- suppressWarnings(fitted(f))
  se <- sqrt(diag(as.matrix(vcov(f, full = TRUE))))
  cat(sprintf("REQ2 c %g conv %d fitted SE NaN %d se %s (se/c intercept %.4g)\n", c0, f$opt$convergence, sum(is.nan(fe[, 2])), paste(signif(se, 4), collapse = " "), se[1] / c0))
}
# what the old steps read: reverted arm, at y * 1000
cat("DONE\n")
