# Lane wcheck2, the asym_laplace refusal: for each term the Laplace route
# integrates in mu, frm()'s verdict, and the Laplace objective (TMB's,
# from frm(dry_run = "objective")) profiled along the term's scale: every
# log sd of its block is moved by t from the optimum of the objective,
# every other outer parameter re-optimized. The proof says the profile
# never falls for the refused terms; the cases the check leaves open
# are profiled too.
#   Rscript dev/wcheck2-ald-profile.R <base|lane|lib> [design ...]
# Data: seed 20261114, 80 rows, groups of 4, as the review's g4 design
# (dev/reviews/wcheck2-review-al.R), plus z and a time index.
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
ns <- asNamespace("frmtmb")
set.seed(20261114)
n <- 80
d <- data.frame(x = rnorm(n), z = runif(n), g = factor(rep(1:20, each = 4)),
                tt = factor(rep(1:4, 20)))
u <- rnorm(20, 0, 0.7)[as.integer(d$g)]
d$y <- 1 + 0.5 * d$x + sin(3 * d$z) + u + rnorm(n)
d$yp <- exp(d$y / 3)
q8 <- quantile(d$y, 0.85)
d$cc <- ifelse(d$y > q8, "right", "none")
d$yc <- pmin(d$y, q8)
AL <- asym_laplace()
designs <- list(
  g = list(bf(y ~ x + (1 | g)), AL),
  slope = list(bf(y ~ x + (1 + x | g)), AL),
  slope_uncor = list(bf(y ~ x + (1 + x || g)), AL),
  s = list(bf(y ~ x + s(z)), AL),
  t2 = list(bf(y ~ x + t2(z, x)), AL),
  gp = list(bf(y ~ x + gp(z)), AL),
  hsgp = list(bf(y ~ x + gp(z, k = 8, c = 5 / 4)), AL),
  ar1 = list(bf(y ~ x + ar1(tt + 0 | g)), AL),
  zial = list(bf(y ~ x + (1 | g)), zero_inflated_asym_laplace()),
  # left open
  log_link = list(bf(yp ~ x + (1 | g)), asym_laplace(link = "log")),
  sigma_re = list(bf(y ~ x, sigma ~ (1 | g)), AL),
  cens = list(bf(yc | cens(cc) ~ x + (1 | g)), AL),
  me = list(bf(y ~ me(x, sdx) + (1 | g)), AL),
  cs = list(bf(y ~ x + cs(1 + x | g)), AL),
  homdiag = list(bf(y ~ x + homdiag(1 + x | g)), AL),
  rr = list(bf(y ~ x + rr(1 + x | g, d = 1)), AL),
  grcov = list(bf(y ~ x + (1 | gr(g, cov = A))), AL))
A <- diag(20) + 0.3
dimnames(A) <- list(levels(d$g), levels(d$g))
d$sdx <- 0.2
pick <- a[-1]
if (!length(pick)) pick <- names(designs)
# the theta positions that hold a block's log sds: a field's natural
# rows of type "sd" (their jacobian), else the first ones of a grouping
# block (one per term, or one for a hom* structure)
sd_thetas <- function(fr) {
  out <- integer()
  for (bk in fr$re_blocks) {
    nr <- tryCatch(ns$cs_natural_rows(bk, rep(0, length(bk$theta_idx))),
                   error = function(e) NULL)
    if (!is.null(nr)) {
      g <- as.matrix(nr$grad)
      rows <- which(nr$type == "sd")
      cols <- which(colSums(abs(g[rows, , drop = FALSE])) > 0)
      out <- c(out, bk$theta_idx[cols])
    } else {
      k <- if (grepl("^hom", bk$covstruct)) 1L else bk$dim
      out <- c(out, bk$theta_idx[seq_len(min(k, length(bk$theta_idx)))])
    }
  }
  out
}
for (nm in pick) {
  ds <- designs[[nm]]
  cat(sprintf("--- %s: %s\n", nm, deparse(ds[[1]]$formula)))
  cat("  frm():", outcome(frm(ds[[1]], family = ds[[2]], data = d,
                               data2 = list(A = A))), "\n")
  fr <- q(frm(ds[[1]], family = ds[[2]], data = d, data2 = list(A = A),
              dry_run = "frame"))
  ob <- q(frm(ds[[1]], family = ds[[2]], data = d, data2 = list(A = A),
              dry_run = "objective"))
  obj <- if (!is.null(ob$obj)) ob$obj else ob
  p <- obj$par
  ith <- which(names(p) == "theta")[sd_thetas(fr)]
  fr0 <- function(v) {
    r <- tryCatch(obj$fn(v), error = function(e) NA_real_)
    if (!is.finite(r)) 1e10 else r
  }
  op <- nlminb(p, fr0, control = list(iter.max = 500, eval.max = 800))
  p <- op$par
  base_sd <- p[ith]
  free <- setdiff(seq_along(p), ith)
  g <- c(0, 1, 2, 4, 8)
  v <- vapply(g, function(t) {
    p[ith] <- base_sd + t
    f <- function(r) { q2 <- p; q2[free] <- r; fr0(q2) }
    o <- nlminb(p[free], f, control = list(iter.max = 500, eval.max = 800))
    if (o$objective < 1e9) p[free] <<- o$par
    -o$objective
  }, 0)
  cat(sprintf("  optimum of the objective: -fn %.4f at log sd %s\n",
              -op$objective, paste(format(base_sd, digits = 4),
                                   collapse = " ")))
  show_profile("  profile, log sd + t", g, v)
  cat(sprintf("  never falls: %s (largest fall %.3g)\n",
              all(diff(v) >= -1e-6 * pmax(1, abs(v[-1]))),
              max(0, -diff(v))))
}
cat("DONE\n")
