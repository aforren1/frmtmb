# Lane alapse, backlog 5449: the structures check_ald_laplace() did not
# profile, profiled as dev/wcheck2-ald-profile.R did: TMB's Laplace
# objective (frm(dry_run = "objective")) along every log sd of the block
# (the structure's own `sd_idx`), moved by t from the objective's
# optimum, every other outer parameter re-optimized. Also checks that
# each structure's covariance is homogeneous of degree 2 in its sds (the
# premise of the proof): vcov(theta with log sds + log c) / vcov(theta)
# = c^2 at a random theta.
#   Rscript dev/alapse-ald-profile.R <base|lane> [design ...]
# Data: seed 20261116, 20 groups x 4 times x 2 rows (160), so no effect
# is observation-level.
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/alapse-common.R")
set.seed(20261116)
n <- 160
d <- data.frame(g = factor(rep(1:20, each = 8)),
                tt = factor(rep(rep(1:4, each = 2), 20)))
d$x <- rnorm(n)
d$tn <- num_factor(c(0, 0.3, 0.9, 1.5)[as.integer(d$tt)])
d$pos <- num_factor(c(0.2, 1.1, 1.7, 3.0)[as.integer(d$tt)])
u <- matrix(rnorm(20 * 4, 0, 0.7), 20, 4)
d$y <- 1 + 0.5 * d$x + u[cbind(as.integer(d$g), as.integer(d$tt))] +
  rnorm(n)
Q <- diag(20) * 2
for (i in 1:19) Q[i, i + 1] <- Q[i + 1, i] <- -0.9
dimnames(Q) <- list(levels(d$g), levels(d$g))
AL <- asym_laplace()
designs <- list(
  ou = bf(y ~ x + ou(tn + 0 | g)),
  toep = bf(y ~ x + toep(tt + 0 | g)),
  hetar1 = bf(y ~ x + hetar1(tt + 0 | g)),
  homcs = bf(y ~ x + homcs(tt + 0 | g)),
  homtoep = bf(y ~ x + homtoep(tt + 0 | g)),
  exp = bf(y ~ x + exp(pos + 0 | g)),
  gau = bf(y ~ x + gau(pos + 0 | g)),
  mat = bf(y ~ x + mat(pos + 0 | g)),
  grprec = bf(y ~ x + (1 | gr(g, prec = Q))))
pick <- a[-1]
if (!length(pick)) pick <- names(designs)
q <- function(expr) suppressWarnings(suppressMessages(expr))
for (nm in pick) {
  fo <- designs[[nm]]
  cat(sprintf("--- %s: %s\n", nm, deparse(fo$formula)))
  r <- cap(frm(fo, family = AL, data = d, data2 = list(Q = Q)))
  cat("  frm():", if (inherits(r$v, "error")) {
    paste("ERROR", class(r$v)[1], substr(conditionMessage(r$v), 1, 160))
  } else {
    sprintf("fit, code %d, logLik %.4f, warnings %s", r$v$opt$convergence,
            as.numeric(logLik(r$v)), paste(unique(r$wc), collapse = ","))
  }, "\n")
  fr <- q(frm(fo, family = AL, data = d, data2 = list(Q = Q),
              dry_run = "frame"))
  ob <- q(frm(fo, family = AL, data = d, data2 = list(Q = Q),
              dry_run = "objective"))
  obj <- if (!is.null(ob$obj)) ob$obj else ob
  bk <- fr$re_blocks[[1]]
  e <- ns$cs_entry(bk$covstruct)
  si <- e$sd_idx(bk)
  # the premise: homogeneous of degree 2 in the sds
  th <- rnorm(length(bk$theta_idx), 0, 0.3)
  V1 <- e$vcov(th, bk)
  th2 <- th
  th2[si] <- th2[si] + log(3)
  V2 <- e$vcov(th2, bk)
  cat(sprintf("  vcov(sds x 3) / vcov: range %.15g to %.15g (9 is degree 2)\n",
              min(as.matrix(V2) / as.matrix(V1)),
              max(as.matrix(V2) / as.matrix(V1))))
  p <- obj$par
  ith <- which(names(p) == "theta")[bk$theta_idx[si]]
  f0 <- function(v) {
    v <- tryCatch(obj$fn(v), error = function(e) NA_real_)
    if (!is.finite(v)) 1e10 else v
  }
  op <- nlminb(p, f0, control = list(iter.max = 500, eval.max = 800))
  p <- op$par
  base_sd <- p[ith]
  free <- setdiff(seq_along(p), ith)
  g <- c(0, 1, 2, 4, 8)
  v <- vapply(g, function(t) {
    p[ith] <- base_sd + t
    f <- function(rr) { q2 <- p; q2[free] <- rr; f0(q2) }
    o <- nlminb(p[free], f, control = list(iter.max = 500, eval.max = 800))
    if (o$objective < 1e9) p[free] <<- o$par
    -o$objective
  }, 0)
  cat(sprintf("  optimum of the objective: -fn %.4f at log sd %s\n",
              -op$objective, paste(format(base_sd, digits = 4),
                                   collapse = " ")))
  cat("  -fn at log sd + 0, 1, 2, 4, 8:", sprintf("%.3f", v), "\n")
  cat(sprintf("  never falls: %s (largest fall %.3g)\n",
              all(diff(v) >= -1e-6 * pmax(1, abs(v[-1]))),
              max(0, -diff(v))))
}
cat("DONE\n")
