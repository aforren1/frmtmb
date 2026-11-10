# Part C: the quantities se_tier3() names the flat parameters with, on
# the retaken Hessian (dev/code8-ridge-exp.R's "retake" variant).
#   CODE8_ARM=before Rscript dev/code8-ridge-tau.R <j...>
js <- as.integer(commandArgs(TRUE))
src <- readLines("dev/code8-ridge-exp.R")
a <- c("retake", "tau")
eval(parse(text = src[grep("^source", src):grep("^set.seed", src)[1] - 1]))
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
y0 <- dn$yn
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
for (j in js) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  f <- quiet(frm(form, data = dn))
  p <- f$opt$par
  h <- fns$fit_outer_hessian(f)
  H <- h$H; E <- h$E
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  Es <- E / outer(D, D)
  nk <- sqrt(colSums((Es %*% e$vectors)^2))
  big <- max(abs(e$values))
  thr <- pmax(fns$se_flat_tol * big, fns$se_dir_mult * nk)
  v <- c(0, 0, p[3], -p[4], 0, 1); tu <- D * v; tu <- tu / sqrt(sum(tu^2))
  cat(sprintf("== j %d a %.2f k %.4f\n", j, p[3], p[4]))
  cat(" ev:", format(e$values, digits = 3), "\n nk:", format(nk, digits = 3),
      "\n thr:", format(thr, digits = 3), "\n")
  rem <- which(e$values <= thr)
  keep <- setdiff(seq_along(e$values), rem)
  gap <- min(e$values[keep]) - max(e$values[rem])
  tau <- max(sqrt(fns$se_flat_tol), fns$se_noise_mult * sqrt(sum(nk[rem]^2)) /
               gap)
  Fr <- e$vectors[, rem, drop = FALSE]
  proj <- sqrt(rowSums(Fr^2))
  cat(" rem:", rem, " gap", format(gap, digits = 3), " tau",
      format(tau, digits = 3), " flat_tol", format(min(tau, 0.5 * max(proj)),
                                                   digits = 3), "\n proj:",
      format(proj, digits = 3), "\n true tangent (unit diag):",
      format(tu, digits = 3), "\n tangent vs removed vec:",
      format(drop(crossprod(Fr, tu)), digits = 6), "\n")
}
cat("DONE\n")
