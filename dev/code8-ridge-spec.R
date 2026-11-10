# Part C: the spectrum of the unit-diagonal outer Hessian on the ridge of
# c0 + exp(a)^k, and where the exact ridge tangent lies in it. The
# invariance (a0, k, theta_1) -> (s a0, k / s, theta_1 + log s) has the
# tangent (0, 0, a0, -k, 0, 1) in (c0, c0_x, a0, k, sigma, theta_1).
#   CODE8_ARM=before|after Rscript dev/code8-ridge-spec.R <j...>
js <- as.integer(commandArgs(TRUE))
source("dev/code8-common.R")
options(width = 200)
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
  nm <- fns$outer_par_names(f)
  h <- fns$fit_outer_hessian(f)
  H <- h$H
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  v <- c(0, 0, p[3], -p[4], 0, 1)
  tu <- D * v; tu <- tu / sqrt(sum(tu^2))
  g <- drop(f$obj$gr(p))
  cat(sprintf("== j %d  a %.3f k %.4f theta %.3f  code %d\n", j, p[3], p[4],
              p[6], f$opt$convergence))
  cat("grad:", format(g, digits = 3), "\n")
  cat("ev / max:", format(e$values / max(e$values), digits = 3), "\n")
  cat("noise ||Es v_k||:", format(sqrt(colSums(((h$E / outer(D, D)) %*%
                                                    e$vectors)^2)), digits = 3),
      "\n")
  cat("tangent loading on eigvecs:", format(drop(crossprod(e$vectors, tu)),
                                            digits = 3), "\n")
  cat("H tangent curvature (unit diag) t'Ht:",
      format(sum(tu * ((H / outer(D, D)) %*% tu)), digits = 3), "\n")
  print(round(e$vectors[, 6:4], 4))
  # the Hessian with the tangent projected out, then c0's variance
  P <- diag(6) - tcrossprod(tu)
  Hp <- P %*% (H / outer(D, D)) %*% P
  ep <- eigen(Hp, symmetric = TRUE)
  cat("projected ev / max:", format(ep$values / max(ep$values), digits = 3),
      "\n")
  Vp <- ep$vectors[, 1:5] %*% (t(ep$vectors[, 1:5]) / ep$values[1:5])
  cat("projected SEs:", format(sqrt(diag(Vp)) / D, digits = 4), "\n")
  sdr <- quiet(fns$sdr_of(f))
  cat("lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
}
cat("DONE\n")
