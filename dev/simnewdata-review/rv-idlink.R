# A partial re_formula on a block shared across dpars by an ID,
# (1 + x | p | g) in mu and (1 | p | g) in sigma: which columns are kept,
# and is the dropped slope drawn from its exact conditional law given
# BOTH kept intercepts?
#   Rscript dev/simnewdata-review/rv-idlink.R > .../log/idlink.txt
source("dev/simnewdata-review/rv-prelude.R")
set.seed(41)
n_g <- 40
d <- data.frame(x = rnorm(n_g * 25), g = factor(rep(seq_len(n_g), 25)))
Vt <- matrix(c(1, .6, .4, .6, .8, .5, .4, .5, .5), 3)
U <- matrix(rnorm(n_g * 3), n_g) %*% chol(Vt)
d$y <- 1 + d$x + U[d$g, 1] + U[d$g, 2] * d$x +
  rnorm(nrow(d), 0, exp(-0.5 + 0.5 * U[d$g, 3]))
fit <- frm(bf(y ~ x + (1 + x | p | g), sigma ~ (1 | p | g)), data = d)
bl <- fit$frame$re_blocks
cat("blocks:", length(bl), " dims:", vapply(bl, `[[`, 0, "dim"), "\n")
for (cp in bl[[1]]$components) {
  cat("  component", cp$lp_key %||% "?", "offset", cp$offset,
      "cols", length(cp$cnms %||% NA), "\n")
}
th <- fit$estimates$theta
V <- as.matrix(frmtmb:::covstruct_registry[[bl[[1]]$covstruct]]$vcov(
  th[bl[[1]]$theta_idx], bl[[1]]))
cat("V:\n"); print(round(V, 3))
vc <- varcorr_matrices(fit)
cat("varcorr_matrices names:", names(vc), "\n")
print(lapply(vc, round, 3))
for (rf in list(~ (1 | g), ~ (1 + x | g))) {
  plan <- tryCatch(frmtmb:::sim_re_plan(fit, rf), error = function(e) e)
  if (inherits(plan, "error")) {
    cat(deparse1(rf), "ERROR:", conditionMessage(plan), "\n"); next
  }
  if (!length(plan$cond)) {
    cat(deparse1(rf), ": redraw", length(plan$redraw), "cond none\n"); next
  }
  keep <- plan$cond[[1]]$keep
  cat(deparse1(rf), ": kept columns", paste(which(keep), collapse = ","),
      "of", length(keep), "\n")
  K <- which(keep); D <- which(!keep)
  B0 <- matrix(fit$estimates$b[bl[[1]]$b_idx], nrow = 3)
  lev <- 5
  A <- V[D, K, drop = FALSE] %*% solve(V[K, K])
  m <- as.vector(A %*% B0[K, lev])
  C <- V[D, D, drop = FALSE] - A %*% V[K, D, drop = FALSE]
  set.seed(1)
  R <- 20000
  dr <- replicate(R, matrix(frmtmb:::sim_draw_b(fit, plan)[bl[[1]]$b_idx],
                            nrow = 3)[D, lev])
  dr <- matrix(dr, nrow = length(D))
  cat(sprintf("   mean z %s  var ratio %s\n",
              paste(format((rowMeans(dr) - m) / sqrt(diag(C) / R), digits = 3),
                    collapse = " "),
              paste(format(apply(dr, 1, var) / diag(C), digits = 4),
                    collapse = " ")))
}
