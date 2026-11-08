# Review of lane nlstart: the units of nl_par_units() beside lane bigx's
# optimizer basis (opt_basis()), on four builds: rellib-r18, nlstart, bigx,
# and a trial merge (rev-nlstart-lib/merge-bigx, R/fit.R by diff3, bigx's
# other R files copied).
# Usage: REV_ARM=<arm or lib> Rscript dev/review-nlstart-merge.R
source("dev/review-nlstart-common.R")
ns <- asNamespace("frmtmb")
cat("has opt_basis:", exists("opt_basis", ns), " has nl_par_units:",
    exists("nl_par_units", ns), "\n")
show <- function(label, f, ref) {
  cat(sprintf("%-46s logLik %.6f  ref %.6f  short %.3g  conv %s  nwarn %d  units %s\n",
              label, llv(f), ref, ref - llv(f),
              if (iserr(f)) NA else f$opt$convergence, length(attr(f, "w")),
              if (iserr(f)) "ERR" else paste(format(f$par_units %||% 1, digits = 3),
                                             collapse = ",")))
  if (!iserr(f) && !is.null(f$obj$env$.frmtmb_opt_basis)) {
    cat("   opt basis positions:", f$obj$env$.frmtmb_opt_basis$pos, "\n")
  }
}
n <- 200
for (seed in c(6201, 6202, 6203)) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), wb = 1e6 + rnorm(n))
  d$y <- 0.5 * d$x + 0.3 * d$x * (d$wb - 1e6) + rnorm(n)
  d$xt <- d$x * 1e-20
  d$xtw <- d$xt * d$wb
  # the linear model y ~ 0 + xt + xtw is the body (a0 + a1 * wb) * xt
  lin <- tt(frm(y ~ 0 + xt + xtw, data = d))
  lin2 <- tt(frm(y ~ 0 + x + x:wb, data = d))
  ref <- max(llv(lin), llv(lin2), na.rm = TRUE)
  cat("== seed", seed, "\n")
  show("linear y ~ 0 + xt + xtw", lin, ref)
  show("linear y ~ 0 + x + x:wb", lin2, ref)
  show("body a * xt, a ~ 1 + wb", tt(frm(bf(y ~ a * xt, a ~ 1 + wb,
                                             nl = TRUE), data = d)), ref)
  show("body a * x, a ~ 1 + wb (no tiny x)", tt(frm(bf(y ~ a * x, a ~ 1 + wb,
                                                       nl = TRUE), data = d)), ref)
  # the lane's own defects on this build
  set.seed(seed)
  d2 <- data.frame(x = rnorm(n)); d2$y <- 0.5 * d2$x + rnorm(n)
  d2$xt <- d2$x * 1e-20
  show("lane defect 2: a * xt", tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE),
                                       data = d2)),
       llv(tt(frm(y ~ 0 + xt, data = d2))))
}
