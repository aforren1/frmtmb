# Is a whole model's tape the same with the per-scale floor as without
# it? One process, two tapes of one frame: the lane's, and one built
# with sd_floored() bound to exp() in the registry's functions. fn and
# gr, joint and Laplace, at random points. Fixture: test-reparam.R's
# rp_cor_case() data, y ~ x + (x | g).
#   Rscript dev/moparity-floor-ident.R
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(31)
ng <- 12L
dd <- data.frame(x = stats::rnorm(120),
                 g = factor(rep(seq_len(ng), length.out = 120)))
u <- matrix(stats::rnorm(2 * ng), 2)
dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + 0.7 * u[1, dd$g] +
                       0.4 * u[2, dd$g] * dd$x, 1)
fit <- frm(bf(y ~ x + (x | g)), family = gaussian(), data = dd)
ns <- asNamespace("frmtmb")
reg_new <- get("covstruct_registry", ns)
reg_old <- reg_new
old_env <- new.env(parent = ns)
old_env$sd_floored <- exp
for (nm in names(reg_old)) {
  f <- reg_old[[nm]]$nll
  if (is.function(f) && identical(environment(f), ns)) {
    environment(f) <- old_env
    reg_old[[nm]]$nll <- f
  }
}
mk <- function(reg, random) {
  unlockBinding("covstruct_registry", ns)
  assign("covstruct_registry", reg, envir = ns)
  lockBinding("covstruct_registry", ns)
  RTMB::MakeADFun(frmtmb:::build_objective(fit$frame), fit$frame$par_template,
                  random = random, silent = TRUE)
}
for (random in list(NULL, "b")) {
  o_new <- mk(reg_new, random)
  o_old <- mk(reg_old, random)
  mk(reg_new, NULL)
  set.seed(1)
  same_fn <- same_gr <- 0
  n <- 50
  for (i in seq_len(n)) {
    p <- o_new$par + stats::rnorm(length(o_new$par), 0, 0.5)
    same_fn <- same_fn + identical(o_new$fn(p), o_old$fn(p))
    same_gr <- same_gr + identical(o_new$gr(p), o_old$gr(p))
  }
  # the control: the two tapes are different tapes where the floor acts
  p <- o_new$par
  p[which(names(p) == "theta")[1L]] <- -1137.64
  cat(sprintf("%s: %d points, fn identical %d, gr identical %d | at log sd -1137.64: new %s, old %s\n",
              if (is.null(random)) "joint" else "Laplace", n, same_fn,
              same_gr, format(o_new$fn(p)), format(o_old$fn(p))))
}
