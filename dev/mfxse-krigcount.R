# Why the kriging row is not counted beside a new level (review m5): the
# extra variance of each part at the reviewer's rows, with and without the
# new level. Data seed 8150 (dev/review-mfxse-vcovwarn.R's gp + (1 | g)).
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
ns <- asNamespace("frmtmb")
set.seed(8150)
mg <- 90
dg <- data.frame(x = runif(mg, 0, 6), g = factor(sample(1:6, mg, TRUE)))
dg$y <- sin(dg$x) + rnorm(6, 0, 0.5)[dg$g] + rnorm(mg, sd = 0.3)
fg <- quiet(frm(y ~ gp(x) + (1 | g), data = dg))
th <- fg$estimates$theta
for (nd in list(data.frame(x = c(dg$x[1], 3.33, 9), g = factor(c("1", "new", "2"))),
                data.frame(x = c(dg$x[1], 3.33, 9), g = factor(c("1", "3", "2"))))) {
  lp <- fg$frame$linpreds[[1]]
  ed <- ns$lp_eta_design(fg, lp, nd, TRUE, TRUE)
  cat("n", ed$n, " sm_parts", length(ed$sm_parts), "\n")
  for (sp in ed$sm_parts) {
    cat(" block", sp$bk$covstruct, " theta_idx", sp$bk$theta_idx,
        " extra_var", format(sp$extra_var, digits = 3),
        " share", format(sp$extra_var / ns$sd_floored(th[sp$bk$theta_idx[1]])^2,
                         digits = 3), "\n")
  }
  ev <- ns$lp_extra_var(fg, ed, TRUE)
  cat(" new-level rows:", unlist(lapply(ev$new_levels, `[[`, "nas")),
      " gp parts:", length(ev$gp), "\n")
}
