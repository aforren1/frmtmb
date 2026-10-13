# Reviewer of lane seflags, claim 1's split: a downward direction keeps
# 0.69.0's zeroed loadings. Where a downward direction is NOT confirmed
# concave (se_line_probe() FALSE, so the likelihood is flat along it),
# the lane names its dominant parameters "flat" but zeroes the loadings
# of the rest, so a parameter with loading just under half the largest
# keeps a finite SE and a prediction along it is called estimable. This
# prints, for every se_tier3() call, each downward direction, whether it
# was confirmed concave, and the largest loading it zeroed relative to
# its largest.
#   Rscript dev/seflags-rev-neg.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
trace("se_tier3", exit = quote({
  if (exists("neg") && length(neg)) {
    for (ii in seq_along(neg)) {
      kk <- neg[ii]
      av <- abs(e$vectors[, kk])
      zero <- mark == ""
      cat(sprintf(paste0("  TIER3 neg ev %.3g thr %.3g concave %s | largest",
                         " %.3f | largest zeroed %.3f (%s) ratio %.2f\n"),
                  ev[kk], thr[kk], conc[ii], max(av),
                  if (any(zero)) max(av[zero]) else 0,
                  if (any(zero)) nm[free][zero][which.max(av[zero])] else "-",
                  if (any(zero)) max(av[zero]) / max(av) else 0))
    }
  }
}), print = FALSE, where = ns)
show <- function(lab, r) {
  f <- r$v
  if (inherits(f, "rev_err")) { cat(lab, "ERROR", f, "\n"); return() }
  invisible(suppressWarnings(ns$sdr_of(f)))
  cat(sprintf("%-40s code %d | lost %s | warn %d\n", lab,
              f$opt$convergence, lost_of(f), length(r$w)))
}
set.seed(8)
n <- 120
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
show("spline by + gp", cap(frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),
                              data = d)))
# mo() seed 12 of lane nanse's study (brms_monotonic's data code)
for (seed in c(12, 71)) {
  set.seed(seed)
  income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(income_options, 100, TRUE),
                   levels = income_options, ordered = TRUE)
  mean_ls <- c(30, 60, 70, 75)
  ls <- mean_ls[income] + rnorm(100, sd = 7)
  age <- rnorm(100, mean = 40, sd = 10)
  dm <- data.frame(income, ls, age)
  show(sprintf("mo(income) * age seed %d", seed),
       cap(frm(ls ~ mo(income) * age, data = dm)))
}
# gr(g, by = f) fixture (lane nanse review)
src <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/setier-rev-a-rev-grby.R"
if (file.exists(src)) cat("(gr-by: see the lane's attack set, 0 diffs)\n")
# c0k sweep seeds: nonlinear with a group effect, FD Hessian
for (seed in 1:12) {
  set.seed(seed)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  show(sprintf("c0k seed %d", seed),
       cap(frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
                  nl = TRUE), data = dn)))
}
