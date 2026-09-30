# Reviewer check (lane ceplot): mm(g1, g2, by = cbind(f1, f2)) on draws,
# against brms 2.23.0 at the same draws (fixed_param, one chain per
# draw) where the two packages read the same new levels, and against
# the law of the new effect (per-draw residual over that draw's by-level
# sd) everywhere. Data seed 45 (the worker's test case), draws 1,
# curves 1 to 3 (brms) and 7 (residuals).
#   Rscript dev/ceplot-rev-mmby-brms.R > dev/ceplot-rev-log/mmby-brms.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
N <- 200
set.seed(45)
d <- data.frame(x = rnorm(300), g1 = factor(sample(1:10, 300, TRUE)),
                g2 = factor(sample(1:10, 300, TRUE)))
fl <- rep(c("a", "b"), each = 5)
d$f1 <- factor(fl[d$g1])
d$f2 <- factor(fl[d$g2])
u <- rnorm(10, 0, 1)
d$y <- rnorm(300, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]), 0.5)
fit <- frm(bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))),
           family = gaussian(), data = d)
for (bk in fit$frame$re_blocks) {
  cat("block", bk$term_label, "theta_idx", bk$theta_idx, "by level",
      bk$by$level, "\n")
}
ds <- hand(fit, N)
M <- ds$draws
rl <- sprintf("r_mmg1g2[%d,Intercept]", 1:10)
stopifnot(all(rl %in% colnames(M)))
th <- exp(M[, grep("^theta_", colnames(M))])
bylev <- vapply(fit$frame$re_blocks, function(bk) bk$by$level, "")
sd_a <- th[, which(bylev %in% c("a", "1"))]
sd_b <- th[, which(bylev %in% c("b", "2"))]
inits <- lapply(seq_len(N), function(i) {
  sdv <- c(sd_a[i], sd_b[i])
  list(b = array(M[i, "b_x"], 1),
       Intercept = M[i, "b_Intercept"] + mean(d$x) * M[i, "b_x"],
       sigma = M[i, "sigma"], sd_1 = matrix(sdv, 1, 2),
       z_1 = matrix(M[i, rl] / sdv[c(rep(1, 5), rep(2, 5))], 1))
})
b <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2))), data = d,
      algorithm = "fixed_param", chains = N, iter = 1, warmup = 0,
      init = inits, refresh = 0, seed = 1, silent = 2)))
bm <- as_draws_matrix(b)
say("brms draws %d; r max diff %.3g; sd names %s", ndraws(b),
    max(abs(bm[, rl] - M[, rl])),
    paste(grep("^sd_", colnames(bm), value = TRUE), collapse = ","))
ce <- function(o, cond, ...) {
  suppressMessages(suppressWarnings(conditional_effects(
    o, "x", resolution = 3, re_formula = NULL, conditions = cond, ...)))$x
}
conds <- list(
  C1_both_unset_a = list(f1 = "a", f2 = "a"),
  C2_g1seen_g2unset_b = list(g1 = "2", f1 = "a", f2 = "b"),
  C3_both_unset_a_b = list(f1 = "a", f2 = "b"),
  C0_both_seen = list(g1 = "2", g2 = "7", f1 = "a", f2 = "b"))
for (nm in names(conds)) {
  for (s in 1:3) {
    set.seed(s)
    a <- tryCatch(ce(b, conds[[nm]], sample_new_levels = "gaussian"), error = function(e) conditionMessage(e))
    if (is.character(a)) {
      say("%s seed %d: brms ERROR %s", nm, s, a)
      next
    }
    f <- ce(ds, conds[[nm]], seed = s)
    say("%s seed %d: brms est %s lo %s | frmtmb est %s lo %s | max diff est %.3g lo %.3g up %.3g",
        nm, s, f4(a$estimate__), f4(a$lower__), f4(f$estimate__),
        f4(f$lower__), max(abs(a$estimate__ - f$estimate__)),
        max(abs(a$lower__ - f$lower__)), max(abs(a$upper__ - f$upper__)))
  }
}
## the law, per draw: residual after the fixed part and the seen members
law <- function(cond, seen, var_fun) {
  e <- ce(ds, cond, seed = 7, spaghetti = TRUE)
  sp <- attr(e, "spaghetti")
  sp <- sp[order(sp$sample__, sp$x), ]
  xs <- sort(unique(sp$x))
  E <- matrix(sp$estimate__, ncol = length(xs), byrow = TRUE)
  known <- M[, "b_Intercept"] + outer(M[, "b_x"], xs)
  for (l in seen) known <- known + 0.5 * M[, l]
  R <- E - known
  z <- R[, 1] / sqrt(var_fun(seq_len(N)))
  c(flat = max(abs(R - R[, 1])), mean = mean(z), var = var(z))
}
say("law C1 (one shared level in block a, weights add to 1): %s",
    f4(law(conds$C1_both_unset_a, character(0), function(i) sd_a[i]^2)))
say("law C2 (member 2 new in block b, weight 0.5): %s",
    f4(law(conds$C2_g1seen_g2unset_b, "r_mmg1g2[2,Intercept]",
           function(i) 0.25 * sd_b[i]^2)))
say("law C3 (one new level in each block): %s",
    f4(law(conds$C3_both_unset_a_b, character(0),
           function(i) 0.25 * (sd_a[i]^2 + sd_b[i]^2))))
say("(var SE at N = %d is %.3f)", N, sqrt(2 / (N - 1)))
say("done")
