# Short check, item 5: gp_block_lp() names the right predictor on every shape, and units follow. Then what happens on the spellings that still reach the live scan.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
set.seed(3); n <- 60; x <- runif(n, 0, 6); z <- runif(n, 0, 6)
d <- data.frame(x, z, f = factor(sample(c("a", "b"), n, TRUE)), y = 1000 * sin(x) + 100 * rnorm(n), y2 = cos(x) + rnorm(n, sd = .1), g = factor(rep(1:6, 10)))
sp <- list(one = y ~ gp(x), two = y ~ gp(x) + gp(z), by = y ~ gp(x, by = f), re = y ~ (1 | g) + gp(x), sig = bf(y ~ gp(x), sigma ~ gp(z)), mv = bf(y ~ gp(x)) + bf(y2 ~ gp(z)) + set_rescor(FALSE), nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE))
for (nm in names(sp)) {
  f <- suppressWarnings(frm(sp[[nm]], d))
  bks <- Filter(function(b) identical(b$covstruct, "gp"), f$frame$re_blocks)
  u <- ns$predictor_units(f$frame)
  cat(sprintf("RKY %-4s keys %s | term %s | units %s\n", nm, paste(vapply(bks, ns$gp_block_lp, "", frame = f$frame), collapse = ","), paste(vapply(bks, function(b) b$dpar %||% "?", ""), collapse = ","), paste(names(u), signif(u, 3), collapse = " ")))
}
cat("DONE\n")
