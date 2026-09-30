.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(44)
d4 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)))
d4$f <- factor(ifelse(as.integer(d4$g) <= 6, "a", "b"))
sdg <- ifelse(1:12 <= 6, 0.3, 1.5)
d4$y <- rnorm(240, 1 + 0.5 * d4$x + rnorm(12, 0, sdg)[d4$g], 0.5)
f4 <- frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(), data = d4)
for (bk in f4$frame$re_blocks) {
  cat("----\n"); str(bk[c("group_name", "term_label", "dim", "n_levels", "levels", "b_idx", "theta_idx", "by")], max.level = 2)
  str(lapply(bk$components, function(cp) cp[intersect(names(cp), c("by", "by_level", "by_var", "levels", "group_name", "label"))]), max.level = 3)
}
