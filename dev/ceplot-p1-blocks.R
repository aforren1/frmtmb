# Lane ceplot punch 1: the grouping-factor fields of the blocks of the
# shapes B1 names, to key brms's per-group "old_levels" choice.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(3)
d <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)), h = factor(rep(1:8, 30)))
d$f <- factor(ifelse(as.integer(d$g) <= 6, "a", "b"))
d$g2 <- factor(sample(1:12, 240, TRUE))
d$y <- rnorm(240, 1 + rnorm(12)[d$g])
d$yl <- exp(d$y / 4)
show <- function(fit) {
  for (bk in fit$frame[["re_blocks"]]) {
    cat(sprintf("  group_name=%s term=%s cov=%s by=%s nlev=%d first=%s\n",
                format(bk[["group_name"]]), format(bk[["term_label"]]),
                bk[["covstruct"]], format(bk[["by"]][["level"]]),
                length(bk[["levels"]]), paste(head(bk[["levels"]], 3), collapse = ",")))
  }
}
cat("mu+sigma:\n"); show(frm(bf(yl ~ x + (1 | g), sigma ~ (1 | g)), family = lognormal(), data = d))
cat("||:\n"); show(frm(bf(y ~ x + (1 + x || g)), family = gaussian(), data = d))
cat("gr by:\n"); show(frm(bf(y ~ x + (1 | gr(g, by = f))), family = gaussian(), data = d))
cat("mm:\n"); show(frm(bf(y ~ x + (1 | mm(g, g2))), family = gaussian(), data = d))
cat("g:h:\n"); show(frm(bf(y ~ x + (1 | g:h)), family = gaussian(), data = d))
