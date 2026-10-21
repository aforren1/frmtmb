# Fix round 2: what a block carries to name its term as written.
.libPaths(c("C:/Users/adf44/source/r/wt-fuzzfix-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1); n <- 60
d <- data.frame(ga = factor(rep(1:12, each = 5)), y = rnorm(n), x = rnorm(n), tt = factor(rep(1:5, 12)))
fr <- frm(bf(y ~ x + (1 | ga) + ar1(tt + 0 | ga), sigma ~ (1 | ga)), data = d, dry_run = "frame")
for (bk in fr$re_blocks) { cat("---", bk$term_label, bk$dpar, "\n"); str(bk$components[[1]][c("written_text", "label", "dpar", "lp_key")]) ; print(frmtmb:::block_written_label(bk)) }
