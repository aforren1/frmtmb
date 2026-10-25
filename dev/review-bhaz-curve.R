# Reviewer: frm_curve() on the two spellings, each after the same seed.
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.survival))
set.seed(2); n <- 300; x <- rnorm(n)
tt <- stats::rexp(n, exp(-0.5 + 0.7 * x)); ct <- stats::rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
f6b <- frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(), data = dd)
f6c <- frm(bf(time | cens(cens) ~ x), family = cox(df = 6), data = dd)
g <- data.frame(x = c(-1, 0, 1))
set.seed(9); a <- frm_curve(f6b, g); set.seed(9); b <- frm_curve(f6c, g)
k <- setdiff(names(a), character())
cat("columns identical:", all(vapply(k, function(j) identical(a[[j]], b[[j]]), NA)), "\n")
