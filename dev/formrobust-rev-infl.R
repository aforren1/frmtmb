# Reviewer: why influence() differs between weights(wt * k) and
# weights(wk). Seed 101 data as formrobust-rev-aterm2.R.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(101)
n <- 60
d <- data.frame(x = rnorm(n), t = runif(n, 0.5, 2), c = rpois(n, 4) + 1L,
                wt = runif(n, 0.5, 2), s = runif(n, 0.2, 0.6),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$y <- rnorm(n, 0.5 * d$x)
k <- 3
d$wk <- d$wt * 3
fa <- frm(bf(y | weights(wt * k) ~ x + (1 | f)), data = d)
fb <- frm(bf(y | weights(wk) ~ x + (1 | f)), data = d)
ia <- influence(fa, groups = "f")
ib <- influence(fb, groups = "f")
str(unclass(ia), max.level = 1)
for (nm in names(unclass(ia))) {
  cat(nm, identical(unclass(ia)[[nm]], unclass(ib)[[nm]]), "\n")
}
print(all.equal(ia, ib))
