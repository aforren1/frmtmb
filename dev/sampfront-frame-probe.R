# Lane sampfront: what a fit's frame holds for the subset blockers.
#   Rscript dev/sampfront-frame-probe.R > dev/sampfront-out/frame-probe.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rnorm(48),
                t = rep(1:8, 6))
d$y <- 1 + 0.5 * d$x + rnorm(48)
d$xo <- d$x + rnorm(48, 0, 0.1)
d$sx <- 0.1
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
f1 <- q(bf(y ~ x + ar(time = t, gr = g)), data = d)
cat("ar: names(frame$autocor) =", names(f1$frame$autocor), "\n")
str(f1$frame$autocor[[1]], max.level = 1)
cat("spec autocor:\n")
str(f1$spec$responses[[1]]$autocor, max.level = 1)
f2 <- q(bf(y ~ me(xo, sx)), data = d)
cat("me: mi_map names", names(f2$frame$mi_map), "\n")
f3 <- q(bf(y ~ x + cosy(gr = g)), data = d)
cat("cosy: autocor", names(f3$frame$autocor), "\n")
str(f3$frame$autocor[[1]], max.level = 1)
f4 <- q(bf(y ~ x), data = d)
cat("plain: autocor", length(f4$frame$autocor), " mi_map",
    length(f4$frame$mi_map), " subset_rows", length(f4$frame$subset_rows),
    "\n")
cat("rownames(f4$data)[1:3]", rownames(f4$data)[1:3], "\n")
x <- poly(d$x, 2)
cat("class poly", class(x), " after [", class(x[1:3, , drop = FALSE]), "\n")
