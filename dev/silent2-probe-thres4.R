# Lane silent2: the thres(4) fit of test-ordinal-mixture.R that the new
# threshold-SE test warns on. Which threshold, and is it placed?
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  g <- factor(sample(c("a", "b"), n, TRUE))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
  data.frame(y, x, z, g)
}
src <- readLines("tests/testthat/test-ordinal-mixture.R")
i <- grep("^omx_data <- function", src)
j <- grep("^}", src)
j <- j[j > i][1]
eval(parse(text = src[i:j]))
d <- omx_data(20261016)
print(table(d$y))
f <- suppressWarnings(frm(bf(y | thres(4) ~ x),
                          family = mixture(cumulative(), sratio()), data = d))
print(suppressWarnings(fixef(f)))
print(frmtmb:::mixture_ord_degeneracy(f, "y"))
print(frmtmb:::mixture_ord_thres_se(f, "y"))
sr <- f$spec$responses$y$family$mix$ord$views[[2]]$thres$unident
cat("unident of comp 2:", sr, "\n")
p <- f$opt$par
print(round(p, 3))
