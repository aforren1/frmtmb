## brms 2.23.0's default_prior() rows of class lncor for two
## logistic_normal() responses, Y then V, with the data of
## tests/testthat/test-simplex-mvbf.R (smv_data(2)): the response
## order of the rows, which that test pins.
##   Rscript dev/brmsnames-1024-lncor.R > dev/brmsnames-1024-log/lncor.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- readLines("tests/testthat/test-simplex-mvbf.R")
i <- grep("^smv_data <- function", src)
j <- i + which(src[(i + 1):length(src)] == "}")[1L]
eval(parse(text = src[i:j]))
dd <- smv_data(2)
p <- suppressMessages(brms::default_prior(
  brms::mvbf(brms::bf(Y ~ x, family = brms::logistic_normal()),
             brms::bf(V ~ x, family = brms::logistic_normal())),
  data = dd))
p <- as.data.frame(p)
print(p[p$class == "lncor", c("class", "resp")])
