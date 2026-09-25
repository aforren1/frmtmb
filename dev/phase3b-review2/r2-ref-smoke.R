source("dev/phase3b-review2/r2-prelude.R"); r2_lib("lane")
suppressPackageStartupMessages(library(Rmpfr))
x <- mpfr(-40, 200); print(pnorm(x)); print(log(pnorm(mpfr(-1e4, 200))))
