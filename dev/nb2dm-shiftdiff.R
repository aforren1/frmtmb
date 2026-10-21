# Would lgamma_shift_diff() make the dirichlet-multinomial density
# accurate at large phi? Its accuracy is documented for shifts up to
# 100; the counts here reach 2000. Against the exact rising-factorial
# sum, on doubles.
#   Rscript dev/nb2dm-shiftdiff.R > dev/nb2dm-out/shiftdiff.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
sd_ <- frmtmb:::lgamma_shift_diff
exact <- function(a, s) if (s == 0) 0 else sum(log(a + seq_len(s) - 1))
g <- expand.grid(a = 10^c(-6, -3, -1, 0, 1, 3, 6, 9, 12, 15),
                 s = c(0, 1, 3, 17, 99, 400, 2000))
g$exact <- mapply(exact, g$a, g$s)
g$shift <- mapply(sd_, g$a, g$s)
g$naive <- lgamma(g$a + g$s) - lgamma(g$a)
rel <- function(v) ifelse(g$exact == 0, abs(v), abs(v - g$exact) /
                            abs(g$exact))
g$rel_shift <- rel(g$shift)
g$rel_naive <- rel(g$naive)
print(aggregate(cbind(rel_shift, rel_naive) ~ a, data = g, FUN = max),
      digits = 3)
print(aggregate(cbind(rel_shift, rel_naive) ~ s, data = g, FUN = max),
      digits = 3)
