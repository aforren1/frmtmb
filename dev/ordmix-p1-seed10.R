# Punch round 1, B2: what the reviewer's label flags on the identified
# designs' seed 10 (reach 18.7 to 41, no collapsed gap), and on
# rev_cum_cum_500 seed 30.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cut4 <- function(lat) 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
set.seed(10); n <- 500
x <- rnorm(n); z <- rnorm(n); cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
d <- data.frame(x, y = cut4(lat))
print(table(d$y))
f <- frm(bf(y ~ x), family = mixture(cumulative(), cumulative()), data = d)
print(fixef(f)); print(f$opt$par); print(frmtmb:::mixture_ord_degeneracy(f, "y"))
print(eigen(f$obj$he(f$opt$par))$values)
