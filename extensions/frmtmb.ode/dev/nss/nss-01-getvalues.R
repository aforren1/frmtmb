source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
library(RTMB)

gv <- RTMB:::getValues
ctx <- RTMB:::ad_context

cat("outside a tape, ad_context():", ctx(), "\n")

f <- MakeTape(function(x) {
  cat("  inside tape, ad_context():", ctx(), "\n")
  v <- try(gv(x), silent = TRUE)
  cat("  getValues(x): ")
  print(v)
  y <- x * 3 + 1
  v2 <- try(gv(y), silent = TRUE)
  cat("  getValues(x*3+1): ")
  print(v2)
  sum(y)
}, c(2, 5))

cat("tape value at (2,5):", f(c(2, 5)), "\n")
cat("tape value at (7,9):", f(c(7, 9)), "\n")
