# Does frmtmb fit brms's nearest spelling to a multivariate normal
# mixture, bf(mvbind(y1, y2) ~ 1) with a mixture family (one class
# variable per response)? The answer goes in ?mixture_mvn.
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 300
d <- data.frame(y1 = c(rnorm(n / 2), rnorm(n / 2, 4)),
                y2 = c(rnorm(n / 2, 1), rnorm(n / 2, -3)))
r <- tryCatch({
  f <- frm(bf(mvbind(y1, y2) ~ 1) + set_rescor(FALSE),
           family = mixture(gaussian(), gaussian()), data = d)
  print(logLik(f))
  f1 <- frm(bf(y1 ~ 1), family = mixture(gaussian(), gaussian()), data = d)
  f2 <- frm(bf(y2 ~ 1), family = mixture(gaussian(), gaussian()), data = d)
  cat("sum of the two univariate fits:",
      as.numeric(logLik(f1)) + as.numeric(logLik(f2)), "\n")
  "OK"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, "\n")
