source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
set.seed(4)
n <- 120
d <- data.frame(x = rnorm(n))
d$y1 <- 1 + 0.5 * d$x + rnorm(n)
d$o <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = c("a","b","c"), ordered_result = FALSE)
f <- frm(mvbf(bf(y1 ~ x) + gaussian(), bf(o ~ x) + categorical()), data = d)
r <- tryCatch(fitted(f), error = function(e) conditionMessage(e)); print(r)
print(dim(fitted(f, resp = "o")))
