source("dev/refits-lib.R")
set.seed(77)
n <- 300
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(d2$fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n); u2 <- runif(n)
d2$yo <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
d2$yo2 <- ifelse(u2 < p1, 1L, ifelse(u2 < p1 + p2, 2L, 3L))
mv2 <- suppressWarnings(
  frm(bf(yo ~ x) + bf(yo2 ~ z), family = sratio(), data = d2))
nd <- data.frame(x = 0, z = 0, fc = factor("a", levels = c("a", "b", "c")))
options(error = NULL)
r <- tryCatch(withCallingHandlers(fitted(mv2, newdata = nd),
  error = function(e) print(sys.calls())), error = function(e) e)
print(r)
