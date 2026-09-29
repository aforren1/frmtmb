.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
set.seed(1)
d <- expand.grid(t = 1:6, g = factor(1:4))
d <- d[-c(3, 11), ]
d$x <- rnorm(nrow(d)); d$y <- rnorm(nrow(d)) + 3
d$w <- runif(nrow(d), .5, 2)
d$cc <- sample(c(0, 0, 0, 1, -1), nrow(d), TRUE)
d$y2 <- rnorm(nrow(d))
try_it <- function(lab, f, ...) {
  r <- tryCatch({ make_stancode(f, data = d, ...); "OK" },
                error = function(e) paste("ERR:", conditionMessage(e)))
  cat(lab, ": ", r, "\n", sep = "")
}
try_it("weights+ar", bf(y | weights(w) ~ x + ar(t, g)), family = gaussian())
try_it("cens+arma", bf(y | cens(cc) ~ x + arma(t, g, p = 1, q = 1)),
       family = gaussian())
try_it("trunc+ma", bf(y | trunc(lb = 0) ~ x + ma(t, g, q = 2)),
       family = gaussian())
try_it("sigma~x+ar2", bf(y ~ x + ar(t, g, p = 2), sigma ~ x),
       family = gaussian())
try_it("student arma22 +(1|g)", bf(y ~ x + arma(t, g, p = 2, q = 2) +
                                     (1 | g)), family = student())
try_it("mv rescor + ar both", bf(y ~ x + ar(t, g)) +
         bf(y2 ~ x + ar(t, g)) + set_rescor(TRUE), family = gaussian())
