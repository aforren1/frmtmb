source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# B1: what the default's pre-fit hides or breaks. Each case is fitted
# under the default and under autoscale = FALSE; messages AND warnings
# are collected. Seed 505.
ctlF <- frmtmb_control(autoscale = FALSE)
cond <- function(expr) {
  w <- character(0); m <- character(0)
  f <- tryCatch(withCallingHandlers(expr,
         warning = function(x) { w <<- c(w, conditionMessage(x))
                                 invokeRestart("muffleWarning") },
         message = function(x) { m <<- c(m, conditionMessage(x))
                                 invokeRestart("muffleMessage") }),
       error = function(e) e)
  list(fit = f, warn = w, msg = m)
}
show <- function(lab, r) {
  cat(sprintf("%-38s %s | warnings=%d messages=%d\n", lab,
              if (inherits(r$fit, "error"))
                paste("ERROR:", conditionMessage(r$fit)) else
                sprintf("logLik %.6f code %d", as.numeric(logLik(r$fit)),
                        r$fit$opt$convergence),
              length(r$warn), length(r$msg)))
  for (x in unique(c(r$warn, r$msg))) cat("     |", substr(x, 1, 150), "\n")
}
set.seed(505)
n <- 200
x <- rnorm(n) * 1e-4

cat("== (a) bounded-link init at its bound: all-ones bernoulli ==\n")
d <- data.frame(y = rep(1L, n), x)
show("default", cond(frm(y ~ x, family = bernoulli(), data = d)))
show("FALSE", cond(frm(y ~ x, family = bernoulli(), data = d,
                       control = ctlF)))

cat("== (b) nonlinear start at the prior location ==\n")
t <- runif(n, 0, 5)
d2 <- data.frame(t, x, y = 3 * exp(-0.7 * t) + 2000 * x + rnorm(n, 0, 0.1))
nlf <- bf(y ~ a * exp(-k * t) + c, a ~ 1, k ~ 1, c ~ 0 + x, nl = TRUE)
pr <- set_prior("normal(0.5, 1)", nlpar = "k")
show("default", cond(frm(nlf, data = d2, prior = pr)))
show("FALSE", cond(frm(nlf, data = d2, prior = pr, control = ctlF)))

cat("== (c) a user start with a bound on b_x ==\n")
d3 <- data.frame(x, y = 1 + 5e3 * x + rnorm(n))
pr3 <- set_prior("", class = "b", coef = "x", lb = 100)
st <- list(beta = c(1, 4000))
show("default", cond(frm(y ~ x, data = d3, prior = pr3, start = st)))
show("FALSE", cond(frm(y ~ x, data = d3, prior = pr3, start = st,
                       control = ctlF)))

cat("== (d) sparse_x with a small slope column ==\n")
g <- factor(rep(1:20, each = 10))
d4 <- data.frame(g, x = rnorm(n) * 0.01)
d4$y <- 1 + 50 * d4$x + rnorm(20)[g] + rnorm(20, 0, 30)[g] * d4$x + rnorm(n)
show("default", cond(frm(y ~ x + (1 + x | g), data = d4,
                         control = frmtmb_control(sparse_x = TRUE))))
show("FALSE", cond(frm(y ~ x + (1 + x | g), data = d4,
                       control = frmtmb_control(sparse_x = TRUE,
                                                autoscale = FALSE))))

cat("== (e) a singular random slope (true sd 0), small column ==\n")
d5 <- data.frame(g, x = rnorm(n) * 0.01)
d5$y <- 1 + 50 * d5$x + rnorm(20)[g] + rnorm(n)
show("default", cond(frm(y ~ x + (1 + x | g), data = d5)))
show("FALSE", cond(frm(y ~ x + (1 + x | g), data = d5, control = ctlF)))

cat("== (f) a single-level grouping factor with a small slope ==\n")
d6 <- data.frame(g1 = factor(rep("a", n)), x = rnorm(n) * 0.01)
d6$y <- 1 + 50 * d6$x + rnorm(n)
show("default", cond(frm(y ~ x + (1 + x | g1), data = d6)))
show("FALSE", cond(frm(y ~ x + (1 + x | g1), data = d6, control = ctlF)))

cat("== (g) an observation-level random slope (olre check) ==\n")
d7 <- data.frame(id = factor(seq_len(n)), x = rnorm(n) * 0.01)
d7$y <- rpois(n, exp(1 + 30 * d7$x))
show("default", cond(frm(y ~ x + (1 + x | id), family = poisson(),
                         data = d7)))
show("FALSE", cond(frm(y ~ x + (1 + x | id), family = poisson(), data = d7,
                       control = ctlF)))
