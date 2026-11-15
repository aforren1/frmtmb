# Lane wcheck2, backlog 4985: an invalid prior on a weighted model with
# an observation-level effect. The data of
# dev/reviews/wcheck-short-attacks.R part couple (seed 20261111, n = 80,
# weights 1 + 8 / n): which error frm() gives first.
#   Rscript dev/wcheck2-4985.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)))
d$y <- 0.5 + 0.4 * d$x + rnorm(n, 0, ifelse(d$fac == "a", 0.6, 1.2))
d$w <- 1 + 8 / n
PS <- function(...) set_prior(..., dpar = "sigma")
fx <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + x)
first <- function(lab, expr) {
  m <- tryCatch(q(expr), error = function(e) {
    paste("ERROR", substr(gsub("\\s+", " ", conditionMessage(e)), 1, 140))
  })
  if (!is.character(m)) m <- "fit"
  cat(sprintf("%-58s %s\n", lab, m))
}
first("normal(0, 1) class b coef Intercept on sigma (invalid)",
      frm(fx, family = gaussian(), data = d,
          prior = PS("normal(0, 1)", class = "b", coef = "Intercept")))
first("the same, unweighted (no refusal to meet)",
      frm(bf(y ~ x + (1 | obs), sigma ~ 1 + x), family = gaussian(),
          data = d, prior = PS("normal(0, 1)", class = "b",
                               coef = "Intercept")))
first("a coef the model does not have (invalid)",
      frm(fx, family = gaussian(), data = d,
          prior = PS("normal(0, 1)", class = "b", coef = "z")))
# the absent case: a valid prior that does not bound the model
first("student_t(3, 0, 2.5) class Intercept on sigma (valid)",
      frm(fx, family = gaussian(), data = d,
          prior = PS("student_t(3, 0, 2.5)", class = "Intercept")))
first("normal(0, 5) class b on mu (valid)",
      frm(fx, family = gaussian(), data = d,
          prior = set_prior("normal(0, 5)", class = "b")))
first("no prior",
      frm(fx, family = gaussian(), data = d))
cat("DONE\n")
