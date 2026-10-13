# Lane optrap: the one-level group with a bound (test-se-check.R m3,
# seed 101) straight through nlminb: which of the start, the scale and
# the bound makes PORT zigzag to its evaluation limit.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(101)
d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
u <- suppressWarnings(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
  data = d3, prior = set_prior("", class = "b", ub = 0.1),
  dry_run = "objective"))
obj <- u$obj
p0 <- obj$par
cat("start", p0, "\n")
up <- c(Inf, 0.1, Inf, Inf)
for (th in c(0, log(2))) for (sc in list(c(1, 1, 1, 1), c(0.5, 0.5, 1, 1))) {
  p <- p0
  p[4] <- th
  r <- nlminb(p, obj$fn, obj$gr, upper = up, scale = sc,
              control = list(eval.max = 1000, iter.max = 1000))
  cat(sprintf("theta0 %.3f scale %s: code %d evals %d obj %.6f theta %.3f | %s\n",
              th, paste(sc, collapse = ","), r$convergence,
              r$evaluations[1], r$objective, r$par[4], r$message))
}
