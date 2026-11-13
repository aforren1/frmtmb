# Review of lane nlstop: the cost of nl_limit_message() on its own, against
# the fit it runs on (min of 5 calls each), and the 450-row scan's bodies
# at seed 101 with the check on and off, interleaved, in one process.
# Usage: Rscript dev/reviews/nlstop-review-nlcost.R
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
lim <- ns$nl_limit_message
tm <- function(expr, k = 5) {
  min(vapply(seq_len(k), function(i) {
    t0 <- proc.time()[["elapsed"]]; force(expr()); proc.time()[["elapsed"]] - t0
  }, 0))
}
rows <- list()
set.seed(1); n <- 2000
d <- data.frame(x = runif(n, 0, 10)); d$y <- 10 * d$x / (1 + d$x) + rnorm(n, 0, 0.5)
rows$mm2000 <- list(f = bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE), d = d, fam = gaussian())
set.seed(2); n <- 2000
d <- data.frame(x = runif(n, 0, 5), g = factor(sample(40, n, TRUE)))
d$y <- (4 + rnorm(40, 0, 0.5)[d$g]) * exp(-0.5 * d$x) + rnorm(n, 0, 0.3)
rows$decay_re <- list(f = bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
                      d = d, fam = gaussian(),
                      pr = c(set_prior("normal(4, 3)", nlpar = "a"),
                             set_prior("normal(0.5, 1)", nlpar = "b")))
set.seed(3); n <- 3000
d <- data.frame(x = runif(n, 0, 10), g = factor(sample(60, n, TRUE)))
d$y <- rpois(n, exp(rnorm(60, 0, 0.3)[d$g]) * 8 * d$x / (1.5 + d$x))
rows$poismm_re <- list(f = bf(y ~ log(Vm * x / (K + x)) + u, Vm + K ~ 1, u ~ 0 + (1 | g), nl = TRUE),
                       d = d, fam = poisson(),
                       pr = c(set_prior("normal(8, 5)", nlpar = "Vm"),
                              set_prior("normal(1.5, 2)", nlpar = "K")))
cat("ROWHEAD row t_fit t_check ratio null\n")
for (rn in names(rows)) {
  rw <- rows[[rn]]
  r <- tt(frm(rw$f, family = rw$fam, data = rw$d, prior = rw$pr))
  if (iserr(r)) { cat("ROW", rn, "ERR", r$msg, "\n"); next }
  tf <- tm(function() suppressWarnings(suppressMessages(
    frm(rw$f, family = rw$fam, data = rw$d, prior = rw$pr))), k = 3)
  tc <- tm(function() lim(r))
  cat(sprintf("ROW %s %.3f %.4f %.4f %s warn %s\n", rn, tf, tc, tc / tf,
              is.null(lim(r)), wsum(r)))
  flush.console()
}
cat("DONE\n")
