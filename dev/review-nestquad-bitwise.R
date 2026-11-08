# Attack 7: one-block quadrature and Laplace fits, bitwise against
# rellib-r18. Run once per arm (ARM=lane|base) pinned to one core type;
# writes dev/review-nestquad-out/bitwise-<arm>-<tag>.rds.
#   ARM=lane Rscript dev/review-nestquad-bitwise.R <tag>
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
tag <- commandArgs(TRUE)[1] %||% "1"
set.seed(77)
n <- 240
d <- data.frame(g = factor(rep(1:24, each = 10)), h = factor(rep(1:8, 30)),
                x = rnorm(n))
d$c <- factor(paste(d$g, rep(1:2, each = 5)))
d$f <- factor(c("a", "b", "c")[(as.integer(d$g) %% 3) + 1])
eta <- 0.2 + 0.5 * d$x + rnorm(24, 0, 0.7)[d$g]
d$yb <- rbinom(n, 1, plogis(eta))
d$yp <- rpois(n, exp(eta))
d$nt <- sample(1:8, n, TRUE); d$yt <- rbinom(n, d$nt, plogis(eta))
d$yg <- eta + rnorm(n)
d$ynb <- rnbinom(n, size = 2, mu = exp(eta))
d$w <- runif(n, 0.5, 1.5)
d$yo <- factor(cut(eta + rlogis(n), c(-Inf, -0.5, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
mods <- list(
  q_bern = list(bf(yb ~ x + (1 | g)) + bernoulli(), TRUE),
  q_pois = list(bf(yp ~ x + (1 | g)) + poisson(), TRUE),
  q_binom = list(bf(yt | trials(nt) ~ x + (1 | g)) + binomial(), TRUE),
  q_gaus = list(bf(yg ~ x + (1 | g)) + gaussian(), TRUE),
  q_bysplit = list(bf(yb ~ x + (1 | gr(g, by = f))) + bernoulli(), TRUE),
  q_nbw = list(bf(ynb | weights(w) ~ x + (1 | g)) + negbinomial(), TRUE),
  l_nested = list(bf(yb ~ x + (1 | g/c)) + bernoulli(), FALSE),
  l_crossed = list(bf(yp ~ x + (1 | g) + (1 | h)) + poisson(), FALSE),
  l_slope = list(bf(yg ~ x + (1 + x | g)) + gaussian(), FALSE),
  l_nb = list(bf(ynb ~ x + (1 | g)) + negbinomial(), FALSE),
  l_cum = list(bf(yo ~ x + (1 | g)) + cumulative(), FALSE),
  l_lm = list(bf(yg ~ x) + gaussian(), FALSE))
out <- lapply(names(mods), function(nm) {
  m <- mods[[nm]]
  f <- try_msg(frm(m[[1]], data = d, quadrature = m[[2]]))
  if (inherits(f, "error")) return(list(err = conditionMessage(f)))
  list(par = f$opt$par, obj = f$opt$objective, ll = as.numeric(logLik(f)),
       vc = tryCatch(as.matrix(vcov(f)), error = function(e) NULL),
       b = f$estimates$b, warns = attr(f, "warns"))
})
names(out) <- names(mods)
saveRDS(out, sprintf("%s/dev/review-nestquad-out/bitwise-%s-%s.rds", wt, arm, tag))
cat("DONE", arm, tag, "\n")
