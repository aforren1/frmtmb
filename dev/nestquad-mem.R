# The memory guard's prediction against the measured peak, per family.
# Lane copy of dev/review-nestquad-mem.R (same data, seed 20261108, 20
# rows per school in 4 classes of 5), with more families and the tape
# bytes the guard predicts before its family factor.
#   Rscript dev/nestquad-mem.R <family> <rows> [ko ki]
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
a <- commandArgs(TRUE)
fam <- a[1]
n <- as.integer(a[2] %||% "3000")
nodes <- if (length(a) >= 4) as.integer(a[3:4]) else c(11L, 15L)
set.seed(20261108)
n_out <- n %/% 20
d <- data.frame(school = factor(rep(seq_len(n_out), each = 20)),
                cls = rep(rep(1:4, each = 5), n_out), x = rnorm(n))
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
eta <- 0.3 + 0.4 * d$x + rnorm(n_out, 0, 0.5)[d$school] +
  rnorm(nlevels(d$cls_uni), 0, 0.5)[d$cls_uni]
re <- "x + (1 | school/cls_uni)"
lhs <- "y"
spec <- switch(fam,
  negbinomial = { d$y <- rnbinom(n, size = 3, mu = exp(eta)); negbinomial() },
  Gamma = { d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta)); Gamma(link = "log") },
  student = { d$y <- eta + rt(n, 4); student() },
  cumulative = { d$y <- factor(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                                   labels = FALSE), ordered = TRUE); cumulative() },
  hurdle_poisson = { d$y <- ifelse(runif(n) < 0.3, 0, rpois(n, exp(eta)) + 1)
                     hurdle_poisson() },
  zi_negbinomial = { d$y <- ifelse(runif(n) < 0.2, 0, rnbinom(n, 3, mu = exp(eta)))
                     zero_inflated_negbinomial() },
  skew_normal = { d$y <- eta + rnorm(n); skew_normal() },
  exgaussian = { d$y <- eta + rnorm(n) + rexp(n, 2); exgaussian() },
  weibull = { d$y <- rweibull(n, 2, exp(eta)); weibull() },
  cens_lognormal = { d$y <- exp(eta + rnorm(n, 0, 0.5))
                     d$ce <- ifelse(d$y > 3, "right", "none"); d$y <- pmin(d$y, 3)
                     lhs <- "y | cens(ce)"; lognormal() },
  cens_beta = { p <- plogis(eta); d$y <- rbeta(n, p * 8, (1 - p) * 8)
                d$ce <- ifelse(d$y > 0.8, "right", "none"); d$y <- pmin(d$y, 0.8)
                lhs <- "y | cens(ce)"; Beta() },
  cens_poisson = { d$y <- rpois(n, exp(eta))
                   d$ce <- ifelse(d$y >= 4, "right", "none")
                   d$y <- pmin(d$y, 4); lhs <- "y | cens(ce)"; poisson() },
  binomial = { d$nt <- 10; d$y <- rbinom(n, 10, plogis(eta))
               lhs <- "y | trials(nt)"; binomial() },
  lognormal = { d$y <- exp(eta + rnorm(n, 0, 0.5)); lognormal() },
  von_mises = { mu <- 2 * atan(0.5 * eta)
                d$y <- ((mu + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi
                von_mises() },
  beta_binomial = { d$nt <- 10; p <- plogis(eta)
                    d$y <- rbinom(n, 10, rbeta(n, p * 8, (1 - p) * 8))
                    lhs <- "y | trials(nt)"; beta_binomial() },
  bernoulli = { d$y <- rbinom(n, 1, plogis(eta)); bernoulli() },
  poisson = { d$y <- rpois(n, exp(eta)); poisson() },
  gaussian = { d$y <- eta + rnorm(n, 0, 0.6); gaussian() },
  beta = { p <- plogis(eta); d$y <- rbeta(n, p * 8, (1 - p) * 8); Beta() })
fo <- bf(as.formula(paste(lhs, "~", re))) + spec
# Laplace first, so the fixed part (package, data, Laplace tape) is in the
# baseline the way it is inside frm() before the nested route starts
fl <- try_msg(frm(fo, data = d))
gc(); gc()
base <- cur_mb()
peak0 <- peak_mb()
msgs <- character(0)
t0 <- proc.time()[[3]]
fq <- withCallingHandlers(
  try_msg(frm(fo, data = d, quadrature = TRUE,
              control = frmtmb_control(verbose = 1, quad_max_gb = 60,
                                       quad_nodes = nodes))),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage") })
el <- proc.time()[[3]] - t0
peak1 <- peak_mb()
pm <- grep("MB predicted", msgs, value = TRUE)[1]
pred_mb <- as.numeric(sub(".*about ([0-9.e+]+) MB predicted.*", "\\1", pm))
tape_mb <- as.numeric(sub(".*[(]tape ([0-9.e+]+) MB.*", "\\1", pm))
key <- sub(".*family ([^)]+)[)].*", "\\1", pm)
vals <- NA
if (!inherits(fq, "error")) {
  ts <- ns$quad_tape_size(fq$obj)
  vals <- ts[2] / (nrow(d) * prod(nodes))
}
delta <- peak1 - base
cat(sprintf(paste("MEM %-15s key %-18s n %5d nodes %dx%d time %6.1fs",
                  "delta %6.0f MB tape %7.2f MB factor_needed %5.2f",
                  "perrow_kB %6.2f pred %6.0f MB pred/delta %.2f",
                  "values/row %.1f %s\n"),
            fam, key, n, nodes[1], nodes[2], el, delta, tape_mb,
            (delta - 128) / tape_mb, (delta - 128) * 1024 / (n * prod(nodes)),
            pred_mb, pred_mb / delta, vals,
            if (inherits(fq, "error")) paste("ERROR", conditionMessage(fq))
            else ""))
