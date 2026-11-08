# Attack 5: the memory guard's prediction against the measured peak, on
# families the lane did not calibrate on. One configuration per process.
#   Rscript dev/review-nestquad-mem.R <family> <rows> [ko ki]
#   Rscript dev/review-nestquad-mem.R noinfo
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
a <- commandArgs(TRUE)
fam <- a[1]
n <- as.integer(a[2] %||% "3000")
nodes <- if (length(a) >= 4) as.integer(a[3:4]) else c(11L, 15L)
seed <- 20261108
set.seed(seed)

if (fam == "noinfo") {
  # TMB:::info() gone: in-memory replacement that errors, as when TMB
  # drops or renames it. The guard should fall back, and the fallback
  # should not refuse a small model the lane's tests fit.
  utils::assignInNamespace("info", function(...) stop("no info()"), "TMB")
  set.seed(20261020)
  ng <- 12; nt <- 6; nn <- ng * nt
  d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = nn)),
                  x = rnorm(nn))
  d$y <- rbinom(nn, 1, 0.5)
  for (lim in c(4, 64)) {
    r <- try_msg(frm(bf(y ~ x + (1 | ga/gb)) + bernoulli(), data = d,
                     quadrature = TRUE,
                     control = frmtmb_control(quad_max_gb = lim)))
    cat(sprintf("NOINFO rows %d limit %g: %s\n", nn, lim,
                if (inherits(r, "error")) conditionMessage(r) else
                  sprintf("fit, logLik %.6f", as.numeric(logLik(r)))))
  }
  quit()
}

n_out <- n %/% 20
d <- data.frame(school = factor(rep(seq_len(n_out), each = 20)),
                cls = rep(rep(1:4, each = 5), n_out), x = rnorm(n))
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
eta <- 0.3 + 0.4 * d$x + rnorm(n_out, 0, 0.5)[d$school] +
  rnorm(nlevels(d$cls_uni), 0, 0.5)[d$cls_uni]
re <- "x + (1 | school/cls_uni)"
spec <- switch(fam,
  negbinomial = { d$y <- rnbinom(n, size = 3, mu = exp(eta)); negbinomial() },
  Gamma = { d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta)); Gamma(link = "log") },
  student = { d$y <- eta + rt(n, 4); student() },
  lognormal = { d$y <- exp(eta + rnorm(n, 0, 0.5)); lognormal() },
  von_mises = { d$y <- atan2(sin(eta / 3 + rnorm(n, 0, 0.5)), cos(eta / 3 + rnorm(n, 0, 0.5)))
                von_mises() },
  binomial = { d$nt <- sample(1:10, n, TRUE); d$y <- rbinom(n, d$nt, plogis(eta))
               re <- "x + (1 | school/cls_uni)"; binomial() },
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
                     re <- paste("x + (1 | school/cls_uni)"); lognormal() },
  bernoulli = { d$y <- rbinom(n, 1, plogis(eta)); bernoulli() })
lhs <- if (fam == "cens_lognormal") "y | cens(ce)" else if (fam == "binomial") "y | trials(nt)" else "y"
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
    cat(sprintf("STAGE peak %6.0f cur %6.0f | %s", peak_mb(), cur_mb(),
                conditionMessage(m)))
    flush(stdout())
    invokeRestart("muffleMessage") })
el <- proc.time()[[3]] - t0
peak1 <- peak_mb()
pred <- regmatches(msgs, regexpr("about [0-9.e+]+ MB predicted", msgs))
pred_mb <- as.numeric(sub("about ([0-9.e+]+) MB.*", "\\1", pred))
vals <- NA
if (!inherits(fq, "error")) {
  ts <- ns$quad_tape_size(fq$obj)
  plan <- ns$quad_check_blocks(fq$frame)
  vals <- ts[2] / (nrow(d) * prod(nodes))
}
cat(sprintf("MEM %-15s n %5d nodes %dx%d time %6.1fs base %6.0f MB peak_before %6.0f peak_after %6.0f delta %6.0f MB pred %6.0f MB ratio pred/delta %.2f values/stacked-row %.1f %s\n",
            fam, n, nodes[1], nodes[2], el, base, peak0, peak1,
            peak1 - base, pred_mb, pred_mb / (peak1 - base), vals,
            if (inherits(fq, "error")) paste("ERROR", conditionMessage(fq)) else ""))
if (length(attr(fq, "warns"))) cat("WARN", attr(fq, "warns"), "\n")
