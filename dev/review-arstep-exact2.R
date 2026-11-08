# Review of lane arstep, short check of the fix round: the AR marginal on
# numbers (projection route, direct route, QR fallback) and the constant-sigma
# marginal against fit$obj$fn, and escaped warnings.
# Usage: Rscript review-arstep-exact2.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
rv_libs(commandArgs(TRUE)[1])
ns <- asNamespace("frmtmb")
# count which route the AR marginal takes
CNT <- new.env(); CNT$qr <- 0L; CNT$proj <- 0L
o_qr <- ns$gp_ar_proj_qr
o_pj <- ns$gp_ar_proj
rv_put("gp_ar_proj_qr", function(s, phi) { CNT$qr <- CNT$qr + 1L; o_qr(s, phi) })
rv_put("gp_ar_proj", function(s, phi) { CNT$proj <- CNT$proj + 1L; o_pj(s, phi) })
marg_at <- function(s, pl) {
  if (!is.null(s$ar)) return(ns$gp_marginal_nll(s)(pl))
  nll <- ns$gp_marginal_nll(s)
  ob <- RTMB::MakeADFun(nll, pl[intersect(names(pl), c("beta", "betad", "theta"))],
                        silent = TRUE)
  ob$fn(ob$par)
}
check <- function(lbl, f, d, phis = list(NULL), npt = 4, ...) {
  r <- rv_quiet(frm(f, data = d, family = gaussian(), ...))
  fit <- r$value
  if (inherits(fit, "error")) {
    cat(sprintf("EX2 %-34s fit error: %s\n", lbl, conditionMessage(fit))); return(invisible())
  }
  s <- ns$gp_marginal_setup(fit$frame)
  used <- !is.null(fit$obj$env$.frmtmb_gp_marginal)
  if (is.null(s)) {
    cat(sprintf("EX2 %-34s gate REFUSED, start used %s, warnings %d, ll %.6f\n",
                lbl, used, length(r$warnings), as.numeric(logLik(fit))))
    return(invisible())
  }
  p0 <- fit$opt$par
  iac <- which(names(p0) == "thetaac")
  worst <- 0
  CNT$qr <- 0L
  for (k in seq_len(max(npt, length(phis)))) {
    p <- p0
    if (k > 1L) p <- p0 + 0.15 * cos(seq_along(p0) * k)
    ph <- phis[[min(k, length(phis))]]
    if (length(iac) && !is.null(ph)) p[iac] <- ph
    v <- fit$obj$fn(p)
    pl <- fit$obj$env$parList(p)
    pl <- pl[intersect(names(pl), c("beta", "betad", "theta", "thetaac"))]
    m <- marg_at(s, pl)
    rel <- abs(m - v) / abs(v)
    worst <- max(worst, rel)
    cat(sprintf("EX2P %-34s pt %d ar %-16s fn %.10f marg %.10f rel %.1e\n", lbl, k,
                if (length(iac)) paste(format(pl$thetaac, digits = 4), collapse = ",") else "-",
                v, m, rel))
  }
  cat(sprintf("EX2 %-34s admitted u %d n %d direct %s used %s warnings %d qr-calls %d worst rel %.1e\n",
              lbl, s$u, s$n, if (is.null(s$ar)) "-" else s$ar$direct, used,
              length(r$warnings), CNT$qr, worst))
}
panel <- function(S, T, seed, sizes = NULL, phi = 0.5, re = FALSE) {
  set.seed(seed)
  sizes <- sizes %||% rep(T, S)
  d <- do.call(rbind, lapply(seq_along(sizes), function(i) {
    data.frame(t = seq_len(sizes[i]), s = i)
  }))
  d$s <- factor(d$s)
  d$x <- d$t / 5
  d$y <- sin(d$x) + (if (re) rnorm(length(sizes), 0, 0.4)[d$s] else 0) +
    unlist(lapply(sizes, function(m) as.numeric(arima.sim(list(ar = phi), m, sd = 0.2))))
  d[sample.int(nrow(d)), ]
}
ph1 <- list(NULL, 0.3, 0.95, 1.05, 1.3, 1.8, -0.9)
check("panel 12x50", y ~ gp(x) + ar(time = t, gr = s), panel(12, 50, 1), ph1)
check("panel unequal 7,20,33 + (1|s)", y ~ gp(x) + (1 | s) + ar(time = t, gr = s),
      panel(3, 0, 2, sizes = c(7, 20, 33), re = TRUE), ph1)
check("panel p2 8x30", y ~ gp(x) + ar(time = t, gr = s, p = 2), panel(8, 30, 3),
      list(NULL, c(0.5, 0.3), c(1.2, -0.1), c(1.5, 0.2)))
set.seed(4); n <- 60; x <- sort(runif(n, 0, 10))
d1 <- data.frame(x = x, t = seq_len(n), y = sin(x) + 2 * (x > 5) +
                   as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.2)))
check("single n60 direct", y ~ gp(x) + ar(time = t), d1, ph1)
check("single n60 shuffled", y ~ gp(x) + ar(time = t), d1[sample.int(n), ], ph1)
# the lane's ninth row, step 160 seed 2 (ga_step / dev/arstep-seedscan.R)
n <- 160; set.seed(74000 + 2 * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
d2 <- data.frame(x = x, t = seq_len(n),
                 y = 2 * (x > 5) - 1 + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03)))
check("step 160 seed 2", y ~ gp(x) + ar(time = t), d2,
      list(NULL, 0.5, 1.0, 1.02, 1.05, 1.1, 0.9), npt = 7)
cat(sprintf("EX2 step 160 seed 2 fit ll %.4f\n",
            as.numeric(logLik(frm(y ~ gp(x) + ar(time = t), data = d2)))))
# constant sigma, with and without se()
set.seed(12); n <- 50
d3 <- data.frame(x = round(runif(n, 0, 10), 1), t = seq_len(n), se = runif(n, 0.05, 0.3))
d3$y <- sin(d3$x) + 2 * (d3$x > 5) + as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.2))
check("sigma = 0.3", bf(y ~ gp(x), sigma = 0.3), d3)
check("se() alone", y | se(se) ~ gp(x), d3)
check("se(sigma = TRUE)", y | se(se, sigma = TRUE) ~ gp(x), d3)
check("se(sigma = TRUE), sigma = 0.3", bf(y | se(se, sigma = TRUE) ~ gp(x), sigma = 0.3), d3)
check("se() alone, sigma = 0.3", bf(y | se(se) ~ gp(x), sigma = 0.3), d3)
check("ar, sigma = 0.3", bf(y ~ gp(x) + ar(time = t), sigma = 0.3), d3)
check("ar, sigma = 0.2", bf(y ~ gp(x) + ar(time = t), sigma = 0.2), d3)
cat("DONE\n")
