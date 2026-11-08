# Review of lane arstep: is the AR marginal exact on every structure the gate
# admits, and does the gate refuse what it should?
# Three readings per point: fit$obj$fn (the tape, base code), the lane's
# gp_marginal_nll(), and rv_ref() below (dense, written here from the model
# statement, no frmtmb code). Rows are SHUFFLED in most cases so data order
# is not time order.
# Usage: Rscript review-arstep-exact.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")

# A: within each group, rows in time order; lag k is the row k places
# earlier in that order (brms ignores the spacing of time). e = A (y - mu).
rv_A <- function(t, g, ar) {
  n <- length(t)
  A <- diag(n)
  for (grp in unique(g)) {
    r <- which(g == grp)
    r <- r[order(t[r])]
    for (k in seq_along(ar)) {
      for (i in seq_along(r)) if (i > k) A[r[i], r[i - k]] <- A[r[i], r[i - k]] - ar[k]
    }
  }
  A
}
rv_ker <- function(x, sdgp, lscale) {
  sdgp^2 * (exp(-outer(x, x, "-")^2 / (2 * lscale^2)) +
              1e-12 * outer(x, x, "=="))
}
# nll of y given everything; Vf the covariance of the random part
rv_ref <- function(y, mu, Vf, sigma, A) {
  n <- length(y)
  C <- A %*% Vf %*% t(A) + sigma^2 * diag(n)
  R <- chol(C)
  e <- backsolve(R, as.vector(A %*% (y - mu)), transpose = TRUE)
  sum(log(diag(R))) + 0.5 * sum(e^2) + 0.5 * n * log(2 * pi)
}

mk <- function(n, seed, p = 1L, groups = NULL, gaps = FALSE, shuffle = TRUE,
               re = FALSE, off = FALSE) {
  set.seed(seed)
  sizes <- if (is.null(groups)) n else groups
  g <- rep(seq_along(sizes), sizes)
  t <- unlist(lapply(sizes, function(m) {
    if (gaps) sort(sample.int(3L * m, m)) else seq_len(m)
  }))
  x <- round(runif(length(g), 0, 10), 3)
  arc <- c(0.5, -0.2, 0.1)[seq_len(p)]
  e <- unlist(lapply(sizes, function(m) {
    as.numeric(stats::arima.sim(list(ar = arc), m, sd = 0.2))
  }))
  gr2 <- factor(sample(letters[1:5], length(g), TRUE))
  ug <- rnorm(5, 0, 0.5)
  o <- if (off) runif(length(g), -0.5, 0.5) else 0
  y <- sin(x) + 2 * (x > 5) + e + o + if (re) ug[as.integer(gr2)] else 0
  d <- data.frame(y = y, x = x, t = t, g = factor(g), h = gr2, o = o,
                  z = rnorm(length(g)))
  if (shuffle) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

check <- function(lbl, f, d, phis, gcol = NULL, re = FALSE, off = FALSE, ...) {
  fit <- tryCatch(frm(f, data = d, family = gaussian(), ...),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("EX %-28s fit error: %s\n", lbl, conditionMessage(fit)))
    return(invisible())
  }
  s <- ns$gp_marginal_setup(fit$frame)
  used <- !is.null(fit$obj$env$.frmtmb_gp_marginal)
  if (is.null(s)) {
    cat(sprintf("EX %-28s gate REFUSED, marginal start used %s\n", lbl, used))
    return(invisible())
  }
  nll <- ns$gp_marginal_nll(s)
  p0 <- fit$opt$par
  iac <- which(names(p0) == "thetaac")
  g <- if (is.null(gcol)) rep(1L, nrow(d)) else d[[gcol]]
  worst_m <- 0; worst_r <- 0
  for (k in seq_along(phis)) {
    p <- p0
    if (k > 1L) p <- p0 + 0.15 * cos(seq_along(p0) * k)
    if (!is.null(phis[[k]])) p[iac] <- phis[[k]]
    v <- fit$obj$fn(p)
    pl <- fit$obj$env$parList(p)
    m <- nll(list(beta = pl$beta, betad = pl$betad, theta = pl$theta,
                  thetaac = pl$thetaac))
    # reference: the gp block and, if present, the (1|h) block
    bks <- fit$frame$re_blocks
    Vf <- 0
    for (bk in bks) {
      th <- pl$theta[bk$theta_idx]
      if (identical(bk$covstruct, "gp")) {
        Vf <- Vf + rv_ker(d$x, exp(th[1]), exp(th[2]))
      } else {
        Zh <- outer(as.integer(d$h), seq_len(nlevels(d$h)), "==") * 1
        Vf <- Vf + exp(2 * th[1]) * tcrossprod(Zh)
      }
    }
    mu <- pl$beta[1] + if (off) d$o else 0
    r <- rv_ref(d$y, mu, Vf, exp(pl$betad[1]), rv_A(d$t, g, pl$thetaac))
    rel_m <- abs(m - v) / abs(v); rel_r <- abs(r - v) / abs(v)
    worst_m <- max(worst_m, rel_m); worst_r <- max(worst_r, rel_r)
    cat(sprintf("EXP %-28s pt %d ar %-22s fn %.10f marg %.10f ref %.10f relm %.1e relr %.1e\n",
                lbl, k, paste(format(pl$thetaac, digits = 4), collapse = ","),
                v, m, r, rel_m, rel_r))
  }
  cat(sprintf("EX %-28s admitted, used %s, opt ar %s, worst rel marg-vs-fn %.1e ref-vs-fn %.1e\n",
              lbl, used, paste(format(p0[iac], digits = 4), collapse = ","),
              worst_m, worst_r))
}

ph1 <- list(NULL, 0.3, 0.95, 1.0, 1.05, 1.3, -0.9)
check("p1 n60 shuffled", y ~ gp(x) + ar(time = t), mk(60, 1), ph1)
check("p1 n60 gaps", y ~ gp(x) + ar(time = t), mk(60, 2, gaps = TRUE), ph1)
check("p2 n60", y ~ gp(x) + ar(time = t, p = 2), mk(60, 3, p = 2),
      list(NULL, c(0.5, 0.3), c(1.2, -0.1), c(0.6, 0.5)))
check("p3 n70", y ~ gp(x) + ar(time = t, p = 3), mk(70, 4, p = 3),
      list(NULL, c(0.5, -0.2, 0.1), c(0.9, 0.3, -0.25), c(1.1, 0.1, 0.05)))
check("p1 gr unequal 5,17,38", y ~ gp(x) + ar(time = t, gr = g),
      mk(60, 5, groups = c(5, 17, 38)), ph1, gcol = "g")
check("p2 gr with size-1 group", y ~ gp(x) + ar(time = t, gr = g, p = 2),
      mk(40, 6, p = 2, groups = c(1, 2, 15, 22)),
      list(NULL, c(0.4, 0.2), c(1.05, 0.02)), gcol = "g")
check("p1 gr gaps unequal", y ~ gp(x) + ar(time = t, gr = g),
      mk(50, 7, groups = c(10, 40), gaps = TRUE), ph1, gcol = "g")
check("p1 + (1|h)", y ~ gp(x) + (1 | h) + ar(time = t),
      mk(60, 8, re = TRUE), ph1, re = TRUE)
check("p1 + offset", y ~ gp(x) + offset(o) + ar(time = t),
      mk(60, 9, off = TRUE), ph1, off = TRUE)
check("p1 no shuffle n120", y ~ gp(x) + ar(time = t), mk(120, 10, shuffle = FALSE),
      ph1)
d11 <- mk(60, 11); d11$x <- round(d11$x)  # repeated positions
check("p1 repeated x", y ~ gp(x) + ar(time = t), d11, ph1)

# Gate edge: each should be REFUSED (or, if admitted, exact)
d <- mk(50, 12)
d$se <- runif(nrow(d), 0.1, 0.3); d$w <- runif(nrow(d), 0.5, 2)
check("edge se()", y | se(se, sigma = TRUE) ~ gp(x) + ar(time = t), d, ph1)
check("edge weights()", y | weights(w) ~ gp(x) + ar(time = t), d, ph1)
check("edge cov=TRUE", y ~ gp(x) + ar(time = t, cov = TRUE), d, ph1)
check("edge ma", y ~ gp(x) + ma(time = t), d, ph1)
check("edge arma", y ~ gp(x) + arma(time = t), d, ph1)
check("edge sigma ~ z", bf(y ~ gp(x) + ar(time = t), sigma ~ z), d, ph1)
check("edge sigma ~ 1 explicit", bf(y ~ gp(x) + ar(time = t), sigma ~ 1), d, ph1)
check("edge sigma constant 0.3", bf(y ~ gp(x) + ar(time = t), sigma = 0.3), d,
      ph1)
check("edge log link", y ~ gp(x) + ar(time = t), d, ph1)  # control
fit <- tryCatch(frm(y ~ gp(x) + ar(time = t), data = transform(d, y = exp(y)),
                    family = gaussian(link = "log")), error = function(e) e)
cat("EX edge log link gate:",
    if (inherits(fit, "error")) conditionMessage(fit) else
      is.null(ns$gp_marginal_setup(fit$frame)), "\n")
for (fam in c("student", "skew_normal")) {
  fit <- tryCatch(frm(y ~ gp(x) + ar(time = t), data = d,
                      family = get(fam, asNamespace("frmtmb"))()),
                  error = function(e) e)
  cat("EX edge", fam, "gate refused:",
      if (inherits(fit, "error")) conditionMessage(fit) else
        is.null(ns$gp_marginal_setup(fit$frame)), "\n")
}
# a prior on the AR coefficient: the objective holds a term the marginal
# does not, so the agreement guard must refuse the start
pr <- tryCatch(frm(y ~ gp(x) + ar(time = t), data = d, family = gaussian(),
                   prior = prior(normal(0, 0.5), class = ar)),
               error = function(e) e)
cat("EX edge prior on ar:",
    if (inherits(pr, "error")) conditionMessage(pr) else
      paste("marginal used", !is.null(pr$obj$env$.frmtmb_gp_marginal),
            "ll", format(as.numeric(logLik(pr)), digits = 10)), "\n")
cat("DONE\n")
