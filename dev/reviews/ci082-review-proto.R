# Review of lane ci082: a structural "fixed row" test, for comparison with
# the lane's 16-eps rule. A row is fixed when it takes the same value
# (to 1e-10 relative) at the estimate, at the structure's start, and at
# the estimate moved by +/-1 in each theta: a row a structure fixes is
# constant everywhere, an estimated row flat at the estimate is not.
#   Rscript dev/reviews/ci082-review-proto.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")

rows_g <- function(bk, tt) {
  V <- ns$cs_entry(bk[["covstruct"]])$vcov(tt, bk)
  out <- log(sqrt(diag(V)))
  if (ns$block_has_cor(bk)) {
    C <- stats::cov2cor(V)
    c0 <- ns$varcorr_cor_clamp
    out <- c(out, atanh(pmin(pmax(C[lower.tri(C)], -c0), c0)))
  }
  out
}
compare <- function(tag, fit) {
  th <- fit$estimates[["theta"]]
  for (bk in fit$frame[["re_blocks"]]) {
    t0 <- th[bk[["theta_idx"]]]
    g0 <- rows_g(bk, t0)
    # the lane's rule
    J <- vapply(seq_along(t0), function(i) {
      h <- 1e-5 * max(abs(t0[i]), 1)
      tp <- t0; tp[i] <- tp[i] + h
      tm <- t0; tm[i] <- tm[i] - h
      (rows_g(bk, tp) - rows_g(bk, tm)) / (2 * h)
    }, numeric(length(g0)))
    J <- matrix(J, nrow = length(g0))
    lane <- (rowSums(abs(sweep(J, 2, 2e-5 * pmax(abs(t0), 1), `*`)) >
                       16 * .Machine$double.eps * pmax(abs(g0), 1)) == 0)
    # the structural probe
    st <- ns$cs_entry(bk[["covstruct"]])$start(bk)
    pts <- list(st)
    for (i in seq_along(t0)) for (s in c(-1, 1)) {
      tt <- t0; tt[i] <- tt[i] + s; pts[[length(pts) + 1L]] <- tt
    }
    moved <- Reduce(`|`, lapply(pts, function(p) {
      gp <- rows_g(bk, p)
      # relative to the larger of the two values, with no floor: a
      # tiny row that is estimated moves by its own size
      !is.finite(gp) | abs(gp - g0) > 1e-10 * pmax(abs(g0), abs(gp))
    }))
    cat(sprintf("RESULT %s block %s rows %d lane_fixed %d probe_fixed %d disagree %d\n",
                tag, bk[["term_label"]], length(g0), sum(lane),
                sum(!moved), sum(lane != !moved)))
  }
}

# fixed correlations 0.5 and 0.999 (registered), 20 perturbations each
for (rr in c(0.5, 0.999)) {
  xc <- function(dd) matrix(rr, dd, dd) + diag(1 - rr, dd)
  frmtmb_register_covstruct(
    "xcs",
    lpdf = function(b, theta, blk) {
      dd <- blk[["dim"]]
      Sigma <- exp(2 * theta[1]) * xc(dd)
      dim(b) <- c(dd, length(b) %/% dd)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * xc(blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L, replace = TRUE)
  set.seed(11L)
  G <- 30L
  d0 <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:3,
                    g = factor(sprintf("g%02d", seq_len(G))))
  u <- t(crossprod(chol(0.8^2 * xc(3L)), matrix(rnorm(3 * G), 3)))
  d0$x <- rnorm(nrow(d0))
  d0$y <- 0.4 + 0.7 * d0$x + u[cbind(as.integer(d0$g), as.integer(d0$t))] +
    rnorm(nrow(d0), 0, 0.6)
  for (k in 0:19) {
    d <- d0
    d$y <- d0$y * (1 + 2e-15 * k)
    fit <- suppressWarnings(frm(bf(y ~ x + xcs(0 + t | g)),
                                family = gaussian(), data = d))
    compare(sprintf("xcs%g k%d", rr, k), fit)
  }
}

# ou() on data with no correlation across times (natural fit)
set.seed(1)
G <- 50L
times <- c(0, 0.5, 1, 2, 8, 20, 40)
nt <- length(times)
b <- rnorm(G * nt)
d <- data.frame(y = rep(b, each = 2L) + rnorm(2L * G * nt, 0, 0.5),
                g = factor(rep(seq_len(G), each = 2L * nt)),
                tim = num_factor(rep(rep(times, each = 2L), G)))
compare("ouind", frm(y ~ 1 + ou(tim + 0 | g), family = gaussian(), data = d))

# homcs started at -30 on its logistic correlation
set.seed(1)
G <- 40L
d <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:4,
                 g = factor(sprintf("g%02d", seq_len(G))))
z <- matrix(rnorm(3 * G), 3)
u <- t(sweep(z, 2, colMeans(z)))
d$y <- 1 + 1.2 * u[cbind(as.integer(d$g), as.integer(d$t))] +
  rnorm(nrow(d), 0, 0.3)
compare("homcs-30", suppressWarnings(frm(y ~ 1 + homcs(0 + t | g),
                                         family = gaussian(), data = d,
                                         start = list(theta = c(0, -30)))))

# an ordinary us() block, interior
set.seed(62)
n_g <- 30
g <- factor(rep(seq_len(n_g), each = 8))
U <- matrix(stats::rnorm(n_g * 2), n_g) %*%
  chol(matrix(c(0.64, 0.2, 0.2, 0.25), 2))
dd <- data.frame(g = g, x = stats::rnorm(n_g * 8))
dd$y <- 1 + 0.5 * dd$x + U[g, 1] + U[g, 2] * dd$x +
  stats::rnorm(n_g * 8, 0, 0.6)
compare("us", frm(bf(y ~ x + (x | g)) + gaussian(), data = dd))
