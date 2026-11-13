# Review of lane arspat: vcov_cluster() over ar() series against an
# independent sandwich: dense per-cluster log-densities, numDeriv scores
# and Hessian, CR0. Output to stdout.
source("dev/reviews/arspat-review-common.R")
suppressPackageStartupMessages(library(numDeriv))

set.seed(17)
G <- 30; K <- 6
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
n <- nrow(d)
d$x <- rnorm(n) + rep(rnorm(G), each = K)
d$h <- factor((as.integer(d$g) - 1L) %/% 3L + 1L)   # 3 series per h
e <- as.vector(vapply(seq_len(G), function(i)
  as.vector(arima.sim(list(ar = 0.5), K, sd = 0.7)), numeric(K)))
d$y <- 1 + 0.6 * d$x + e + rnorm(nlevels(d$h), 0, 0.5)[d$h] * 0
d$yre <- 1 + 0.6 * d$x + e + rnorm(nlevels(d$h), 0, 0.6)[d$h]
# a cluster variable that splits series
d$split <- factor(ifelse(d$t <= 3, paste0(d$g, "a"), paste0(d$g, "b")))
d$one <- factor(1)

# dense ar(1) density of cluster rows i at theta = (b0, b1, log sigma,
# atanh phi[, log tau])
cl_ll <- function(th, i, yv, re = FALSE) {
  sig <- exp(th[3]); phi <- tanh(th[4])
  mu <- th[1] + th[2] * d$x[i]
  tt <- d$t[i]; gg <- d$g[i]
  S <- sig^2 / (1 - phi^2) * phi^abs(outer(tt, tt, "-")) *
    outer(gg, gg, "==")
  if (re) S <- S + exp(2 * th[5]) * outer(d$h[i], d$h[i], "==")
  mvtnorm::dmvnorm(yv[i], mu, S, log = TRUE)
}
indep <- function(fit, yv, cl, re = FALSE) {
  th0 <- c(0, 0, 0, 0, if (re) 0)
  tot <- function(th) sum(vapply(split(seq_len(n), d$g), function(i)
    cl_ll(th, i, yv, FALSE), 0))
  if (re) tot <- function(th) sum(vapply(split(seq_len(n), d$h), function(i)
    cl_ll(th, i, yv, TRUE), 0))
  o <- optim(th0, function(th) -tot(th), method = "BFGS",
             control = list(reltol = 1e-14, maxit = 1000))
  th <- o$par
  H <- hessian(tot, th)
  S <- t(vapply(split(seq_len(n), cl), function(i)
    grad(function(t) cl_ll(t, i, yv, re), th), numeric(length(th))))
  B <- solve(-H)
  V <- B %*% crossprod(S) %*% B
  list(th = th, V = V[1:2, 1:2], S = S[, 1:2])
}

cmp <- function(label, fit, yv, clvar, re = FALSE) {
  chk(label, {
    ref <- indep(fit, yv, d[[clvar]], re)
    V <- vcov_cluster(fit, stats::as.formula(paste("~", clvar)),
                      type = "CR0")
    cat("beta frmtmb:", format(unname(fixef(fit)[, 1]), digits = 10),
        " indep:", format(ref$th[1:2], digits = 10), "\n")
    cat("V frmtmb:", format(as.vector(V), digits = 8), "\n")
    cat("V indep :", format(as.vector(ref$V), digits = 8), "\n")
    cat(sprintf("max rel diff V: %.3e\n",
                max(abs(as.vector(V) - as.vector(ref$V)) /
                      abs(as.vector(ref$V)))))
    sc <- cluster_scores(fit, stats::as.formula(paste("~", clvar)))
    cat("cluster_scores dim", dim(sc), "\n")
    sb <- as.matrix(sc)[, 1:2]
    cat(sprintf("max |score - indep| / max|indep|: %.3e\n",
                max(abs(sb - ref$S)) / max(abs(ref$S))))
    NULL
  })
}

f1 <- frm(y ~ x + ar(t, g, cov = TRUE), data = d)
cmp("cov=TRUE, cluster = g", f1, d$y, "g")
cmp("cov=TRUE, cluster = h (3 series each)", f1, d$y, "h")
chk("cov=TRUE, cluster splits series", vcov_cluster(f1, ~ split))
chk("cov=TRUE, one cluster", vcov_cluster(f1, ~ one))
chk("cov=TRUE, cluster = t (crosses every series)", vcov_cluster(f1, ~ t))

f2 <- frm(yre ~ x + (1 | h) + ar(t, g, cov = TRUE), data = d)
cmp("cov=TRUE + (1 | h), cluster = h", f2, d$yre, "h", re = TRUE)
chk("cov=TRUE + (1 | h), cluster = g (re crosses)", vcov_cluster(f2, ~ g))

# cov = FALSE, brms's conditional form: mu_t + phi (y_{t-1} - mu_{t-1})
f3 <- frm(y ~ x + ar(t, g), data = d)
chk("cov=FALSE, cluster = g, against conditional likelihood", {
  llc <- function(th, i) {
    sig <- exp(th[3]); phi <- tanh(th[4])
    mu <- th[1] + th[2] * d$x[i]
    r <- d$y[i] - mu
    m <- mu + phi * c(0, r[-length(r)])
    sum(dnorm(d$y[i], m, sig, log = TRUE))
  }
  tot <- function(th) sum(vapply(split(seq_len(n), d$g), function(i)
    llc(th, i), 0))
  o <- optim(c(0, 0, 0, 0), function(th) -tot(th), method = "BFGS",
             control = list(reltol = 1e-14, maxit = 1000))
  H <- hessian(tot, o$par)
  S <- t(vapply(split(seq_len(n), d$g), function(i)
    grad(function(t) llc(t, i), o$par), numeric(4)))
  B <- solve(-H)
  Vr <- (B %*% crossprod(S) %*% B)[1:2, 1:2]
  V <- vcov_cluster(f3, ~ g, type = "CR0")
  cat("V frmtmb:", format(as.vector(V), digits = 8), "\n")
  cat("V indep :", format(as.vector(Vr), digits = 8), "\n")
  cat(sprintf("max rel diff V: %.3e\n",
              max(abs(as.vector(V) - as.vector(Vr)) / abs(as.vector(Vr)))))
  NULL
})
chk("cov=FALSE, cluster splits series", vcov_cluster(f3, ~ split))

chk("marginaleffects vcov = ~ g", {
  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    a <- marginaleffects::avg_slopes(f1, variables = "x", vcov = ~ g)
    V <- vcov_cluster(f1, ~ g, type = "CR1")
    cat("avg_slopes SE", a$std.error, " sqrt(V_CR1[x,x])", sqrt(V[2, 2]),
        " sqrt(V_CR0)", sqrt(vcov_cluster(f1, ~g, type = "CR0")[2, 2]), "\n")
  }
  NULL
})
chk("HAC refusal text", marginaleffects::avg_slopes(f1, variables = "x",
                                                    vcov = "HAC"))
chk("fcor refused", {
  M <- 0.5^abs(outer(1:n, 1:n, "-"))
  ff <- frm(y ~ x + fcor(M), data = d, data2 = list(M = M))
  vcov_cluster(ff, ~ g)
})
