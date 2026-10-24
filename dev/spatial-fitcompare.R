# Part C: fits must be unchanged by the move. Fits the car() and spde()
# models of the moved tests (test-car-spde.R, test-predict-re-
# uncertainty.R, test-tmb-examples.R) and records logLik, every
# estimate, every standard error, fitted() at the data, fitted() and
# predict(re_formula = NA) on newdata with an unseen location,
# ranef() conditional SDs, confint_varcorr() and simulate().
#
#   Rscript dev/spatial-fitcompare.R fit <out.rds>   (SPATIAL_LIBS,
#                                                     SPATIAL_PKGS)
#   Rscript dev/spatial-fitcompare.R compare <before.rds> <after.rds>
#
# Run each arm in ONE process; dev/rtmb-pitfalls.md item 21 says a
# fresh process can differ by a few ulps by core type, so the compare
# reports the largest relative difference and does not assume bits.
a <- commandArgs(trailingOnly = TRUE)
if (a[1] == "compare") {
  b <- readRDS(a[2]); f <- readRDS(a[3])
  cat("before:", b$lib, b$version, "| after:", f$lib, f$version,
      paste(f$pkgs, collapse = ","), "\n")
  for (q in names(f$pkg_libs)) cat("lib:", q, f$pkg_libs[[q]], "\n")
  rel <- function(x, y) {
    x <- unlist(x); y <- unlist(y)
    if (length(x) != length(y)) return(Inf)
    ok <- is.finite(x) & is.finite(y)
    if (!identical(is.finite(x), is.finite(y))) return(Inf)
    if (!any(ok)) return(0)
    max(abs(x[ok] - y[ok]) / pmax(abs(y[ok]), 1e-300))
  }
  for (m in names(b$fits)) {
    x <- b$fits[[m]]; y <- f$fits[[m]]
    cat(sprintf("%-14s", m))
    for (k in names(x)) {
      r <- rel(y[[k]], x[[k]])
      ident <- identical(unname(unlist(y[[k]])), unname(unlist(x[[k]])))
      cat(sprintf(" %s=%s", k, if (ident) "identical" else
        format(r, digits = 2)))
    }
    cat("\n")
  }
  quit(save = "no")
}
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pk <- strsplit(Sys.getenv("SPATIAL_PKGS"), ",", fixed = TRUE)[[1]]
for (p in pk) suppressMessages(library(p, character.only = TRUE))

lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c))
  n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
  }
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
car_data <- function(seed, r = 4, c = 4, per = 6) {
  set.seed(seed)
  W <- lattice_W(r, c)
  n <- nrow(W)
  K <- diag(rowSums(W)) - W + matrix(1 / (1e-3 * n)^2, n, n)
  phi <- 1.2 * drop(crossprod(chol(solve(K)), stats::rnorm(n)))
  loc <- factor(rep(rownames(W), each = per), levels = rownames(W))
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + stats::rnorm(nrow(d), 0, 0.5)
  d$cnt <- stats::rpois(nrow(d), exp(0.2 + 0.3 * d$x + 0.5 * phi[d$loc]))
  list(d = d, W = W)
}
record <- function(fit, nd, re_new = TRUE) {
  est <- unlist(fit$estimates)
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  fi <- fitted(fit)
  out <- list(logLik = as.numeric(logLik(fit)), est = est, se = unname(se),
              fitted = unname(fi[, c("Estimate", "Est.Error")]),
              fitted_new = unname(fitted(fit, newdata = nd,
                                         allow_new_levels = TRUE)[
                                           , c("Estimate", "Est.Error")]),
              pop_new = unname(fitted(fit, newdata = nd,
                                      re_formula = NA)[, "Estimate"]),
              condsd = unname(attr(ranef(fit, condVar = TRUE)[[1]],
                                   "condSD")),
              vc = confint_varcorr(fit)[, c("estimate", "lwr", "upr")],
              sim = unname(as.matrix(simulate(fit, nsim = 3, seed = 7))))
  out
}
fits <- list()
s <- car_data(1)
nd <- s$d[c(1, 30, 50), ]
nd$loc <- factor(c("L1", "L9", "nowhere"))
for (ty in c("escar", "esicar", "icar", "bym2")) {
  fo <- stats::as.formula(paste0("y ~ x + car(W, gr = loc, type = '", ty,
                                 "')"))
  f <- suppressWarnings(frm(bf(fo), family = gaussian(), data = s$d,
                            data2 = list(W = s$W)))
  fits[[paste0("gauss_", ty)]] <- record(f, nd)
}
f <- suppressWarnings(frm(bf(cnt ~ x + car(W, gr = loc, type = "icar")),
                          family = poisson(), data = s$d,
                          data2 = list(W = s$W)))
fits$pois_icar <- record(f, nd)
# a disconnected graph: two 2 x 2 blocks and a singleton
W2 <- matrix(0, 9, 9, dimnames = list(paste0("A", 1:9), paste0("A", 1:9)))
for (blk in list(1:4, 5:8)) {
  Wb <- lattice_W(2, 2)
  W2[blk, blk] <- unname(Wb)
}
set.seed(3)
d2 <- data.frame(loc = factor(rep(rownames(W2), each = 5),
                              levels = rownames(W2)), x = rnorm(45))
d2$y <- 1 + d2$x + rnorm(9)[d2$loc] + rnorm(45, 0, 0.4)
nd2 <- d2[c(1, 40), ]
for (ty in c("esicar", "icar")) {
  f <- suppressWarnings(frm(bf(stats::as.formula(paste0(
    "y ~ x + car(W2, gr = loc, type = '", ty, "')"))),
    family = gaussian(), data = d2, data2 = list(W2 = W2)))
  fits[[paste0("disc_", ty)]] <- record(f, nd2)
}
# spde over a 1-D chain
set.seed(5)
nn <- 20; h <- 0.5
C0 <- diag(rep(h, nn)); C0[1, 1] <- h / 2; C0[nn, nn] <- h / 2
G <- matrix(0, nn, nn)
for (i in seq_len(nn - 1)) {
  G[i, i] <- G[i, i] + 1 / h; G[i + 1, i + 1] <- G[i + 1, i + 1] + 1 / h
  G[i, i + 1] <- -1 / h; G[i + 1, i] <- -1 / h
}
fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
k2 <- 0.8^2
u <- drop(crossprod(chol(solve(k2 * k2 * C0 + 2 * k2 * G + fem$g2)),
                    rnorm(nn)))
node <- factor(rep(seq_len(nn), each = 6))
d3 <- data.frame(node = node, x = rnorm(nn * 6))
d3$y <- 0.3 + 0.4 * d3$x + u[as.integer(d3$node)] + rnorm(nrow(d3), 0, 0.3)
nd3 <- d3[c(2, 40, 100), ]
f <- frm(bf(y ~ x + spde(fem, gr = node)), family = gaussian(), data = d3)
fits$spde_chain <- record(f, nd3)
# spde on an fmesher mesh, Weibull with right censoring
if (requireNamespace("fmesher", quietly = TRUE)) {
  set.seed(126)
  n <- 200
  loc <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  mesh <- fmesher::fm_mesh_2d(loc = loc, max.edge = c(0.2, 0.5),
                              cutoff = 0.05)
  fm <- fmesher::fm_fem(mesh)
  fm <- list(c0 = fm$c0, g1 = fm$g1, g2 = fm$g2)
  idx <- mesh$idx$loc
  Q_t <- exp(-1)^2 * (exp(1.5)^4 * fm$c0 + 2 * exp(1.5)^2 * fm$g1 + fm$g2)
  x_t <- as.vector(backsolve(chol(as.matrix(Q_t)), rnorm(nrow(fm$c0))))
  sexv <- rbinom(n, 1, 0.5); agev <- rnorm(n)
  eta <- -0.5 + 0.3 * sexv + 0.2 * agev + x_t[idx]
  time <- rweibull(n, shape = 1.2, scale = exp(-eta))
  dd <- data.frame(time = pmin(time, 2), cens = as.integer(time > 2),
                   sex = sexv, age = agev,
                   node = factor(idx, levels = seq_len(mesh$n)))
  f <- frm(bf(time | cens(cens) ~ sex + age + spde(fm, gr = node)),
           family = weibull(), data = dd)
  fits$spde_weibull <- record(f, dd[c(1, 2, 3), ])
}
out <- a[2]
saveRDS(list(fits = fits, lib = dirname(find.package("frmtmb")),
             pkg_libs = vapply(pk, function(q) dirname(find.package(q)), ""),
             version = as.character(packageVersion("frmtmb")), pkgs = pk),
        out)
cat("wrote", out, "models", length(fits), "\n")
