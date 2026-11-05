# The mean on a new row against its closed form within one fit, on the routes beside a
# latent residual: exp(eta + s^2 / 2) on a log link with s the node's sd, times 1 - zi, or
# the hurdle's own mean. Est.Error against the draw's sd alone, which is its lower bound.
# Not the comparison across fits of dev/rel079-probe-routes.R, whose no-latent intercept
# has already absorbed s^2 / 2, so its "ref" over-reads the mean.
#   Rscript dev/rel079-probe-routes3.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r17", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spatial); library(frmtmb.latent)})
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")),
    format(file.mtime(system.file("DESCRIPTION", package = "frmtmb"))), "\n")
ns <- asNamespace("frmtmb")
q <- function(x) signif(x, 4)
go <- function(lab, form, d, nd, family = poisson(), data2 = NULL, scale_m = function(fit, nd) 1, ...) {
  f <- tryCatch(suppressMessages(suppressWarnings(frm(form, data = d, family = family, data2 = data2))),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-22s fit REFUSED: %s\n", lab, substr(gsub("\n", " ", conditionMessage(f)), 1, 130)))
    return(invisible())
  }
  a <- tryCatch(suppressWarnings(suppressMessages(fitted(f, newdata = nd, ...))), error = function(e) e)
  if (inherits(a, "error")) {
    cat(sprintf("%-22s fitted REFUSED: %s\n", lab, substr(gsub("\n", " ", conditionMessage(a)), 1, 130)))
    return(invisible())
  }
  rn <- names(f$spec$responses)[1L]
  rs <- f$spec$responses[[rn]]
  L <- ns$me_row_loads(f, rs, nd, TRUE, FALSE, NULL)
  s2 <- L$s2
  eta <- as.numeric(suppressWarnings(frm_linpred(f, newdata = nd, type = "link", dpar = "mu")))
  m <- exp(eta + s2 / 2) * scale_m(f, nd)
  draw <- m * sqrt(exp(s2) - 1)
  cat(sprintf("%-22s node sd %s  mean ratio %s  se / draw %s  width %s\n", lab, q(sqrt(s2[1L])),
              paste(q(a[, "Estimate"] / m), collapse = " "), paste(q(a[, "Est.Error"] / draw), collapse = " "),
              paste(q(a[, "Q97.5"] - a[, "Q2.5"]), collapse = " ")))
}
ar <- " + ar(time = t, gr = s, cov = TRUE)"
fm <- function(x) as.formula(paste0(x, ar))
for (seed in c(7405, 7406)) {
  cat("\n===== seed", seed, "\n")
  set.seed(seed)
  d <- expand.grid(t = 1:25, s = factor(1:8))
  n <- nrow(d)
  d$x <- rnorm(n); d$z <- rnorm(n)
  d$e <- ave(rnorm(n), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
  re <- rnorm(8, 0, 0.5)
  d$y <- rpois(n, exp(2 + 0.3 * d$x + re[d$s] + 0.3 * d$e))
  lat <- rnorm(8)
  d$xg <- (lat + rnorm(8, 0, 0.3))[d$s]; d$sg <- 0.3
  d$y2 <- rpois(n, exp(1.5 + 0.5 * lat[d$s] + 0.3 * d$e))
  d$zo <- d$z + rnorm(n, 0, 0.3); d$sz <- 0.3
  d$y3 <- rpois(n, exp(1.5 + 0.5 * d$z + 0.3 * d$e))
  d$h <- d$y * rbinom(n, 1, 0.7)
  d$s2 <- factor(sample(1:8, n, TRUE), levels = 1:8)
  d$node <- factor(sample(1:20, n, TRUE), levels = 1:20)
  nd <- data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), s2 = factor(c(3, 4), levels = 1:8), x = c(0, 1),
                   z = 0, node = factor(c(1, 2), levels = 1:20))
  go("(1 | s)", fm("y ~ x + (1 | s)"), d, nd)
  go("me(gr =), seen", fm("y2 ~ me(xg, sg, gr = s)"), d, data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), xg = c(0, 1), sg = 0.3))
  go("me()", fm("y3 ~ me(zo, sz)"), d, data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), zo = c(0, 1), sz = 0.3))
  go("gp(k =)", fm("y ~ gp(x, k = 10)"), d, nd)
  go("exact gp", fm("y ~ gp(x)"), d, nd)
  go("s()", fm("y ~ s(x)"), d, nd)
  W <- matrix(0, 8, 8, dimnames = list(1:8, 1:8))
  W[cbind(1:7, 2:8)] <- 1
  W <- W + t(W)
  go("car()", fm("y ~ x + car(W, gr = s)"), d, nd, data2 = list(W = W))
  A <- diag(8); dimnames(A) <- list(1:8, 1:8)
  go("gr(cov = A)", fm("y ~ x + (1 | gr(s, cov = A))"), d, nd, data2 = list(A = A))
  go("mm()", fm("y ~ x + (1 | mm(s, s2))"), d, nd)
  C <- diag(20)
  G <- diag(c(1, rep(2, 18), 1)); G[cbind(1:19, 2:20)] <- -1; G[cbind(2:20, 1:19)] <- -1
  fem <- list(c0 = C, g1 = G, g2 = G %*% G)
  go("spde()", fm("y ~ x + spde(fem, gr = node)"), d, nd, data2 = list(fem = fem))
  go("rr", fm("y ~ x + (0 + x + z | s, rr = 1)"), d, nd)
  zi <- function(f, nd) 1 - stats::plogis(unname(fixef(f)["zi_Intercept", "Estimate"]))
  go("zero_inflated_poisson", fm("h ~ x"), d, nd, family = zero_inflated_poisson(),
     scale_m = function(f, nd) 1 - as.numeric(frm_linpred(f, newdata = nd, dpar = "zi", type = "response")))
  go("hmm()", fm("y ~ x"), d, nd, family = hmm(K = 2, poisson(), time = t, group = s))
  go("lca family", fm("y ~ x"), d, nd, family = mixture(poisson(), poisson()))
}
cat("PROBE END\n")
