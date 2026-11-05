# More routes beside a latent residual: rr(), spde(), sar() with ar(), mixture_mvn(),
# cs() and me() together, zi and hurdle families. The same reference as
# dev/rel079-probe-routes.R, or a refusal as it stands.
#   Rscript dev/rel079-probe-routes2.R <library>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spatial); library(frmtmb.latent)})
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")),
    format(file.mtime(system.file("DESCRIPTION", package = "frmtmb"))), "\n")
ns <- asNamespace("frmtmb")
q <- function(x) signif(x, 4)
try_fit <- function(lab, e) {
  r <- tryCatch(suppressMessages(suppressWarnings(e)), error = function(err) err)
  if (inherits(r, "error")) {
    cat(sprintf("%-26s REFUSED: %s\n", lab, substr(gsub("\n", " ", conditionMessage(r)), 1, 170)))
    return(NULL)
  }
  r
}
report <- function(lab, f, f0, nd, ...) {
  if (is.null(f)) return(invisible())
  a <- tryCatch(suppressWarnings(suppressMessages(fitted(f, newdata = nd, ...))), error = function(e) e)
  if (inherits(a, "error")) {
    cat(sprintf("%-26s fitted REFUSED: %s\n", lab, substr(gsub("\n", " ", conditionMessage(a)), 1, 170)))
    return(invisible())
  }
  b <- if (!is.null(f0)) suppressWarnings(suppressMessages(fitted(f0, newdata = nd, ...)))
  rn <- names(f$spec$responses)[1L]
  sl <- ns$latres_row_sd(f, rn, NULL)
  cat(sprintf("%-26s sl %s\n", lab, format(sl, digits = 3)))
  cat("   with latent  est", q(a[, "Estimate"]), " se", q(a[, "Est.Error"]), " Q", q(a[, "Q2.5"]), "/", q(a[, "Q97.5"]), "\n")
  if (!is.null(b)) {
    m <- b[, "Estimate"] * exp(sl^2 / 2)
    v <- sqrt(b[, "Est.Error"]^2 * exp(sl^2) + (b[, "Estimate"] * exp(sl^2 / 2) * sqrt(exp(sl^2) - 1))^2)
    cat("   no latent    est", q(b[, "Estimate"]), " se", q(b[, "Est.Error"]), "\n")
    cat("   log-link ref est", q(m), " se", q(v), " est ratio", q(a[, "Estimate"] / m), " se ratio", q(a[, "Est.Error"] / v), "\n")
  }
}
ar <- "ar(time = t, gr = s, cov = TRUE)"
for (seed in c(7403, 7404)) {
  cat("\n===== seed", seed, "\n")
  set.seed(seed)
  d <- expand.grid(t = 1:25, s = factor(1:8))
  n <- nrow(d)
  d$x <- rnorm(n); d$z <- rnorm(n)
  d$e <- ave(rnorm(n), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
  d$re <- rnorm(8, 0, 0.5)[d$s]
  d$y <- rpois(n, exp(2 + 0.3 * d$x + d$re + 0.3 * d$e))
  nd <- data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), x = c(0, 1), z = 0)

  # rr()
  f <- try_fit("rr()", frm(as.formula(paste("y ~ x + (0 + x + z | rr(s, d = 1)) +", ar)), data = d, family = poisson()))
  f0 <- try_fit("rr(), no latent", frm(y ~ x + (0 + x + z | rr(s, d = 1)), data = d, family = poisson()))
  report("rr()", f, f0, nd)
  # spde(): positions on a mesh of the groups
  f <- try_fit("spde()", {
    set.seed(seed + 1L)
    d$node <- factor(sample(1:20, n, TRUE), levels = 1:20)
    C <- diag(20)
    G <- diag(c(1, rep(2, 18), 1))
    G[cbind(1:19, 2:20)] <- -1
    G[cbind(2:20, 1:19)] <- -1
    fem <- list(c0 = C, g1 = G, g2 = G %*% G)
    frm(as.formula(paste("y ~ x + spde(fem, gr = node) +", ar)), data = d, family = poisson(), data2 = list(fem = fem))
  })
  # sar() with ar(): two residual terms
  f <- try_fit("sar() + ar()", {
    W <- matrix(0, n, n)
    W[cbind(1:(n - 1), 2:n)] <- 1
    W <- W + t(W)
    frm(as.formula(paste("y ~ x + sar(W) +", ar)), data = d, family = poisson(), data2 = list(W = W))
  })
  if (!is.null(f)) report("sar() + ar()", f, NULL, nd)
  # mixture_mvn(): the D columns under a latent
  f <- try_fit("mixture_mvn()", {
    d$Y <- cbind(d$x + d$e + rnorm(n), d$x - d$e + rnorm(n))
    frm(as.formula(paste("bf(Y ~ 1 +", ar, ")")), data = d, family = mixture_mvn(K = 2, D = 2))
  })
  f <- try_fit("mixture_mvn(), bf", {
    d$Y <- cbind(d$x + d$e + rnorm(n), d$x - d$e + rnorm(n))
    frm(bf(as.formula(paste("Y ~ 1 +", ar))) + mixture_mvn(K = 2, D = 2), data = d)
  })
  if (!is.null(f)) report("mixture_mvn()", f, NULL, nd)
  # hurdle and zero-inflated rows, closed form on the count
  d$h <- d$y * rbinom(n, 1, 0.7)
  f <- try_fit("zero_inflated_poisson", frm(as.formula(paste("h ~ x +", ar)), data = d, family = zero_inflated_poisson()))
  f0 <- try_fit("zip, no latent", frm(h ~ x, data = d, family = zero_inflated_poisson()))
  report("zero_inflated_poisson", f, f0, nd)
  f <- try_fit("hurdle_poisson", frm(as.formula(paste("h ~ x +", ar)), data = d, family = hurdle_poisson()))
  f0 <- try_fit("hurdle, no latent", frm(h ~ x, data = d, family = hurdle_poisson()))
  report("hurdle_poisson", f, f0, nd)
}
cat("PROBE END\n")
