# The routes the release review left unchecked, beside a latent residual: me(gr =), hsgp,
# spatial and latent covstructs, mixtures. Each reads fitted(newdata) against a reference
# with no frmtmb reader of the term: the same fit with the other block out, and the
# latent's draw by its marginal sd. A refusal is reported as it stands.
#   Rscript dev/rel079-probe-routes.R <library>
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
    cat(sprintf("%-26s REFUSED: %s\n", lab, substr(gsub("\n", " ", conditionMessage(r)), 1, 150)))
    return(NULL)
  }
  r
}
base <- function(seed, sd_e = 0.3) {
  set.seed(seed)
  d <- expand.grid(t = 1:25, s = factor(1:8))
  d$x <- rnorm(nrow(d))
  d$zo <- rnorm(nrow(d))
  d$sz <- 0.3
  d$e <- ave(rnorm(nrow(d)), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
  d$re <- rnorm(8, 0, 0.5)[d$s]
  d$y <- rpois(nrow(d), exp(2 + 0.3 * d$x + d$re + sd_e * d$e))
  d
}
ar <- "ar(time = t, gr = s, cov = TRUE)"
report <- function(lab, f, f0, nd, ...) {
  if (is.null(f)) return(invisible())
  a <- tryCatch(suppressWarnings(suppressMessages(fitted(f, newdata = nd, ...))),
                error = function(e) e)
  if (inherits(a, "error")) {
    cat(sprintf("%-26s fitted REFUSED: %s\n", lab, substr(gsub("\n", " ", conditionMessage(a)), 1, 150)))
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
for (seed in c(7401, 7402)) {
  cat("\n===== seed", seed, "\n")
  d <- base(seed)
  # me(gr =): a latent value per group level
  d$lat <- rnorm(8)[d$s]
  d$xg <- d$lat + rnorm(8, 0, 0.3)[d$s]
  d$sg <- 0.3
  d$y2 <- rpois(nrow(d), exp(1.5 + 0.5 * d$lat + 0.3 * d$e))
  f <- try_fit("me(gr =)", frm(as.formula(paste("y2 ~ me(xg, sg, gr = s) +", ar)), data = d, family = poisson()))
  f0 <- try_fit("me(gr =), no latent", frm(y2 ~ me(xg, sg, gr = s), data = d, family = poisson()))
  nd <- data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), xg = c(0, 1), sg = 0.3)
  report("me(gr =), seen level", f, f0, nd)
  nd2 <- data.frame(t = 26, s = factor("new"), xg = c(0, 1), sg = 0.3)
  report("me(gr =), new level", f, f0, nd2, allow_new_levels = TRUE)

  # hsgp
  f <- try_fit("hsgp", frm(as.formula(paste("y ~ gp(x, k = 10) +", ar)), data = d, family = poisson()))
  f0 <- try_fit("hsgp, no latent", frm(y ~ gp(x, k = 10), data = d, family = poisson()))
  nd <- data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), x = c(0, 1))
  report("gp(k =)", f, f0, nd)

  # smooth
  f <- try_fit("s()", frm(as.formula(paste("y ~ s(x) +", ar)), data = d, family = poisson()))
  f0 <- try_fit("s(), no latent", frm(y ~ s(x), data = d, family = poisson()))
  report("s()", f, f0, nd)

  # a spatial covstruct: car() over the groups
  W <- matrix(0, 8, 8, dimnames = list(1:8, 1:8))
  W[cbind(1:7, 2:8)] <- 1
  W <- W + t(W)
  f <- try_fit("car", frm(as.formula(paste("y ~ x + car(W, gr = s) +", ar)), data = d, family = poisson(), data2 = list(W = W)))
  f0 <- try_fit("car, no latent", frm(y ~ x + car(W, gr = s), data = d, family = poisson(), data2 = list(W = W)))
  nd3 <- data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), x = c(0, 1))
  report("car()", f, f0, nd3)

  # mixture_mvn from frmtmb.latent and a mixture of families beside a latent
  f <- try_fit("mixture(poisson, poisson)", frm(as.formula(paste("y ~ x +", ar)), data = d, family = mixture(poisson(), poisson())))
  report("mixture", f, NULL, nd3)
  d$ym <- as.numeric(scale(d$x)) + 0.3 * d$e + rnorm(nrow(d))
  f <- try_fit("mixture_mvn", frm(as.formula(paste("mvbind(ym, y) ~ x +", ar)), data = d, family = mixture_mvn(2)))
  if (!is.null(f)) report("mixture_mvn", f, NULL, nd3)
  # a latent covstruct on a latent family
  f <- try_fit("gr(cov = A)", {
    A <- diag(8); dimnames(A) <- list(1:8, 1:8)
    frm(as.formula(paste("y ~ x + (1 | gr(s, cov = A)) +", ar)), data = d, family = poisson(), data2 = list(A = A))
  })
  f0 <- try_fit("gr(cov = A), no latent", {
    A <- diag(8); dimnames(A) <- list(1:8, 1:8)
    frm(y ~ x + (1 | gr(s, cov = A)), data = d, family = poisson(), data2 = list(A = A))
  })
  report("gr(cov = A)", f, f0, nd3)
  f <- try_fit("mm()", {
    d$s2 <- factor(sample(1:8, nrow(d), TRUE), levels = 1:8)
    frm(as.formula(paste("y ~ x + (1 | mm(s, s2)) +", ar)), data = d, family = poisson())
  })
  f0 <- try_fit("mm(), no latent", frm(y ~ x + (1 | mm(s, s2)), data = d, family = poisson()))
  nd4 <- data.frame(t = 26, s = factor(c(1, 2), levels = 1:8), s2 = factor(c(3, 4), levels = 1:8), x = c(0, 1))
  report("mm()", f, f0, nd4)
}
cat("PROBE END\n")
