# Reviewer of lane illse, round 2, items 1 to 3. (a) Guard absent: on
# well-conditioned fits vcov_cluster(), hypothesis(vcov = vcov_cluster()),
# insight::get_varcov() and marginaleffects must be 0.70.0's, not NaN.
# Saves them for a bitwise comparison of the two arms. (b) On a
# degree-8 fit held as a factor: a factorless matrix passed to
# hypothesis(vcov =), a single coefficient with it, vcov_cluster(full =
# FALSE), and summary(vcov = vcov_cluster()).
#   REV_ARM=lane|base Rscript dev/illse-rev2-guard.R <out.rds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
out <- commandArgs(TRUE)[1]
res <- list()
for (s in 1:3) {
  set.seed(s)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                  g = factor(rep(1:20, length.out = n)))
  u <- stats::rnorm(20, 0, 0.5)[d$g]
  d$y <- 1 + 0.5 * d$x + sin(3 * d$z) + u + stats::rnorm(n, 0, 0.4)
  d$yb <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$x + u))
  specs <- list(lmm = list(y ~ x + z + (1 | g), gaussian()),
                glmm = list(yb ~ x + z + (1 | g), bernoulli()),
                glm = list(yb ~ x + z, bernoulli()),
                poly3 = list(y ~ poly(x, 3) + z, gaussian()))
  nd <- d[c(1, 50, 120), ]
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    f <- suppressWarnings(frm(sp[[1]], family = sp[[2]], data = d))
    vc <- cap(vcov_cluster(f, ~ g, full = TRUE))
    hs <- "Intercept + x - 2 * z = 0"
    hc <- if (nm == "poly3") NULL else cap(hypothesis(f, hs, vcov = vc$v))
    iv <- cap(insight::get_varcov(f))
    mf <- cap(marginaleffects::predictions(f, newdata = nd))
    res[[paste(nm, s)]] <- list(
      root = !is.null(ns$sdr_of(f)$cov_fixed_root),
      vc = unclass(vc$v)[, ], vc_w = vc$w,
      hc = if (!is.null(hc)) hc$v$hypothesis$Est.Error, hc_w = hc$w,
      iv = iv$v, iv_w = iv$w,
      mf = if (!inherits(mf$v, "rev_err")) mf$v$std.error, mf_w = mf$w)
    cat(sprintf("%-6s s%d | root %s | vcov_cluster finite %s warn %s | hyp(vcov) %s | insight finite %s warn %s | mfx finite %s warn %s\n",
                nm, s, res[[paste(nm, s)]]$root, all(is.finite(vc$v)), wtag(vc$w),
                if (is.null(hc)) "-" else format(hc$v$hypothesis$Est.Error, digits = 6),
                all(is.finite(iv$v)), wtag(iv$w),
                all(is.finite(res[[paste(nm, s)]]$mf)), wtag(mf$w)))
  }
}
saveRDS(res, out)
# (b) degree 8 held as a factor
set.seed(1)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, length.out = 200)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(200, 0, 0.3)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + (1 | g)
f <- suppressWarnings(frm(fo, data = d))
cat("deg 8 root", !is.null(ns$sdr_of(f)$cov_fixed_root), "\n")
V <- suppressWarnings(vcov(f, full = TRUE))
hs <- paste0("Intercept + ", paste0(1.5^(1:8), " * x", 1:8, collapse = " + "), " = 0")
h1 <- cap(hypothesis(f, hs, vcov = V))
h2 <- cap(hypothesis(f, "x8 = 0", vcov = V))
h3 <- cap(hypothesis(f, "x8 = 0"))
cat(sprintf("factorless V: combination SE %s warn %s | x8 alone SE %s (own %s) warn %s\n",
            format(h1$v$hypothesis$Est.Error), wtag(h1$w),
            format(h2$v$hypothesis$Est.Error), format(h3$v$hypothesis$Est.Error), wtag(h2$w)))
cat("  combination warning:", substr(paste(h1$w, collapse = " / "), 1, 300), "\n")
vf <- cap(vcov_cluster(f, ~ g))
cat("vcov_cluster(full = FALSE): dim", dim(vf$v), "root attr", !is.null(attr(vf$v, "root")),
    "finite", all(is.finite(vf$v)), "\n")
hf <- cap(hypothesis(f, hs, vcov = vf$v))
cat("  hypothesis(vcov = full = FALSE matrix):", if (inherits(hf$v, "rev_err")) substr(hf$v, 1, 120) else format(hf$v$hypothesis$Est.Error), wtag(hf$w), "\n")
vt <- cap(vcov_cluster(f, ~ g, full = TRUE))
vs <- vt$v
vs2 <- vs[, ]
cat("  subsetting [, ] keeps root:", !is.null(attr(vs2, "root")), "\n")
hs2 <- cap(hypothesis(f, hs, vcov = vs2))
cat("  hypothesis(vcov = V[, ]):", format(hs2$v$hypothesis$Est.Error), wtag(hs2$w), "\n")
hs3 <- cap(hypothesis(f, hs, vcov = vs))
cat("  hypothesis(vcov = vcov_cluster(full = TRUE)):", format(hs3$v$hypothesis$Est.Error), wtag(hs3$w), "\n")
hs4 <- cap(hypothesis(f, hs, vcov = function(x) vcov_cluster(x, ~ g, full = TRUE)))
cat("  hypothesis(vcov = function):", if (inherits(hs4$v, "rev_err")) substr(hs4$v, 1, 100) else format(hs4$v$hypothesis$Est.Error), wtag(hs4$w), "\n")
sm <- cap(summary(f, vcov = vs))
cat("  summary(vcov = vcov_cluster(full = TRUE)):", if (inherits(sm$v, "rev_err")) substr(sm$v, 1, 100) else
  paste(signif(sm$v$fixed[, "Est.Error"][1:3], 5), collapse = " "), wtag(sm$w), "| sqrt diag",
  paste(signif(sqrt(diag(vs))[1:3], 5), collapse = " "), "\n")
