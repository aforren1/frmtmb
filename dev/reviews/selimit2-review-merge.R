# Review of lane selimit2, item 7: the merge with quadseam2 (and
# slopescale). Loads a scratch tree where R/fit.R and R/se-check.R are
# git merge-file of selimit2, quadseam2 and slopescale over ec8178d3, and
# asks whether se_limit_ranges() (which calls obj$fn, under one-block
# quadrature the tape frozen at the estimate) measures the level range on
# the objective the limit verdict used (se_probe_fn(), a rule taped at
# each point).
#   Rscript selimit2-review-merge.R <merged tree>
tree <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(tree, export_all = FALSE, quiet = TRUE,
                                   helpers = FALSE))
ns <- asNamespace("frmtmb")
cat("loaded from", getNamespaceInfo(ns, "path"), "\n")
cat("se_limit_ranges uses obj$fn:",
    any(grepl("obj$fn", deparse(ns$se_limit_ranges), fixed = TRUE)), "\n")
cat("se_edge_limit uses se_probe_fn:",
    any(grepl("se_probe_fn", deparse(ns$se_edge_limit), fixed = TRUE)), "\n")
cat("slopescale se_sd_gain_up0 present:", exists("se_sd_gain_up0", ns), "\n")
cat("selimit2 se_limit_joint_dirs present:",
    exists("se_limit_joint_dirs", ns), "\n")
cat("quadseam2 quad_reanchored_fn present:",
    exists("quad_reanchored_fn", ns), "\n")
quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    cat("   W:", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("   M:", substr(conditionMessage(m), 1, 200), "\n")
    invokeRestart("muffleMessage")
  })
}
# the level range per direction with a given objective, as
# se_limit_ranges() does it
ranges_with <- function(fit, names, fnq) {
  li <- ns$se_limit_info(fit)
  om <- ns$outer_par_map(fit)
  p <- fit$opt$par
  u <- ns$se_limit_units(fit, om)
  f0 <- fnq(p)
  out <- list()
  for (n in names) {
    w <- numeric(length(p))
    w[match(n, om$names)] <- 1
    w <- w / sqrt(sum((w / u)^2))
    ext <- vapply(c(-1, 1), function(s) {
      t <- 0
      for (k in 1:12) {
        tk <- 0.5 * 2^(k - 1L)
        d <- fnq(p + s * tk * w) - f0
        if (!is.finite(d) || d > 1e-6) break
        t <- tk
      }
      t
    }, 0)
    out[[n]] <- ext
  }
  if (is.matrix(li$jw)) {
    for (k in seq_len(ncol(li$jw))) {
      w <- li$jw[, k]
      w <- w / sqrt(sum((w / u)^2))
      ext <- vapply(c(-1, 1), function(s) {
        t <- 0
        for (kk in 1:12) {
          tk <- 0.5 * 2^(kk - 1L)
          d <- fnq(p + s * tk * w) - f0
          if (!is.finite(d) || d > 1e-6) break
          t <- tk
        }
        t
      }, 0)
      out[[paste0("joint", k)]] <- ext
    }
  }
  out
}
probe <- function(lbl, fit, nd) {
  sl <- ns$sdr_of(fit)$se_lost
  cat(lbl, "| code", fit$opt$convergence, "| quadrature", fit$quadrature,
      "| lost", paste(names(sl), sl), "\n")
  ln <- names(sl)[sl == "limit"]
  if (!length(ln)) return(invisible())
  fz <- function(q) tryCatch(fit$obj$fn(q), error = function(e) NA_real_)
  pf <- ns$se_probe_fn(fit)
  fp <- function(q) tryCatch(pf(q), error = function(e) NA_real_)
  t0 <- proc.time()[["elapsed"]]
  rz <- ranges_with(fit, ln, fz)
  t1 <- proc.time()[["elapsed"]]
  rp <- ranges_with(fit, ln, fp)
  t2 <- proc.time()[["elapsed"]]
  for (n in names(rz)) {
    cat(sprintf("   %s level range: frozen tape [-%g, %g] (%.1f s) | taped at each point [-%g, %g] (%.1f s)\n",
                n, rz[[n]][1], rz[[n]][2], t1 - t0, rp[[n]][1], rp[[n]][2],
                t2 - t1))
  }
  p <- fit$opt$par
  om <- ns$outer_par_map(fit)
  for (n in ln) {
    j <- match(n, om$names)
    for (t in c(-512, -64, -8, 8, 64, 512)) {
      q <- p
      q[j] <- p[j] + t * ns$se_limit_units(fit, om)[j]
      cat(sprintf("     %s %+5g units: frozen %.3g | re-taped %.3g\n", n, t,
                  fz(q) - fz(p), fp(q) - fp(p)))
    }
  }
  if (!is.null(nd)) {
    f1 <- quiet(stats::fitted(fit, newdata = nd))
    print(cbind(nd, f1[, c("Estimate", "Est.Error")]))
  }
}
# A: poisson, a level all zero, (1 | g), x; quadseam2's runoff data
set.seed(5)
d <- data.frame(g = factor(rep(1:20, each = 9)),
                f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
d$y[d$f == "c"] <- 0L
for (quad in c(FALSE, TRUE)) {
  fit <- quiet(frm(y ~ f + x + (1 | g), family = poisson(), data = d,
                   quadrature = quad))
  probe(paste("A runoff fc, quad", quad), fit,
        data.frame(f = c("a", "c"), x = c(0, 0), g = d$g[1]))
}
# B: zero-inflated poisson whose zi runs to 0, zi ~ x, (1 | g) in mu;
# a new row far outside x's range is asked by the range test
set.seed(9)
d2 <- data.frame(g = factor(rep(1:25, each = 8)), x = rnorm(200))
d2$y <- rpois(200, exp(0.8 + 0.3 * d2$x + rnorm(25, 0, 0.4)[d2$g]))
for (quad in c(FALSE, TRUE)) {
  fit <- quiet(frm(bf(y ~ x + (1 | g), zi ~ x),
                   family = zero_inflated_poisson(), data = d2,
                   quadrature = quad))
  probe(paste("B zi ~ x to 0, quad", quad), fit,
        data.frame(x = c(0, 4, -12, 12), g = d2$g[1]))
}
# C: poisson y ~ x + (1 | g), one count at the largest x (a joint end)
set.seed(4)
d3 <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
d3$y <- 0L
d3$y[which.max(d3$x)] <- 1L
for (quad in c(FALSE, TRUE)) {
  fit <- quiet(frm(y ~ x + (1 | g), family = poisson(), data = d3,
                   quadrature = quad))
  probe(paste("C one count at the top, quad", quad), fit,
        data.frame(x = c(0, max(d3$x) + 1), g = d3$g[1]))
}
cat("END\n")
