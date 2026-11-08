# Review of lane latperf, attack 4: refusals and limits, the routes that still
# skip, the factor-held covariance (dense route), and the coef_pos trim.
# Usage: Rscript dev/review-latperf-refuse.R lane|base   (seed 11)
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-latperf-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat(sprintf("lib: %s frmtmb %s\n", dirname(find.package("frmtmb")),
            format(packageVersion("frmtmb"))))
ns <- asNamespace("frmtmb")
sparse_of <- function(fit) {
  f <- get0("jc_is_sparse", envir = ns)
  if (is.null(f)) NA else f(frmtmb:::get_joint_cov(fit))
}
catch <- function(expr) {
  msgs <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    message = function(m) {
      msgs <<- c(msgs, paste0("[", class(m)[1], "] ", conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) {
      msgs <<- c(msgs, paste0("ERROR [", paste(class(e), collapse = ","),
                              "] ", conditionMessage(e)))
      NULL
    })
  list(v = v, msgs = msgs)
}
pr <- function(lbl, x) {
  cat(" ", lbl, "\n")
  for (m in x$msgs) cat("    ", substr(gsub("\n", " ", m), 1, 400), "\n")
}
set.seed(11)
G <- 30
K <- 12
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$u <- runif(n, 0, 10)
d$s <- factor(sample.int(6, n, TRUE))
e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), K, sd = 0.6)))
eta <- 0.6 + 0.4 * d$x + e
d$cnt <- rpois(n, exp(eta))
d$cntp <- rpois(n, exp(0.3 + 0.1 * d$u - 0.008 * d$u^2 + e))
d$ord <- cut(eta + rlogis(n), c(-Inf, 0, 1, 2, Inf), ordered_result = TRUE)
nd <- data.frame(t = 1:3, g = factor(rep("new", 3)), x = 0,
                 s = factor(1:3, levels = levels(d$s)))
q <- function(e) suppressWarnings(suppressMessages(e))

cat("\n== R1 refusals of a dense block (frmtmb.latres_dense_rows = 10)\n")
fit <- q(frm(cnt ~ x + (1 | s) + ar(t, g, cov = TRUE), poisson(), data = d))
cat("  sparse:", sparse_of(fit), "\n")
local({
  op <- options(frmtmb.latres_dense_rows = 10L)
  on.exit(options(op))
  pr("frm_joint_cov(fit)", catch(frm_joint_cov(fit)))
  pr("frm_lp_basis(fit)", catch(frm_lp_basis(fit)))
  r <- catch(frm_lp_basis(fit, newdata = nd))
  pr(sprintf("frm_lp_basis(fit, newdata) ok %s", !is.null(r$v)), r)
  r <- catch(conditional_effects(fit, "x"))
  pr(sprintf("conditional_effects(fit, 'x') ok %s", !is.null(r$v)), r)
})

cat("\n== R2 coef_pos of frm_lp_basis(newdata) with x = 0 on every row\n")
fit0 <- q(frm(cnt ~ x + (1 | s), poisson(), data = d))
for (ff in list(list("latent", fit), list("no latent", fit0))) {
  lb <- frm_lp_basis(ff[[2]], newdata = nd)
  cat(sprintf("  %-9s ncol(A) %d coef_names: %s\n", ff[[1]], ncol(lb$A),
              paste(lb$coef_names, collapse = ", ")))
}

cat("\n== R3 ordinal fitted rows (frmtmb.latres_se_rows = 10)\n")
fo <- q(frm(ord ~ x + ar(t, g, cov = TRUE), cumulative(), data = d))
local({
  op <- options(frmtmb.latres_se_rows = 10L)
  on.exit(options(op))
  r <- catch(fitted(fo))
  pr(sprintf("fitted(): Est.Error NA %d of %d", sum(is.na(r$v[, "Est.Error", ])),
             length(r$v[, "Est.Error", ])), r)
  r <- catch(frm_linpred(fo, se.fit = TRUE))
  pr(sprintf("frm_linpred(se.fit): NA %d of %d", sum(is.na(r$v$se.fit)),
             length(r$v$se.fit)), r)
})

cat("\n== R4 cs() layers\n")
fc <- catch(frm(ord ~ cs(x) + ar(t, g, cov = TRUE), acat(), data = d))
pr(sprintf("acat cs(x) + ar fit ok %s", !is.null(fc$v)), fc)
if (!is.null(fc$v)) {
  fc <- fc$v
  a <- catch(frm_linpred(fc, se.fit = TRUE))
  cat("  default limit: se.fit dims", paste(dim(a$v$se.fit), collapse = "x"),
      "NA", sum(is.na(a$v$se.fit)), "\n")
  local({
    op <- options(frmtmb.latres_se_rows = 10L)
    on.exit(options(op))
    r <- catch(frm_linpred(fc, se.fit = TRUE))
    pr(sprintf("se_rows = 10: se.fit NA %d of %d", sum(is.na(r$v$se.fit)),
               length(r$v$se.fit)), r)
    r <- catch(fitted(fc))
    pr(sprintf("se_rows = 10: fitted Est.Error NA %d of %d",
               sum(is.na(r$v[, "Est.Error", ])), length(r$v[, "Est.Error", ])), r)
  })
  cat("  se.fit first row:", sprintf("%.12g", head(a$v$se.fit[1, ], 3)), "\n")
}

cat("\n== R5 a covariance held as a factor: the dense route\n")
fr <- q(frm(cntp ~ poly(u, 6, raw = TRUE) + ar(t, g, cov = TRUE), poisson(),
            data = d))
sdr <- frmtmb:::sdr_of(fr)
cat("  cov_fixed_root", !is.null(sdr$cov_fixed_root), "sparse", sparse_of(fr),
    "\n")
r <- catch(fitted(fr))
cat("  fitted Est.Error:", sprintf("%.12g", r$v[c(1, 100, 300), "Est.Error"]),
    "NA", sum(is.na(r$v[, "Est.Error"])), "\n")
r <- catch(frm_linpred(fr, se.fit = TRUE))
cat("  linpred se:", sprintf("%.12g", r$v$se.fit[c(1, 100, 300)]), "\n")
local({
  op <- options(frmtmb.latres_se_rows = 10L, frmtmb.latres_dense_rows = 10L)
  on.exit(options(op))
  fr2 <- fr
  fr2$cache <- list2env(mget(setdiff(ls(fr$cache), "Vjoint"), envir = fr$cache))
  r <- catch(fitted(fr2))
  pr(sprintf("limits 10: fitted Est.Error NA %d of %d",
             sum(is.na(r$v[, "Est.Error"])), nrow(d)), r)
  r <- catch(fitted(fr2, newdata = nd))
  pr(sprintf("limits 10: fitted(newdata) Est.Error NA %d",
             sum(is.na(r$v[, "Est.Error"]))), r)
})

cat("\n== R6 frm_joint_cov() transient memory against the 8 p^2 it returns\n")
for (GG in c(100, 250)) {
  set.seed(12)
  d2 <- expand.grid(t = seq_len(K), g = factor(seq_len(GG)))
  d2$x <- rnorm(nrow(d2))
  e2 <- as.vector(replicate(GG, arima.sim(list(ar = 0.6), K, sd = 0.6)))
  d2$cnt <- rpois(nrow(d2), exp(0.5 + 0.4 * d2$x + e2))
  f2 <- q(frm(cnt ~ x + ar(t, g, cov = TRUE), poisson(), data = d2))
  invisible(frmtmb:::get_joint_cov(f2))
  invisible(gc(reset = TRUE))
  g0 <- sum(gc()[, 2])
  jc <- frm_joint_cov(f2)
  g1 <- sum(gc()[, 6])
  p <- nrow(jc$V)
  cat(sprintf("  p %d: returned V %.1f MB, peak heap above start %.1f MB (%.1fx)\n",
              p, 8 * p^2 / 2^20, g1 - g0, (g1 - g0) / (8 * p^2 / 2^20)))
}
