# Short check, item 6: fits without mi() on a registered term, lane against
# base, bitwise. Run each arm pinned to one core (dev/rtmb-pitfalls.md 21):
#   Rscript dev/reviews/arspat2-short-ident.R <lane|base>
# The base arm compares against the lane's rds.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-arspat2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
outdir <- "dev/reviews/arspat2-short-out"
rook <- function(k) {
  g <- expand.grid(r = seq_len(k), c = seq_len(k))
  n <- nrow(g)
  outer(seq_len(n), seq_len(n), function(i, j) {
    as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
  })
}
capture <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(w) {
      msgs <<- c(msgs, paste("W:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, paste("M:", conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  list(val = val, msgs = msgs)
}
rec <- function(fit, msgs) {
  p <- fit$opt$par
  pts <- lapply(0:2, function(k) p + 0.05 * k * seq_along(p) / length(p))
  list(par = p, ll = as.numeric(logLik(fit)),
       cov = frmtmb:::sdr_of(fit)$cov.fixed,
       fn = vapply(pts, function(q) fit$obj$fn(q), 0),
       gr = lapply(pts, function(q) as.vector(fit$obj$gr(q))),
       fixef = fixef(fit), msgs = msgs,
       cv = tryCatch(suppressWarnings(confint_varcorr(fit)),
                     error = function(e) conditionMessage(e)))
}
n <- 60
M <- 0.5^abs(outer(1:n, 1:n, "-"))
set.seed(4)
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)), z = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$h] +
  as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
W <- rook(7)
W <- W / rowSums(W)
ns_ <- nrow(W)
set.seed(21)
ds <- data.frame(x = rnorm(ns_), h = factor(rep(1:7, length.out = ns_)))
ds$y <- as.vector(solve(diag(ns_) - 0.6 * W, 1 + 0.5 * ds$x +
                          rnorm(7)[ds$h] + rnorm(ns_)))
dm <- d
dm$y[c(3, 17, 40)] <- NA
dx <- d
dx$z[c(5, 12, 33)] <- NA
R <- list()
c1 <- capture(frm(y ~ x + fcor(M), data = d, data2 = list(M = M)))
R$fcor_gauss <- rec(c1$val, c1$msgs)
c2 <- capture(frm(y ~ x + (1 | h) + fcor(M), data = d, data2 = list(M = M),
                  family = student()))
R$fcor_student_re <- rec(c2$val, c2$msgs)
c3 <- capture(frm(y ~ x + sar(W), data = ds, data2 = list(W = W)))
R$sar_lag_gauss <- rec(c3$val, c3$msgs)
c4 <- capture(frm(y ~ x + (1 | h) + sar(W, type = "error"), data = ds,
                  data2 = list(W = W), family = student()))
R$sar_error_student_re <- rec(c4$val, c4$msgs)
c5 <- capture(frm(y | mi() ~ x + (1 | h), data = dm, REML = TRUE))
R$reml_gauss_mi <- rec(c5$val, c5$msgs)
c6 <- capture(frm(bf(y ~ x + fcor(M)) + bf(z | mi() ~ x) + set_rescor(FALSE),
                  data = dx, data2 = list(M = M)))
R$mv_fcor_plus_mi_other <- rec(c6$val, c6$msgs)
saveRDS(R, file.path(outdir, sprintf("ident-%s.rds", arm)))
if (arm == "base") {
  L <- readRDS(file.path(outdir, "ident-lane.rds"))
  for (k in names(R)) {
    a <- L[[k]]
    b <- R[[k]]
    same <- vapply(c("par", "ll", "cov", "fn", "gr", "fixef", "msgs"),
                   function(f) identical(a[[f]], b[[f]]), NA)
    maxd <- max(abs(unlist(a[c("par", "ll", "fn")]) -
                      unlist(b[c("par", "ll", "fn")])))
    cvs <- identical(a$cv, b$cv)
    cat(sprintf("%-24s identical: %s | max |diff| par/ll/fn %.1e | confint_varcorr identical %s%s | msgs %d\n",
                k, paste(names(same), same, sep = "=", collapse = " "), maxd,
                cvs, if (!cvs && is.data.frame(a$cv)) sprintf(
                  " (lane type %s vs %s; lwr/upr %s / %s vs %s / %s)",
                  a$cv$type[1], b$cv$type[1], format(a$cv$lwr[1]),
                  format(a$cv$upr[1]), format(b$cv$lwr[1]),
                  format(b$cv$upr[1])) else "", length(a$msgs)))
  }
}
