# Reviewer round 2: extra_cov's three forms (empty sparse, sparse, dense)
# through every consumer. For each form: the class frm_lp_basis()
# returns, and frm_curve() (simultaneous), a difference, a derivative, a
# feature, emmeans and fitted() run without error and agree with a dense
# reference built here (as.matrix of the returned extra_cov, or the same
# unseen rows asked alone, where the dense path is taken).
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
set.seed(5)
d <- data.frame(x = round(runif(60, 0, 6), 1),
                f = factor(rep(c("a", "b", "c"), 20)))
d$y <- sin(d$x) + as.numeric(d$f) * 0.3 + rnorm(60, 0, 0.3)
d$yo <- factor(cut(d$y, 3, labels = FALSE), ordered = TRUE)
obs <- sort(unique(d$x))
fit <- frm(bf(y ~ f + gp(x)), data = d)
fs <- frm(bf(y ~ f + s(x, k = 6)), data = d)
fc <- frm(bf(y ~ gp(x, by = f, cmc = FALSE)), data = d)
grids <- list(
  empty = list(fs, data.frame(x = seq(0.5, 5.5, length.out = 40),
                              f = factor("a", levels = levels(d$f)))),
  sparse = list(fit, data.frame(x = c(obs[1:36], 6.3, 6.6, 7, 7.4),
                                f = factor("a", levels = levels(d$f)))),
  dense = list(fit, data.frame(x = seq(6.1, 8, length.out = 40),
                               f = factor("a", levels = levels(d$f)))),
  sparse_cmcF = list(fc, data.frame(x = c(obs[1:36], 6.3, 6.6, 7, 7.4),
                                    f = factor(rep(c("a", "b"), 20),
                                               levels = levels(d$f)))))
try_msg <- function(expr) tryCatch({force(expr); "ok"},
                                   error = function(e) conditionMessage(e))
for (nm in names(grids)) {
  m <- grids[[nm]][[1]]; nd <- grids[[nm]][[2]]
  lb <- frm_lp_basis(m, newdata = nd, extra_cov = TRUE)
  E <- lb$extra_cov
  cls <- class(E)[1]
  A <- as.matrix(lb$A)
  Sref <- A %*% lb$V %*% t(A) + as.matrix(E)
  cv <- frm_curve(m, newdata = nd, nsim = 2000, seed = 1)
  cat(sprintf(paste0("== %-11s extra_cov class %s, nnz %d of %d | ",
                     "frm_curve Sigma vs dense ref max rel %.2e | ",
                     "crit %.4f\n"),
              nm, cls, if (inherits(E, "Matrix")) length(E@x) else
                sum(E != 0), length(E),
              max(abs(attr(cv, "Sigma") - Sref)) / max(abs(Sref)),
              cv$.crit_sim[1]))
  unseen <- which(!nd$x %in% obs)
  if (length(unseen) && length(unseen) < nrow(nd)) {
    lb2 <- frm_lp_basis(m, newdata = nd[unseen, , drop = FALSE],
                        extra_cov = TRUE)
    cat(sprintf("   unseen rows asked alone (%s): max |E[u,u] - E_alone| %.2e\n",
                class(lb2$extra_cov)[1],
                max(abs(as.matrix(E)[unseen, unseen] -
                          as.matrix(lb2$extra_cov)))))
  }
  nd2 <- nd
  nd2$f <- factor("b", levels = levels(d$f))
  cat("   difference:", try_msg(frm_curve(m, newdata = nd, contrast = nd2,
                                          nsim = 1000, seed = 1)),
      "| deriv:", try_msg(frm_curve_deriv(m, var = "x", newdata = nd,
                                          nsim = 1000, seed = 1)),
      "| deriv o2:", try_msg(frm_curve_deriv(m, var = "x", order = 2,
                                             newdata = nd, nsim = 1000,
                                             seed = 1)),
      "| feature:", try_msg(frm_curve_feature(m, var = "x",
                                              type = "extremum",
                                              newdata = nd)),
      "\n")
}
# emmeans and fitted() on the sparse and dense forms
em_at <- function(xv) {
  s <- as.data.frame(emmeans::emmeans(fit, "f", at = list(x = xv)))
  p <- frm_linpred(fit, newdata = data.frame(f = factor(levels(d$f),
                                                        levels = levels(d$f)),
                                             x = xv), se.fit = TRUE)
  max(abs(s$SE / as.numeric(p$se.fit) - 1))
}
cat(sprintf("emmeans SE vs frm_linpred se.fit, x observed %.2e, x unseen %.2e\n",
            em_at(obs[3]), em_at(6.7)))
em2 <- tryCatch(as.data.frame(emmeans::emmeans(fit, ~ x, at = list(
  x = c(obs[1:30], 6.4, 7)))), error = function(e) conditionMessage(e))
cat("emmeans over a mostly observed x grid:",
    if (is.character(em2)) em2 else sprintf("ok, %d rows, SE finite %s",
                                            nrow(em2), all(is.finite(em2$SE))),
    "\n")
fo <- suppressWarnings(frm(bf(yo ~ gp(x)), data = d, family = cumulative()))
ndo <- data.frame(x = c(obs[1:10], 6.5, 7))
fv <- tryCatch(fitted(fo, newdata = ndo)[, "Est.Error", ],
               error = function(e) conditionMessage(e))
cat("fitted() cumulative Est.Error on a mixed grid:",
    if (is.character(fv)) fv else sprintf("finite %s, rows past data > inside %s",
                                          all(is.finite(fv)),
                                          all(fv[11:12, 2] > max(fv[1:10, 2]))),
    "\n")
cat("DONE\n")
