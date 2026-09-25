source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# Downstream of an engaged fit at a benign slope scale (sd 0.03, where
# autoscale = FALSE converges): profile confint, predict se.fit,
# bootstrap refits, update(), and a multivariate ID-linked slope.
# Seed 31.
set.seed(31)
ng <- 20; per <- 15; n <- ng * per
g <- factor(rep(seq_len(ng), each = per))
xs <- rnorm(n)
y <- 1 + 0.7 * xs + rnorm(ng, 0, 0.8)[g] + rnorm(ng, 0, 0.5)[g] * xs + rnorm(n)
y2 <- -1 + 0.3 * xs + rnorm(ng, 0, 0.6)[g] * xs + rnorm(n)
d <- data.frame(y, y2, x = xs * 0.03, g)
ctlF <- frmtmb_control(autoscale = FALSE)
a <- frm(y ~ x + (1 + x | g), data = d)
b <- frm(y ~ x + (1 + x | g), data = d, control = ctlF)
cat("engaged:", !is.null(a$par_units), "/", !is.null(b$par_units), "\n")
cat(sprintf("logLik diff %.3e\n", as.numeric(logLik(a) - logLik(b))))
ca <- tryCatch(confint(a, method = "profile"), error = function(e) e)
cb <- tryCatch(confint(b, method = "profile"), error = function(e) e)
if (inherits(ca, "error") || inherits(cb, "error")) {
  cat("profile confint ERROR:", if (inherits(ca, "error")) conditionMessage(ca),
      "|", if (inherits(cb, "error")) conditionMessage(cb), "\n")
} else {
  cat("profile confint default:\n"); print(signif(ca, 6))
  cat(sprintf("profile confint max rel diff default vs FALSE: %.2e\n", rel(ca, cb)))
}
nd <- data.frame(x = c(-0.05, 0, 0.05), g = factor(c("1", "2", "zz")))
pa <- frm_linpred(a, newdata = nd, se.fit = TRUE, allow_new_levels = TRUE)
pb <- frm_linpred(b, newdata = nd, se.fit = TRUE, allow_new_levels = TRUE)
cat(sprintf("frm_linpred fit rel %.2e se.fit rel %.2e\n",
            rel(pa$fit, pb$fit), rel(pa$se.fit, pb$se.fit)))
ua <- update(a, . ~ . + 0 + 1)
cat("update() keeps engagement:", !is.null(ua$par_units),
    sprintf("logLik diff to a %.2e\n", as.numeric(logLik(ua) - logLik(a))))
bt <- tryCatch(frm_bootstrap(a, nsim = 5, seed = 2), error = function(e) e)
bf_ <- tryCatch(frm_bootstrap(b, nsim = 5, seed = 2), error = function(e) e)
if (inherits(bt, "error") || inherits(bf_, "error")) {
  cat("bootstrap ERROR:", if (inherits(bt, "error")) conditionMessage(bt), "|",
      if (inherits(bf_, "error")) conditionMessage(bf_), "\n")
} else {
  ta <- bt$t %||% bt[["t"]]; tb <- bf_$t %||% bf_[["t"]]
  cat(sprintf("bootstrap replicate estimates max rel diff: %.2e\n", rel(ta, tb)))
}

cat("\n-- multivariate, slope ID-linked across two responses, x at 1e-3 --\n")
d3 <- d; d3$x <- xs * 1e-3
dref <- d; dref$x <- xs
mvf <- mvbf(bf(y ~ x + (1 + x | p | g)), bf(y2 ~ x + (1 + x | p | g)),
            rescor = FALSE)
r_ref <- frm(mvf, data = dref, control = ctlF)
r_def <- fitw(mvf, data = d3)
r_F <- frm(mvf, data = d3, control = ctlF)
cat(sprintf("ref %.6f c%d | default %.6f c%d engaged %s | FALSE %.6f c%d\n",
            as.numeric(logLik(r_ref)), r_ref$opt$convergence,
            as.numeric(logLik(r_def$fit)), r_def$fit$opt$convergence,
            !is.null(r_def$tpl), as.numeric(logLik(r_F)), r_F$opt$convergence))
pl <- frmtmb:::autoscale_plan(r_def$fit$frame)
for (k in names(pl)) cat("  ", k, "z thetas:",
                         vapply(pl[[k]]$z, function(e) e$theta, 1L), "\n")
cat("  template move:", tpl_move(r_def), "\n")
