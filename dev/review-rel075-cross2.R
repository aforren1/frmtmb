# Review: more cross-lane probes (Dirichlet class sbhaz off a cox model;
# sar() with a group-level exp() term; cox + bhaz prior table order with
# group-level terms against brms).
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spatial); library(frmtmb.survival)})
show <- function(label, expr) {
  cat("\n=====", label, "\n")
  r <- tryCatch(expr, error = function(e) {
    cat("ERROR [", class(e)[1], "]:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) print(r)
  invisible(r)
}
ptab <- function(p) as.data.frame(p)[, c("prior", "class", "coef", "group", "dpar", "source")]
set.seed(3)
d <- data.frame(x = rnorm(40), y = rnorm(40))
show("dirichlet on sbhaz in a gaussian model (frmtmb)",
     validate_prior(set_prior("dirichlet(2)", class = "sbhaz"), y ~ x, data = d))
show("dirichlet on sbhaz in a gaussian model (brms)",
     ptab(brms::validate_prior(brms::set_prior("dirichlet(2)", class = "sbhaz"), y ~ x, data = d)))
show("normal on sbhaz in a gaussian model (frmtmb)",
     validate_prior(set_prior("normal(0, 1)", class = "sbhaz"), y ~ x, data = d))
k <- 6
g <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A); n <- nrow(W)
set.seed(20261024)
d2 <- data.frame(x = runif(n, 0, 1.5), z = rnorm(n), g = factor(rep(1:6, 6)))
d2$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * exp(d2$x) + rnorm(6)[d2$g] + rnorm(n, 0, 0.5)))
f <- y ~ exp(x) + (1 + exp(x) | g) + sar(W)
show("default_prior exp + RE exp + sar (frmtmb)", ptab(default_prior(f, data = d2, data2 = list(W = W))))
show("default_prior exp + RE exp + sar (brms)", ptab(brms::default_prior(f, data = d2, data2 = list(W = W))))
fit <- show("frm", frm(f, data = d2, data2 = list(W = W)))
if (!is.null(fit)) {
  show("variables", variables(fit))
  show("fitted re_formula = ~ (1 + exp(x) | g)", head(fitted(fit, re_formula = ~ (1 + exp(x) | g))))
}
set.seed(20261024)
ds <- data.frame(x = runif(150, 0, 1.5), z = rnorm(150), g = factor(rep(1:10, 15)))
ds$time <- rexp(150, exp(-1 + 0.5 * exp(ds$x) + rnorm(10, 0, .3)[ds$g]))
ds$cens <- as.numeric(runif(150) < 0.2)
fs <- time | cens(cens) + bhaz(df = 4) ~ exp(x) + z + (1 | g)
show("cox + bhaz + RE default_prior (frmtmb)", ptab(default_prior(fs, data = ds, family = cox())))
show("cox + bhaz + RE default_prior (brms)", ptab(brms::default_prior(fs, data = ds, family = brms::cox())))
