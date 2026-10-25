# Fix round 1 probe: the map check, class priors, same-row newdata,
# the complete-graph warning, the autocor= refusal. Seed 2.
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat("---", lab, "\n"); print(r)
}
set.seed(2)
g <- expand.grid(r = 1:7, c = 1:7); n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
f <- frm(y ~ x + sar(W), data = d, data2 = list(W = W))
show("pars", frmtmb_autocor_pars(f))
show("default_prior", default_prior(y ~ x + sar(W), data = d,
                                    data2 = list(W = W)))
fp <- frm(y ~ x + sar(W), data = d, data2 = list(W = W),
          prior = set_prior("normal(0, 0.2)", class = "lagsar"))
show("lagsar with normal(0, 0.2) prior", fp$estimates$thetaac)
show("lb on lagsar", frm(y ~ x + sar(W), data = d, data2 = list(W = W),
                         prior = set_prior("", class = "lagsar", lb = 0.3,
                                           ub = 0.9))$estimates)
nd <- d; nd$x <- nd$x + 1
show("fitted on same rows newdata, minus A^-1 eta_new", {
  rho <- unname(frmtmb:::autocor_natural(f$estimates$thetaac,
                                         f$frame$autocor[[1]]))
  max(abs(fitted(f, newdata = nd)[, 1] -
            solve(diag(n) - rho * W, frm_linpred(f, newdata = nd))))
})
show("simulate same rows", dim(simulate(f, newdata = nd, nsim = 2)))
show("predict same rows", dim(predict(f, newdata = nd)))
show("fitted other rows", fitted(f, newdata = nd[1:5, ]))
K <- (matrix(1, 9, 9) - diag(9)) / 8
dk <- data.frame(x = rnorm(9), y = rnorm(9))
show("complete graph", frm(y ~ x + sar(K), data = dk, data2 = list(K = K)))
