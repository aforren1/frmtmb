# How rellib-r8 treats a matrix response (multinomial) and categorical on
# every post-fit route. Output: dev/famsimplex-explore-multi.log
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(2)
n <- 200
dd <- data.frame(x = rnorm(n), g = gl(10, 20))
eta <- cbind(0, 0.3 + 0.5 * dd$x, -0.2 - 0.4 * dd$x)
P <- exp(eta) / rowSums(exp(eta))
Y <- t(apply(P, 1, function(p) rmultinom(1, 20, p)))
colnames(Y) <- c("a", "b", "c")
dd$Y <- Y
dd$tr <- 20
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1) cat(r, "\n") else {
    cat("class:", class(r), " dim:", dim(r), "\n")
    if (!is.null(dimnames(r))) str(dimnames(r))
  }
  invisible(r)
}
fit <- frm(bf(Y | trials(tr) ~ x + (1 | g)), family = multinomial(K = 3),
           data = dd)
print(fit)
show("fitted", fitted(fit))
show("fitted summary=FALSE", fitted(fit, summary = FALSE))
show("predict", predict(fit))
show("posterior_epred", posterior_epred(fit, ndraws = 5))
show("posterior_predict", posterior_predict(fit, ndraws = 5))
show("posterior_linpred", posterior_linpred(fit, ndraws = 5))
show("log_lik", log_lik(fit, ndraws = 5))
show("residuals", residuals(fit))
show("ce", conditional_effects(fit, "x"))
show("ce cat", conditional_effects(fit, "x", categorical = TRUE))
show("default_prior", default_prior(fit))
show("simulate", simulate(fit, nsim = 2))
show("pp_check", pp_check(fit, ndraws = 5))
dd$f <- factor(sample(c("a", "b", "c"), n, TRUE))
fit2 <- frm(bf(f ~ x + (1 | ID | g)), family = categorical(), data = dd)
print(fit2)
show("cat fitted", fitted(fit2))
show("cat posterior_epred", posterior_epred(fit2, ndraws = 5))
show("cat ce", conditional_effects(fit2, "x", categorical = TRUE))
show("cat predict", predict(fit2))
show("cat posterior_predict", posterior_predict(fit2, ndraws = 5))
show("cat log_lik", log_lik(fit2, ndraws = 5))
print(VarCorr(fit2))
