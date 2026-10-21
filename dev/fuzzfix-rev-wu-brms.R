# Reviewer: brms 2.23.0 on gaussian + weights + (1 | obs) with
# sum(w) > n (data of dev/fuzzfix-rev-wu.R, seed 20261021), and the
# sigma ~ 0 + xs model that check_weights_unbounded() refuses.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
cat(sprintf("sum(w) - n = %.4f\n", sum(d$w) - n))
show <- function(lbl, f) {
  s <- summary(f)
  rh <- max(c(s$fixed$Rhat, s$spec_pars$Rhat,
              unlist(lapply(s$random, function(r) r$Rhat))), na.rm = TRUE)
  dr <- as.matrix(f)
  cn <- grep("^sigma$|^b_sigma_xs$", colnames(dr), value = TRUE)[1]
  q <- quantile(dr[, cn], c(0.05, 0.5, 0.95))
  cat(sprintf("%-26s %s q05/50/95 %.4g %.4g %.4g max Rhat %.3f divergent %d\n",
              lbl, cn, q[1], q[2], q[3], rh,
              sum(nuts_params(f)$Value[nuts_params(f)$Parameter ==
                                         "divergent__"])))
}
f1 <- brm(y | weights(w) ~ x + (1 | obs), data = d, chains = 2, iter = 1000,
          seed = 1, refresh = 0, silent = 2)
show("sum(w) > n", f1)
f2 <- brm(bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + xs), data = d,
          chains = 2, iter = 1000, seed = 1, refresh = 0, silent = 2)
show("sigma ~ 0 + xs", f2)
print(fixef(f2))
