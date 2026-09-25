# How brms 2.23.0 counts lags of ar(time, gr, cov = TRUE) at newdata,
# against frmtmb on the same data. Worker's fixture (seed 8, 40 x 6).
#   Rscript dev/simnewdata-review/rv-brms-time.R > .../log/brms-time.txt
source("dev/simnewdata-review/rv-prelude.R")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(8)
n_g <- 40
d <- expand.grid(time = 1:6, g = factor(seq_len(n_g)))
d$y <- 1 + unlist(lapply(seq_len(n_g), function(i) {
  as.vector(stats::arima.sim(list(ar = 0.7), 6))
}))
fit <- brm(y ~ 1 + ar(time = time, gr = g, cov = TRUE), data = d,
           chains = 2, iter = 2000, refresh = 0, seed = 1, silent = 2)
phi <- mean(as_draws_matrix(fit)[, "ar[1]"])
cat(sprintf("brms posterior mean ar[1] %.3f; phi^2 %.3f\n", phi, phi^2))
ff <- frm(frmtmb::bf(y ~ 1 + ar(time, gr = g, cov = TRUE)), data = d)
Rf <- autocor_matrix(ff)
cat(sprintf("frmtmb R[1,2] %.3f R[1,3] %.3f\n", Rf[1, 2], Rf[1, 3]))
pc <- function(nd) {
  pp <- posterior_predict(fit, newdata = cbind(nd, y = 0),
                          allow_new_levels = TRUE)
  round(stats::cor(pp), 3)
}
fc <- function(nd) {
  r <- tryCatch({
    s <- as.matrix(simulate(ff, nsim = 4000, seed = 1, newdata = nd))
    round(stats::cor(t(s)), 3)
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
  r
}
show <- function(label, nd) {
  cat("\n==", label, ": time =", paste(nd$time, collapse = ", "),
      " g =", paste(nd$g, collapse = ", "), "\n-- brms\n")
  print(pc(nd))
  cat("-- frmtmb\n"); print(fc(nd))
}
show("A seen times 2, 3", data.frame(time = c(2, 3), g = factor("a")))
show("B seen times 2, 4 (a gap)", data.frame(time = c(2, 4), g = factor("a")))
show("C unseen 9 after 2", data.frame(time = c(2, 9), g = factor("a")))
show("D rows reversed 4, 2", data.frame(time = c(4, 2), g = factor("a")))
show("E 2, 9, 3 unsorted", data.frame(time = c(2, 9, 3), g = factor("a")))
show("F 1, 6 far apart", data.frame(time = c(1, 6), g = factor("a")))
