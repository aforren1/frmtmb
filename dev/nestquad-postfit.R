# Post-fit methods on a nested quadrature fit: profile intervals (which
# refit with a map), update(), predict() on new data, simulate().
# Test design, seed 20261020.
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n), w = runif(n, 0.5, 2))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$y <- rbinom(n, 1, plogis(0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
                             rnorm(nlevels(ab), 0, 1)[ab]))
fq <- frm(bf(y ~ 1 + x + (1 | ga/gb)), family = bernoulli(), data = d,
          quadrature = TRUE)
run <- function(tag, expr) {
  r <- tryCatch(expr, error = function(e) e)
  cat("POSTFIT", tag, if (inherits(r, "error")) paste("ERROR",
      conditionMessage(r)) else "ok", "\n")
  r
}
w <- run("wald", confint(fq))
print(w)
p <- run("profile", confint(fq, method = "profile", parm = "x"))
print(p)
u <- run("update", update(fq, control = frmtmb_control(quad_nodes = 21)))
if (!inherits(u, "error")) {
  cat(sprintf("POSTFIT loglik default %.8f nodes21 %.8f\n",
              as.numeric(logLik(fq)), as.numeric(logLik(u))))
}
pr <- run("predict", predict(fq, newdata = d[1:5, ]))
print(pr)
s <- run("simulate", simulate(fq, nsim = 2, seed = 1))
r <- run("ranef", ranef(fq))
