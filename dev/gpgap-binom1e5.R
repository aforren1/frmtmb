# The review's M2 row it could not read: binomial at 1e5 trials, seed 5601, frm()
# 0.85 nats ABOVE the non-centered reference at brms's line. Same RNG path as
# dev/review-gpgap-silent.R. Is the reference short, or is the model above it?
#   Rscript dev/gpgap-binom1e5.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
put <- function(val) {
  unlockBinding("gp_corr", ns)
  assign("gp_corr", val, envir = ns)
  lockBinding("gp_corr", ns)
}
at <- function(nu, abs) {
  function(theta, blk) {
    Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
    I <- diag(1, nrow(Q))
    exp(-Q) + I * (if (abs) nu * exp(-2 * theta[1]) else nu)
  }
}
set.seed(5601)
n <- 60
x <- round(runif(n, 0, 6), 3)
e <- rnorm(n)
d <- data.frame(x = x, z = round(runif(n, 0, 6), 3),
                g = factor(rep(1:6, length.out = n)))
for (lm in c(4, 8, 12, 16)) d$cnt <- rpois(n, exp(lm + sin(2 * d$x)))
for (tr in c(1e3, 1e5)) {
  d$tr <- tr
  d$suc <- rbinom(n, d$tr, plogis(sin(2 * d$x)))
}
fr <- frm(bf(suc | trials(tr) ~ gp(x)), family = binomial(), data = d)
bk <- fr$frame$re_blocks[[1]]
th <- fr$estimates$theta[bk$theta_idx]
cat(sprintf("frm(): sdgp %.4f rho %.4f ll %.4f\n", exp(th[1]), exp(th[2]),
            -fr$opt$objective))
frn <- fr$frame
frn$ncp_blocks <- 1L
run <- function(lab, tpl) {
  fn <- frmtmb::build_objective(frn)
  o <- RTMB::MakeADFun(fn, tpl, random = "b", map = frn$map, silent = TRUE)
  op <- nlminb(o$par, o$fn, o$gr, control = list(iter.max = 2000, eval.max = 3000))
  tt <- op$par[names(op$par) == "theta"]
  cat(sprintf("%-34s ll %.4f sdgp %.4f rho %.4f conv %d gr %.1e\n", lab, -op$objective,
              exp(tt[1]), exp(tt[2]), op$convergence, max(abs(o$gr(op$par)))))
  o
}
tp <- frn$par_template
tp$b[] <- 0
for (nu in list(list(1e-12, TRUE), list(1e-12, FALSE), list(1e-9, FALSE),
                list(1e-6, FALSE))) {
  put(at(nu[[1]], nu[[2]]))
  run(sprintf("default start, %g %s", nu[[1]], if (nu[[2]]) "abs" else "rel"), tp)
  t2 <- tp
  t2$theta <- fr$estimates$theta
  t2$beta <- fr$estimates$beta
  run(sprintf("from frm(), %g %s", nu[[1]], if (nu[[2]]) "abs" else "rel"), t2)
  t3 <- t2
  t3$theta[bk$theta_idx] <- log(c(2.0, 1.0))
  run(sprintf("sd 2 rho 1, %g %s", nu[[1]], if (nu[[2]]) "abs" else "rel"), t3)
}
put(orig)
# one-row-per-position residual on the link: how far the data sit from the smooth
# mean, against the binomial noise, tells whether a white noise term pays
p <- d$suc / d$tr
eta <- qlogis(p)
w <- d$tr * p * (1 - p)
f <- as.numeric(fr$estimates$b[bk$b_idx])[match(d$x, sort(unique(d$x)))] +
  fr$estimates$beta[1]
cat(sprintf("z of the residual at frm(): sd %.3f (1 if binomial alone)\n",
            sd((eta - f) * sqrt(w))))
cat(sprintf("nugget variance %.3g, median working noise %.3g\n",
            1e-6 * exp(2 * th[1]), median(1 / w)))
