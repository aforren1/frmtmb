# Short check, item 2: check_laplace() on a MAP fit, whose priors make the block
# qualify for the non-centered plan. Which nugget does the diagnostic run read, against
# the 1e-6 of the fit it checks? Seed 6601.
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib","C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(6601); n <- 30; x <- round(runif(n,0,6),3); e0 <- rnorm(n)
for (s in c(0.1, 1e-3)) {
  d <- data.frame(x=x,y=sin(2*x)+s*e0)
  f <- suppressWarnings(frm(bf(y~gp(x)),family=gaussian(),data=d,
    prior=c(set_prior("student_t(3,0,2.5)",class="sdgp"),set_prior("inv_gamma(2,1)",class="lscale"),set_prior("student_t(3,0,2.5)",class="sigma"))))
  ds <- suppressWarnings(frm_sample(f,chains=1,iter=30,warmup=15,refresh=0,seed=5,.diagnostic=TRUE))
  bk <- Filter(function(b) identical(b$covstruct,"gp"), ds$fit$frame$re_blocks)[[1]]
  cat(sprintf("DIAG2 noise %g fit sigma %.3g | block nugget on the diagnostic run: %s\n", s, as.numeric(sigma(f)), format(bk$gp_nugget %||% "NULL (1e-6)")))
  t0 <- proc.time()[["elapsed"]]
  cl <- suppressWarnings(check_laplace(f, chains=2, iter=600, refresh=0, seed=5))
  print(as.data.frame(cl)[,1:6], digits=3)
  cat(sprintf("DIAG2 check_laplace %.0f s\n", proc.time()[["elapsed"]]-t0))
}
