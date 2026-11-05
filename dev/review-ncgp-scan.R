# Item 1 and 2: what gp_start_scan() reads. Unit sd and the start sigma, whatever y's
# units; the objective at each trial rho, on y and on y * 1000. Seeds below.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
show <- function(lab, x, y) {
  d <- data.frame(y = y, x = x)
  f <- suppressWarnings(frm(y ~ gp(x), d))
  bk <- f$frame$re_blocks[[1]]
  dm <- sqrt(max(Reduce(`+`, bk$aux_D2)))
  tp <- f$frame$par_template
  o <- f$obj
  cat(sprintf("RSC %s: start theta %s betad %.3f beta %.3f; fit sd %.4g rho %.4g sigma %.4g ll %.4f conv %d\n", lab, paste(round(tp$theta, 3), collapse = ","), tp$betad[1], tp$beta[1],
              exp(f$estimates$theta[1]), exp(f$estimates$theta[2]), exp(f$estimates$betad[1]), logLik(f), f$opt$convergence))
  tv <- function(th) { t2 <- tp; t2$theta <- th; p <- ns$outer_from_template(t2, o, f$frame); o$fn(p) }
  base <- tv(tp$theta)
  cat("RSC   base start objective", base, "\n")
  for (r in dm / 4^(0:3)) cat(sprintf("RSC   rho %.4g objective %.5f\n", r, tv(c(tp$theta[1], log(r)))))
  sds <- sd(y)
  for (r in dm / 4^(0:3)) cat(sprintf("RSC   rho %.4g sd=sd(y) objective %.5f\n", r, tv(c(log(sds), log(r)))))
}
set.seed(1080); x <- sort(round(runif(80, 0, 10), 4)); e <- rnorm(80)
show("mix s=1e-2", x, 2 * sin(3 * x) + x + 1e-2 * e)
set.seed(4140); x <- sort(runif(40, 0, 10)); e <- rnorm(40)
show("sin6 y", x, sin(6 * x) + 1e-3 * e)
show("sin6 y*1000", x, 1000 * (sin(6 * x) + 1e-3 * e))
show("sin6 y/1000", x, 1e-3 * (sin(6 * x) + 1e-3 * e))
cat("DONE\n")
