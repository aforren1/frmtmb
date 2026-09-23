# Why sn_sigma_noint and sn_mu_noint still move: the residual START of
# an INTERCEPT dpar (the lane's purpose), or a placement into an
# intercept-less design (the blocker)? The starts decide it.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
set.seed(23)
n <- 250
g <- factor(rep(1:4, length.out = n))
x <- rnorm(n)
dd <- data.frame(x = x, g = g,
                 ysn = 0.4 + 0.8 * x + 1.2 * (0.9486833 * abs(rnorm(n)) +
                                                0.3162278 * rnorm(n)))
show <- function(lab, fo) {
  u <- frm(fo, family = skew_normal(), data = dd, dry_run = "objective")
  cat(sprintf("%-18s beta %s\n%-18s betad %s\n", lab,
              paste(format(u$estimates$beta, digits = 5), collapse = " "),
              "", paste(paste0(names(u$estimates$betad), "=",
                               format(u$estimates$betad, digits = 5)),
                        collapse = " ")))
}
show("sn_alpha_noint", bf(ysn ~ x, sigma ~ 1, alpha ~ 0 + g))
show("sn_sigma_noint", bf(ysn ~ x, sigma ~ 0 + g, alpha ~ 1))
show("sn_mu_noint", bf(ysn ~ 0 + g, sigma ~ 1, alpha ~ 1))
cat("\nA sigma ~ 0 + g start of all zeros means sigma was SKIPPED,",
    "\nwhich is the pre-lane rule; only alpha may be placed.\n")
