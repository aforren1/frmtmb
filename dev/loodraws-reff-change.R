# The size of the r_eff change on loo(ndraws =): 0.20.0 used r_eff = 1
# on a subset of draws, this build uses relative_eff() on one chain, as
# brms does. Same draws, same rows, both rules. Data seed 1 (60 rows);
# frm_sample(chains = 4, iter = 2000, warmup = 1000, seed = 1).
#   Rscript dev/loodraws-reff-change.R -> dev/loodraws-log/reff-change.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 6)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.7)[dd$g])
for (f in list(bf(y ~ x), bf(y ~ x + (1 | g)))) {
  ds <- q(frm_sample(f, family = gaussian(), data = dd, chains = 4,
                     iter = 2000, warmup = 1000, seed = 1, refresh = 0))
  for (nd in c(500, 1000, 3000)) {
    new <- q(loo(ds, ndraws = nd))
    ll <- log_lik(ds, ndraws = nd)
    old <- q(loo::loo(ll, r_eff = 1))
    one <- q(loo::relative_eff(exp(ll), chain_id = rep(1L, nrow(ll))))
    cat(sprintf(paste0("%-18s ndraws %4d: elpd_loo old %.6f new %.6f ",
                       "diff %.3g (se %.3g); p_loo diff %.3g; ",
                       "r_eff range %.3f %.3f; new == one-chain: %s\n"),
                deparse(f[[1]]), nd, old$estimates[1, 1],
                new$estimates[1, 1],
                new$estimates[1, 1] - old$estimates[1, 1],
                new$estimates[1, 2],
                new$estimates[2, 1] - old$estimates[2, 1],
                min(one), max(one),
                identical(new$estimates,
                          q(loo::loo(ll, r_eff = one))$estimates)))
  }
}
