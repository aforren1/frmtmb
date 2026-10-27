# Reviewer: frm_sample()'s ar posterior (after) against brms's OWN
# posterior under brms's own defaults, unweighted: the distance that
# remains from the sigma parameterization (innovation vs marginal sd).
# Reads dev/acprior-out/agree-{after,brms}.rds (seed 3, 4 x 4000).
# Log: dev/review-acprior-log/unweighted.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
o <- "C:/Users/adf44/source/r/frmtmb-wt-acprior/dev/acprior-out"
A <- readRDS(file.path(o, "agree-after.rds"))$ar$par
B <- readRDS(file.path(o, "agree-brms.rds"))$ar$par
m <- function(x) c(mean(x), posterior::mcse_mean(matrix(x, ncol = 4)))
a <- m(A); b <- m(B)
cat(sprintf("after %.4f (mcse %.5f); brms unweighted %.4f (mcse %.5f); z %.2f\n",
            a[1], a[2], b[1], b[2], (a[1] - b[1]) / sqrt(a[2]^2 + b[2]^2)))
