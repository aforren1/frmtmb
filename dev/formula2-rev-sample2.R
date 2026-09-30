# Reviewer: does sigma1 carry sigma2's value on every frm_sample() draw?
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
q <- function(expr) suppressWarnings(suppressMessages(expr))
fe <- q(frm(bf(y ~ x, sigma1 = "sigma2"),
            family = mixture(gaussian(), gaussian()), data = d))
s <- q(frm_sample(fe, chains = 2, iter = 300, warmup = 150, refresh = 0,
                  seed = 3))
dr <- posterior::as_draws_matrix(s)
a <- as.numeric(dr[, "sigma1"]); b <- as.numeric(dr[, "sigma2"])
cat("n draws", length(a), "all ==", all(a == b), "identical numeric",
    identical(a, b), "\n")
str(dr[1:2, c("sigma1", "sigma2")])
da <- posterior::as_draws_array(s)
cat("array: identical per chain",
    identical(as.numeric(da[, , "sigma1"]), as.numeric(da[, , "sigma2"])),
    "\n")
dl <- posterior::as_draws_df(s)
cat("df:", identical(dl$sigma1, dl$sigma2), "\n")
# the inverse direction (natural -> sampler scale) must drop sigma1
nc <- frmtmb.sample:::draws_natural_cols(fe)
str(nc$equated); print(nc$extra)
m <- as.matrix(dr)[1:3, ]
inv <- frmtmb.sample:::draws_to_natural(m, fe, inverse = TRUE)
cat("inverse cols:", colnames(inv), "\n")
