# Punch round 1: where the hand closed form and frm_extra_cov_deriv()
# part by 5e-6 relative.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(5)
x <- sort(stats::runif(50, 0, 4))
dg <- data.frame(x = x, y = sin(1.3 * x) + stats::rnorm(50, 0, 0.15))
fg <- frm(bf(y ~ gp(x)), family = gaussian(), data = dg)
th <- fg$estimates$theta
gi <- fg$frame$linpreds[["y.mu"]]$gps[[1]]
bk <- fg$frame$re_blocks[[gi$block_id]]
str(gi[setdiff(names(gi), "positions")], max.level = 1)
pos <- gi$positions[, 1]
cat("npos", length(pos), "range", range(pos), "theta", th, "\n")
l <- exp(th[2])
K <- exp(-outer(pos, pos, "-")^2 / (2 * l^2)) + diag(1e-6, length(pos))
thc <- th; thc[1] <- 0
K2 <- unname(ns$covstruct_registry[["gp"]]$vcov(thc, bk))
cat("max |K - K2|", max(abs(K - K2)), "\n")
print(K2[1:3, 1:3] - K[1:3, 1:3])
xs <- c(3.0, 4.2, 4.8, 5.5)
r <- outer(xs, pos, "-")
s2 <- exp(2 * th[1])
k0 <- exp(-r^2 / (2 * l^2))
k1 <- -r / l^2 * k0
v1 <- s2 * (1 / l^2 - rowSums((k1 %*% solve(K)) * k1))
nd <- data.frame(x = xs)
e1 <- diag(frm_extra_cov_deriv(fg, nd, var = "x", order = 1))
e1a <- diag(frm_extra_cov_deriv(fg, nd, var = "x", order = 1,
                                re_formula = NA))
env <- fg$spec$responses[[1]]$formula_env
e1b <- diag(ns$gp_krig_deriv_cov(fg, gi, nd, env, "x", 1L))
print(rbind(v1, e1, e1a, e1b, e1 / v1 - 1))
