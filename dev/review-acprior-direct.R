# Reviewer: replace the lane's importance weighting by a direct fit.
# brms 2.23.0's own Stan code for ar(cov = TRUE) on the data of
# dev/acprior-agree.R, with ONE edit: the student_t(3, 0, 2.5) prior
# is put on the marginal sd sigma_m = sigma / sqrt(1 - ar^2) (frmtmb's
# parameter) with the change-of-variables term -log(sqrt(1 - ar^2)),
# instead of on brms's innovation sd. Then brms's chains target exactly
# frmtmb.sample's model and no weights are needed.
# Usage: Rscript dev/review-acprior-direct.R
# Log: dev/review-acprior-log/direct.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "rstan",
    format(packageVersion("rstan")), "\n")
set.seed(20261027)
ng <- 8; k <- 6
da <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                 x = rnorm(ng * k))
e <- numeric(nrow(da))
for (i in seq_len(nrow(da))) {
  e[i] <- if (da$t[i] == 1) rnorm(1) else 0.9 * e[i - 1] +
    rnorm(1, 0, sqrt(1 - 0.81))
}
da$y <- 1 + 0.5 * da$x + e
fo <- y ~ x + ar(time = t, gr = g, cov = TRUE)
print(as.data.frame(default_prior(fo, data = da))[, c("prior", "class",
                                                       "lb", "ub")])
code <- stancode(fo, data = da)
old <- "lprior += student_t_lpdf(sigma | 3, 0, 2.5)"
stopifnot(sum(grepl(old, strsplit(code, "\n")[[1]], fixed = TRUE)) == 1L)
new <- paste0("lprior += student_t_lpdf(sigma / sqrt(1 - square(ar[1])) ",
              "| 3, 0, 2.5) - 0.5 * log1m(square(ar[1]))")
code2 <- sub(old, new, code, fixed = TRUE)
cat(grep("lprior \\+=", strsplit(code2, "\n")[[1]], value = TRUE), sep = "\n")
sdat <- standata(fo, data = da)
mod <- rstan::stan_model(model_code = code2)
sf <- rstan::sampling(mod, data = sdat, chains = 4, iter = 4000, seed = 3,
                      refresh = 0)
a <- as.matrix(sf, pars = c("ar[1]", "sigma"))
ess <- posterior::ess_bulk(matrix(a[, 1], ncol = 4))
saveRDS(list(par = a[, 1], sigma = a[, 2]),
        file.path(wt, "dev", "review-acprior-log", "direct-brms.rds"))
A <- readRDS(file.path(wt, "dev", "acprior-out", "agree-after.rds"))$ar$par
B <- readRDS(file.path(wt, "dev", "acprior-out", "agree-before.rds"))$ar$par
R <- readRDS(file.path(wt, "dev", "acprior-out", "agree-brms.rds"))$ar
essA <- posterior::ess_bulk(matrix(A, ncol = 4))
essB <- posterior::ess_bulk(matrix(B, ncol = 4))
row <- function(lab, x, es) {
  cat(sprintf("%-28s mean %.4f sd %.4f q05 %.4f q50 %.4f q95 %.4f ESS %.0f mcse %.5f\n",
              lab, mean(x), sd(x), quantile(x, 0.05), median(x),
              quantile(x, 0.95), es, sd(x) / sqrt(es)))
  c(mean(x), sd(x) / sqrt(es))
}
d0 <- row("brms, prior on sigma_m", a[, 1], ess)
a0 <- row("frm_sample after (lane)", A, essA)
b0 <- row("frm_sample before (lane)", B, essB)
z <- function(u, v) (u[1] - v[1]) / sqrt(u[2]^2 + v[2]^2)
cat(sprintf("z(after - direct) %.2f; z(before - direct) %.2f\n",
            z(a0, d0), z(b0, d0)))
# the lane's weights, for their Kish ESS
pt <- function(s) dt(s / 2.5, 3) / 2.5
s <- sqrt(1 - R$par^2)
w <- pt(R$sigma / s) / (pt(R$sigma) * s)
wn <- w / sum(w)
cat(sprintf("lane weights: Kish ESS %.0f of %d; max weight share %.5f\n",
            1 / sum(wn^2), length(w), max(wn)))
cat(sprintf("weighted brms mean %.4f; direct %.4f; diff / mcse(direct) %.2f\n",
            sum(wn * R$par), d0[1], (sum(wn * R$par) - d0[1]) / d0[2]))
cat(sprintf("KS after vs direct: D %.4f p %.3g (draws autocorrelated)\n",
            suppressWarnings(ks.test(A, a[, 1]))$statistic,
            suppressWarnings(ks.test(A, a[, 1]))$p.value))
