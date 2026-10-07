# Reviewer of lane surface, claim 1: pp_mixture() Est.Error on a fit.
# The lane's construction (dev/surface-ppmix-mc.R: data seed 4, 300
# rows, MC seed 1, R = 4000 draws from N(estimate, V)), extended to ALL
# rows, binned by the estimate, and against one alternative: the SD of
# plogis(N(logit p, lse^2)), the law the reported interval assumes.
#   Rscript dev/surface-rev-ppmix.R > dev/surface-rev-out/ppmix.txt
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages(library(frmtmb))
set.seed(4)
dm <- data.frame(y = c(rnorm(150, -1.5), rnorm(150, 1.5)))
fmx <- frm(y ~ 1, family = mixture(gaussian(), gaussian()), data = dm)
pf <- pp_mixture(fmx)
ds <- frmtmb:::fit_draw_space(fmx)
v0 <- frmtmb:::fit_outer_vector(fmx, ds$map)
set.seed(1)
R <- 4000
Z <- MASS::mvrnorm(R, v0, ds$V)
P <- vapply(seq_len(R), function(r) {
  mixture_probs(frmtmb:::fit_set_outer(fmx, Z[r, ], ds$map))[, 1]
}, numeric(nrow(dm)))
mc_sd <- apply(P, 1, sd)
p <- pf[, "Estimate", 1]
se <- pf[, "Est.Error", 1]
lse <- se / (p * (1 - p))
# SD of the logit-normal the interval is built from, by Gauss-Hermite
gh <- statmod_free_gh <- function(n) {
  # Golub-Welsch for probabilists' Hermite nodes
  i <- seq_len(n - 1)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- sqrt(i)
  J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}
g <- gh(40)
ln_sd <- vapply(seq_along(p), function(i) {
  v <- plogis(qlogis(p[i]) + lse[i] * g$x)
  m <- sum(g$w * v)
  sqrt(max(sum(g$w * (v - m)^2), 0))
}, 0)
bins <- cut(pmin(p, 1 - p), c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5),
            right = FALSE)
q <- apply(P, 1, quantile, c(0.025, 0.975))
tab <- do.call(rbind, lapply(split(seq_along(p), bins), function(ix) {
  if (!length(ix)) return(NULL)
  data.frame(n = length(ix),
             med_delta_over_mc = median(se[ix] / mc_sd[ix]),
             min_delta_over_mc = min(se[ix] / mc_sd[ix]),
             med_logitnormal_over_mc = median(ln_sd[ix] / mc_sd[ix]),
             med_lo_err_logit = median(abs(qlogis(pf[ix, "Q2.5", 1]) -
                                             qlogis(q[1, ix]))),
             med_hi_err_logit = median(abs(qlogis(pf[ix, "Q97.5", 1]) -
                                             qlogis(q[2, ix]))),
             covers_mc_mean = mean(pf[ix, "Q2.5", 1] <= rowMeans(P[ix, ,
               drop = FALSE]) & pf[ix, "Q97.5", 1] >=
               rowMeans(P[ix, , drop = FALSE])))
}))
cat("rows binned by min(p, 1 - p) of component 1 (300 rows)\n")
print(cbind(bin = rownames(tab), signif(tab, 3)), row.names = FALSE)
cat("\nall rows: delta/MC median", signif(median(se / mc_sd), 4),
    "| logit-normal/MC median", signif(median(ln_sd / mc_sd), 4), "\n")
cat("rows with MC SD > 0.01 (lane's filter): delta/MC median",
    signif(median((se / mc_sd)[mc_sd > 0.01]), 4),
    "| logit-normal/MC median",
    signif(median((ln_sd / mc_sd)[mc_sd > 0.01]), 4), "\n")
cat("rows with MC SD <= 0.01:", sum(mc_sd <= 0.01), "\n")
cat("logit SE vs MC SD of logit(p): median ratio",
    signif(median(lse / apply(qlogis(P), 1, sd)), 4), "\n")
saveRDS(list(pf = pf, mc_sd = mc_sd, ln_sd = ln_sd, q = q,
             mc_mean = rowMeans(P)),
        "dev/surface-rev-out/ppmix.rds")
