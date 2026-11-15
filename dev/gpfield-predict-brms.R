# Lane gpfield, item 1919: predict()'s spread on a gp() field against brms
# 2.23.0's posterior_predict(), inside and past the positions. Data as
# test-gp-predict-draw.R (set.seed(9200), 30 rows, 2 sin(x) + N(0, 0.2)).
# frmtmb arms in one process: "lane" and "base" (gp_drawn_b() mocked to
# return nothing, which is 0.82.0's draw). brms: default priors, 2 chains of
# 2000, seed 1. Args: none
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
if (nzchar(Sys.getenv("GF_PATCH"))) gf_patch("R/predict-brms.R", c("predict_b_drawer", "gp_drawn_b", "gp_field_draw_plan", "gp_field_draw"))
suppressPackageStartupMessages(library(brms))
set.seed(9200)
n <- 30
d <- data.frame(x = sort(stats::runif(n, 0, 10)))
d$y <- 2 * sin(d$x) + stats::rnorm(n, 0, 0.2)
nd <- data.frame(x = c(5.5, 10.5, 12, 15, 60))
fit <- frmtmb::frm(frmtmb::bf(y ~ gp(x)), data = d, family = gaussian())
cat("frm sdgp", exp(fit$estimates$theta[1]), "lscale(data units)",
    exp(fit$estimates$theta[2]), "sigma", exp(fit$estimates$betad[1]), "\n")
real <- get("gp_drawn_b", asNamespace("frmtmb"))
out <- list()
for (arm in c("lane", "base")) {
  gf_put("gp_drawn_b", if (arm == "base") function(fit) integer(0) else real)
  set.seed(1)
  p <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000)
  out[[arm]] <- p
}
gf_put("gp_drawn_b", real)
bfit <- brms::brm(brms::bf(y ~ gp(x)), data = d, family = gaussian(),
                chains = 2, iter = 2000, seed = 1, refresh = 0,
                backend = "rstan")
pb <- brms::posterior_predict(bfit, newdata = nd)
s <- summary(bfit)
cat("brms sdgp", s$gp$Estimate[1], "lscale", s$gp$Estimate[2], "sigma",
    s$spec_pars$Estimate[1], "\n")
w <- function(v) diff(stats::quantile(v, c(0.025, 0.975)))
for (j in seq_len(nrow(nd))) {
  cat(sprintf(paste0("PRED x %5.1f | sd lane %.3f base %.3f brms %.3f | ",
                     "95%% width lane %.3f base %.3f brms %.3f\n"),
              nd$x[j], sd(out$lane[, j]), sd(out$base[, j]), sd(pb[, j]),
              w(out$lane[, j]), w(out$base[, j]), w(pb[, j])))
}
cat("DONE\n")
