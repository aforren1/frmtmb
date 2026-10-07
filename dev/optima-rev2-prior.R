# Reviewer of lane optima, re-check (b): the sampler's simplex density
# alone is uniform on the simplex. (1) frmtmb.sample's mo_simplex_nlp()
# taped by itself and sampled with tmbstan, D = 3 and D = 5; (2) the
# whole frm_sample() pipeline on a model whose likelihood barely reads
# the simplex (b held near 0 by a normal(0, 1e-4) prior). Test: on
# draws thinned to near independence, Kolmogorov-Smirnov of w_1
# against Beta(1, D - 1) (the Dirichlet(1) marginal) and of
# w_1 / (w_1 + w_2) against U(0, 1). Control: the same density without
# the Jacobian (flat in the softmax coordinates) must be rejected.
#   Rscript dev/optima-rev2-prior.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
softmax <- function(z) {
  e <- exp(c(0, z))
  e / sum(e)
}
test_w <- function(label, W, D) {
  k1 <- stats::ks.test(W[, 1], "pbeta", 1, D - 1)
  r <- W[, 1] / (W[, 1] + W[, 2])
  k2 <- stats::ks.test(r, "punif")
  cat(sprintf(paste0("%-34s n %5d | mean w %s | KS w1~Beta(1,%d) ",
                     "p %.3g | KS w1/(w1+w2)~U p %.3g | min w q01 %.2g\n"),
              label, nrow(W), paste(format(colMeans(W), digits = 3),
                                    collapse = ","), D - 1, k1$p.value,
              k2$p.value, stats::quantile(apply(W, 1, min), 0.01)))
}
for (D in c(3, 5)) {
  set.seed(D)
  d <- data.frame(x = sample(0:D, 200, TRUE))
  d$y <- 0.3 * d$x + rnorm(200)
  fit <- frm(bf(y ~ mo(x)), data = d, family = gaussian())
  snlp <- frmtmb.sample:::mo_simplex_nlp(fit)
  zn <- frmtmb::mo_frame_terms(fit)[[1]]$zeta
  for (arm in c("lane", "no Jacobian")) {
    fn <- if (arm == "lane") {
      function(p) snlp(p)
    } else {
      # the control: a flat density on the softmax coordinates, kept
      # proper by a wide normal so tmbstan has a target
      function(p) sum(p[[zn]]^2) / (2 * 15^2)
    }
    pl <- stats::setNames(list(numeric(D - 1)), zn)
    obj <- RTMB::MakeADFun(fn, pl, silent = TRUE)
    sf <- suppressWarnings(tmbstan::tmbstan(obj, chains = 4, iter = 6000,
                                            warmup = 1000, seed = 3,
                                            refresh = 0))
    z <- as.matrix(sf)[, seq_len(D - 1), drop = FALSE]
    W <- t(apply(z, 1, softmax))
    sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
    div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    # thin to near independence: every 10th of 20000
    Wt <- W[seq(1, nrow(W), by = 10), , drop = FALSE]
    test_w(sprintf("D=%d tape alone, %s (div %d)", D, arm, div),
           Wt, D)
  }
}
# (2) the whole pipeline, D = 3 (x in 0:3)
set.seed(9)
d <- data.frame(x = sample(0:3, 200, TRUE))
d$y <- rnorm(200)
fit <- frm(bf(y ~ mo(x)), data = d, family = gaussian())
pr <- set_prior("normal(0, 0.0001)", class = "b", coef = "mox")
s <- suppressWarnings(suppressMessages(
  frm_sample(fit, prior = pr, chains = 4, iter = 6000, warmup = 1000,
             seed = 3, cores = 1, refresh = 0)))
m <- as.matrix(s)
W <- m[, grep("^simo_", colnames(m)), drop = FALSE]
cat("pipeline columns:", colnames(W), "| max |row sum - 1|",
    format(max(abs(rowSums(W) - 1)), digits = 3), "| b sd",
    format(stats::sd(m[, grep("^bsp_", colnames(m))]), digits = 3), "\n")
test_w("D=3 frm_sample, b ~ normal(0, 1e-4)",
       W[seq(1, nrow(W), by = 10), , drop = FALSE], ncol(W))
