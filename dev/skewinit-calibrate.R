# Calibration data for the alpha = 0 stall detector, measured on the
# BASE build. One row per fit: what a stall looks like and what a
# genuinely symmetric fit looks like, in every signal a detector could
# read.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

# Genuinely symmetric: alpha = 0, same covariate shape as the stall
# design so the marginal response is still skewed by the covariate.
make_sym <- function(seed, n = 200) {
  set.seed(1000 + seed)
  xs <- -abs(rnorm(n)) * 3
  data.frame(y = xs + rnorm(n, 0, 1.5), xs = xs)
}
# Mild skew: alpha small but nonzero, the hardest case for any
# near-zero threshold.
make_mild <- function(seed, alpha = 1, n = 200) {
  set.seed(2000 + seed)
  xs <- -abs(rnorm(n)) * 3
  d <- alpha / sqrt(1 + alpha^2)
  e <- 1.5 * (d * abs(rnorm(n)) + sqrt(1 - d^2) * rnorm(n))
  data.frame(y = xs + e, xs = xs)
}

alpha_idx <- function(f) grep("^alpha_", names(f$estimates$betad))
outer_alpha <- function(f) {
  pn <- names(f$opt$par)
  which(pn == "betad")[alpha_idx(f)]
}

one <- function(dd, label, seed) {
  f <- frm_sn(dd)
  ia <- outer_alpha(f)
  p <- f$opt$par
  a_hat <- p[[ia]]
  f0 <- f$obj$fn(p)
  # probe: is the optimum actually a maximum along alpha? the joint
  # objective with alpha alone displaced is an upper bound on the
  # profile, so a drop here is conclusive
  probe <- vapply(c(0.25, 0.5, 1, 2), function(dl) {
    q1 <- p; q1[ia] <- a_hat + dl
    q2 <- p; q2[ia] <- a_hat - dl
    f0 - min(f$obj$fn(q1), f$obj$fn(q2))
  }, 0)
  se_a <- sqrt(diag(vcov(f))[["alpha_Intercept"]])
  fp <- frm_sn(dd, start = list(betad = c(log(sd(dd$y)), 2)))
  fm <- frm_sn(dd, start = list(betad = c(log(sd(dd$y)), -2)))
  best <- max(as.numeric(logLik(fp)), as.numeric(logLik(fm)))
  s <- sn_ll(dd)
  mu <- as.numeric(fitted(f))
  c(seed = seed, ll = as.numeric(logLik(f)), ll_sn = s[["ll"]],
    alpha_sn = s[["alpha"]], alpha = a_hat, se_alpha = se_a,
    ll_refit = best,
    skew_lm = skew(residuals(lm(y ~ xs, data = dd))),
    skew_fit = skew(dd$y - mu),
    probe_0.25 = probe[1], probe_0.5 = probe[2],
    probe_1 = probe[3], probe_2 = probe[4],
    n = nrow(dd))
}

run <- function(label, gen, seeds = 1:40) {
  rows <- lapply(seeds, function(s) {
    r <- try(one(gen(s), label, s), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    r
  })
  rows <- Filter(Negate(is.null), rows)
  d <- as.data.frame(do.call(rbind, rows))
  d$arm <- label
  d
}

out <- rbind(
  run("stall_stated", function(s) make_data(s, FALSE)),
  run("stall_dead",   function(s) make_data(s, TRUE)),
  run("sym_n200",     function(s) make_sym(s, 200)),
  run("sym_n50",      function(s) make_sym(s, 50)),
  run("mild_a1",      function(s) make_mild(s, 1, 200)),
  run("mild_a0.5",    function(s) make_mild(s, 0.5, 200))
)
p <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/calib.csv"
utils::write.csv(out, p, row.names = FALSE)
cat("rows", nrow(out), "->", p, "\n")
print(round(out[, c("seed", "ll", "ll_sn", "alpha", "se_alpha", "ll_refit",
                    "skew_lm", "skew_fit", "probe_0.5", "probe_2")], 4))
