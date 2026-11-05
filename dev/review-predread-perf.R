# Item 3: per-reader time on a large fit, one build per process. Seeds fixed; blocks
# grown past 1.2 s; minimum of rounds; a fixed arithmetic control in each process.
#   REV_LIB=before Rscript dev/review-predread-perf.R <out.tsv>
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
outf <- commandArgs(TRUE)[1]
set.seed(8001)
n <- 200000
d <- data.frame(x = rnorm(n), z = rnorm(n),
                f = factor(sample(letters[1:6], n, TRUE)),
                g = factor(sample(sprintf("g%03d", 1:200), n, TRUE)))
re <- rnorm(200, 0, 0.5)[d$g]
d$y <- 0.5 + d$x - 0.4 * d$z + re + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.3 * d$x + 0.3 * re))
d$ord <- as.integer(cut(d$x + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
sm <- d[seq_len(5000), ]
fit <- quiet(frm(bf(y ~ x * f + z + (1 | g)), family = gaussian(), data = d))
fit2 <- quiet(frm(bf(y ~ s(x, k = 10) + z), family = gaussian(), data = sm))
fit3 <- quiet(frm(bf(ord ~ x + z), family = cumulative(), data = sm))
fitp <- quiet(frm(bf(cnt ~ x + z), family = poisson(), data = d))
ds <- rev_draws_at(fit, 40L)
nd <- d[seq_len(50000), ]
gr <- data.frame(x = seq(-3, 3, length.out = 400), z = 0, f = factor("a", levels = levels(d$f)),
                 g = factor("g001", levels = levels(d$g)))
readers <- list(
  linpred_nd = function() frm_linpred(fit, newdata = nd),
  linpred_se = function() frm_linpred(fit, newdata = nd[1:5000, ], se.fit = TRUE),
  linpred_in = function() frm_linpred(fit),
  fitted_nd = function() fitted(fitp, newdata = nd),
  basis = function() frm_lp_basis(fit, newdata = gr, extra_cov = TRUE),
  curve = function() frm_curve(fit, newdata = gr, seed = 1),
  curve_deriv = function() frm_curve_deriv(fit, "x", newdata = gr, seed = 1),
  curve_smooth = function() frm_curve(fit2, newdata = data.frame(x = seq(-2, 2, length.out = 300), z = 0), seed = 1),
  predict_nd = function() predict(fit, newdata = nd[1:2000, ], ndraws = 20),
  simulate_nd = function() simulate(fit, nsim = 5, seed = 1, newdata = nd[1:20000, ]),
  emm = function() as.data.frame(emmeans::emmeans(fit, ~ f)),
  ce = function() conditional_effects(fit, "x", resolution = 50),
  get_predict = function() marginaleffects::get_predict(fit, newdata = nd),
  ordinal_fitted = function() fitted(fit3, newdata = sm[1:2000, ]),
  draws_epred = function() frmtmb.sample::posterior_epred(ds, newdata = nd[1:5000, ], draw_ids = 1:40),
  draws_linpred = function() frmtmb.sample::posterior_linpred(ds, newdata = nd[1:5000, ], draw_ids = 1:40),
  draws_pp = function() frmtmb.sample::posterior_predict(ds, newdata = nd[1:5000, ], draw_ids = 1:40),
  control = function() {
    s <- 0
    for (i in 1:3000000) s <- s + sqrt(i)
    s
  })
out <- list()
for (r in names(readers)) {
  invisible(quiet(readers[[r]]()))
  reps <- 1L
  repeat {
    t <- system.time(for (i in seq_len(reps)) quiet(readers[[r]]()))[["elapsed"]]
    if (t >= 1.2 || reps >= 64L) break
    reps <- reps * 2L
  }
  ts <- t / reps
  for (k in 1:3) {
    t2 <- system.time(for (i in seq_len(reps)) quiet(readers[[r]]()))[["elapsed"]] / reps
    ts <- c(ts, t2)
  }
  out[[r]] <- data.frame(reader = r, reps = reps, tmin = min(ts), tmed = median(ts))
  cat(r, reps, round(min(ts), 4), "\n")
}
write.table(do.call(rbind, out), outf, sep = "\t", row.names = FALSE, quote = FALSE)
