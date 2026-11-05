# Short check 2, item 4, R4: the exact change, tried in-process. A guard at
# get_joint_cov() with the limit at 1000 rows, on 2000 rows: what does each route do?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
ns <- asNamespace("frmtmb")
orig <- get("get_joint_cov", ns)
g <- function(fit) {
  if (length(frame_latres(fit$frame))) {
    m <- fit$frame[["n_obs"]]
    lim <- getOption("frmtmb.latres_dense_rows", 20000L)
    if (m > lim) {
      frm_stop("This standard error reads the whole joint covariance of a latent residual as a dense matrix, about ",
               format(16 * m^2 / 1e9, digits = 2), " GB at ", m, " rows. options(frmtmb.latres_dense_rows = ) moves the limit of ",
               lim, " rows; frm_linpred(), fitted() on new rows, emmeans() and predict() do not read it", call. = FALSE)
    }
  }
  orig(fit)
}
environment(g) <- ns
assignInNamespace("get_joint_cov", g, ns = "frmtmb")
options(frmtmb.latres_dense_rows = 1000L)
set.seed(1)
N <- 2000
d <- data.frame(t = seq_len(N), g = factor(1), x = rnorm(N))
d$y <- rpois(N, exp(.5 + .3 * d$x + as.numeric(arima.sim(list(ar = .5), N, sd = .5))))
f <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
nd <- d[1:20, ]
for (lab in c("fitted", "residuals", "fitted_new", "linpred_new", "ce", "curve", "basis", "basis_new", "predictions", "slopes", "emmeans", "predict", "predict_new", "confint", "summary", "hypothesis", "simulate")) {
  r <- msg_of(switch(lab,
    fitted = fitted(f), residuals = residuals(f), fitted_new = fitted(f, newdata = nd),
    linpred_new = frm_linpred(f, newdata = nd, se.fit = TRUE), ce = conditional_effects(f, "x"),
    curve = frm_curve(f, newdata = nd), basis = frm_lp_basis(f), basis_new = frm_lp_basis(f, newdata = nd),
    predictions = marginaleffects::predictions(f, newdata = nd), slopes = marginaleffects::avg_slopes(f, variables = "x"),
    emmeans = emmeans::emmeans(f, ~ x, at = list(x = 0)), predict = predict(f, ndraws = 10), predict_new = predict(f, newdata = nd, ndraws = 10),
    confint = confint(f), summary = summary(f), hypothesis = hypothesis(f, "x > 0"), simulate = simulate(f, nsim = 1)))
  cat(sprintf("%-12s %s\n", lab, if (is.null(r$err)) "ok" else paste("REFUSED:", substr(r$err, 1, 110))))
}
