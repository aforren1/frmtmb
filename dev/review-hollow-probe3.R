# Reviewer probe 3: the cens false accept end to end, make_conditions()
# on draws, and the size the stored $data adds to a saved draws object
.libPaths(c("C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb.sample"), "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
draws_at <- function(fit, k = 5) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  est <- coef(fit)
  m <- matrix(0, k, length(lab), dimnames = list(NULL, lab))
  # parameters at the ML estimate, perturbed per draw, on frmtmb's chart
  set.seed(11)
  m[, "b_Intercept"] <- est[["b_Intercept"]] %||% 1
  structure(list(stanfit = NULL, draws = m, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}

# 1. interval cens: newdata whose upper bound differs at an interval row
set.seed(2)
n <- 40
dd <- data.frame(x = rnorm(n), cc = rep(c(0, 0, 0, 2), 10))
dd$y <- 1 + 0.5 * dd$x + rnorm(n)
dd$y2 <- dd$y + 1
fc <- q(frm(y | cens(cc, y2) ~ x, data = dd))
lab <- c(frmtmb::brms_par_labels(fc), "lp__")
cat("labels:", paste(lab, collapse = " "), "\n")
m <- matrix(0, 3, length(lab), dimnames = list(NULL, lab))
m[, "b_Intercept"] <- c(0.9, 1, 1.1)
m[, "b_x"] <- 0.5
m[, "sigma"] <- 1
ds <- structure(list(stanfit = NULL, draws = m, fit = fc, data = fc$data),
                class = "frmtmb_draws")
ll0 <- log_lik(ds)
d1 <- dd; d1$y2[4] <- d1$y2[4] + 3
ll1 <- tryCatch(log_lik(ds, newdata = d1), error = function(e) e)
if (inherits(ll1, "error")) {
  cat("cens changed bound: refused:", conditionMessage(ll1), "\n")
} else {
  cat("cens changed bound: accepted; identical to log_lik(ds):",
      identical(ll1, ll0), "\n")
  mu <- 1 + 0.5 * dd$x[4]
  cat(sprintf("row 4 (interval [%.4f, %.4f] then [%.4f, %.4f]): log_lik %.6f; ",
              dd$y[4], dd$y2[4], d1$y[4], d1$y2[4], ll1[2, 4]),
      sprintf("density on the changed interval %.6f, on the fitted %.6f\n",
              log(pnorm(d1$y2[4], mu) - pnorm(d1$y[4], mu)),
              log(pnorm(dd$y2[4], mu) - pnorm(dd$y[4], mu))))
}

# 2. make_conditions() on draws against on the fit, poly() model
fp <- q(frm(y ~ poly(x, 2), data = dd))
dp <- structure(list(stanfit = NULL, draws = matrix(0, 1, 1), fit = fp,
                     data = fp$data), class = "frmtmb_draws")
r_fit <- tryCatch(make_conditions(fp, "x"), error = function(e) e)
r_drw <- tryCatch(make_conditions(dp, "x"), error = function(e) e)
cat("make_conditions(fit, 'x'):",
    if (inherits(r_fit, "error")) conditionMessage(r_fit) else
      paste(format(unlist(r_fit$x), digits = 4), collapse = " "), "\n")
cat("make_conditions(draws, 'x'):",
    if (inherits(r_drw, "error")) conditionMessage(r_drw) else
      paste(format(unlist(r_drw$x), digits = 4), collapse = " "), "\n")

# 3. size: serialized bytes of the frame against the fit and a draws
# object with 4000 draws of the fixed parameters and the group effects
for (nn in c(1e4, 1e5)) {
  set.seed(3)
  big <- data.frame(x = rnorm(nn), z = rnorm(nn),
                    g = factor(sample(1:200, nn, TRUE)))
  big$y <- 1 + big$x + rnorm(200)[big$g] + rnorm(nn)
  fb <- q(frm(y ~ x + z + (1 | g), data = big))
  labb <- c(frmtmb::brms_par_labels(fb), "lp__")
  mm <- matrix(rnorm(4000 * length(labb)), 4000, length(labb))
  without <- list(stanfit = NULL, draws = mm, fit = fb)
  with <- c(without, list(data = fb$data))
  sz <- function(o) length(serialize(o, NULL))
  cat(sprintf(paste("n %g: frame %.2f MB; fit %.2f MB; draws object",
                    "without $data %.2f MB, with %.2f MB (+%.1f%%)\n"),
              nn, sz(fb$data) / 2^20, sz(fb) / 2^20, sz(without) / 2^20,
              sz(with) / 2^20, 100 * (sz(with) / sz(without) - 1)))
}
