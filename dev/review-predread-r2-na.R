# Short check of fix round 2, item 5: NA in x or sdx across readers and families the lane's
# script does not cover. Seed 7421.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7421)
n <- 240
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4, g = factor(rep(sprintf("g%02d", 1:24), each = 10)))
d$xo <- d$lat + rnorm(n, 0, 0.4)
gl <- rnorm(24)[d$g]
d$sg <- 0.4
d$xg <- gl + rnorm(24, 0, 0.4)[d$g]
eta <- 0.4 + 0.7 * d$lat
d$y <- eta + rnorm(n, 0, 0.5)
d$yg <- 1 + gl + rnorm(n, 0, 0.5)
d$cnt <- rpois(n, exp(eta))
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
d$zc <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$pos <- exp(eta + rnorm(n, 0, 0.5))
nd <- data.frame(xo = c(1, NA, 0, 1), sx = c(0.4, 0.4, NA, 0.4), z = c(0, 0, 0, 1))
sh <- function(lab, e) {
  w <- character()
  cl <- character()
  r <- tryCatch(withCallingHandlers(e, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    cl <<- c(cl, class(c)[1L])
    invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")), error = function(err) structure("ERR", m = conditionMessage(err)))
  txt <- if (identical(c(r), "ERR")) paste("ERROR:", substr(attr(r, "m"), 1, 110)) else {
    v <- if (is.list(r) && !is.data.frame(r)) r[[1L]] else r
    paste(paste(dim(v) %||% length(v), collapse = "x"), " : ", paste(signif(as.numeric(as.matrix(v))[seq_len(min(8, length(unlist(v))))], 3), collapse = " "))
  }
  cat(sprintf("%-30s warn %d %s | %s\n", lab, length(w), paste(unique(cl), collapse = ","), substr(txt, 1, 110)))
}
for (sp in list(list("gaussian", bf(y ~ me(xo, sx) + z), gaussian()),
                list("poisson", bf(cnt ~ me(xo, sx) + z), poisson()),
                list("cumulative", bf(ord ~ me(xo, sx) + z), cumulative()),
                list("categorical", bf(cat ~ me(xo, sx) + z), categorical()),
                list("zip", bf(zc ~ me(xo, sx), zi ~ me(xo, sx)), zero_inflated_poisson()),
                list("lognormal", bf(pos ~ me(xo, sx)), lognormal()),
                list("inter", bf(y ~ me(xo, sx) + me(xo, sx):z), gaussian()))) {
  cat("\n===", sp[[1]], "\n")
  fit <- quiet(frm(sp[[2]], family = sp[[3]], data = d))
  ds <- rev_draws_at(fit, 20L)
  sh("frm_linpred", frm_linpred(fit, newdata = nd))
  sh("frm_linpred se", frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit)
  sh("fitted", fitted(fit, newdata = nd))
  sh("predict", predict(fit, newdata = nd, ndraws = 20))
  sh("simulate", as.matrix(simulate(fit, nsim = 2, seed = 1, newdata = nd)))
  sh("frm_lp_basis extra", frm_lp_basis(fit, newdata = nd)$extra_var)
  sh("frm_curve sim", frm_curve(fit, newdata = nd, seed = 1)$.se)
  sh("frm_curve contrast", frm_curve(fit, newdata = nd, contrast = transform(nd, z = z + 1), simultaneous = FALSE)$.se)
  sh("frm_curve_deriv z", frm_curve_deriv(fit, "z", newdata = nd, simultaneous = FALSE)$.se)
  sh("get_predict", marginaleffects::get_predict(fit, newdata = nd)$estimate)
  sh("predictions", marginaleffects::predictions(fit, newdata = nd)$std.error)
  sh("insight", as.numeric(insight::get_predicted(fit, data = nd)))
  sh("posterior_epred", frmtmb.sample::posterior_epred(ds, newdata = nd))
  sh("posterior_predict", frmtmb.sample::posterior_predict(ds, newdata = nd))
}
cat("\n=== me(gr = )\n")
fg <- quiet(frm(bf(yg ~ me(xg, sg, gr = g) + z), family = gaussian(), data = d))
ng <- data.frame(xg = c(0.2, NA, 0.2, 1), sg = c(0.4, 0.4, 0.4, 0.4), z = 0, g = c("g01", "g01", "g01", "g02"))
sh("frm_linpred", frm_linpred(fg, newdata = ng))
sh("frm_linpred se", frm_linpred(fg, newdata = ng, se.fit = TRUE)$se.fit)
sh("fitted", fitted(fg, newdata = ng))
L <- frm_lp_basis(fg, newdata = ng, extra_cov = TRUE)
print(round(as.matrix(L$extra_cov), 3))
ds <- rev_draws_at(fg, 200L)
set.seed(3)
e <- quiet(frmtmb.sample::posterior_epred(ds, newdata = ng))
cat("draws: identical rows 1, 3:", identical(e[, 1], e[, 3]), " NA row 2:", all(is.na(e[, 2])), "\n")
cat("\n=== columns\n")
fit <- quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = d))
nc <- nd
nc$sx <- NULL
sh("absent sdx", frm_linpred(fit, newdata = nc))
nc <- nd
nc$sx <- "a"
sh("character sdx", frm_linpred(fit, newdata = nc))
nc <- nd
nc$xo <- NULL
sh("absent x", frm_linpred(fit, newdata = nc))
nc <- nd
nc$xo <- rep(NA, 4)
sh("all NA logical x", frm_linpred(fit, newdata = nc))
nc$xo <- rep(NA_real_, 4)
sh("all NA numeric x", frm_linpred(fit, newdata = nc))
nb <- rev_pair("bf(y ~ me(xo, sx) + z)", "gaussian()", d, ndraws = 20)
sh("brms all NA logical x", posterior_epred(nb$fx, newdata = transform(nd, xo = rep(NA, 4))))
sh("brms one NA", posterior_epred(nb$fx, newdata = nd))
dn <- d
dn$xo[c(5, 9)] <- NA
fn <- quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = dn))
sh("fit's own data, fitted", fitted(fn, newdata = dn))
cat("   NA rows:", which(is.na(quiet(fitted(fn, newdata = dn))[, 1])), "\n")
sh("conditional_effects", conditional_effects(fn, "xo", resolution = 3)[[1]]$estimate__)
