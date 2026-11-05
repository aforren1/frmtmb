# Item 8: semantic conflicts on the trial merge of predread, latres and nlfix.
# Run with REV_LIB=C:/Users/adf44/source/r/rev-predread-lib/merge. Data seed 8101.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(8101)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(sprintf("g%02d", 1:30), each = 10)),
                lat = rnorm(n), sx = 0.4)
d$t <- rep(1:10, 30)
d$xo <- d$lat + rnorm(n, 0, 0.4)
d$cnt <- rpois(n, exp(0.5 + 0.7 * d$lat))
d$y <- 1 + d$lat + rnorm(n, 0, 0.5)
d$ord <- as.integer(cut(d$x + ifelse(d$z > 0, 1.2, -1.2) * rlogis(n), c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
msg <- function(e) {
  v <- tryCatch(quiet(e), error = function(err) structure("ERR", m = conditionMessage(err)))
  if (identical(c(v), "ERR")) return(paste("ERROR:", substr(attr(v, "m"), 1, 200)))
  if (is.list(v) && !is.data.frame(v)) return(paste("list", paste(names(v), collapse = ",")))
  if (is.data.frame(v)) return(paste("df", nrow(v), paste(names(v), collapse = ",")))
  paste(paste(dim(v), collapse = "x"), paste(signif(as.numeric(v)[seq_len(min(8, length(v)))], 5), collapse = " "))
}
nd <- data.frame(x = c(-1, 0, 1), z = c(0.5, -0.5, 1), xo = c(-1, 0.5, 1), sx = c(0.3, 0.6, 0.9), t = 1:3,
                 g = factor("g01", levels = levels(d$g)), y = 1)
cat("=== 1. latres fitted mean x me(): poisson, ar(cov = TRUE) latent residuals, me()\n")
f1 <- tryCatch(quiet(frm(bf(cnt ~ me(xo, sx) + ar(time = t, gr = g, cov = TRUE)), family = poisson(), data = d)), error = function(e) e)
if (inherits(f1, "error")) {
  cat("frm ERR:", conditionMessage(f1), "\n")
} else {
  cat("frame_latres:", length(frmtmb:::frame_latres(f1$frame)), "\n")
  cat("fitted nd:", msg(fitted(f1, newdata = nd)), "\n")
  cat("frm_linpred response:", msg(frm_linpred(f1, newdata = nd, type = "response")), "\n")
  cat("frm_linpred link se:", msg(frm_linpred(f1, newdata = nd, se.fit = TRUE)$se.fit), "\n")
  cat("frm_lp_basis extra:", msg(frm_lp_basis(f1, newdata = nd)$extra_var), "\n")
  cat("predict:", msg(predict(f1, newdata = nd, ndraws = 50)), "\n")
  cat("simulate:", msg(as.matrix(simulate(f1, nsim = 2, seed = 1, newdata = nd))), "\n")
  cat("frm_curve:", msg(frm_curve(f1, newdata = nd, simultaneous = FALSE)$.se), "\n")
  cat("CE:", msg(conditional_effects(f1, "xo", resolution = 3)[[1]]$estimate__), "\n")
  ds <- rev_draws_at(f1, 4L)
  cat("draws epred:", msg(frmtmb.sample::posterior_epred(ds, newdata = nd)), "\n")
}
cat("\n=== 2. latres x cs(): frm_linpred(se.fit) with several layers under latres_skip_se\n")
f2 <- tryCatch(quiet(frm(bf(ord ~ x + cs(z) + ar(time = t, gr = g, cov = TRUE)), family = sratio(), data = d)), error = function(e) e)
if (inherits(f2, "error")) {
  cat("frm ERR:", substr(conditionMessage(f2), 1, 200), "\n")
} else {
  v <- frm_linpred(f2, newdata = NULL)
  cat("in-sample layers dim:", dim(v), "\n")
  cat("se.fit:", msg(frm_linpred(f2, se.fit = TRUE)$se.fit), "\n")
  cat("newdata layers:", msg(frm_linpred(f2, newdata = nd)), "\n")
  cat("frm_curve:", msg(frm_curve(f2, newdata = nd, simultaneous = FALSE)$.se), "\n")
  cat("fitted linear:", msg(dim(fitted(f2, newdata = nd, scale = "linear"))), "\n")
}
cat("\n=== 3. curve under a latent term: curve_parts' cov = FALSE check against latres frames\n")
f3 <- tryCatch(quiet(frm(bf(cnt ~ x + ar(time = t, gr = g, cov = FALSE)), family = poisson(), data = d)), error = function(e) e)
if (inherits(f3, "error")) {
  cat("frm ERR:", substr(conditionMessage(f3), 1, 240), "\n")
} else {
  cat("frm_curve:", msg(frm_curve(f3, newdata = nd, simultaneous = FALSE)$.se), "\n")
  cat("frm_linpred nd:", msg(frm_linpred(f3, newdata = nd)), "\n")
  cat("predict nd:", msg(predict(f3, newdata = nd, ndraws = 20)), "\n")
  cat("simulate nd:", msg(as.matrix(simulate(f3, nsim = 2, seed = 1, newdata = nd))), "\n")
}
f3b <- tryCatch(quiet(frm(bf(cnt ~ x + ar(time = t, gr = g, cov = TRUE)), family = poisson(), data = d)), error = function(e) e)
if (!inherits(f3b, "error")) {
  cat("cov TRUE latent, frm_curve:", msg(frm_curve(f3b, newdata = nd, simultaneous = FALSE)$.se), "\n")
  cat("cov TRUE latent, frm_lp_basis:", msg(frm_lp_basis(f3b, newdata = nd)$extra_var), "\n")
  cat("cov TRUE latent, fitted nd:", msg(fitted(f3b, newdata = nd)), "\n")
  cat("cov TRUE latent, frm_linpred response:", msg(frm_linpred(f3b, newdata = nd, type = "response")), "\n")
  cat("cov TRUE latent, predict nd:", msg(predict(f3b, newdata = nd, ndraws = 20)), "\n")
}
cat("\n=== 4. nl x ordinal (nlfix D1) through the layered readers and dpar_response_link\n")
f4 <- tryCatch(quiet(frm(bf(ord ~ a + b * x, a ~ z, b ~ 1, nl = TRUE), family = cumulative(), data = d)), error = function(e) e)
if (inherits(f4, "error")) {
  cat("frm ERR:", substr(conditionMessage(f4), 1, 240), "\n")
} else {
  cat("frm_linpred link:", msg(frm_linpred(f4, newdata = nd)), "\n")
  cat("frm_linpred dpar mu response:", msg(frm_linpred(f4, newdata = nd, dpar = "mu", type = "response")), "\n")
  cat("plogis of link:", paste(signif(plogis(as.numeric(frm_linpred(f4, newdata = nd))), 5), collapse = " "), "\n")
  cat("frm_linpred response:", msg(frm_linpred(f4, newdata = nd, type = "response")), "\n")
  cat("fitted dpar mu:", msg(fitted(f4, newdata = nd, dpar = "mu")), "\n")
  cat("fitted:", msg(fitted(f4, newdata = nd)), "\n")
  cat("frm_curve:", msg(frm_curve(f4, newdata = nd, simultaneous = FALSE)$.estimate), "\n")
  cat("frm_curve transform dpar mu:", msg(frm_curve(f4, newdata = nd, dpar = "mu", transform = TRUE, simultaneous = FALSE)$.estimate), "\n")
  cat("CE dpar mu:", msg(conditional_effects(f4, "x", dpar = "mu", resolution = 3)[[1]]$estimate__), "\n")
  ds <- rev_draws_at(f4, 2L)
  cat("draws linpred:", msg(frmtmb.sample::posterior_linpred(ds, newdata = nd)), "\n")
  cat("draws epred dpar mu:", msg(frmtmb.sample::posterior_epred(ds, newdata = nd, dpar = "mu")), "\n")
}
cat("\n=== 5. me() inside a body, with nlfix's nl_body_value()\n")
f5 <- quiet(frm(bf(y ~ a + b * z, a ~ me(xo, sx), b ~ 1, nl = TRUE), family = gaussian(), data = d))
ds5 <- rev_draws_at(f5, 2000L)
set.seed(1)
e5 <- frmtmb.sample::posterior_epred(ds5, newdata = nd)
cat("draws epred sd:", signif(apply(e5, 2, sd), 4), " expected about", signif(abs(f5$estimates$beta[2]) * nd$sx, 4), "\n")
cat("predict sd:", signif(apply(predict(f5, newdata = nd, summary = FALSE, ndraws = 2000, propagate_error = FALSE), 2, sd), 4), "\n")
f5b <- tryCatch(quiet(frm(bf(ord ~ a + b * z, a ~ me(xo, sx), b ~ 1, nl = TRUE), family = cumulative(), data = d)), error = function(e) e)
if (!inherits(f5b, "error")) {
  cat("ordinal body with me(): frm_linpred:", msg(frm_linpred(f5b, newdata = nd)), "\n")
  cat("ordinal body with me(): fitted Est.Error:", msg(fitted(f5b, newdata = nd)[, "Est.Error", 1]), "\n")
  dsb <- rev_draws_at(f5b, 500L)
  set.seed(2)
  cat("draws epred cat 1 sd:", signif(apply(frmtmb.sample::posterior_epred(dsb, newdata = nd)[, , 1], 2, sd), 4), "\n")
  cat("predict:", msg(predict(f5b, newdata = nd, ndraws = 20)), "\n")
} else cat("ordinal body with me(): frm ERR", substr(conditionMessage(f5b), 1, 200), "\n")
cat("\n=== 6. cov = FALSE on a body (nlfix), the new frm_curve refusal, and simulate\n")
d$yar <- d$y
f6 <- tryCatch(quiet(frm(bf(yar ~ a + b * z, a ~ 1 + xo, b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), family = gaussian(), data = d)), error = function(e) e)
if (inherits(f6, "error")) {
  cat("frm ERR:", substr(conditionMessage(f6), 1, 240), "\n")
} else {
  ndy <- transform(nd, yar = c(1, 2, 3))
  cat("frm_linpred:", msg(frm_linpred(f6, newdata = ndy)), "\n")
  cat("frm_curve:", msg(frm_curve(f6, newdata = ndy, simultaneous = FALSE)$.se), "\n")
  s6 <- t(vapply(simulate(f6, nsim = 2000, seed = 2, newdata = ndy), as.numeric, numeric(3)))
  cat("simulate mean", signif(colMeans(s6), 4), " sd ", signif(apply(s6, 2, sd), 4), "\n")
  p6 <- predict(f6, newdata = ndy, summary = FALSE, ndraws = 2000, propagate_error = FALSE)
  cat("predict mean", signif(colMeans(p6), 4), " sd ", signif(apply(p6, 2, sd), 4), "\n")
}
