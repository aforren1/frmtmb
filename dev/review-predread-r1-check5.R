# Short check, items 5, 6: the multivariate shape change and its callers, the constant
# dpar on the linear scale against a predicted one, NA handling, and simulate()'s cost.
# Seed 8301. REV_LIB=before gives the before arm of the shapes and the timing.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
before <- identical(Sys.getenv("REV_LIB"), "before")
set.seed(8301)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
d$y <- 1 + d$x + rnorm(n)
d$y2 <- -d$x + 0.5 * d$z + rnorm(n)
d$pos <- exp(0.3 + 0.5 * d$x + rnorm(n, 0, 0.4))
d$ym <- 1 + d$lat + rnorm(n, 0, 0.5)
msg <- function(e) {
  r <- tryCatch(quiet(e), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 130)))
  if (is.character(r) && length(r) == 1L) return(r)
  if (is.list(r) && !is.data.frame(r)) return(paste("list:", paste(names(r), vapply(r, function(v) paste(dim(v) %||% length(v), collapse = "x"), ""), collapse = " ")))
  if (is.data.frame(r)) return(paste("df", nrow(r), "x", ncol(r)))
  paste(paste(dim(r) %||% length(r), collapse = "x"), "|", paste(colnames(r) %||% "", collapse = ","))
}
nd <- data.frame(x = c(-1, 0, 1), z = 0.2)
mv <- quiet(frm(frmtmb::mvbf(bf(y ~ x), bf(y2 ~ x + z), rescor = FALSE), family = gaussian(), data = d))
cat("=== mv with no resp, every reader that could assume one column\n")
R <- list(
  frm_linpred = function() frm_linpred(mv, newdata = nd),
  frm_linpred_response = function() frm_linpred(mv, newdata = nd, type = "response"),
  frm_linpred_se = function() frm_linpred(mv, newdata = nd, se.fit = TRUE),
  frm_linpred_in = function() frm_linpred(mv),
  frm_lp_basis = function() frm_lp_basis(mv, newdata = nd),
  frm_curve = function() frm_curve(mv, newdata = nd, simultaneous = FALSE),
  frm_curve_deriv = function() frm_curve_deriv(mv, "x", newdata = nd),
  frm_curve_feature = function() frm_curve_feature(mv, "x", type = "crossing", at = 0, newdata = data.frame(x = seq(-3, 3, by = 0.5), z = 0)),
  fitted = function() fitted(mv, newdata = nd),
  fitted_linear = function() fitted(mv, newdata = nd, scale = "linear"),
  predict = function() predict(mv, newdata = nd, ndraws = 5),
  residuals = function() residuals(mv),
  get_predict = function() marginaleffects::get_predict(mv, newdata = nd),
  get_predict_link = function() marginaleffects::get_predict(mv, newdata = nd, type = "link"),
  predictions = function() marginaleffects::predictions(mv, newdata = nd),
  predictions_resp = function() marginaleffects::predictions(mv, newdata = nd, resp = "y2"),
  insight = function() as.numeric(insight::get_predicted(mv, data = nd)),
  insight_link = function() as.numeric(insight::get_predicted(mv, data = nd, predict = "link")),
  emmeans = function() as.data.frame(emmeans::emmeans(mv, ~ x, at = list(x = 1))),
  emmeans_rep = function() as.data.frame(emmeans::emmeans(mv, ~ x | rep.meas, at = list(x = 1))),
  CE = function() names(conditional_effects(mv, "x", resolution = 3)),
  simulate = function() simulate(mv, nsim = 1, seed = 1, newdata = nd),
  loo_like = function() logLik(mv))
for (nm in names(R)) cat(sprintf("%-22s %s\n", nm, substr(msg(R[[nm]]()), 1, 150)))
if (!before) {
  bm <- rev_pair("mvbf(bf(y ~ x), bf(y2 ~ x + z), rescor = FALSE)", "gaussian()", d)
  pl <- posterior_linpred(bm$fx, newdata = nd)
  cat("brms posterior_linpred dim", dim(pl), " dn3 ", dimnames(pl)[[3]], "\n")
  cat("rel of frm_linpred to brms:", rev_E$eq_rel(frm_linpred(mv, newdata = nd), pl[1, , ]), " colnames ", colnames(frm_linpred(mv, newdata = nd)), "\n")
  fb <- fitted(bm$fx, newdata = nd)
  cat("brms fitted dim", dim(fb), " dn3 ", dimnames(fb)[[3]], "\n")
  cat("brms emmeans:", msg(as.data.frame(emmeans::emmeans(bm$fx, ~ x, at = list(x = 1)))), "\n")
  cat("brms predictions:", msg(marginaleffects::predictions(bm$fx, newdata = nd)), "\n")
  cat("frm_lp_basis no resp:", msg(frm_lp_basis(mv, newdata = nd)), "\n")
}
cat("\n=== constant dpar on the linear scale, and the same dpar predicted\n")
if (!before) {
  for (sp in list(c("bf(pos ~ x)", "Gamma('log')", "shape"), c("bf(pos ~ x, shape ~ 1)", "Gamma('log')", "shape"))) {
    pr <- rev_pair(sp[1], sp[2], d)
    g <- data.frame(x = c(-1, 1))
    cat(sp[1], ": brms linear", signif(fitted(pr$fx, newdata = g, dpar = sp[3], scale = "linear", summary = FALSE)[1, ], 5),
        " response", signif(fitted(pr$fx, newdata = g, dpar = sp[3], summary = FALSE)[1, ], 5),
        " | frmtmb linear", signif(fitted(pr$fit, newdata = g, dpar = sp[3], scale = "linear")[, "Estimate"], 5),
        " response", signif(fitted(pr$fit, newdata = g, dpar = sp[3])[, "Estimate"], 5), "\n")
  }
  rd <- tools::Rd_db("frmtmb")
  tx <- capture.output(tools::Rd2txt(rd[["frm_linpred.Rd"]], fragment = TRUE))
  cat(grep("no predictor|constant|fixed dpar", tx, value = TRUE), sep = "\n")
}
cat("\n=== NA in new data: a plain predictor against me()\n")
pl <- quiet(frm(bf(ym ~ x + z), family = gaussian(), data = d))
fm <- quiet(frm(bf(ym ~ me(xo, sx) + z), family = gaussian(), data = d))
na1 <- data.frame(x = c(1, NA), z = 0)
na2 <- data.frame(xo = c(1, NA), sx = 0.4, z = 0)
na3 <- data.frame(xo = c(1, 1), sx = c(0.4, NA), z = 0)
na4 <- data.frame(xo = c(1, 1), sx = 0.4, z = c(0, NA))
cat("plain, NA in x: frm_linpred", msg(frm_linpred(pl, newdata = na1)), " values ", tryCatch(paste(signif(frm_linpred(pl, newdata = na1), 4), collapse = " "), error = function(e) "refused"), "\n")
cat("plain, NA in x: fitted", tryCatch(paste(signif(fitted(pl, newdata = na1)[, 1], 4), collapse = " "), error = function(e) substr(conditionMessage(e), 1, 100)), "\n")
cat("me, NA in x: ", msg(frm_linpred(fm, newdata = na2)), "\n")
cat("me, NA in sdx: ", msg(frm_linpred(fm, newdata = na3)), "\n")
cat("me, NA in z:  ", tryCatch(paste(signif(frm_linpred(fm, newdata = na4), 4), collapse = " "), error = function(e) substr(conditionMessage(e), 1, 100)), "\n")
dn <- d
dn$xo[c(3, 9)] <- NA
fna <- quiet(frm(bf(ym ~ me(xo, sx) + z), family = gaussian(), data = dn))
cat("fit with NA in x, rows kept:", nobs(fna), " of", nrow(dn), "\n")
cat("fitted() in sample:", msg(fitted(fna)), "\n")
cat("the fit's own data as newdata:", msg(fitted(fna, newdata = dn)), "\n")
cat("predict(the fit's own data):", msg(predict(fna, newdata = dn, ndraws = 5)), "\n")
dn2 <- d
dn2$sx[c(3, 9)] <- NA
fna2 <- tryCatch(quiet(frm(bf(ym ~ me(xo, sx) + z), family = gaussian(), data = dn2)), error = function(e) e)
cat("fit with NA in sdx:", if (inherits(fna2, "error")) substr(conditionMessage(fna2), 1, 130) else paste("rows", nobs(fna2)), "\n")
if (!before) {
  bfit <- rev_pair("bf(ym ~ me(xo, sx) + z)", "gaussian()", d, ndraws = 50)
  cat("brms on the row with NA:", tryCatch(paste(signif(colMeans(suppressWarnings(posterior_epred(bfit$fx, newdata = na2))), 4), collapse = " "), error = function(e) conditionMessage(e)), "\n")
  cat("brms fit with NA in x:", tryCatch({
    b0 <- quiet(brms::brm(rev_ev("bf(ym ~ me(xo, sx) + z)", "brms"), data = dn, empty = TRUE))
    paste("rows", nrow(b0$data))
  }, error = function(e) substr(conditionMessage(e), 1, 120)), "\n")
}
cat("\n=== simulate(newdata) cost on a plain me() model and a plain model\n")
set.seed(8302)
N <- 100000
big <- data.frame(z = rnorm(N), lat = rnorm(N), sx = 0.4, x = rnorm(N),
                  g = factor(sample(sprintf("g%03d", 1:100), N, TRUE)))
big$xo <- big$lat + rnorm(N, 0, 0.4)
big$y <- 1 + big$lat + 0.3 * big$z + rnorm(100, 0, 0.3)[big$g] + rnorm(N)
fb1 <- quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = big[1:20000, ]))
fb2 <- quiet(frm(bf(y ~ x + z + (1 | g)), family = gaussian(), data = big[1:20000, ]))
ctrl <- function() {
  s <- 0
  for (i in 1:3000000) s <- s + sqrt(i)
  s
}
tm <- function(lab, e) {
  ts <- numeric(3)
  for (k in 1:3) {
    invisible(gc())
    ts[k] <- system.time(quiet(eval(e)))[["elapsed"]]
  }
  cat(sprintf("%-46s min %.3f s (%s)\n", lab, min(ts), paste(round(ts, 3), collapse = " ")))
}
tm("control", quote(ctrl()))
tm("simulate me(), nsim 20, 100000 new rows", quote(simulate(fb1, nsim = 20, seed = 1, newdata = big)))
tm("simulate plain + (1 | g), nsim 20, 100000 rows", quote(simulate(fb2, nsim = 20, seed = 1, newdata = big)))
tm("predict me(), 20 draws, 20000 new rows", quote(predict(fb1, newdata = big[1:20000, ], ndraws = 20)))
