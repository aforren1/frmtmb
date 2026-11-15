# Lane silent2: reproducers of items 5235, 5246, 5353, 5362 and the
# print order of draws. Usage: Rscript dev/silent2-repro.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"),
    packageDescription("frmtmb")$Version, "\n")
try_s <- function(e) {
  ws <- character()
  r <- tryCatch(withCallingHandlers({force(e); "ran"},
                  warning = function(w) {
                    ws <<- c(ws, conditionMessage(w))
                    invokeRestart("muffleWarning")
                  }),
                error = function(err) paste("ERROR:",
                  substr(conditionMessage(err), 1, 200)))
  paste0(r, if (length(ws)) paste0(" [", length(ws), " warnings: ",
                                   substr(ws[1], 1, 120), "]"))
}

## item 5235: an addition term of length n by coincidence
scaled <- frmtmb_family(
  "scaled_normal", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- dpars[["mu"]] * aterms[["vreal1"]]
    -0.5 * ((y - m) / dpars[["sigma"]])^2 - log(dpars[["sigma"]]) -
      0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y / aterms[["vreal1"]]),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vreal1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]])
  })
set.seed(3)
n <- 120
d <- data.frame(x = rnorm(n))
d$s <- rep_len(c(1, 3), n)
d$y <- (1 + 0.5 * d$x) * d$s + rnorm(n, 0, 0.3)
f_expr <- frm(bf(y | vreal(rep_len(c(1, 3), 120)) ~ x), family = scaled,
              data = d)
cat("5235 CE predict, resolution 120:",
    try_s(conditional_effects(f_expr, effects = "x", method = "predict",
                              resolution = 120)), "\n")
cat("5235 CE mean, resolution 120:",
    try_s(conditional_effects(f_expr, effects = "x", resolution = 120)),
    "\n")
nd <- data.frame(x = seq(-2, 2, length.out = 120))
cat("5235 predict(newdata = 120 rows):",
    try_s(fitted(f_expr, newdata = nd)), "\n")
p <- tryCatch(fitted(f_expr, newdata = nd)[, 1], error = function(e) NULL)
if (!is.null(p)) cat("  rows 1:2 ratio", p[2] / p[1], "(x nearly equal)\n")
cat("5235 in-sample fitted():", try_s(fitted(f_expr)), "\n")
cat("5235 predict(newdata = 5 rows):",
    try_s(fitted(f_expr, newdata = nd[1:5, , drop = FALSE])), "\n")

## item 5246: anova() on fits of different rows of the same count
set.seed(5)
d6 <- data.frame(x = rnorm(80))
d6$y1 <- d6$x + rnorm(80)
d6$y2 <- d6$x + rnorm(80)
d6$s1 <- seq_len(80) <= 72
d6$s9 <- seq_len(80) > 8
f <- frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s1) ~ x) +
           set_rescor(FALSE), data = d6)
f9 <- frm(bf(y1 | subset(s9) ~ x) + bf(y2 | subset(s9) ~ 1) +
            set_rescor(FALSE), data = d6)
cat("5246 anova(rows 1..72, rows 9..80):", try_s(anova(f, f9)), "\n")
f1b <- frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s1) ~ 1) +
             set_rescor(FALSE), data = d6)
cat("5246 control anova(rows 1..72, rows 1..72):", try_s(anova(f1b, f)),
    "\n")
# univariate: NA rows of different predictors, same count
set.seed(6)
du <- data.frame(x = rnorm(60), z = rnorm(60))
du$y <- du$x + rnorm(60)
du$x[1:4] <- NA
du$z[5:8] <- NA
ux <- frm(y ~ x, data = du)
uz <- frm(y ~ z, data = du)
cat("5246 anova(y ~ x, y ~ z), NA in different rows:",
    try_s(anova(ux, uz)), "\n")
# same rows, different data frames of equal size
du2 <- du; du2$y <- rev(du2$y)
cat("5246 anova(y ~ x on d, y ~ x + z on rev(y)):",
    try_s(anova(frm(y ~ x, data = du[9:60, ]),
                frm(y ~ x + z, data = du2[9:60, ]))), "\n")
cat("5246 control anova(y ~ x, y ~ x + z) on one frame:",
    try_s(anova(frm(y ~ x, data = du[9:60, ]),
                frm(y ~ x + z, data = du[9:60, ]))), "\n")

## item 5353: group-level design under the session's contrasts
set.seed(7)
dc <- data.frame(f = factor(sample(c("a", "b", "c"), 90, TRUE)),
                 g = factor(rep(1:10, each = 9)))
dc$y <- as.numeric(dc$f) + rnorm(10)[dc$g] + rnorm(90)
old <- options(contrasts = c("contr.sum", "contr.poly"))
fc <- suppressWarnings(frm(y ~ f + (1 + f | g), data = dc))
p_sum <- fitted(fc, newdata = dc[1:6, ])[, 1]
options(old)
cat("5353 predict after options() reset:",
    try_s(fitted(fc, newdata = dc[1:6, ])), "\n")
p_tr <- tryCatch(fitted(fc, newdata = dc[1:6, ])[, 1],
                 error = function(e) rep(NA, 6))
cat("  max |diff| against the prediction under the fit's contrasts",
    format(max(abs(p_tr - p_sum)), digits = 4), "\n")
cat("  in-sample under the fit's contrasts vs newdata:",
    format(max(abs(fitted(fc)[1:6, 1] - p_sum)), digits = 4), "\n")

## item 5362: laplace draws, CE posterior_predict with sigma ~ (1 | h)
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  list(full = full, lap = lap)
}
set.seed(5)
n <- 60
dl <- data.frame(x = rnorm(n), h = factor(rep(1:5, each = 12)))
dl$y <- 0.3 * dl$x + rnorm(n, 0, exp(rnorm(5, 0, 0.5))[dl$h])
fl <- suppressWarnings(frm(bf(y ~ x, sigma ~ (1 | h)), data = dl))
pr <- lap_pair(fl)
ce_nd <- data.frame(h = factor(2, levels = 1:5))
r <- try_s(rr <- conditional_effects(pr$lap, effects = "x", re_formula = NULL,
                                     method = "posterior_predict",
                                     conditions = ce_nd))
cat("5362 CE predict, laplace, sigma ~ (1 | h):", r, "\n")
if (exists("rr")) cat("  NA cells", sum(is.na(rr[[1]]$estimate__)), "of",
                      nrow(rr[[1]]), "\n")
cat("5362 posterior_predict on the same draws:",
    try_s(posterior_predict(pr$lap, newdata = dl[1:3, ], re_formula = NULL)),
    "\n")
cat("5362 control CE predict, laplace, re_formula = NA:",
    try_s(conditional_effects(pr$lap, effects = "x", re_formula = NA,
                              method = "posterior_predict")), "\n")
cat("5362 control CE predict, full draws, re_formula = NULL:",
    try_s(conditional_effects(pr$full, effects = "x", re_formula = NULL,
                              method = "posterior_predict",
                              conditions = ce_nd)), "\n")
fln <- suppressWarnings(frm(bf(y ~ x + (1 | h)), data = dl))
pn <- lap_pair(fln)
cat("5362 control CE predict, laplace, mu ~ (1 | h):",
    try_s(conditional_effects(pn$lap, effects = "x", re_formula = NULL,
                              method = "posterior_predict",
                              conditions = ce_nd)), "\n")

## unfiled: print() of draws lists gp columns in stored order
set.seed(9)
dg <- data.frame(x = runif(60), f = factor(rep(c("a", "b"), 30)))
dg$y <- sin(3 * dg$x) + rnorm(60, 0, 0.3)
fg <- suppressWarnings(frm(y ~ gp(x, by = f, k = 5), data = dg))
pg <- lap_pair(fg)
cat("print order: variables()",
    paste(grep("gp", variables(pg$full), value = TRUE), collapse = " "), "\n")
out <- capture.output(print(pg$full))
rn <- sub(" .*$", "", out[-(1:3)])
cat("print order: print()     ",
    paste(grep("gp", rn, value = TRUE), collapse = " "), "\n")
