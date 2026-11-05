# Short check of fix round 2, item 2: when does frmtmb_interop_se_part fire, what does it
# say, and does re_formula = NA agree as the wording claims? Against brms sampled, seed 7411.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
.libPaths(c(.libPaths()[1L], "C:/Users/adf44/source/r/rev-predread-lib/extra", .libPaths()[-1L]))
set.seed(7411)
J <- 12
m <- 15
n <- J * m
d <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(J)), each = m)), x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)))
u <- rnorm(J, 0, 0.7)[d$g]
d$cnt <- rpois(n, exp(0.4 + 0.5 * d$x + 0.4 * (d$f == "b") + u))
d$y <- 0.4 + 0.5 * d$x + 0.4 * (d$f == "b") + u + rnorm(n, 0, 0.8)
mix <- quiet(frm(bf(cnt ~ x + f + (1 | g)), family = poisson(), data = d))
plain <- quiet(frm(bf(cnt ~ x + f), family = poisson(), data = d))
nd <- d[c(1, 20, 40), ]
count <- function(lab, e) {
  w <- character()
  r <- tryCatch(withCallingHandlers(e, warning = function(c) {
    if (inherits(c, "frmtmb_interop_se_part")) w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")), error = function(err) paste("ERROR", substr(conditionMessage(err), 1, 70)))
  cat(sprintf("%-44s warned %d  %s\n", lab, length(w), if (is.character(r) && length(r) == 1L) r else ""))
  invisible(w)
}
grDevices::pdf(NULL)
cat("=== mixed poisson fit\n")
w1 <- count("predictions() default", marginaleffects::predictions(mix, newdata = nd))
count("predictions(vcov = FALSE)", marginaleffects::predictions(mix, newdata = nd, vcov = FALSE))
count("predictions(re_formula = NA)", marginaleffects::predictions(mix, newdata = nd, re_formula = NA))
count("predictions(re_formula = ~ (1 | g))", marginaleffects::predictions(mix, newdata = nd, re_formula = ~ (1 | g)))
count("avg_predictions()", marginaleffects::avg_predictions(mix))
count("avg_slopes()", marginaleffects::avg_slopes(mix))
count("avg_slopes(vcov = FALSE)", marginaleffects::avg_slopes(mix, vcov = FALSE))
count("avg_comparisons(variables = 'f')", marginaleffects::avg_comparisons(mix, variables = "f"))
count("hypotheses('b2 = 0')", marginaleffects::hypotheses(mix, hypothesis = "b2 = 0"))
count("hypotheses(vcov = FALSE)", marginaleffects::hypotheses(mix, hypothesis = "b2 = 0", vcov = FALSE))
count("get_predict()", marginaleffects::get_predict(mix, newdata = nd))
count("get_coef()", marginaleffects::get_coef(mix))
count("set_coef()", marginaleffects::set_coef(mix, marginaleffects::get_coef(mix)))
count("get_vcov() direct", marginaleffects::get_vcov(mix))
count("datagrid()", marginaleffects::datagrid(model = mix, x = 1))
count("plot_predictions(condition = 'x')", marginaleffects::plot_predictions(mix, condition = "x"))
count("plot_predictions(draw = FALSE, vcov = FALSE)", marginaleffects::plot_predictions(mix, condition = "x", draw = FALSE, vcov = FALSE))
count("plot_slopes(variables = 'x', condition = 'f')", marginaleffects::plot_slopes(mix, variables = "x", condition = "f"))
count("predictions(vcov = matrix)", marginaleffects::predictions(mix, newdata = nd, vcov = as.matrix(marginaleffects::get_vcov(plain) * 1)))
count("emmeans()", as.data.frame(emmeans::emmeans(mix, ~ f)))
count("insight::get_predicted()", insight::get_predicted(mix, data = nd))
count("insight::get_varcov()", insight::get_varcov(mix))
count("fitted(newdata)", fitted(mix, newdata = nd))
count("summary()", summary(mix))
cat("=== plain poisson fit\n")
count("predictions() default", marginaleffects::predictions(plain, newdata = nd))
count("avg_slopes()", marginaleffects::avg_slopes(plain))
cat("\nwording:\n", paste(strwrap(w1[1], 90), collapse = "\n "), "\n")
cat("\n=== the wording's claim: re_formula = NA agrees\n")
pna <- quiet(marginaleffects::predictions(mix, newdata = nd, re_formula = NA))
lna <- frm_linpred(mix, newdata = nd, re_formula = NA, type = "link", se.fit = TRUE)
fna <- fitted(mix, newdata = nd, re_formula = NA)
cat("predictions se ", signif(pna$std.error, 4), "\nfitted Est.Error", signif(fna[, "Est.Error"], 4), "\n")
gnd <- d[c(1, 20, 40), ]
gna <- quiet(frm(bf(y ~ x + f + (1 | g)), family = gaussian(), data = d))
cat("gaussian, predictions(re_formula = NA) se", signif(quiet(marginaleffects::predictions(gna, newdata = nd, re_formula = NA))$std.error, 4),
    "\n          frm_linpred(re_formula = NA) se ", signif(frm_linpred(gna, newdata = nd, re_formula = NA, se.fit = TRUE)$se.fit, 4), "\n")
cat("\n=== sampled brms, flat default priors, the same data: is the warning's direction true?\n")
b <- quiet(brms::brm(brms::bf(cnt ~ x + f + (1 | g)), family = stats::poisson(), data = d, chains = 2, iter = 1500, refresh = 0, seed = 3, silent = 2))
be <- posterior_epred(b, newdata = nd)
cat("level at seen group: brms sd", signif(apply(be, 2, sd), 4), "\n")
cat("  frmtmb fitted se  ", signif(fitted(mix, newdata = nd)[, "Est.Error"], 4), "\n")
cat("  marginaleffects se", signif(quiet(marginaleffects::predictions(mix, newdata = nd))$std.error, 4), "\n")
ds <- function(f) {
  dh <- d
  dl <- d
  dh$x <- d$x + 1e-3
  dl$x <- d$x - 1e-3
  rowMeans((posterior_epred(b, newdata = dh) - posterior_epred(b, newdata = dl)) / 2e-3)
}
sb <- ds()
da <- d
db <- d
da$f <- factor("a", levels = levels(d$f))
db$f <- factor("b", levels = levels(d$f))
cb <- rowMeans(posterior_epred(b, newdata = db) - posterior_epred(b, newdata = da))
ap <- rowMeans(posterior_epred(b))
cat("avg_predictions: brms sd", signif(sd(ap), 4), " marginaleffects", signif(quiet(marginaleffects::avg_predictions(mix))$std.error, 4), "\n")
cat("avg_slopes(x):   brms sd", signif(sd(sb), 4), " marginaleffects", signif(quiet(marginaleffects::avg_slopes(mix, variables = "x"))$std.error, 4), "\n")
cat("avg_comparisons(f): brms sd", signif(sd(cb), 4), " marginaleffects", signif(quiet(marginaleffects::avg_comparisons(mix, variables = "f"))$std.error, 4), "\n")
cat("brms's own marginaleffects avg_slopes(x): CI", signif(unlist(quiet(marginaleffects::avg_slopes(b, variables = "x"))[c("conf.low", "conf.high")]), 4), "\n")
cat("frmtmb coefficient of x, link scale, se:", signif(summary(mix)$fixed["x", 1:2], 4), "\n")
