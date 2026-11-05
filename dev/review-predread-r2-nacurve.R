# Short check of fix round 2, item 5: frm_curve() with a row of NA at its default
# simultaneous = TRUE, plain predictor against me(), both builds. Seed 7422.
#   REV_LIB=before Rscript dev/review-predread-r2-nacurve.R
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7422)
n <- 200
d <- data.frame(z = rnorm(n), x = rnorm(n), sx = 0.4)
d$xo <- d$x + rnorm(n, 0, 0.4)
d$y <- 1 + d$x + rnorm(n)
d$ord <- as.integer(cut(d$x + ifelse(d$z > 0, 1, -1) * rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
sh <- function(lab, e) {
  r <- tryCatch(suppressWarnings(suppressMessages(e)), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 100)))
  cat(sprintf("%-36s %s\n", lab, if (is.character(r)) r else paste(signif(r, 3), collapse = " ")))
}
pl <- quiet(frm(bf(y ~ x + z), family = gaussian(), data = d))
nd <- data.frame(x = c(1, NA, 0), z = 0)
sh("plain simultaneous = FALSE", frm_curve(pl, newdata = nd, simultaneous = FALSE)$.se)
sh("plain default", frm_curve(pl, newdata = nd, seed = 1)$.se)
sh("plain default crit_sim", frm_curve(pl, newdata = nd, seed = 1)$.crit_sim)
cs <- quiet(frm(bf(ord ~ x + cs(z)), family = sratio(), data = d))
sh("cs() contrast with NA in x", frm_curve(cs, newdata = nd, contrast = transform(nd, z = 1), simultaneous = FALSE)$.se)
ct <- quiet(frm(bf(ord ~ x), family = categorical(), data = d))
sh("categorical contrast with NA in x", frm_curve(ct, newdata = nd, contrast = transform(nd, x = x + 1), simultaneous = FALSE)$.se)
if (!identical(Sys.getenv("REV_LIB"), "before")) {
  me <- quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = d))
  nm <- data.frame(xo = c(1, NA, 0), sx = 0.4, z = 0)
  sh("me() simultaneous = FALSE", frm_curve(me, newdata = nm, simultaneous = FALSE)$.se)
  sh("me() default", frm_curve(me, newdata = nm, seed = 1)$.se)
  mc <- quiet(frm(bf(ord ~ me(xo, sx)), family = categorical(), data = d))
  sh("categorical me() contrast", frm_curve(mc, newdata = transform(nm, z = NULL), contrast = transform(nm, z = NULL, xo = xo + 1), simultaneous = FALSE)$.se)
}
