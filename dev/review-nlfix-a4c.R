# Area 4: simulate(newdata) and frm_curve under cov = FALSE, nl against linear.
source("dev/review-nlfix-common.R")
set.seed(4401)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)), t = rep(1:20, 10))
d$y <- 1 + 0.5 * d$x + ave(rnorm(n, 0, 0.5), d$g, FUN = function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
st <- list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.5))
ter <- "ar(time = t, gr = g, cov = FALSE)"
fn <- suppressWarnings(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d, start = st))
fl <- suppressWarnings(frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d))
b <- fn$estimates$beta
nd <- d[1:30, ]
body_nd <- b[[1]] + b[[2]] * nd$x
sim_nl <- as.matrix(simulate(fn, nsim = 4000, newdata = nd, seed = 2))
sim_li <- as.matrix(simulate(fl, nsim = 4000, newdata = nd, seed = 2))
sim_nd2 <- nd; sim_nd2$y <- NULL
cat("simulate(newdata = rows with y): row means minus body\n")
rs <- cbind(t = nd$t, nl = rowMeans(sim_nl) - body_nd, lin = rowMeans(sim_li) - body_nd, shift_fitted = fitted(fn, newdata = nd)[, "Estimate"] - body_nd)
print(round(rs[c(1:6, 21:24), ], 4))
cat("max |nl - lin| of row means:", max(abs(rs[, "nl"] - rs[, "lin"])), " MC se about", 0.5 / sqrt(4000), "\n")
cat("slope of nl's row-mean excess on the observed shift:", coef(lm(rs[, "nl"] ~ rs[, "shift_fitted"]))[2], "\n")
cat("slope of lin's row-mean excess on the observed shift:", coef(lm(rs[, "lin"] ~ rs[, "shift_fitted"]))[2], "\n")
x <- tt(as.matrix(simulate(fn, nsim = 3, newdata = sim_nd2, seed = 2))); show("simulate(newdata without y), nl", if (iserr(x)) x else dim(x))
x <- tt(as.matrix(simulate(fl, nsim = 3, newdata = sim_nd2, seed = 2))); show("simulate(newdata without y), lin", if (iserr(x)) x else dim(x))
cat("\nfrm_curve\n")
cn <- tt(as.data.frame(frm_curve(fn, newdata = nd)))
cl <- tt(as.data.frame(frm_curve(fl, newdata = nd)))
show("curve nl", if (iserr(cn)) cn else 1); show("curve lin", if (iserr(cl)) cl else 1)
if (!iserr(cn)) {
  cat("curve nl .estimate - body: max", max(abs(cn$.estimate - body_nd)), "\n")
  cat("curve nl .estimate - fitted (one-step): max", max(abs(cn$.estimate - fitted(fn, newdata = nd)[, "Estimate"])), "\n")
  print(head(cn[, c(".estimate", ".se")]))
}
lp <- tt(frm_linpred(fn, newdata = nd, se.fit = TRUE)); show("linpred se.fit nl", if (iserr(lp)) lp else 1)
lpl <- tt(frm_linpred(fl, newdata = nd, se.fit = TRUE)); show("linpred se.fit lin", if (iserr(lpl)) lpl else 1)
if (!iserr(lp)) print(head(as.data.frame(lp)))
f0 <- tt(fitted(fn)); f1 <- tt(fitted(fl))
print(head(f0, 3)); print(head(f1, 3))
cat("Est.Error NA in nl:", sum(is.na(f0[, "Est.Error"])), " in lin:", sum(is.na(f1[, "Est.Error"])), "\n")
r0 <- tt(residuals(fn)); r1 <- tt(residuals(fl)); print(head(r0, 3)); print(head(r1, 3))
fc <- suppressWarnings(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)), data = d, start = st))
flc <- suppressWarnings(frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = TRUE)), data = d))
print(head(tt(fitted(fc)), 3)); print(head(tt(fitted(flc)), 3))
cat("\nthe same fit with one more check: frm_linpred(dpar = 'mu') vs. simulate's mean reader\n")
ns <- asNamespace("frmtmb")
cat("\nordinal categorical: curve on a body\n")
print(tt(frm_compat("ar", "nl")))
print(tt(frm_compat("group:autocor", "nl")))
s <- tt(frm_compat()); if (!iserr(s)) { s <- as.data.frame(s); print(s[s$feature_a == "nl" & grepl("ar|autocor", s$feature_b) | (s$feature_b == "nl" & s$feature_a %in% c("ar", "ma", "arma", "cosy", "unstr", "fcor")), c("feature_a", "feature_b", "status")]) }
cat("\nbrms: + acformula twice\n")
.pick <- function(txt) { r <- brms::stancode(eval(parse(text = txt), envir = asNamespace("brms")), data = d); ln <- strsplit(as.character(r), "\n")[[1]]; paste(grep("cosy|cholesky_cor_ar1|cholesky_cor_cosy", ln, value = TRUE)[1:2], collapse = " | ") }
cat("nl:", .pick("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE) + acformula(~ ar(time = t, gr = g, cov = TRUE)) + acformula(~ cosy(time = t, gr = g))"), "\n")
cat("lin:", .pick("bf(y ~ x) + acformula(~ ar(time = t, gr = g, cov = TRUE)) + acformula(~ cosy(time = t, gr = g))"), "\n")
cat("bf(bf(), autocor =):", .pick("bf(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = TRUE)), autocor = ~ cosy(time = t, gr = g))"), "\n")
cat("nlf(mu) with acformula:", substr(paste(grep("cholesky|chol_cor =|ar", strsplit(as.character(brms::stancode(brms::bf(y ~ 1) + brms::nlf(mu ~ a + b * x) + brms::lf(a ~ 1) + brms::lf(b ~ 1) + brms::acformula(~ ar(time = t, gr = g, cov = TRUE)), data = d)), "\n")[[1]], value = TRUE)[1:4], collapse = " | "), 1, 300), "\n")
