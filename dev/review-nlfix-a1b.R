# Area 1: b1 ~ 1 + z, b0 + b1 x ridge, nonlinear ridge, thres(gr), mv.
source("dev/review-nlfix-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(4102)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:15, each = 20)))
lat <- 0.8 * d$x + 0.5 * d$x * d$z + rlogis(n)
d$y <- cut(lat, c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
d$gg <- factor(rep(c("a", "b"), length.out = n))
nd <- d[1:7, ]; nd$x <- nd$x + 1
cat("\n--- b1 ~ 1 + z vs x + x:z\n")
fn <- tt(frm(bf(y ~ b1 * x, b1 ~ 1 + z, nl = TRUE), family = cumulative(), data = d))
fl <- tt(frm(y ~ 0 + x + x:z, family = cumulative(), data = d))
cat("logLik", as.numeric(logLik(fn)), as.numeric(logLik(fl)), "\n")
print(unlist(fn$estimates)); print(unlist(fl$estimates))
for (nm in c("fitted", "predict")) {
  a <- tt(get(nm)(fn, newdata = nd)); b <- tt(get(nm)(fl, newdata = nd))
  cat(nm, "max diff", max(abs(as.numeric(unlist(a)) - as.numeric(unlist(b)))), "\n")
}
a <- tt(frm_linpred(fn, newdata = nd)); b <- tt(frm_linpred(fl, newdata = nd)); cat("linpred", max(abs(a - b)), "\n")
set.seed(3); a <- tt(residuals(fn, type = "osa")); set.seed(3); b <- tt(residuals(fl, type = "osa")); print(if (iserr(a)) a$msg else max(abs(as.numeric(a) - as.numeric(b))))
a <- as.data.frame(frm_curve(fn, newdata = nd)); b <- as.data.frame(frm_curve(fl, newdata = nd)); print(names(a)); print(a[1:2, ]); print(b[1:2, ])
cat("\n--- body with its own intercept: b0 + b1 * x\n")
fb <- tt(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = cumulative(), data = d))
show("fb", if (iserr(fb)) fb else structure(as.numeric(logLik(fb)), w = attr(fb, "w"), m = attr(fb, "m")))
if (!iserr(fb)) {
  print(unlist(fb$estimates)); print(tt(summary(fb)))
  print(fb$opt$convergence); print(fb$opt$message)
  print(tt(sqrt(diag(vcov(fb)))))
}
for (fam in c("sratio", "cratio", "acat", "hurdle_cumulative")) {
  dd <- d
  if (fam == "hurdle_cumulative") dd$y <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
  r <- tt(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = get(fam)(), data = dd))
  show(paste("ridge", fam), if (iserr(r)) r else structure(as.numeric(logLik(r)), w = attr(r, "w"), m = attr(r, "m")))
}
cat("\n--- nonlinear ridge: b1 * exp(b2 * x)\n")
d2 <- d; set.seed(4935)
d2$x <- rnorm(n); d2$z <- rnorm(n)
d2$y <- cut(0.8 * d2$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
for (fam in c("cumulative", "sratio")) {
  r <- tt(frm(bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE), family = get(fam)(), data = d2,
              start = list(beta = c("b1_(Intercept)" = 0.5, "b2_(Intercept)" = 0.5))))
  if (iserr(r)) show(fam, r) else {
    cat(fam, "ll", as.numeric(logLik(r)), "conv", r$opt$convergence, r$opt$message, "\n")
    print(unlist(r$estimates))
    cat(" warn:", substr(attr(r, "w"), 1, 300), sep = "\n  "); cat(" msg:", substr(attr(r, "m"), 1, 300), sep = "\n  ")
    cat("\n")
    print(tt(summary(r)$fixed))
  }
  r <- tt(frm(bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE), family = get(fam)(), data = d2))
  if (iserr(r)) show(paste(fam, "defstart"), r) else {
    cat(fam, "default start: ll", as.numeric(logLik(r)), "conv", r$opt$convergence, "\n"); print(unlist(r$estimates))
    cat(" warn:", substr(attr(r, "w"), 1, 300), sep = "\n  "); cat("\n")
  }
}
cat("\n--- thres gr\n")
d3 <- d
r <- tt(frm(bf(y | thres(gr = gg) ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d3))
l <- tt(frm(y | thres(gr = gg) ~ 0 + x, family = cumulative(), data = d3))
show("nl gr", if (iserr(r)) r else as.numeric(logLik(r))); show("lin gr", if (iserr(l)) l else as.numeric(logLik(l)))
if (!iserr(r) && !iserr(l)) {
  a <- tt(fitted(r, newdata = d3[1:9, ])); b <- tt(fitted(l, newdata = d3[1:9, ]))
  if (iserr(a) || iserr(b)) { print(a); print(b) } else cat("gr fitted diff", max(abs(a - b)), "\n")
  cat("gr linpred diff", max(abs(frm_linpred(r, newdata = d3[1:9, ]) - frm_linpred(l, newdata = d3[1:9, ]))), "\n")
}
cat("\n--- thres equidistant\n")
r <- tt(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(threshold = "equidistant"), data = d))
l <- tt(frm(y ~ 0 + x, family = cumulative(threshold = "equidistant"), data = d))
cat(as.numeric(logLik(r)), as.numeric(logLik(l)), "\n")
cat("\n--- RE in nlpar\n")
r <- tt(frm(bf(y ~ b1 * x, b1 ~ 1 + (1 | g), nl = TRUE), family = cumulative(), data = d))
l <- tt(frm(y ~ 0 + x + (0 + x | g), family = cumulative(), data = d))
cat(as.numeric(logLik(r)), as.numeric(logLik(l)), "\n")
print(tt(max(abs(fitted(r, newdata = nd) - fitted(l, newdata = nd)))))
print(tt(max(abs(fitted(r, newdata = nd, re_formula = NA) - fitted(l, newdata = nd, re_formula = NA)))))
ndn <- nd; ndn$g <- factor("new")
x <- tt(fitted(r, newdata = ndn, allow_new_levels = TRUE)); if (iserr(x)) print(x$msg) else print(dim(x))
x <- tt(ranef(r)); print(if (iserr(x)) x$msg else "ranef ok")
cat("\n--- offset as a data column in the body\n")
d$o <- rep(c(0.5, -0.5), length.out = n)
nd <- d[1:7, ]; nd$x <- nd$x + 1
r <- tt(frm(bf(y ~ b1 * x + o, b1 ~ 1, nl = TRUE), family = cumulative(), data = d))
l <- tt(frm(y ~ 0 + x + offset(o), family = cumulative(), data = d))
cat(as.numeric(logLik(r)), as.numeric(logLik(l)), "\n")
print(tt(max(abs(fitted(r, newdata = nd) - fitted(l, newdata = nd)))))
cat("\n--- integer y\n")
d$yi <- as.integer(d$y)
r <- tt(frm(bf(yi ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d))
l <- tt(frm(yi ~ 0 + x, family = cumulative(), data = d))
show("int y", if (iserr(r)) r else as.numeric(logLik(r))); show("int y lin", if (iserr(l)) l else as.numeric(logLik(l)))
cat("\n--- b1 ~ 0\n")
r <- tt(frm(bf(y ~ b1 * x, b1 ~ 0, nl = TRUE), family = cumulative(), data = d))
show("b ~ 0", if (iserr(r)) r else as.numeric(logLik(r)))
l <- tt(frm(y ~ 1, family = cumulative(), data = d)); cat("int only:", as.numeric(logLik(l)), "\n")
cat("\n--- mv with an ordinal body and a gaussian body\n")
d$y2 <- rnorm(n) + d$x
r <- tt(frm(mvbf(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE, family = cumulative()),
                 bf(y2 ~ a + c * x, a + c ~ 1, nl = TRUE, family = gaussian()), rescor = FALSE), data = d))
show("mv", if (iserr(r)) r else as.numeric(logLik(r)))
if (!iserr(r)) {
  s <- sum(as.numeric(logLik(frm(y ~ 0 + x, family = cumulative(), data = d))), as.numeric(logLik(frm(y2 ~ x, data = d))))
  cat("sum separate", s, "\n")
  x <- tt(fitted(r, resp = "y", newdata = nd)); print(if (iserr(x)) x$msg else dim(x))
}
