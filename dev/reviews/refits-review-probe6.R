# Review of lane refits, 3188: which fits get marginaleffects_modeldata,
# and predictions() on cbind() binomial and mvbind() fits, before/after.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
set.seed(3)
n <- 80
d <- data.frame(x = rnorm(n), N = 12L)
d$s <- rbinom(n, 12, plogis(0.2 + 0.5 * d$x)); d$f <- 12L - d$s
d$y1 <- rnorm(n, d$x); d$y2 <- rnorm(n, -d$x)
fb <- frm(cbind(s, f) ~ x, data = d, family = binomial())
fmv <- frm(bf(mvbind(y1, y2) ~ x), data = d)
for (nm in c("fb", "fmv")) {
  fit <- get(nm)
  md <- attr(fit, "marginaleffects_modeldata")
  go(paste(nm, "attribute columns"), if (is.null(md)) "NULL" else names(md))
  go(paste(nm, "predictions() head"), {
    p <- as.data.frame(marginaleffects::predictions(fit))
    list(dim = dim(p), head = head(p[, intersect(c("group", "estimate", "std.error"), names(p))], 3))
  })
  go(paste(nm, "avg_predictions(by = NULL)"),
     as.data.frame(marginaleffects::avg_predictions(fit))[, intersect(c("group", "estimate", "std.error"), names(as.data.frame(marginaleffects::avg_predictions(fit))))])
}
go("glm reference for cbind binomial avg prediction",
   { g <- glm(cbind(s, f) ~ x, data = d, family = binomial())
     as.data.frame(marginaleffects::avg_predictions(g))[, c("estimate", "std.error")] })
