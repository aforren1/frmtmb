source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# Item 6: vcov() on REML / profile fits with no outer parameter, which
# errored on 0.62.0. Reference 1: the ML fit's vcov() on the same data
# (the ML outer Hessian is the same observed information). Reference 2:
# my own numerical Hessian of a hand-written log-likelihood. Seed 707.
set.seed(707)
n <- 300
f <- factor(sample(letters[1:4], n, TRUE))
x <- rnorm(n)
eta <- c(a = 0.2, b = 0.8, c = -0.3, d = 1.1)[as.character(f)] + 0.4 * x
d <- data.frame(f, x, yg = rgeom(n, 1 / (1 + exp(eta))),
                yp = rpois(n, exp(eta)))
X <- model.matrix(~ 0 + f + x, d)
nll_geo <- function(b) {
  mu <- exp(drop(X %*% b)); -sum(dgeom(d$yg, 1 / (1 + mu), log = TRUE))
}
nll_poi <- function(b) -sum(dpois(d$yp, exp(drop(X %*% b)), log = TRUE))
numhess <- function(fn, b, h = 1e-4) {
  k <- length(b); H <- matrix(0, k, k)
  for (i in 1:k) for (j in 1:k) {
    ei <- replace(numeric(k), i, h); ej <- replace(numeric(k), j, h)
    H[i, j] <- (fn(b + ei + ej) - fn(b + ei - ej) - fn(b - ei + ej) +
                  fn(b - ei - ej)) / (4 * h * h)
  }
  H
}
for (fam in c("geometric", "poisson")) {
  y <- if (fam == "geometric") "yg" else "yp"
  fo <- as.formula(paste(y, "~ 0 + f + x"))
  famo <- if (fam == "geometric") geometric() else poisson()
  ml <- frm(fo, family = famo, data = d)
  Vn <- solve(numhess(if (fam == "geometric") nll_geo else nll_poi,
                      unname(fixef(ml)[, "Estimate"])))
  for (mode in c("REML", "profile")) {
    fit <- if (mode == "REML") frm(fo, family = famo, data = d, REML = TRUE)
           else frm(fo, family = famo, data = d,
                    control = frmtmb_control(profile = TRUE))
    V <- tryCatch(vcov(fit), error = function(e) e)
    if (inherits(V, "error")) { cat(fam, mode, "ERROR", conditionMessage(V), "\n"); next }
    cat(sprintf("%-10s %-8s dim %s names %s\n", fam, mode,
                paste(dim(V), collapse = "x"),
                paste(rownames(V), collapse = ",")))
    cat(sprintf("   vs ML vcov: max rel %.2e | vs numerical Hessian: max rel %.2e | fixef vs ML rel %.2e\n",
                rel(V, vcov(ml)), rel(V, Vn),
                rel(fixef(fit)[, "Estimate"], fixef(ml)[, "Estimate"])))
    s <- tryCatch(summary(fit), error = function(e) e)
    cat("   summary():", if (inherits(s, "error")) conditionMessage(s) else "ok",
        " se from fixef vs sqrt(diag(vcov)):",
        sprintf("%.2e", rel(fixef(fit)[, "Est.Error"], sqrt(diag(V)))), "\n")
  }
}
cat("\n-- the naming guard's other side: a model WITH an outer parameter --\n")
fitg <- frm(yp ~ 0 + f + x + (1 | f), family = poisson(), data = d, REML = TRUE)
print(dim(vcov(fitg)))
