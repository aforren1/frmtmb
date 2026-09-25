source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# cs() in predict(): my own design. Four categories, a cs() effect that
# is not proportional, a group intercept, and a second cs() response in
# a multivariate fit. Seed 404. z is in units of the binomial Monte Carlo
# error; propagate_error = FALSE so the parameter is held at the MLE.
set.seed(404)
n <- 500
ng <- 25
g <- factor(rep(seq_len(ng), each = n / ng))
x <- rnorm(n); w <- runif(n, 1, 3)
u <- rnorm(ng, 0, 0.3)[g]
sim_ord <- function(cuts, slopes, xx, uu) {
  # sequential draw: category k with p_k = plogis(cut_k + slope_k*x + u)
  K <- length(cuts) + 1
  out <- rep(K, length(xx))
  alive <- rep(TRUE, length(xx))
  for (k in seq_along(cuts)) {
    stop_k <- alive & runif(length(xx)) < plogis(cuts[k] + slopes[k] * xx + uu)
    out[stop_k] <- k
    alive <- alive & !stop_k
  }
  out
}
yo <- sim_ord(c(-0.5, 0.3, 0.2), c(1.4, -1.0, 0.3), x, u)
yo2 <- sim_ord(c(0.2, -0.4), c(-0.8, 0.9), log(w), 0)
y2 <- 0.5 * x + rnorm(n)
d <- data.frame(x, w, g, yo, yo2, y2)
ND <- 20000L
maxz <- function(prop, p) {
  max(abs(prop - p) / sqrt(pmax(p * (1 - p), 1e-12) / ND))
}
# predict() returns category frequencies over the draws for an ordinal
# response (brms's shape); take them against fitted()'s exact probability
cmp <- function(lab, fit, newdata = NULL, ...) {
  set.seed(9)
  pr <- tryCatch(predict(fit, newdata = newdata, ndraws = ND,
                         propagate_error = FALSE, ...),
                 error = function(e) e)
  fi <- tryCatch(fitted(fit, newdata = newdata, ...), error = function(e) e)
  if (inherits(pr, "error") || inherits(fi, "error")) {
    cat(sprintf("%-44s ERROR predict=%s | fitted=%s\n", lab,
                if (inherits(pr, "error")) conditionMessage(pr) else "ok",
                if (inherits(fi, "error")) conditionMessage(fi) else "ok"))
    return(invisible())
  }
  p <- if (length(dim(fi)) == 3L) fi[, "Estimate", ] else fi
  if (is.null(dim(p))) p <- matrix(p, 1)
  if (is.null(dim(pr))) pr <- matrix(pr, 1)
  cat(sprintf("%-44s rows=%d dim(pred)=%s max|z|=%.2f max|diff|=%.4f\n",
              lab, nrow(p), paste(dim(pr), collapse = "x"),
              maxz(pr, p), max(abs(pr - p))))
}
nd <- data.frame(x = c(2.5, -2.5, 0, 1), w = c(1, 3, 2, 1.5),
                 g = factor(c("1", "2", "3", "4"), levels = levels(g)))
nd1 <- nd[1, , drop = FALSE]
ndnew <- data.frame(x = c(2.5, -2.5), w = c(1, 3), g = factor(c("zz", "zz")))
for (fam in list(sratio(), acat(), cratio())) {
  lab <- fam$family
  f <- frm(bf(yo ~ cs(x) + (1 | g)), family = fam, data = d)
  cmp(paste(lab, "in sample"), f)
  cmp(paste(lab, "newdata 4 rows"), f, nd)
  cmp(paste(lab, "newdata 1 row"), f, nd1)
  cmp(paste(lab, "re_formula = NA, 4 rows"), f, nd, re_formula = NA)
  cmp(paste(lab, "re_formula = NA, in sample"), f, NULL, re_formula = NA)
  # an unseen level: predict() draws the level's effect, fitted() holds
  # it at the population value, so these differ by the group variance
  # only; dropping cs() would move the x = 2.5 row by about 0.5
  cmp(paste(lab, "unseen level, allow_new_levels"), f, ndnew,
      allow_new_levels = TRUE)
  f0 <- frm(bf(yo ~ x + (1 | g)), family = fam, data = d)
  cat(sprintf("   same row without cs(): fitted P(Y=1) at x=2.5: cs %.3f, no-cs %.3f\n",
              fitted(f, newdata = nd1)[1, "Estimate", 1],
              fitted(f0, newdata = nd1)[1, "Estimate", 1]))
}

cat("\n-- multivariate, two cs() responses and a gaussian --\n")
mf <- frm(mvbf(bf(yo ~ cs(x) + (1 | g), family = sratio()),
               bf(yo2 ~ cs(log(w)), family = acat()),
               bf(y2 ~ x, family = gaussian()), rescor = FALSE), data = d)
for (r in c("yo", "yo2")) {
  cmp(paste("mv", r, "in sample"), mf, NULL, resp = r)
  cmp(paste("mv", r, "newdata 4 rows"), mf, nd, resp = r)
  cmp(paste("mv", r, "newdata 1 row"), mf, nd1, resp = r)
  cmp(paste("mv", r, "re_formula = NA"), mf, nd, resp = r, re_formula = NA)
}
# the same responses fitted alone must give the same probabilities
s2 <- frm(bf(yo2 ~ cs(log(w))), family = acat(), data = d)
cat(sprintf("mv yo2 vs univariate fitted: max|diff| %.2e\n",
            max(abs(fitted(mf, resp = "yo2")[, "Estimate", ] -
                    fitted(s2)[, "Estimate", ]))))

cat("\n-- cumulative refuses cs() --\n")
print(tryCatch(frm(bf(yo ~ cs(x)), family = cumulative(), data = d),
               error = function(e) conditionMessage(e)))

cat("\n-- cs() on a factor: single-row newdata --\n")
d$fc <- factor(sample(c("a", "b", "c"), n, TRUE))
ff <- tryCatch(frm(bf(yo ~ cs(fc)), family = sratio(), data = d),
               error = function(e) e)
if (inherits(ff, "error")) {
  cat("refused:", conditionMessage(ff), "\n")
} else {
  print(coef(ff) %||% fixef(ff))
  nd3 <- data.frame(fc = factor(c("a", "b", "c")))
  nd3c <- data.frame(fc = factor("c"))
  nd3b <- data.frame(fc = "b")
  cat("fitted at 3 rows (levels a,b,c):\n")
  print(round(fitted(ff, newdata = nd3)[, "Estimate", ], 4))
  cat("fitted at the single row fc = 'c' (factor with one level):\n")
  print(round(fitted(ff, newdata = nd3c)[, "Estimate", ], 4))
  cat("fitted at the single row fc = 'b' (character):\n")
  print(tryCatch(round(fitted(ff, newdata = nd3b)[, "Estimate", ], 4),
                 error = function(e) conditionMessage(e)))
  set.seed(3)
  cat("predict at single row fc = 'c':\n")
  print(round(predict(ff, newdata = nd3c, ndraws = ND,
                      propagate_error = FALSE), 4))
}
