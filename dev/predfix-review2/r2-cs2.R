source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# cs() follow-up: cumulative refusal, and cs() on a factor at a
# single-row newdata. Seed 405.
set.seed(405)
n <- 500
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(x, fc, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
ND <- 20000L

cat("-- cumulative refuses cs() --\n")
print(tryCatch(frm(bf(yo ~ cs(x)), family = cumulative(), data = d),
               error = function(e) conditionMessage(e)))
cat("-- multivariate with an ordinal response --\n")
d$y2 <- rnorm(n)
print(tryCatch(frm(mvbf(bf(yo ~ cs(x), family = sratio()),
                        bf(y2 ~ x, family = gaussian()), rescor = FALSE),
                   data = d), error = function(e) conditionMessage(e)))

cat("\n-- cs() on a factor --\n")
ff <- tryCatch(frm(bf(yo ~ cs(fc)), family = sratio(), data = d),
               error = function(e) e)
if (inherits(ff, "error")) {
  cat("refused:", conditionMessage(ff), "\n")
} else {
  print(fixef(ff)[, 1:2])
  nd3 <- data.frame(fc = factor(c("a", "b", "c")))
  cat("fitted at 3 rows, levels a, b, c:\n")
  print(round(fitted(ff, newdata = nd3)[, "Estimate", ], 4))
  cat("fitted at the single row fc = factor('c'):\n")
  print(round(fitted(ff, newdata = data.frame(fc = factor("c")))[, "Estimate", ], 4))
  cat("fitted at the single row fc = factor('c', levels = c('a','b','c')):\n")
  print(round(fitted(ff, newdata = data.frame(
    fc = factor("c", levels = c("a", "b", "c"))))[, "Estimate", ], 4))
  cat("predict at the single row fc = factor('c'):\n")
  set.seed(3)
  print(round(predict(ff, newdata = data.frame(fc = factor("c")), ndraws = ND,
                      propagate_error = FALSE), 4))
  cat("empirical category shares by level in the data:\n")
  print(round(prop.table(table(d$fc, d$yo), 1), 4))
  cat("same model with dummy columns, cs(fb) + cs(fcc):\n")
  d$fb <- as.numeric(d$fc == "b"); d$fcc <- as.numeric(d$fc == "c")
  fd <- frm(bf(yo ~ cs(fb) + cs(fcc)), family = sratio(), data = d)
  cat(sprintf("logLik cs(fc) %.4f | dummies %.4f\n",
              as.numeric(logLik(ff)), as.numeric(logLik(fd))))
}
