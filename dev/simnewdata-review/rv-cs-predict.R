# predict() at newdata with an ordinal cs() term: does it use the new
# rows' covariate? Measured on whichever build RV_LIB selects.
#   RV_LIB=base Rscript dev/simnewdata-review/rv-cs-predict.R
# The model's only x effect is cs(x), so a predict() that took the cs()
# offsets from the FITTED rows would ignore newdata's x entirely.
source("dev/simnewdata-review/rv-prelude.R")
set.seed(5)
n <- 300
d <- data.frame(x = rnorm(n))
d$yo <- factor(cut(2 * d$x + rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                   labels = FALSE), ordered = TRUE)
nd <- data.frame(x = c(3, 3, -3))
cat("fitted rows 1..3 have x =", format(d$x[1:3], digits = 3), "\n")
for (fam in c("sratio", "acat")) {
  fo <- frm(bf(yo ~ cs(x)), family = get(fam, asNamespace("frmtmb"))(),
            data = d)
  nd_draws <- 20000
  pnd <- unclass(predict(fo, newdata = nd, ndraws = nd_draws))
  pin <- unclass(predict(fo, ndraws = nd_draws))[1:3, ]
  fnd <- fitted(fo, newdata = nd)
  ex <- sapply(seq_len(dim(fnd)[3]), function(k) fnd[, "Estimate", k])
  cat("\nfamily", fam, "\n")
  cat("  exact P(Y = k) at newdata x = 3, 3, -3 (fitted(newdata)):\n")
  print(round(ex, 4))
  cat("  predict(newdata = x 3, 3, -3), ", nd_draws, " draws:\n", sep = "")
  print(round(pnd, 4))
  cat("  predict() on fitted rows 1..3:\n")
  print(round(pin, 4))
  se <- sqrt(ex * (1 - ex) / nd_draws)
  cat(sprintf("  max |z| predict(newdata) vs exact: %.1f; vs fitted rows 1..3: %.1f\n",
              max(abs(pnd - ex) / se), max(abs(pnd - pin) / se)))
  big <- data.frame(x = rep(2, n + 5))
  cat("  newdata with n + 5 rows:",
      try_msg(predict(fo, newdata = big, ndraws = 10)), "\n")
  if (!rv_base) {
    s <- as.matrix(simulate(fo, nsim = nd_draws, seed = 1, newdata = nd))
    fr <- t(apply(s, 1, function(r) {
      table(factor(r, levels = levels(d$yo))) / nd_draws
    }))
    cat(sprintf("  lane simulate(newdata) vs exact: max |z| %.1f\n",
                max(abs(fr - ex) / se)))
  }
}
