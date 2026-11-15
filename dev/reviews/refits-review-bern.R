# Review of lane refits, item 1687: how refit(), update(newdata =),
# influence() and frm_bootstrap() read a bernoulli response under
# several codings. Seed 9, n = 40.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
set.seed(9)
n <- 40
x <- rnorm(n)
s <- rbinom(n, 1, plogis(0.3 + x))     # successes, 0/1
codings <- list(
  zero_one = s,
  logical = s == 1,
  one_two = s + 1,                    # failure 1, success 2
  minus_one_one = 2 * s - 1,          # failure -1, success 1
  minus_one_zero = s - 1,             # failure -1, success 0
  factor_noyes = factor(ifelse(s == 1, "yes", "no")),
  char_one_two = as.character(s + 1)
)
cf <- function(f) if (is.character(f)) f else round(fixef(f, flatten = TRUE), 6)
for (nm in names(codings)) {
  cat("\n==========", nm, "==========\n")
  d <- data.frame(x = x)
  d$y <- codings[[nm]]
  fit <- go(paste(nm, "fit"), frm(y ~ x, data = d, family = bernoulli()))
  if (is.character(fit)) next
  lv <- fit$spec$responses[[1L]]$family[["bin_levels"]]
  cat("bin_levels:", lv, "\n")
  go("fit coef", cf(fit))
  go("refit(own values)", cf(refit(fit, d$y)))
  go("refit(codes s)", cf(refit(fit, s)))
  # flip every response: as values and as codes
  flip_val <- codings[[nm]][c(n:1)]
  go("refit(own values reversed) vs refit(codes reversed)",
     list(cf(refit(fit, flip_val)), cf(refit(fit, s[n:1]))))
  # a vector of one value: every row the same number
  for (v in c(0, 1)) {
    go(paste0("refit(rep(", v, ", n)) intercept sign"),
       { r <- refit(fit, rep(v, n)); sign(fixef(r, flatten = TRUE)[1]) })
    d2 <- d
    d2$y <- if (is.factor(d$y)) d$y else rep(v, n)
    if (!is.factor(d$y)) {
      go(paste0("update(newdata y = ", v, ") intercept sign"),
         { r <- update(fit, newdata = d2); sign(fixef(r, flatten = TRUE)[1]) })
    }
  }
  # mixtures of codes and values
  if (nm %in% c("one_two", "minus_one_one", "minus_one_zero")) {
    mix <- ifelse(seq_len(n) %% 2 == 0, s, codings[[nm]])
    go("refit(mix codes/values)", cf(refit(fit, mix)))
  }
}

cat("\n========== bootstrap, near-all-success data ==========\n")
set.seed(21)
n2 <- 12
x2 <- rnorm(n2)
s2 <- c(rep(1, 11), 0)
d01 <- data.frame(x = x2, y = s2)
d12 <- data.frame(x = x2, y = s2 + 1)
f01 <- frm(y ~ x, data = d01, family = bernoulli())
f12 <- frm(y ~ x, data = d12, family = bernoulli())
sims <- simulate(f12, nsim = 40, seed = 5)
allone <- vapply(sims, function(v) all(v == 1), NA)
cat("simulate() draws all 1:", sum(allone), "of 40; values in draws:",
    paste(sort(unique(unlist(sims))), collapse = ","), "\n")
b01 <- suppressWarnings(frm_bootstrap(f01, nsim = 40, seed = 5))
b12 <- suppressWarnings(frm_bootstrap(f12, nsim = 40, seed = 5))
cat("bootstrap t identical 0/1 vs 1/2:",
    identical(b01$t, b12$t), " max|diff|:",
    max(abs(b01$t - b12$t), na.rm = TRUE), "\n")

cat("\n========== influence(), 1/2 against 0/1 ==========\n")
set.seed(9)
d01 <- data.frame(x = x, y = s)
d12 <- data.frame(x = x, y = s + 1)
i01 <- influence(frm(y ~ x, data = d01, family = bernoulli()))
i12 <- influence(frm(y ~ x, data = d12, family = bernoulli()))
cat("influence fixed tables max|diff|:",
    max(abs(unlist(i01[["fixed.effects[-case]"]]) -
            unlist(i12[["fixed.effects[-case]"]]))), "\n")
cat("names(i01):", paste(names(i01), collapse = ", "), "\n")
