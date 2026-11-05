# Small checks left: prior_summary's source column after update(), brms's source column
# of a table's edited row, and predict(newdata) under cov = TRUE on the merged tree.
# Data seed 7601.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7601)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 20)), t = rep(1:20, 10))
d$y <- 10 + 2 * d$x - d$z + rnorm(n)
d$t <- ave(seq_len(n), d$g, FUN = seq_along)
sf <- frmtmb::set_prior
d0 <- quiet(frm_sample(bf(y ~ x + z), family = gaussian(), data = d, chains = 1, iter = 200, refresh = 0,
                       seed = 5, prior = sf("normal(10, 0.1)", class = "b", coef = "x")))
cat("prior_summary, first call:\n")
print(as.data.frame(frmtmb::prior_summary(d0))[, c("prior", "class", "coef", "source")])
d1 <- quiet(update(d0, prior = sf("normal(5, 0.1)", class = "sigma")))
cat("prior_summary, after update():\n")
print(as.data.frame(frmtmb::prior_summary(d1))[, c("prior", "class", "coef", "source")])
tb <- brms::default_prior(brms::bf(y ~ x + z), data = d)
tb$prior[tb$class == "b" & tb$coef == "x"] <- "normal(10, 0.1)"
cat("brms table, source of the edited row:", tb$source[tb$class == "b" & tb$coef == "x"], "\n")
bq <- quiet(brms::brm(brms::bf(y ~ x + z), data = d, chains = 1, iter = 200, refresh = 0, seed = 5, silent = 2, prior = tb))
print(as.data.frame(bq$prior)[, c("prior", "class", "coef", "source")])
fit <- quiet(frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), family = gaussian(), data = d))
nd <- data.frame(x = c(-1, 0, 1), t = 1:3, g = factor(1, levels = levels(d$g)))
r <- tryCatch(predict(fit, newdata = nd, ndraws = 20), error = function(e) substr(conditionMessage(e), 1, 170))
cat("predict(newdata), cov = TRUE, this build:", if (is.character(r)) r else "answers", "\n")
