# Lane parity2: the data every probe of the lane reads, so the brms and
# frmtmb arms see the same rows. Sourced by dev/parity2-brms-probe.R and
# dev/parity2-frm-probe.R. Data seed 1 (the check2.R data of the
# triage, triage-082/scratch-B/check2.R), plus columns for the other
# items.
set.seed(1)
n <- 150
d <- data.frame(x = runif(n, 0.1, 2), z = rnorm(n), k = sample(1:4, n, TRUE),
                w = rnorm(n), g = factor(sample(letters[1:8], n, TRUE)))
d$y <- 1 + d$x * d$z + 0.3 * d$k + rnorm(n)
d$yo <- cut(d$y, quantile(d$y, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
            labels = FALSE)
d$y2 <- 0.5 * d$x + rnorm(n)
d$x2 <- d$x + rnorm(n, sd = 0.3)
d$x2[seq(3, n, 10)] <- NA

# A random intercept model with a few levels for the new-level band
# (item 1497) and the loo subset (item 5331). Seed 11.
set.seed(11)
dg <- data.frame(x = rnorm(48), g = factor(rep(1:6, each = 8)))
dg$y <- 1 + 0.5 * dg$x + rnorm(6, sd = 0.8)[dg$g] + rnorm(48, sd = 0.5)

# drop_unused_levels = FALSE (item 5337): g of 7 levels, 6 used; f of 3,
# 2 used. Seed 3.
set.seed(3)
du <- data.frame(g = factor(sample(1:6, 60, TRUE), levels = 1:7),
                 f = factor(sample(c("a", "b"), 60, TRUE),
                            levels = c("a", "b", "c")),
                 x = rnorm(60))
du$y <- rnorm(60, 1 + 0.4 * du$x)

# Multinomial and ordinal responses (items 5323, 5325). Seed 5.
set.seed(5)
dm <- data.frame(x = rnorm(60))
pm <- cbind(1, exp(0.5 * dm$x), exp(-0.3 * dm$x))
dm$Y <- t(apply(pm / rowSums(pm), 1, function(p) rmultinom(1, 10, p)))
colnames(dm$Y) <- c("a", "b", "c")
dm$size <- 10L
dm$yo <- cut(dm$x + rnorm(60), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
