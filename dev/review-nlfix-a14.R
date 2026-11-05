# Prior tables of an ordinal body and an autocor body against brms's default_prior().
source("dev/review-nlfix-common.R")
set.seed(4401)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)), t = rep(1:20, 10))
d$y <- 1 + 0.5 * d$x + ave(rnorm(n, 0, 0.5), d$g, FUN = function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
d$gg <- factor(rep(c("p", "q"), length.out = n))
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
cmp <- function(lab, ff, fb, fam, famb) {
  a <- tt(as.data.frame(get_prior(ff, data = d, family = fam)))
  if (iserr(a)) { cat(lab, "ERR", a$msg, "\n"); return(invisible()) }
  b <- as.data.frame(brms::default_prior(fb, data = d, family = famb))
  key <- function(z) sort(paste(z$class, z$coef, z$group, z$dpar, z$nlpar, sep = "/"))
  cat(sprintf("%-12s rows frm %d brms %d | only frm: %s | only brms: %s\n", lab, nrow(a), nrow(b), paste(setdiff(key(a), key(b)), collapse = " "), paste(setdiff(key(b), key(a)), collapse = " ")))
}
cmp("plain", bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), brms::bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), cumulative(), brms::cumulative())
cmp("gr", bf(yo | thres(gr = gg) ~ b1 * x, b1 ~ 1, nl = TRUE), brms::bf(yo | thres(gr = gg) ~ b1 * x, b1 ~ 1, nl = TRUE), cumulative(), brms::cumulative())
cmp("equidistant", bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), brms::bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), cumulative(threshold = "equidistant"), brms::cumulative(threshold = "equidistant"))
cmp("disc", bf(yo ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE), brms::bf(yo ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE), cumulative(), brms::cumulative())
cmp("mixture", bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), brms::bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), mixture(cumulative(), cumulative()), brms::mixture(brms::cumulative(), brms::cumulative()))
cmp("ar", bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)), brms::bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)), gaussian(), brms::brmsfamily("gaussian"))
f0 <- frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d)
f1 <- tt(frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d, prior = set_prior("normal(5, 0.1)", class = "Intercept", coef = "1")))
cat("prior on threshold 1: flat", round(summary(f0)$fixed[1, 1], 3), " with normal(5, 0.1):", if (iserr(f1)) f1$msg else round(summary(f1)$fixed[1, 1], 3), "\n")
