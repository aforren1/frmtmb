# Review of lane nlstop2, attack 3 continued: huber() models that fit
# under Laplace on 0.82.0 and that the route sends to quadrature.
# REV_ARM=before|after Rscript dev/reviews/nlstop2-review-huber2.R
source("dev/reviews/nlstop2-review-common.R")
hd <- function(seed, n = 100, ng = 20) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), each = n / ng)),
                  x = rnorm(n))
  b <- rnorm(ng)[d$g]
  d$y <- 1 + 0.5 * d$x + 0.5 * b + rt(n, 3)
  d$h <- factor(rep(1:10, length.out = n))  # crossed with g
  d$z <- rnorm(n)
  d
}
one <- function(tag, expr) {
  t0 <- proc.time()[[3]]
  r <- tt(expr)
  el <- proc.time()[[3]] - t0
  if (inherits(r$v, "tterr")) {
    cat("ROW", tag, "ERROR", substr(r$v$msg, 1, 160), "\n")
  } else {
    cat("ROW", tag, "ok quad", isTRUE(r$v$quadrature), "ll",
        format(llv(r$v), digits = 10), "code", r$v$opt$convergence, "warn",
        if (length(r$w)) paste(unique(r$w), collapse = "+") else "-",
        "t", round(el, 2), "\n")
  }
  invisible(r$v)
}
for (s in 1:10) {
  d <- hd(s)
  one(paste("sigma~(1|g) seed", s),
      frm(bf(y ~ x, sigma ~ (1 | g)), family = huber(), data = d))
}
for (s in 1:10) {
  d <- hd(s)
  one(paste("sigma~(1|g) scale1000 seed", s),
      frm(bf(y ~ x, sigma ~ (1 | g)), family = huber(),
          data = transform(d, y = y * 1000)))
}
for (s in 1:5) {
  d <- hd(s)
  one(paste("crossed seed", s),
      frm(y ~ x + (1 | g) + (1 | h), family = huber(), data = d))
}
for (s in 1:5) {
  d <- hd(s)
  d$a <- factor((as.integer(d$g) - 1) %/% 4)
  one(paste("nested a/g seed", s),
      frm(y ~ x + (1 | a / g), family = huber(), data = d))
}
for (s in 1:5) {
  d <- hd(s)
  one(paste("slope (0+x|g) seed", s),
      frm(y ~ x + (0 + x | g), family = huber(), data = d))
}
for (s in 1:5) {
  d <- hd(s)
  d$pos <- exp(d$y / 4)
  one(paste("log link seed", s),
      frm(pos ~ x + (1 | g), family = huber(link = "log"), data = d))
}
# methods that read the joint covariance of b
d <- hd(4)
fit <- frm(y ~ x + (1 | g), family = huber(), data = d)
cat("quad", isTRUE(fit$quadrature), "\n")
show("vcov_cluster", tt(vcov_cluster(fit, cluster = d$g)))
show("frm_linpred se", tt(frm_linpred(fit, se.fit = TRUE)))
show("predict", tt(predict(fit)))
show("fitted", tt(fitted(fit)))
show("residuals", tt(residuals(fit)))
show("conditional_effects", tt(conditional_effects(fit)))
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  show("avg_slopes", tt(marginaleffects::avg_slopes(fit)))
}
if (requireNamespace("emmeans", quietly = TRUE)) {
  show("emmeans", tt(emmeans::emmeans(fit, ~ x)))
}
cat("DONE\n")
