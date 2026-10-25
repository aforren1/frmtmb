## Reviewer's cases for the brmsnames lane (2026-10-24), seed 2410.
## Sourced by dev/review-brmsnames-brms.R and dev/review-brmsnames-frm.R.
rv_data <- function() {
  set.seed(2410)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  x1 = stats::rnorm(n),
                  xp = stats::runif(n, 0.5, 2),
                  g = factor(rep(1:12, length.out = n)),
                  t = rep(1:4, length.out = n))
  d$fg <- factor(ifelse(as.integer(d$g) %% 2 == 0, "a", "b"))
  d$logxp <- stats::rnorm(n)
  d$IxM1 <- stats::rnorm(n)
  d$y <- 1 + 0.4 * exp(d$x / 2) + 0.3 * d$z +
    stats::rnorm(12, 0, 0.3)[d$g] + stats::rnorm(n, 0, exp(0.1 * d$z))
  d$y2 <- 0.5 * exp(d$x / 3) + stats::rnorm(n)
  d$yn <- 2 * exp(0.3 * d$x) + stats::rnorm(n, 0, 0.2)
  eta2 <- 0.3 * exp(d$x / 2)
  p <- cbind(1, exp(eta2), exp(0.2))
  p <- p / rowSums(p)
  d$ycat <- apply(p, 1, function(pr) sample(1:3, 1, prob = pr))
  d
}

rv_cases <- function(bf, mvbf) {
  g <- "gaussian"
  list(
    list(id = "log1p", form = bf(y ~ log1p(xp)), family = g),
    list(id = "abs", form = bf(y ~ abs(x)), family = g),
    list(id = "exp:z", form = bf(y ~ exp(x):z), family = g),
    list(id = "I(exp)", form = bf(y ~ I(exp(x))), family = g),
    list(id = "(exp|g)", form = bf(y ~ x + (exp(x) | g)), family = g),
    list(id = "gr by", form = bf(y ~ exp(x) + (1 + exp(x) | gr(g, by = fg))),
         family = g),
    list(id = "mv exp exp", form = mvbf(bf(y ~ exp(x)), bf(y2 ~ exp(x)),
                                        rescor = FALSE), family = g),
    list(id = "nl sigma~exp", form = bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1,
                                        sigma ~ exp(z), nl = TRUE),
         family = g),
    list(id = "cat mu2~exp", form = bf(ycat ~ 1, mu2 ~ exp(x)),
         family = "categorical"),
    list(id = "clash log", form = bf(y ~ log(xp) + logxp), family = g),
    list(id = "clash I", form = bf(y ~ I(x - 1) + IxM1), family = g),
    list(id = "clash poly", form = bf(y ~ poly(x1, 2) + poly(x, 12)),
         family = g),
    list(id = "mv nl", form = mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + exp(z),
                                      b ~ 1, nl = TRUE),
                                   bf(y ~ x + exp(z)), rescor = FALSE),
         family = g),
    list(id = "nl re 0+", form = bf(y ~ a + b * x,
                                    a ~ 1 + exp(z) + (1 | g), b ~ 0 + z,
                                    nl = TRUE), family = g),
    list(id = "mv nl nl sigma re",
         form = mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + z + (1 | g), b ~ 1,
                        sigma ~ z, nl = TRUE),
                     bf(y ~ c + d * x, c ~ exp(z), d ~ 1 + (1 | g),
                        sigma ~ x, nl = TRUE), rescor = FALSE),
         family = g)
  )
}
