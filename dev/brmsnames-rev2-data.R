## Reviewer recheck: hostile-name data and formulas shared by the brms and
## frmtmb arms. Data seed 52. Sourced, not run.
suppressMessages(library(splines))
rev2_data <- function() {
  set.seed(52)
  n <- 240
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$f <- factor(sample(c("a b", "c-d", "e:f"), n, TRUE),
                levels = c("a b", "c-d", "e:f"))
  d$f2 <- factor(sample(c("u", "v"), n, TRUE))
  d$g <- factor(paste("lvl", rep(1:12, length.out = n)))
  d$h <- factor(sample(c("p", "q"), n, TRUE))
  d$g1 <- factor(sample(1:10, n, TRUE))
  d$g2 <- factor(sample(1:10, n, TRUE))
  d$xo <- factor(sample(1:4, n, TRUE), ordered = TRUE)
  ug <- rnorm(12, 0, 0.8)
  d$y <- 1 + 0.4 * d$x - 0.3 * d$z + 0.5 * (d$f == "c-d") +
    0.6 * d$x * (d$f == "e:f") + ug[as.integer(d$g)] + rnorm(n, 0, 0.8)
  d$yo <- cut(d$y, quantile(d$y, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
              labels = c("lo", "mid", "hi"), ordered_result = TRUE)
  d$y.a_b <- d$y
  d$y2 <- 0.5 * d$x + ug[as.integer(d$g)] + rnorm(n, 0, 2)
  d$ya <- d$y2 + rnorm(n)
  d$y_a <- d$y + rnorm(n)
  d$yn <- 2 * exp(0.3 * d$x) + ug[as.integer(d$g)] / 4 + rnorm(n, 0, 0.3)
  # an interaction group whose parts carry brms's join character
  d$gi <- factor(sample(c("1_2", "1"), n, TRUE))
  d$hi <- factor(ifelse(d$gi == "1", "2_3", "3"))
  d$hi[sample(n, 40)] <- "2_3"
  d$fab <- factor(sample(c("a b", "ab"), n, TRUE))
  d$IxE2 <- rnorm(n)
  d$sigma_Intercept <- rnorm(n)
  # gr(by =) needs a by-variable constant within each group; appended
  # last so every column above keeps its draws
  d$hg <- factor(ifelse(as.integer(d$g) %% 2 == 0, "p", "q"))
  d
}

rev2_models <- function() {
  list(
    B1 = list(f = y ~ x * z * f + I(x^2) + poly(z, 2, raw = TRUE) +
                (1 + x | g), fam = "gaussian"),
    B2 = list(f = y ~ ns(x, df = 3) + bs(z, df = 4) + mo(xo) +
                (1 | gr(g, by = hg)), fam = "gaussian"),
    B3 = list(f = yo ~ cs(x) + (1 | g), fam = "acat"),
    B4 = list(f = y ~ x + (1 | mm(g1, g2)), fam = "gaussian"),
    B5 = list(f = "nl", fam = "gaussian"),
    B6 = list(f = "mv", fam = "gaussian"),
    B7 = list(f = y ~ s(x, by = f2) + t2(x, z) + gp(z, k = 5, c = 5 / 4),
              fam = "gaussian"),
    B8 = list(f = y ~ x + (1 | gi:hi), fam = "gaussian"),
    B10 = list(f = "sigmacol", fam = "gaussian")
  )
}

## the bf() objects need the package whose bf() builds them
rev2_formula <- function(key, M, bf, mvbf) {
  if (identical(M$f, "nl")) {
    return(bf(yn ~ a1 * exp(b2 * x), a1 ~ 1 + (1 | g), b2 ~ 1, nl = TRUE))
  }
  if (identical(M$f, "mv")) {
    return(mvbf(bf(y.a_b ~ x + (1 | p | g), sigma ~ z),
                bf(y2 ~ x + (1 | p | g)), rescor = FALSE))
  }
  if (identical(M$f, "sigmacol")) {
    return(bf(y ~ x + (1 + sigma_Intercept | g), sigma ~ (1 | g)))
  }
  bf(M$f)
}
