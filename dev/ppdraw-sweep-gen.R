# Data generators and model specs of lane ppdraw's sweep (dev/ppdraw-sweep.R).
# Sourced; expects n in the caller's environment.
gen <- function(fam, case, seed = 101) {
  set.seed(seed)
  x <- stats::rnorm(n)
  inner <- case == "interior"
  zi <- if (inner) 0.3 else 0
  d <- data.frame(x = x)
  zero <- stats::runif(n) < zi
  if (fam == "zero_inflated_poisson") {
    d$y <- ifelse(zero, 0L, stats::rpois(n, exp(0.5 + 0.4 * x)))
  } else if (fam == "zero_inflated_negbinomial") {
    d$y <- ifelse(zero, 0L, stats::rnbinom(n, mu = exp(0.8 + 0.4 * x),
                                           size = 3))
  } else if (fam %in% c("zero_inflated_binomial",
                        "zero_inflated_beta_binomial")) {
    d$N <- sample(5:15, n, TRUE)
    p <- stats::plogis(0.3 + 0.5 * x)
    if (fam == "zero_inflated_beta_binomial") {
      p <- stats::rbeta(n, p * 10, (1 - p) * 10)
    }
    d$y <- ifelse(zero, 0L, stats::rbinom(n, d$N, p))
  } else if (fam == "zero_inflated_beta") {
    m <- stats::plogis(0.2 + 0.5 * x)
    d$y <- ifelse(zero, 0, stats::rbeta(n, m * 10, (1 - m) * 10))
  } else if (fam == "zero_one_inflated_beta") {
    m <- stats::plogis(0.2 + 0.5 * x)
    zoi <- if (inner) 0.2 else 0
    one <- stats::runif(n) < 0.5
    d$y <- ifelse(stats::runif(n) < zoi, as.numeric(one),
                  stats::rbeta(n, m * 10, (1 - m) * 10))
  } else if (fam == "zero_inflated_asym_laplace") {
    d$y <- ifelse(zero, 0, 2 + 0.5 * x + stats::rnorm(n))
  } else if (fam %in% c("hurdle_poisson", "hurdle_negbinomial")) {
    mu <- exp(if (inner) 0.5 + 0.4 * x else 2 + 0.2 * x)
    pos <- vapply(mu, function(m) {
      repeat {
        v <- if (fam == "hurdle_poisson") stats::rpois(1, m) else {
          stats::rnbinom(1, mu = m, size = 3)
        }
        if (v > 0) return(as.integer(v))
      }
    }, 1L)
    d$y <- ifelse(zero, 0L, pos)
  } else if (fam == "hurdle_gamma") {
    d$y <- ifelse(zero, 0, stats::rgamma(n, 3, 3 / exp(0.5 + 0.4 * x)))
  } else if (fam == "hurdle_lognormal") {
    d$y <- ifelse(zero, 0, stats::rlnorm(n, 0.5 + 0.4 * x, 0.5))
  } else if (fam == "hurdle_cumulative") {
    k <- cut(x + stats::rlogis(n), c(-Inf, -0.5, 0.8, Inf), labels = FALSE)
    d$y <- ifelse(zero, 0L, k)
  } else if (fam == "mixture_gaussian") {
    if (inner) {
      c2 <- stats::runif(n) < 0.4
      d$y <- ifelse(c2, 3, -1) + 0.5 * x + stats::rnorm(n)
    } else {
      d$y <- 0.5 * x + stats::rnorm(n)
    }
  } else if (fam == "mixture_poisson") {
    if (inner) {
      c2 <- stats::runif(n) < 0.4
      d$y <- stats::rpois(n, exp(ifelse(c2, 2.3, 0.7) + 0.3 * x))
    } else {
      d$y <- stats::rpois(n, exp(1 + 0.3 * x))
    }
  } else if (fam == "mixture_cumulative") {
    if (inner) {
      c2 <- stats::runif(n) < 0.4
      d$y <- cut(x + ifelse(c2, 3, 0) + stats::rlogis(n),
                 c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
    } else {
      d$y <- cut(x + stats::rlogis(n), c(-Inf, -1, 0.5, 2, Inf),
                 labels = FALSE)
    }
  }
  d
}

spec <- function(fam) {
  switch(fam,
    zero_inflated_binomial = list(bf(y | trials(N) ~ x),
                                  zero_inflated_binomial()),
    zero_inflated_beta_binomial = list(bf(y | trials(N) ~ x),
                                       zero_inflated_beta_binomial()),
    mixture_gaussian = list(bf(y ~ x), mixture(gaussian(), gaussian())),
    mixture_poisson = list(bf(y ~ x), mixture(poisson(), poisson())),
    mixture_cumulative = list(bf(y ~ x),
                              mixture(cumulative(), cumulative())),
    list(bf(y ~ x), get(fam, asNamespace("frmtmb"))()))
}

fams <- c("zero_inflated_poisson", "zero_inflated_negbinomial",
          "zero_inflated_binomial", "zero_inflated_beta_binomial",
          "zero_inflated_beta", "zero_one_inflated_beta",
          "zero_inflated_asym_laplace", "hurdle_poisson",
          "hurdle_negbinomial", "hurdle_gamma", "hurdle_lognormal",
          "hurdle_cumulative", "mixture_gaussian", "mixture_poisson",
          "mixture_cumulative")
