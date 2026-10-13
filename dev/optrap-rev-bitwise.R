# Reviewer of lane optrap, item 1 ("fits as before" at unit 1): a sample
# of fits with no unit, the objective and every outer parameter printed
# in hex, to compare two builds bit for bit. Pin both runs to one core
# type (dev/rtmb-pitfalls.md item 21).
#   Rscript dev/optrap-rev-bitwise.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
mk <- function(s) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100),
                  z = runif(100))
  d$b <- rnorm(20, 0, 0.5)[d$g]
  d$eta <- 0.3 + 0.5 * d$x + d$b
  d
}
lev <- c("a", "b", "c", "d")
models <- list(
  gauss_u1 = function(d) {
    d$y <- d$eta + rnorm(100)
    frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  },
  gauss_slope_u1 = function(d) {
    d$y <- d$eta + d$x * rnorm(20, 0, 0.3)[d$g] + rnorm(100)
    frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)
  },
  student_u1 = function(d) {
    d$y <- d$eta + rt(100, 4)
    frm(bf(y ~ x + (1 | g)), family = student(), data = d)
  },
  poisson = function(d) {
    d$y <- rpois(100, exp(d$eta))
    frm(bf(y ~ x + (1 | g)), family = poisson(), data = d)
  },
  bernoulli = function(d) {
    d$y <- rbinom(100, 1, plogis(d$eta))
    frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d)
  },
  gamma = function(d) {
    d$y <- rgamma(100, 2, 2 / exp(d$eta))
    frm(bf(y ~ x + (1 | g)), family = Gamma(link = "log"), data = d)
  },
  negbin = function(d) {
    d$y <- rnbinom(100, mu = exp(d$eta), size = 2)
    frm(bf(y ~ x + (1 | g)), family = negbinomial(), data = d)
  },
  beta = function(d) {
    m <- plogis(d$eta)
    d$y <- rbeta(100, m * 5, (1 - m) * 5)
    frm(bf(y ~ x + (1 | g)), family = Beta(), data = d)
  },
  cumulative = function(d) {
    d$y <- cut(d$eta + rlogis(100), c(-Inf, -0.5, 0.5, 1.5, Inf),
               labels = lev, ordered_result = TRUE)
    frm(bf(y ~ x + (1 | g)), family = cumulative(), data = d)
  },
  mixture = function(d) {
    d$y <- 0.5 * d$x + ifelse(rbinom(100, 1, 0.4) == 1, rnorm(100, 3),
                              rnorm(100))
    frm(bf(y ~ x), family = mixture(gaussian, gaussian), data = d)
  },
  mo_u1 = function(d) {
    d$o <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
    d$y <- c(0, 0.6, 0.8, 0.9)[d$o] + rnorm(100) * 0.9
    frm(bf(y ~ mo(o) + x), family = gaussian(), data = d)
  },
  zip = function(d) {
    d$y <- rpois(100, exp(d$eta)) * rbinom(100, 1, 0.8)
    frm(bf(y ~ x + (1 | g)), family = zero_inflated_poisson(), data = d)
  },
  sigma_x_u1 = function(d) {
    d$y <- d$eta + exp(0.2 * d$x) * rnorm(100)
    frm(bf(y ~ x + (1 | g), sigma ~ x), family = gaussian(), data = d)
  },
  reml_u1 = function(d) {
    d$y <- d$eta + rnorm(100)
    frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d, REML = TRUE)
  }
)
for (nm in names(models)) for (s in 1:3) {
  d <- mk(s)
  f <- tryCatch(suppressMessages(suppressWarnings(models[[nm]](d))),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%s %d ERROR %s\n", nm, s, conditionMessage(f)))
    next
  }
  un <- tryCatch(length(get("predictor_units",
                            asNamespace("frmtmb"))(f$frame)),
                 error = function(e) NA)
  cat(sprintf("%s %d units %s code %d obj %a par %s\n", nm, s, un,
              f$opt$convergence, f$opt$objective,
              paste(sprintf("%a", f$opt$par), collapse = ",")))
}
