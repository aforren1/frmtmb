# The three filed defects, with what the optimizer did. Seeds as named.
# Usage: NLSTART_ARM=before|after Rscript dev/nlstart-repro.R
source("dev/nlstart-common.R")
diag_fit <- function(label, r) {
  cat(label, ":", ll(r), "\n")
  if (iserr(r)) return(invisible())
  o <- r$opt
  cat("   iters", o$iterations %||% NA, " evals",
      paste(o$evaluations, collapse = "/"), " msg:", o$message, "\n")
  g <- tryCatch(r$obj$gr(o$par), error = function(e) NA)
  cat("   grad", format(g, digits = 3), "\n")
  if (length(attr(r, "w"))) cat("   W:", substr(attr(r, "w"), 1, 160),
                                sep = "\n      ")
  if (length(attr(r, "m"))) cat("   M:", substr(attr(r, "m"), 1, 160),
                                sep = "\n      ")
  cat("\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## 1: b1 * (1 - exp(-b2 * x)) under gaussian, the 450-row scan's data
gen <- function(seed, slope) {
  set.seed(seed); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z
  d$y <- lat + rnorm(n, 0, 0.5)
  d
}
fm <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 ~ 1, b2 ~ 1, nl = TRUE)
for (slope in c(0.6, -0.6)) for (s in 1:5) {
  d <- gen(s * 100 + 1, slope)
  r <- tt(frm(fm, family = gaussian(), data = d))
  diag_fit(sprintf("sat seed %d slope %.1f default", s * 100 + 1, slope), r)
  r2 <- tt(frm(fm, family = gaussian(), data = d,
               start = list(beta = c("b1_(Intercept)" = -2,
                                     "b2_(Intercept)" = -0.5))))
  cat("   start (-2,-0.5):", ll(r2), "\n")
  r3 <- tt(frm(fm, family = gaussian(), data = d,
               start = list(beta = c("b1_(Intercept)" = 2,
                                     "b2_(Intercept)" = 0.5))))
  cat("   start (2,0.5):", ll(r3), "\n")
}

## 2: a * x with x about 1e-20 (dev/review-nlfix-c9.R, seed 6201)
set.seed(6201); n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 20)),
                t = rep(1:20, 10), z = rnorm(n))
d$y <- 0.5 * d$x + rnorm(n)
d$xt <- d$x * 1e-20
diag_fit("tiny x linear", tt(frm(y ~ 0 + xt, data = d)))
diag_fit("tiny x body", tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d)))

## 3: a + b * x with y about 1e6 (dev/review-nlfix-a10.R, seed 9901)
set.seed(9901)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)))
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4,
            ordered_result = TRUE)
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$xs <- d$x * 1e-10
d$ybig <- 1e6 + 5e4 * d$x + rnorm(n, 0, 1e3)
diag_fit("big y linear", tt(frm(ybig ~ x, data = d)))
diag_fit("big y body", tt(frm(bf(ybig ~ a + b * x, a + b ~ 1, nl = TRUE),
                              data = d)))
