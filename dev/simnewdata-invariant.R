# The identity simulate(newdata = <the fitted rows>) must meet: at the
# same seed it draws what simulate() draws on the fitted rows, because
# the two differ only in where the distributional parameters come from
# (the stored design against a design rebuilt on newdata). Run across a
# family grid under re_formula = NULL and NA.
#   Rscript dev/simnewdata-invariant.R > dev/simnewdata-log/invariant.txt
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
set.seed(101)
n <- 240
d <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                g = factor(rep(seq_len(12), each = 20)),
                h = factor(rep(seq_len(8), 30)))
u <- stats::rnorm(12, 0, 0.5)[d$g]
eta <- 0.3 + 0.5 * d$x + u
d$yg <- stats::rnorm(n, eta, 0.7)
d$yp <- stats::rpois(n, exp(eta))
d$nt <- sample(5:15, n, replace = TRUE)
d$yb <- stats::rbinom(n, d$nt, stats::plogis(eta))
d$yz <- ifelse(stats::runif(n) < 0.3, 0L, stats::rpois(n, exp(eta)))
d$yo <- factor(cut(eta + stats::rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                   labels = FALSE), ordered = TRUE)
d$yk <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
d$yt <- stats::rnorm(n, eta, 0.7)
d <- d[d$yt > -0.5, ]
d$ys <- sin(3 * d$z) + stats::rnorm(nrow(d), 0, 0.3)
d$ym <- ifelse(stats::runif(nrow(d)) < 0.4, stats::rnorm(nrow(d), -2),
               stats::rnorm(nrow(d), 2))
d$ystu <- eta[seq_len(nrow(d))] + stats::rt(nrow(d), 4)
d$time <- rep(1:20, length.out = nrow(d))

designs <- list(
  gaussian = list(bf(yg ~ x + (1 + x | g)), gaussian()),
  poisson = list(bf(yp ~ x + (1 | g) + (1 | h)), poisson()),
  binomial = list(bf(yb | trials(nt) ~ x + (1 | g)), binomial()),
  zip = list(bf(yz ~ x + (1 | g)), zero_inflated_poisson()),
  cumulative = list(bf(yo ~ x + (1 | g)), cumulative()),
  cumulative_cs = list(bf(yo ~ cs(x)), sratio()),
  categorical = list(bf(yk ~ x), categorical()),
  trunc = list(bf(yt | trunc(lb = -0.5) ~ x + (1 | g)), gaussian()),
  smooth = list(bf(ys ~ s(z)), gaussian()),
  distributional = list(bf(yg ~ x + (1 | g), sigma ~ x), gaussian()),
  mixture = list(bf(ym ~ 1), mixture(gaussian(), gaussian())),
  student_ar = list(bf(ystu ~ x + ar(time, gr = g, cov = TRUE)),
                    student()),
  offset = list(bf(yp ~ x + offset(z) + (1 | g)), poisson())
)
cmp <- function(a, b) {
  a <- as.matrix(as.data.frame(lapply(a, as.numeric)))
  b <- as.matrix(as.data.frame(lapply(b, as.numeric)))
  if (!identical(dim(a), dim(b))) return("DIM")
  if (identical(unname(a), unname(b))) return("identical")
  r <- max(abs(a - b) / pmax(abs(b), 1))
  sprintf("max rel diff %.3g, %d of %d cells differ", r, sum(a != b),
          length(a))
}
for (nm in names(designs)) {
  ds <- designs[[nm]]
  fit <- tryCatch(suppressWarnings(frm(ds[[1]], family = ds[[2]], data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-15s FIT ERROR %s\n", nm, conditionMessage(fit)))
    next
  }
  for (rf in list(NULL, NA)) {
    lab <- if (is.null(rf)) "NULL" else "NA"
    a <- tryCatch(simulate(fit, nsim = 20, seed = 5, re_formula = rf),
                  error = function(e) e)
    b <- tryCatch(simulate(fit, nsim = 20, seed = 5, re_formula = rf,
                           newdata = d), error = function(e) e)
    res <- if (inherits(a, "error")) {
      paste("in-sample ERROR:", conditionMessage(a))
    } else if (inherits(b, "error")) {
      paste("newdata ERROR:", conditionMessage(b))
    } else {
      cmp(a, b)
    }
    cat(sprintf("%-15s re_formula = %-4s %s\n", nm, lab, res))
  }
}
