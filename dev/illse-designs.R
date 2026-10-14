# Lane illse: the design families every measurement of the lane uses.
# Each entry builds list(formula, data, family, lab) from a seed; X is
# the fixed design of mu when the design has a closed-form reference.
pw <- function(d) paste0(" + I(x^", 2:d, ")", collapse = "")
des_poly <- function(deg, re = FALSE, fam = "gaussian", lo = 1, hi = 2,
                     n = 200) {
  force(deg); force(re); force(fam); force(lo); force(hi); force(n)
  function(seed) {
    set.seed(seed)
    d <- data.frame(x = stats::runif(n, lo, hi),
                    g = factor(rep(seq_len(20), length.out = n)))
    eta <- sin(3 * d$x)
    if (re) eta <- eta + stats::rnorm(20, 0, 0.5)[d$g]
    d$y <- switch(fam,
      gaussian = eta + stats::rnorm(n, 0, 0.3),
      bernoulli = stats::rbinom(n, 1, stats::plogis(2 * eta)),
      poisson = stats::rpois(n, exp(1 + eta)),
      cumulative = cut(2 * eta + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
                       labels = FALSE))
    if (fam == "cumulative") d$y <- factor(d$y, ordered = TRUE)
    rhs <- paste0("x", if (deg > 1) pw(deg) else "", if (re) " + (1 | g)")
    lhs <- "y"
    list(formula = stats::as.formula(paste(lhs, "~", rhs)), data = d,
         family = switch(fam, gaussian = gaussian(),
                         bernoulli = bernoulli(), poisson = poisson(),
                         cumulative = cumulative()),
         lab = sprintf("%s deg %d%s", fam, deg, if (re) " (1|g)" else ""),
         poly_formula = stats::as.formula(paste0(
           "y ~ poly(x, ", deg, ")", if (re) " + (1 | g)")))
  }
}
des_offset <- function(off, re = FALSE) {
  force(off); force(re)
  function(seed) {
    set.seed(seed)
    d <- data.frame(x = off + stats::runif(200),
                    g = factor(rep(seq_len(20), length.out = 200)))
    d$y <- 1 + 2 * (d$x - off) + stats::rnorm(200, 0, 0.3) +
      if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
    list(formula = stats::as.formula(paste0("y ~ x", if (re) " + (1 | g)")),
         data = d, family = gaussian(),
         lab = sprintf("gaussian x + %g%s", off, if (re) " (1|g)" else ""))
  }
}
des_collin <- function(dl, re = FALSE) {
  force(dl); force(re)
  function(seed) {
    set.seed(seed)
    d <- data.frame(x1 = stats::rnorm(200),
                    g = factor(rep(seq_len(20), length.out = 200)))
    d$x2 <- d$x1 + dl * stats::rnorm(200)
    d$y <- 1 + d$x1 + d$x2 + stats::rnorm(200, 0, 0.3) +
      if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
    list(formula = stats::as.formula(paste0("y ~ x1 + x2",
                                            if (re) " + (1 | g)")),
         data = d, family = gaussian(),
         lab = sprintf("gaussian x2 = x1 + %g e%s", dl,
                       if (re) " (1|g)" else ""))
  }
}
des_scales <- function(re = FALSE) {
  force(re)
  function(seed) {
    set.seed(seed)
    d <- data.frame(a = stats::rnorm(200, 0, 1e-6),
                    b = stats::rnorm(200, 0, 1e6),
                    c = 3e5 + stats::rnorm(200, 0, 1),
                    g = factor(rep(seq_len(20), length.out = 200)))
    d$y <- 1 + 1e6 * d$a + 1e-6 * d$b + (d$c - 3e5) +
      stats::rnorm(200, 0, 0.3) +
      if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
    list(formula = stats::as.formula(paste0("y ~ a + b + c",
                                            if (re) " + (1 | g)")),
         data = d, family = gaussian(),
         lab = sprintf("gaussian mixed scales%s", if (re) " (1|g)" else ""))
  }
}
des_b5 <- function(re = FALSE) {
  force(re)
  function(seed) {
    set.seed(seed)
    d <- data.frame(x = stats::runif(200, 1, 2), z = stats::rnorm(200),
                    w = factor(sample(c("a", "b"), 200, TRUE)),
                    g = factor(rep(seq_len(20), length.out = 200)))
    d$y <- sin(3 * d$x) + 0.3 * d$z + 0.2 * (d$w == "b") +
      stats::rnorm(200, 0, 0.3) +
      if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
    list(formula = stats::as.formula(paste0("y ~ x", pw(8), " + z + w",
                                            if (re) " + (1 | g)")),
         data = d, family = gaussian(),
         lab = sprintf("gaussian B5 deg 8 + z + w%s",
                       if (re) " (1|g)" else ""),
         poly_formula = stats::as.formula(paste0(
           "y ~ poly(x, 8) + z + w", if (re) " + (1 | g)")))
  }
}
illse_designs <- function() {
  out <- list()
  for (deg in 3:9) for (re in c(FALSE, TRUE)) {
    out[[length(out) + 1]] <- des_poly(deg, re)
  }
  for (off in c(1e5, 1e6, 3e6, 1e7)) for (re in c(FALSE, TRUE)) {
    out[[length(out) + 1]] <- des_offset(off, re)
  }
  for (dl in c(1e-5, 1e-6, 1e-7)) for (re in c(FALSE, TRUE)) {
    out[[length(out) + 1]] <- des_collin(dl, re)
  }
  for (re in c(FALSE, TRUE)) out[[length(out) + 1]] <- des_scales(re)
  for (re in c(FALSE, TRUE)) out[[length(out) + 1]] <- des_b5(re)
  for (fam in c("bernoulli", "poisson", "cumulative")) {
    for (deg in 6:8) for (re in c(FALSE, TRUE)) {
      out[[length(out) + 1]] <- des_poly(deg, re, fam)
    }
  }
  out
}
