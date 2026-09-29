# The false-alarm rate of the gradient check over a design set the field
# produces. One process per design group, so three run at once.
#
#   Rscript dev/gradcheck-03-design.R base g1
#
# groups: g1 = gaussian poisson binomial
#         g2 = cumulative sratio nearzero
#         g3 = rescor bounded nonlinear
#
# Writes one RDS per (library, group). Every row carries its seed.

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "base"
group <- if (length(args) > 1L) args[[2L]] else "g1"
nrep <- if (length(args) > 2L) as.integer(args[[3L]]) else 20L
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    " group:", group, " reps:", nrep, "\n")

# Every design is a function of (n, seed) returning a fit or an error.
# The truth is generated at known coefficients, so a fit that lands on
# them is correct whatever the gradient reads.
designs <- list(
  gaussian = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    d$y <- rnorm(n, 1 + 0.8 * d$x - 0.4 * d$z, 1.5)
    frm(bf(y ~ x + z), family = gaussian(), data = d)
  },
  poisson = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    d$y <- rpois(n, exp(0.4 + 0.5 * d$x - 0.3 * d$z))
    frm(bf(y ~ x + z), family = poisson(), data = d)
  },
  binomial = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    d$y <- rbinom(n, 1, plogis(0.2 + 0.9 * d$x - 0.5 * d$z))
    frm(bf(y ~ x + z), family = bernoulli(), data = d)
  },
  cumulative = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
    e <- 0.8 * d$x1 - 0.5 * d$x2
    d$yo <- cut(e + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
                labels = FALSE)
    frm(bf(yo ~ x1 + x2), family = cumulative(), data = d)
  },
  sratio = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
    e <- 0.8 * d$x1 - 0.5 * d$x2
    d$yo <- cut(e + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
                labels = FALSE)
    frm(bf(yo ~ x1 + x2), family = sratio(), data = d)
  },
  nearzero = function(n, seed) {
    # a mixed model whose variance component is truly zero: the fit sits
    # on the boundary of the parameter space, where a one-sided
    # likelihood makes the gradient check hardest
    set.seed(seed)
    q <- max(8L, n %/% 25L)
    d <- data.frame(x = rnorm(n), g = factor(rep_len(seq_len(q), n)))
    d$y <- rnorm(n, 1 + 0.6 * d$x, 1)
    frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  },
  rescor = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n))
    e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.6, 0.6, 1), 2))
    d$y1 <- 1 + 0.7 * d$x + e[, 1]
    d$y2 <- -0.5 + 0.3 * d$x + e[, 2]
    frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE), family = gaussian(),
        data = d)
  },
  bounded = function(n, seed) {
    # the bound BINDS: the unconstrained slope is 0.8, the cap is 0.1
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    d$y <- rnorm(n, 1 + 0.8 * d$x - 0.4 * d$z, 1.5)
    frm(bf(y ~ x + z), family = gaussian(), data = d,
        prior = set_prior("", class = "b", ub = 0.1))
  },
  nonlinear = function(n, seed) {
    set.seed(seed)
    d <- data.frame(x = runif(n, 0, 5))
    d$y <- rnorm(n, 3 * exp(-0.8 * d$x), 0.2)
    frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE),
        family = gaussian(), data = d,
        start = list(beta = c(a_Intercept = 1, b_Intercept = 1)))
  }
)

groups <- list(g1 = c("gaussian", "poisson", "binomial"),
               g2 = c("cumulative", "sratio", "nearzero"),
               g3 = c("rescor", "bounded", "nonlinear"))
sizes <- c(200L, 1000L, 5000L, 20000L)
out <- list()
t0 <- proc.time()[["elapsed"]]
for (dn in groups[[group]]) {
  for (n in sizes) {
    for (r in seq_len(nrep)) {
      # one seed per (design, n, replicate), independent across cells so
      # no two cells share a draw
      seed <- 7000L + 997L * match(dn, names(designs)) + 31L * r +
        match(n, sizes) * 100003L
      cp <- gc_catch(designs[[dn]](n, seed))
      f <- cp$value
      if (inherits(f, "gc_error")) {
        out[[length(out) + 1L]] <- list(
          design = dn, n = n, rep = r, seed = seed, err = as.character(f),
          warned = NA, np = NA_integer_, gmax = NA_real_,
          gmax_par = NA_character_, gproj = NA_real_, nactive = NA_integer_,
          decr = NA_real_, gscaled = NA_real_, grel = NA_real_,
          obj = NA_real_, conv = NA_integer_)
        next
      }
      pr <- gc_probe(f)
      out[[length(out) + 1L]] <- c(
        list(design = dn, n = n, rep = r, seed = seed, err = NA_character_,
             warned = gc_warned_grad(cp$warnings)), pr)
    }
    cat(sprintf("%-11s n %6d  done %5.1f s\n", dn, n,
                proc.time()[["elapsed"]] - t0))
    utils::flush.console()
  }
}
f <- file.path("dev", paste0("gradcheck-03-", which_lib, "-", group, ".rds"))
saveRDS(out, f)
cat("wrote", f, " rows", length(out), "\n")
