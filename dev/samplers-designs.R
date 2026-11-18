# Designs shared by the samplers lane's scripts. Sourced, not run.
#
# samplers_lib(): the lane's library paths and environment.
# samplers_eq_env(): the equivalence survey's instrument, loaded as
#   dev/equiv-measure.R loads it, so a survey pair's spec is read from
#   the survey itself rather than restated here.
# samplers_design(id): a list(bf_frm, bf_brms, family_frm, family_brms,
#   data, data2, prior_frm, prior_brms) for a survey pair "pNNN" or for
#   one of the lane's large designs.

samplers_lib <- function() {
  lib <- "C:/Users/adf44/source/r/wt-samplers-lib"
  Sys.setenv(R_USER_CACHE_DIR = file.path(lib, "_cache"),
             R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  .libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  invisible(lib)
}

samplers_eq_env <- function() {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(testthat)
  })
  attach(list(me = brms::me), name = "equiv_brms_me",
         warn.conflicts = FALSE)
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  e
}

samplers_big <- function(id) {
  switch(id,
    # 2000 groups of 5 rows: the large poisson GLMM
    big_pois = {
      set.seed(20261117)
      ng <- 2000L
      d <- data.frame(g = factor(rep(seq_len(ng), each = 5L)),
                      x = rnorm(5L * ng))
      u <- rnorm(ng, 0, 0.5)
      d$y <- rpois(nrow(d), exp(0.5 + 0.3 * d$x + u[as.integer(d$g)]))
      list(bf = "bf(y ~ x + (1 | g))", family = "poisson()", data = d)
    },
    # 300 x 40 crossed bernoulli, 6000 rows
    big_crossed = {
      set.seed(20261118)
      d <- data.frame(a = factor(rep(seq_len(300L), each = 20L)))
      d$b <- factor(sample(seq_len(40L), nrow(d), TRUE), levels = 1:40)
      d$x <- rnorm(nrow(d))
      ua <- rnorm(300, 0, 0.7)
      ub <- rnorm(40, 0, 0.5)
      eta <- -0.3 + 0.5 * d$x + ua[as.integer(d$a)] + ub[as.integer(d$b)]
      d$y <- rbinom(nrow(d), 1, plogis(eta))
      list(bf = "bf(y ~ x + (1 | a) + (1 | b))", family = "bernoulli()",
           data = d)
    },
    # 500 groups of 8, correlated random slope, bernoulli
    big_slope = {
      set.seed(20261119)
      ng <- 500L
      d <- data.frame(g = factor(rep(seq_len(ng), each = 8L)),
                      x = rnorm(8L * ng))
      z1 <- rnorm(ng)
      z2 <- 0.3 * z1 + sqrt(1 - 0.09) * rnorm(ng)
      u0 <- 0.8 * z1
      u1 <- 0.4 * z2
      i <- as.integer(d$g)
      d$y <- rbinom(nrow(d), 1, plogis(-0.2 + 0.6 * d$x + u0[i] +
                                         u1[i] * d$x))
      list(bf = "bf(y ~ x + (1 + x | g))", family = "bernoulli()", data = d)
    },
    stop("no design ", id))
}

# The sampler comparison's own designs (question 3).
samplers_q3 <- function(id) {
  switch(id,
    small_pois = {
      set.seed(301)
      d <- data.frame(g = factor(rep(1:30, each = 10L)), x = rnorm(300))
      u <- rnorm(30, 0, 0.6)
      d$y <- rpois(300, exp(0.4 + 0.3 * d$x + u[as.integer(d$g)]))
      list(bf = "bf(y ~ x + (1 | g))", family = "poisson()", data = d)
    },
    big_pois = samplers_big("big_pois"),
    big_crossed = samplers_big("big_crossed"),
    slope = {
      set.seed(304)
      ng <- 40L
      d <- data.frame(g = factor(rep(seq_len(ng), each = 10L)),
                      x = rnorm(10L * ng))
      z1 <- rnorm(ng)
      z2 <- 0.5 * z1 + sqrt(0.75) * rnorm(ng)
      i <- as.integer(d$g)
      d$y <- 1 + 0.5 * d$x + 0.8 * z1[i] + 0.4 * z2[i] * d$x +
        rnorm(nrow(d), 0, 0.7)
      list(bf = "bf(y ~ x + (1 + x | g))", family = "gaussian()", data = d)
    },
    gp = {
      set.seed(305)
      d <- data.frame(x = sort(runif(60, 0, 10)))
      d$y <- sin(d$x) + rnorm(60, 0, 0.3)
      list(bf = "bf(y ~ gp(x))", family = "gaussian()", data = d)
    },
    ald = {
      set.seed(306)
      d <- data.frame(g = factor(rep(1:20, each = 10L)), x = rnorm(200))
      u <- rnorm(20, 0, 0.5)
      d$y <- 1 + 0.5 * d$x + u[as.integer(d$g)] + rt(200, 4) * 0.5
      list(bf = "bf(y ~ x + (1 | g))", family = "asym_laplace()", data = d)
    },
    weak = {
      # two rows per group: the Laplace approximation is weak
      set.seed(307)
      d <- data.frame(g = factor(rep(1:150, each = 2L)), x = rnorm(300))
      u <- rnorm(150, 0, 1)
      d$y <- rbinom(300, 1, plogis(0.2 + 0.5 * d$x + u[as.integer(d$g)]))
      list(bf = "bf(y ~ x + (1 | g))", family = "bernoulli()", data = d)
    },
    stop("no design ", id))
}

samplers_design <- function(id, eq_env = NULL) {
  if (grepl("^p[0-9]{3}$", id)) {
    smp <- read.delim("dev/equiv-data/sample.tsv")
    row <- smp[smp$id == id, ]
    key <- paste(row$feature_a, row$feature_b, sep = " x ")
    spec <- eq_env$eq_specs[[key]]()
    spec$key <- key
    if (length(spec$data2)) spec$brms_args$data2 <- spec$data2
    return(spec)
  }
  s <- if (startsWith(id, "big_")) samplers_big(id) else samplers_q3(id)
  s$key <- id
  s
}
