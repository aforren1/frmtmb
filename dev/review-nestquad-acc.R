# Attack 1: the nested rule against an independent reference (adaptive
# Gauss-Hermite per integral, each centered on its own mode by
# optimize(), written from the model alone in review-nestquad-common.R).
#   Rscript dev/review-nestquad-acc.R <case>
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
case <- commandArgs(TRUE)[1]
seed <- 20261106
set.seed(seed)

mk <- function(ns_, ncl, nr, sd_a, sd_i, gen) {
  d <- do.call(rbind, lapply(seq_len(ns_), function(s) {
    do.call(rbind, lapply(seq_len(ncl[(s - 1) %% length(ncl) + 1]),
                          function(c) {
      k <- nr[sample.int(length(nr), 1)]
      data.frame(school = s, cls = c, x = rnorm(k))
    }))
  }))
  d$school <- factor(d$school)
  d$cls_uni <- factor(paste0(d$school, "_", d$cls))
  d$cls_rep <- factor(d$cls)
  ua <- rnorm(nlevels(d$school), 0, sd_a)
  ui <- rnorm(nlevels(d$cls_uni), 0, sd_i)
  d$ua <- ua[d$school]; d$ui <- ui[d$cls_uni]
  gen(d)
}
lp_bern <- function(d) function(eta, r) sum(d$w[r] * dbinom(d$y[r], 1, plogis(eta), log = TRUE))
lp_binom <- function(d) function(eta, r) sum(dbinom(d$y[r], d$n[r], plogis(eta), log = TRUE))
lp_pois <- function(d) function(eta, r) sum(dpois(d$y[r], exp(eta), log = TRUE))
lp_nb <- function(d, shape) function(eta, r) sum(dnbinom(d$y[r], size = shape, mu = exp(eta), log = TRUE))

cfg <- switch(case,
  unbal = {
    # unbalanced: 1 class in some schools, 10 in others, singleton rows;
    # weights and an offset
    d <- mk(16, c(1, 10, 2, 1, 5, 3), c(1, 1, 2, 4, 6), 0.8, 1.0, function(d) {
      d$off <- runif(nrow(d), -0.5, 0.5); d$w <- runif(nrow(d), 0.5, 2)
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.6 * d$x + d$off + d$ua + d$ui)); d })
    list(d = d, fo = bf(y | weights(w) ~ x + offset(off) + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = d$off)
  },
  slope = {
    d <- mk(12, c(3, 4), c(3, 5, 8), 0.7, 0.9, function(d) {
      d$w <- 1
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.4 * d$x + d$x * (d$ua + d$ui))); d })
    list(d = d, fo = bf(y ~ x + (0 + x | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = d$x, zi = d$x, off = 0)
  },
  intslope = {
    # (1 | school) + (0 + x | school): one level of each block per school
    d <- mk(25, 1, c(4, 6, 8), 0.8, 0.8, function(d) {
      d$w <- 1
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.4 * d$x + d$ua + d$x * d$ui)); d })
    d$cls_uni <- d$school
    list(d = d, fo = bf(y ~ x + (1 | school) + (0 + x | school)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = d$x, off = 0)
  },
  dummy = {
    # rows with z = 0 read the outer block only (r0 rows of the plan)
    d <- mk(12, c(2, 3), c(3, 5), 0.8, 1.0, function(d) {
      d$w <- 1; d$z <- rbinom(nrow(d), 1, 0.5)
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + d$ua + d$z * d$ui)); d })
    list(d = d, fo = bf(y ~ x + (1 | school) + (0 + z | school:cls_rep)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = d$z, off = 0)
  },
  dummy2 = {
    # rows with z = 0 read neither block
    d <- mk(12, c(2, 3), c(3, 5), 0.8, 1.0, function(d) {
      d$w <- 1; d$z <- rbinom(nrow(d), 1, 0.6)
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + d$z * (d$ua + d$ui))); d })
    list(d = d, fo = bf(y ~ x + (0 + z | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = d$z, zi = d$z, off = 0)
  },
  flat = {
    # three outer levels, no outer effect: a flat outer marginal
    d <- mk(3, 8, c(3, 4), 0, 0.8, function(d) {
      d$w <- 1
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + d$ui)); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  sdzero = {
    # no inner effect in the truth: the inner sd goes to 0
    d <- mk(15, c(3, 4), c(3, 5), 1, 0, function(d) {
      d$w <- 1
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + d$ua)); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  deg_bern = {
    # every response 1 but two: the log-likelihood is near 0
    d <- mk(8, c(2, 3), c(3, 4), 0.5, 0.5, function(d) {
      d$w <- 1; d$y <- 1L; d$y[c(1, 7)] <- 0L; d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  deg_pois = {
    d <- mk(8, c(2, 3), c(3, 4), 0.5, 0.5, function(d) {
      d$y <- 0L; d$y[3] <- 1L; d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + poisson(),
         lp = lp_pois(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  csep = {
    # complete separation on x
    d <- mk(8, c(2, 3), c(3, 4), 0.5, 0.5, function(d) {
      d$w <- 1; d$y <- as.integer(d$x > 0); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  binom = {
    d <- mk(10, c(2, 4), c(2, 4), 0.6, 0.8, function(d) {
      d$n <- sample(1:12, nrow(d), TRUE)
      d$y <- rbinom(nrow(d), d$n, plogis(-0.3 + 0.5 * d$x + d$ua + d$ui)); d })
    list(d = d, fo = bf(y | trials(n) ~ x + (1 | school) + (1 | school:cls_rep)) + binomial(),
         lp = lp_binom(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  negbin = {
    d <- mk(10, c(2, 3), c(2, 3, 6), 0.6, 0.7, function(d) {
      d$y <- rnbinom(nrow(d), size = 2, mu = exp(0.5 + 0.4 * d$x + d$ua + d$ui)); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + negbinomial(),
         lp = NULL, za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  },
  bigout = {
    d <- mk(10, c(3), c(3, 5), 3, 0.1, function(d) {
      d$w <- 1
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + d$ua + d$ui)); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0,
         sd_at = c(3, 0.1))
  },
  bigin = {
    d <- mk(10, c(3), c(3, 5), 0.1, 3, function(d) {
      d$w <- 1
      d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + d$ua + d$ui)); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0,
         sd_at = c(0.1, 3))
  },
  sep = {
    # near-separable: a steep slope on x and few rows per class
    d <- mk(10, c(2, 3), c(2, 3), 0.5, 0.5, function(d) {
      d$w <- 1
      d$y <- as.integer(4 * d$x + d$ua + d$ui + rlogis(nrow(d), scale = 0.3) > 0); d })
    list(d = d, fo = bf(y ~ x + (1 | school/cls_uni)) + bernoulli(),
         lp = lp_bern(d), za = rep(1, nrow(d)), zi = rep(1, nrow(d)), off = 0)
  })
d <- cfg$d
cat(sprintf("CASE %s seed %d rows %d outer %d inner %d\n", case, seed,
            nrow(d), nlevels(d$school), nlevels(droplevels(d$cls_uni))))
mode <- Sys.getenv("MODE", "acc")
t0 <- proc.time()[[3]]
if (mode == "rounds20") {
  # in-memory only: lift the re-anchor cap to 20 rounds to see whether
  # the loop settles or cycles
  f <- ns$quad_nest_fit
  formals(f)$rounds <- 20L
  utils::assignInNamespace("quad_nest_fit", f, "frmtmb")
  mode <- "anchor"
}
if (mode == "anchor") {
  # attack 3: the re-anchor loop's rounds, the final tape's gradient,
  # and whether logLik is the integral re-anchored at the estimate
  fq <- try_msg(frm(cfg$fo, data = d, quadrature = TRUE,
                    control = frmtmb_control(verbose = 2)))
  for (w in attr(fq, "warns")) cat("WARN", w, "\n")
  par <- fq$opt$par
  g <- fq$obj$gr(par)
  cat(sprintf("ANCHOR final-tape max|grad| %.3e (grad_tol %g)\n",
              max(abs(g)), fq$control$grad_tol %||% NA))
  na <- nest_at(fq, par, fq$control$quad_nodes)
  cat(sprintf("ANCHOR logLik %.12f re-anchored-at-par %.12f diff %.3e\n",
              as.numeric(logLik(fq)), na, as.numeric(logLik(fq)) - na))
  cat(sprintf("ANCHOR obj$fn(par) %.12f\n", -fq$obj$fn(par)))
  cat(sprintf("ANCHOR convergence %s message %s\n", fq$opt$convergence,
              fq$opt$message %||% ""))
  quit()
}
if (mode == "intref") {
  fq <- try_msg(frm(cfg$fo, data = d, quadrature = TRUE))
  par <- fq$opt$par
  pl <- ns$imp_par_list(fq$frame$par_template, par, fq$frame$map)
  eta0 <- pl$beta[1] + pl$beta[2] * d$x + cfg$off
  s <- sds_of(fq, par)
  ga <- as.integer(d$school); gi <- as.integer(factor(d$cls_uni))
  # nested stats::integrate(), each integrand scaled by its value at the
  # mode so the tolerance is relative to the integral
  # log int exp(f) over a finite range that holds all but e^-60 of it,
  # split at the mode
  int_log <- function(f, lo0, hi0) {
    m <- optimize(f, c(lo0, hi0), maximum = TRUE, tol = 1e-12)
    fm <- m$objective; x0 <- m$maximum
    h <- 1e-3
    c2 <- -(f(x0 + h) - 2 * fm + f(x0 - h)) / h^2
    sc <- if (is.finite(c2) && c2 > 0) 1 / sqrt(c2) else 1
    lo <- x0 - sc; while (f(lo) > fm - 60) lo <- lo - sc
    hi <- x0 + sc; while (f(hi) > fm - 60) hi <- hi + sc
    g <- function(v) exp(vapply(v, f, 1) - fm)
    v1 <- integrate(g, lo, x0, rel.tol = 1e-10, subdivisions = 5000L)$value
    v2 <- integrate(g, x0, hi, rel.tol = 1e-10, subdivisions = 5000L)$value
    fm + log(v1 + v2)
  }
  tot <- 0
  for (a in unique(ga)) {
    ra <- which(ga == a)
    inner_log <- function(u) {
      sum(vapply(unique(gi[ra]), function(c) {
        rc <- ra[gi[ra] == c]
        fc <- function(v) dnorm(v, 0, s[[2]], log = TRUE) +
          cfg$lp(eta0[rc] + cfg$za[rc] * u + cfg$zi[rc] * v, rc)
        int_log(fc, -60, 60)
      }, 1)) + dnorm(u, 0, s[[1]], log = TRUE)
    }
    tot <- tot + int_log(inner_log, -30, 30)
  }
  cat(sprintf("INTREF %s %.10f logLik %.10f err %.3e\n", case, tot,
              as.numeric(logLik(fq)), abs(tot - as.numeric(logLik(fq)))))
  for (k in list(c(11, 15), c(21, 21), c(31, 31), c(41, 41), c(50, 50))) {
    v <- nest_at(fq, par, k)
    cat(sprintf("INTREF NODES %2d x %2d err %.3e\n", k[1], k[2], abs(v - tot)))
  }
  for (k in c(40, 60, 80, 120)) {
    cat(sprintf("INTREF myGH k %d err %.3e\n", k,
                abs(nest_ref(eta0, cfg$za, cfg$zi, ga, gi, s[[1]], s[[2]],
                             cfg$lp, k) - tot)))
  }
  quit()
}
fq <- try_msg(frm(cfg$fo, data = d, quadrature = TRUE))
cat("fit time", proc.time()[[3]] - t0, "\n")
if (length(attr(fq, "warns"))) cat("WARN", attr(fq, "warns"), sep = "\n  ")
if (inherits(fq, "error")) {
  cat("ERROR", paste(class(fq), collapse = "/"), conditionMessage(fq), "\n")
  quit()
}
cat("SETTLED", format(fq$quad_nested$settled %||% NA), "rounds",
    format(fq$quad_nested$rounds %||% NA), "\n")
par <- fq$opt$par
print(par)
ref_at <- function(par, k) {
  pl <- ns$imp_par_list(fq$frame$par_template, par, fq$frame$map)
  beta <- pl$beta
  eta0 <- beta[1] + beta[2] * d$x + cfg$off
  s <- sds_of(fq, par)
  lp <- cfg$lp
  if (is.null(lp)) lp <- lp_nb(d, exp(pl$betad[1]))
  nest_ref(eta0, cfg$za, cfg$zi, as.integer(d$school),
           as.integer(factor(d$cls_uni)), s[[1]], s[[2]], lp, k)
}
report <- function(tag, par, ll_fit = NULL) {
  r40 <- ref_at(par, 40); r60 <- ref_at(par, 60)
  la <- lap_at(fq, par)
  cat(sprintf("%s sds %s ref40 %.10f ref60 %.10f refconv %.2e laplace_err %.3e\n",
              tag, paste(format(sds_of(fq, par), digits = 4), collapse = "/"),
              r40, r60, abs(r40 - r60), abs(la - r60)))
  if (!is.null(ll_fit)) {
    cat(sprintf("%s logLik(fit) %.10f err %.3e  -opt$objective %.10f\n", tag,
                ll_fit, abs(ll_fit - r60), -fq$opt$objective))
  }
  for (k in list(c(3, 3), c(5, 5), c(7, 7), c(11, 11), c(11, 15), c(15, 15),
                 c(21, 21), c(31, 31))) {
    v <- try(nest_at(fq, par, k), silent = TRUE)
    cat(sprintf("%s NODES %2d x %2d err %s\n", tag, k[1], k[2],
                if (inherits(v, "try-error")) paste("ERROR", v) else
                  format(abs(v - r60), digits = 4)))
  }
}
report("FIT", par, as.numeric(logLik(fq)))
if (!is.null(cfg$sd_at)) {
  p2 <- par
  p2[names(p2) == "theta"] <- log(cfg$sd_at)
  report("AT", p2)
}
cat("DONE\n")
