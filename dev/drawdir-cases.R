# Lane drawdir: the fits the per-direction draw check is measured on.
# Sourced from the worktree root. Extends dev/predraw2-cases.R (groups
# target, sweep, zip, wb, mix) with:
#   mix40   the three mixture designs of predraw2 section 3 (sweep
#           generator): mixture_gaussian:boundary, mixture_poisson:boundary,
#           mixture_cumulative:interior, seeds 1 to 40
#   partR   the short check's part R: the 18 fired mixture fits' data with
#           K = 0, 5, 10, 20 noise covariates in both components
#           (dev/reviews/predraw2-short-a.R)
#   kd      the K designs of dev/predraw2-kdesign.R (y ~ 0 + f)
#   gglm    gaussian GLMs: n 20, 50, 200; k 1, 5, 15 covariates; seeds 1:12;
#           and sigma ~ x at n 30, 100, seeds 1:10
#   gglmm   gaussian GLMMs: (1 | g) with G 5, 10, 30 groups of m 4 or 10
#           rows, sd 0.3 or 1, seeds 1:8; (1 + x | g) with G 8, 20, 40 of 8
#           rows, slope sd 0.2 or 0.5, seeds 1:8
# dc_ids(group) lists ids; dc_fit(id) returns list(fit, d, truth_mean).
source("dev/predraw2-cases.R")

dc_fired <- list(mixture_gaussian = c(20, 23, 25, 29, 36, 40),
                 mixture_poisson = c(5, 11, 13, 18, 20, 24, 25, 26, 31, 33,
                                     36, 39))

dc_kdesigns <- list(
  pois_K100_m5_mu1 = list(K = 100, m = 5, mu = 1, fam = "poisson"),
  pois_K100_m5_mu1.5 = list(K = 100, m = 5, mu = 1.5, fam = "poisson"),
  pois_K200_m5_mu2 = list(K = 200, m = 5, mu = 2, fam = "poisson"),
  pois_K500_m5_mu2 = list(K = 500, m = 5, mu = 2, fam = "poisson"),
  pois_K50_m5_mu2 = list(K = 50, m = 5, mu = 2, fam = "poisson"),
  pois_K100_m10_mu1.5 = list(K = 100, m = 10, mu = 1.5, fam = "poisson"),
  pois_K50_m3_mu1 = list(K = 50, m = 3, mu = 1, fam = "poisson"),
  bern_K60_m20_p03 = list(K = 60, m = 20, mu = 0.3, fam = "bernoulli"),
  bern_K20_m10_p03 = list(K = 20, m = 10, mu = 0.3, fam = "bernoulli"),
  gaus_K10_m3 = list(K = 10, m = 3, mu = 1, fam = "gaussian"),
  gaus_K50_m3 = list(K = 50, m = 3, mu = 1, fam = "gaussian"),
  gaus_K200_m3 = list(K = 200, m = 3, mu = 1, fam = "gaussian"),
  gaus_K500_m3 = list(K = 500, m = 3, mu = 1, fam = "gaussian"),
  gamma_K10_m3 = list(K = 10, m = 3, mu = 1, fam = "gamma"),
  gamma_K50_m3 = list(K = 50, m = 3, mu = 1, fam = "gamma"),
  gamma_K200_m3 = list(K = 200, m = 3, mu = 1, fam = "gamma"),
  gamma_K500_m3 = list(K = 500, m = 3, mu = 1, fam = "gamma"),
  gamma_K20_m2 = list(K = 20, m = 2, mu = 1, fam = "gamma"),
  gamma_K50_m2 = list(K = 50, m = 2, mu = 1, fam = "gamma"))

dc_ids <- function(group) {
  switch(group,
    mix40 = c(paste0("mix40:mixture_gaussian:boundary:", 1:40),
              paste0("mix40:mixture_poisson:boundary:", 1:40),
              paste0("mix40:mixture_cumulative:interior:", 1:40)),
    partR = unlist(lapply(names(dc_fired), function(f) {
      as.vector(outer(dc_fired[[f]], c(0, 5, 10, 20), function(s, k) {
        paste0("partR:", f, ":", k, ":", s)
      }))
    })),
    kd = as.vector(outer(names(dc_kdesigns), 1:10, function(a, s) {
      paste0("kd:", a, ":", s)
    })),
    gglm = c(as.vector(outer(c("20:1", "20:5", "50:1", "50:5", "50:15",
                               "200:1", "200:5", "200:15", "20:15"), 1:12,
                             function(a, s) paste0("gglm:", a, ":", s))),
             as.vector(outer(c("30", "100"), 1:10, function(a, s) {
               paste0("gglms:", a, ":", s)
             }))),
    gglmm = c(as.vector(outer(c("5:4:0.3", "5:4:1", "5:10:0.3", "5:10:1",
                                "10:4:0.3", "10:4:1", "10:10:0.3", "10:10:1",
                                "30:4:0.3", "30:4:1", "30:10:0.3",
                                "30:10:1"), 1:8,
                              function(a, s) paste0("gri:", a, ":", s))),
              as.vector(outer(c("8:0.2", "8:0.5", "20:0.2", "20:0.5",
                                "40:0.2", "40:0.5"), 1:8,
                              function(a, s) paste0("grs:", a, ":", s)))),
    case_ids(group))
}

dc_fit <- function(id) {
  p <- strsplit(id, ":")[[1]]
  kind <- p[1]
  seed <- as.integer(p[length(p)])
  if (kind == "mix40") {
    d <- gen(p[2], p[3], seed)
    s <- spec(p[2])
    return(list(fit = fq(frm(s[[1]], family = s[[2]], data = d)), d = d))
  }
  if (kind == "partR") {
    fam <- p[2]
    K <- as.integer(p[3])
    d <- gen(fam, "boundary", seed)
    set.seed(5000 + seed)
    zs <- character(0)
    for (j in seq_len(K)) {
      zs <- c(zs, paste0("z", j))
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    fo <- stats::as.formula(paste("y ~", paste(c("x", zs), collapse = " + ")))
    tm <- if (fam == "mixture_gaussian") 0.5 * d$x else exp(1 + 0.3 * d$x)
    return(list(fit = fq(frm(bf(fo), family = spec(fam)[[2]], data = d)),
                d = d, truth_mean = tm))
  }
  if (kind == "kd") {
    ds <- dc_kdesigns[[p[2]]]
    set.seed(seed)
    d <- data.frame(f = factor(rep(seq_len(ds$K), each = ds$m)))
    if (ds$fam == "poisson") {
      d$y <- stats::rpois(nrow(d), exp(log(ds$mu) +
                                         stats::rnorm(ds$K, 0, 0.3))[d$f])
      fam <- poisson()
    } else if (ds$fam == "bernoulli") {
      d$y <- stats::rbinom(nrow(d), 1, stats::plogis(
        stats::qlogis(ds$mu) + stats::rnorm(ds$K, 0, 0.3))[d$f])
      fam <- bernoulli()
    } else if (ds$fam == "gaussian") {
      d$y <- stats::rnorm(nrow(d), (ds$mu + stats::rnorm(ds$K, 0, 1))[d$f], 1)
      fam <- gaussian()
    } else {
      mu <- exp(log(ds$mu) + stats::rnorm(ds$K, 0, 0.5))[d$f]
      d$y <- stats::rgamma(nrow(d), 2, 2 / mu)
      fam <- Gamma(link = "log")
    }
    return(list(fit = fq(frm(bf(y ~ 0 + f), family = fam, data = d)), d = d))
  }
  if (kind == "gglm") {
    nn <- as.integer(p[2])
    k <- as.integer(p[3])
    set.seed(seed)
    X <- matrix(stats::rnorm(nn * k), nn, k)
    colnames(X) <- paste0("x", seq_len(k))
    d <- as.data.frame(X)
    d$y <- 1 + drop(X %*% rep(0.3, k)) + stats::rnorm(nn)
    fo <- stats::as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
    return(list(fit = fq(frm(bf(fo), data = d)), d = d))
  }
  if (kind == "gglms") {
    nn <- as.integer(p[2])
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(nn))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(nn, 0, exp(0.3 * d$x))
    return(list(fit = fq(frm(bf(y ~ x, sigma ~ x), data = d)), d = d))
  }
  if (kind == "gri") {
    G <- as.integer(p[2])
    m <- as.integer(p[3])
    sg <- as.numeric(p[4])
    set.seed(seed)
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, sg)[d$g] + stats::rnorm(G * m)
    return(list(fit = fq(frm(bf(y ~ x + (1 | g)), data = d)), d = d))
  }
  if (kind == "grs") {
    G <- as.integer(p[2])
    ss <- as.numeric(p[3])
    set.seed(seed)
    d <- data.frame(g = factor(rep(seq_len(G), each = 8)),
                    x = stats::rnorm(G * 8))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, 0.8)[d$g] +
      stats::rnorm(G, 0, ss)[d$g] * d$x + stats::rnorm(G * 8)
    return(list(fit = fq(frm(bf(y ~ x + (1 + x | g)), data = d)), d = d))
  }
  r <- case_fit(id)
  list(fit = r$fit, d = NULL)
}
