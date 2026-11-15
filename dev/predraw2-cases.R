# Lane predraw2: the fits the detector and the reweighting are measured on.
# Sourced. case_ids(group) lists ids; case_fit(id) fits one and returns
# list(fit, nd, truth), where nd is newdata (or NULL) and truth, where the
# data are simulated from a known model, the true predictive law of nd's
# rows (a function that gives the coverage of an interval) or NULL.
#
# Groups:
#   target  the ppdraw cases: 3269 zib, A1 zi ~ x, A4 hu ~ x, the sweep's
#           three weakly identified mixtures
#   sweep   the 30 fits of dev/ppdraw-sweep.R (seed 101, 300 rows)
#   zip     zero_inflated_poisson on poisson data, seeds 1 to 30
#   wb      well-behaved fits: GLMs, GLMMs, the short check's R and S
#           designs, gaussian (1 | g) 6 groups seeds 1 to 40
#   mix     mixtures on two-component data, seeds 1 to 30 per family
fq <- function(expr) suppressMessages(suppressWarnings(expr))
n <- 300L
source("dev/ppdraw-sweep-gen.R")

case_ids <- function(group) {
  switch(group,
    target = c("zib3269", "A1", "A4", "sweep:mixture_gaussian:boundary",
               "sweep:mixture_poisson:boundary",
               "sweep:mixture_cumulative:interior"),
    sweep = as.vector(outer(fams, c("interior", "boundary"),
                            function(a, b) paste0("sweep:", a, ":", b))),
    zip = paste0("zip:", 1:30),
    wb = c(paste0("glm_gaussian:", 1:10), paste0("glm_poisson:", 1:10),
           paste0("glm_logit:", 1:10), paste0("cumul:", 1:10),
           paste0("negbin:", 1:5), paste0("gamma:", 1:5),
           paste0("beta:", 1:5), paste0("student:", 1:5),
           paste0("G6:", 1:40), paste0("R0.1:", 1:50), paste0("R0.3:", 1:50),
           paste0("S0.1:", 1:25), paste0("S0.2:", 1:25),
           paste0("pglmm:", 1:30), paste0("bglmm:", 1:30)),
    mix = c(paste0("mixg:", 1:30), paste0("mixp:", 1:30)))
}

case_fit <- function(id) {
  parts <- strsplit(id, ":")[[1]]
  kind <- parts[1]
  seed <- suppressWarnings(as.integer(parts[length(parts)]))
  nd <- NULL
  truth <- NULL
  if (kind == "sweep") {
    d <- gen(parts[2], parts[3])
    s <- spec(parts[2])
    fit <- fq(frm(s[[1]], family = s[[2]], data = d))
  } else if (kind %in% c("zib3269", "A1")) {
    set.seed(7)
    dz <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
    dz$y <- stats::rbinom(200, dz$N, stats::plogis(0.3 + 0.5 * dz$x))
    f <- if (kind == "A1") bf(y | trials(N) ~ x, zi ~ x) else
      bf(y | trials(N) ~ x)
    fit <- fq(frm(f, family = zero_inflated_binomial(), data = dz))
  } else if (kind == "A4") {
    set.seed(12)
    dh <- data.frame(x = stats::rnorm(300))
    dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
    fit <- fq(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh))
  } else if (kind == "zip") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(200))
    d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
    fit <- fq(frm(bf(y ~ x), family = zero_inflated_poisson(), data = d))
  } else if (kind == "G6") {
    set.seed(seed)
    d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
    fit <- fq(frm(bf(y ~ x + (1 | g)), data = d))
    nd <- data.frame(x = 0, g = factor("new"))
    truth <- list(mean = 1, sd = sqrt(0.1^2 + 1))
  } else if (kind %in% c("R0.1", "R0.3")) {
    sdg <- as.numeric(substring(kind, 2))
    set.seed(1000 + seed)
    d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, sdg)[d$g] + stats::rnorm(120)
    fit <- fq(frm(bf(y ~ x + (1 | g)), data = d))
    nd <- data.frame(x = 0, g = factor("new"))
    truth <- list(mean = 1, sd = sqrt(sdg^2 + 1))
  } else if (kind %in% c("S0.1", "S0.2")) {
    sds <- as.numeric(substring(kind, 2))
    set.seed(2000 + seed)
    G <- 15
    m <- 8
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    ug <- stats::rnorm(G, 0, 0.8)
    sg <- stats::rnorm(G, 0, sds)
    d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + stats::rnorm(G * m)
    fit <- fq(frm(bf(y ~ x + (1 + x | g)), data = d))
    nd <- data.frame(x = 2, g = factor("new"))
    truth <- list(mean = 2, sd = sqrt(0.8^2 + 4 * sds^2 + 1))
  } else if (kind == "glm_gaussian") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(100))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
    fit <- fq(frm(bf(y ~ x), data = d))
    nd <- data.frame(x = 2)
    truth <- list(mean = 2, sd = 1)
  } else if (kind == "glm_poisson") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(100))
    d$y <- stats::rpois(100, exp(0.5 + 0.4 * d$x))
    fit <- fq(frm(bf(y ~ x), family = poisson(), data = d))
  } else if (kind == "glm_logit") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(100))
    d$y <- stats::rbinom(100, 1, stats::plogis(-0.3 + 0.8 * d$x))
    fit <- fq(frm(bf(y ~ x), family = bernoulli(), data = d))
  } else if (kind == "cumul") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(150))
    d$y <- cut(d$x + stats::rlogis(150), c(-Inf, -1, 0.5, 2, Inf),
               labels = FALSE)
    fit <- fq(frm(bf(y ~ x), family = cumulative(), data = d))
  } else if (kind == "negbin") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(150))
    d$y <- stats::rnbinom(150, mu = exp(1 + 0.3 * d$x), size = 2)
    fit <- fq(frm(bf(y ~ x), family = negbinomial(), data = d))
  } else if (kind == "gamma") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(150))
    d$y <- stats::rgamma(150, 3, 3 / exp(0.5 + 0.3 * d$x))
    fit <- fq(frm(bf(y ~ x), family = Gamma(link = "log"), data = d))
  } else if (kind == "beta") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(150))
    m <- stats::plogis(0.2 + 0.5 * d$x)
    d$y <- stats::rbeta(150, m * 10, (1 - m) * 10)
    fit <- fq(frm(bf(y ~ x), family = Beta(), data = d))
  } else if (kind == "student") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(150))
    d$y <- 1 + 0.5 * d$x + stats::rt(150, 5)
    fit <- fq(frm(bf(y ~ x), family = student(), data = d))
  } else if (kind == "pglmm") {
    set.seed(seed)
    d <- data.frame(g = factor(rep(1:30, each = 10)), x = stats::rnorm(300))
    d$y <- stats::rpois(300, exp(0.5 + 0.3 * d$x +
                                   stats::rnorm(30, 0, 0.5)[d$g]))
    fit <- fq(frm(bf(y ~ x + (1 | g)), family = poisson(), data = d))
  } else if (kind == "bglmm") {
    set.seed(seed)
    d <- data.frame(g = factor(rep(1:20, each = 8)), x = stats::rnorm(160))
    d$y <- stats::rbinom(160, 1, stats::plogis(-0.2 + 0.7 * d$x +
                                                 stats::rnorm(20, 0, 0.8)[d$g]))
    fit <- fq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d))
  } else if (kind == "mixg") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(300))
    c2 <- stats::runif(300) < 0.4
    d$y <- ifelse(c2, 3, -1) + 0.5 * d$x + stats::rnorm(300)
    fit <- fq(frm(bf(y ~ x), family = mixture(gaussian(), gaussian()),
                  data = d))
    nd <- data.frame(x = 0)
    truth <- list(mix = list(w = c(0.6, 0.4), mean = c(-1, 3), sd = c(1, 1)))
  } else if (kind == "mixp") {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(300))
    c2 <- stats::runif(300) < 0.4
    d$y <- stats::rpois(300, exp(ifelse(c2, 2.3, 0.7) + 0.3 * d$x))
    fit <- fq(frm(bf(y ~ x), family = mixture(poisson(), poisson()),
                  data = d))
  } else {
    stop("unknown case ", id)
  }
  list(fit = fit, nd = nd, truth = truth)
}
