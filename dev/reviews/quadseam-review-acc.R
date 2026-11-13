# Review of lane quadseam: one-block quadrature accuracy and settling
# under stress. For each design and seed, fit quadrature = TRUE on one
# build and record: reported logLik, the exact marginal log-likelihood at
# the fit's own estimate (integrate() per group, centred and scaled at
# the group's mode, rel.tol 1e-12), the estimates to 17 digits,
# opt$quad_settle, the gradient of the final tape, warnings and
# messages by class. Arm "glmer" fits glmer(nAGQ = 25) instead and
# records the exact value at its estimate.
#   Rscript quadseam-review-acc.R <base|lane|glmer> <design,...> <seeds>
a <- commandArgs(TRUE)
arm <- a[1]
designs <- strsplit(a[2], ",")[[1]]
seeds <- eval(parse(text = a[3]))
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib:", arm, format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")

gen <- function(des, s) {
  set.seed(s)
  mk <- function(G, m) data.frame(g = factor(rep(seq_len(G), each = m)),
                                  x = rnorm(G * m))
  switch(des,
    bern5 = { d <- mk(5, 50); u <- rnorm(5, 0, 1.5)
      d$y <- rbinom(nrow(d), 1, plogis(-1 + 2 * d$x + u[d$g]))
      list(d = d, fam = "bern") },
    bern500 = { d <- mk(500, 4); u <- rnorm(500, 0, 1.2)
      d$y <- rbinom(nrow(d), 1, plogis(-0.5 + 0.7 * d$x + u[d$g]))
      list(d = d, fam = "bern") },
    bern1000 = { d <- mk(1000, 4); u <- rnorm(1000, 0, 1.2)
      d$y <- rbinom(nrow(d), 1, plogis(-0.5 + 0.7 * d$x + u[d$g]))
      list(d = d, fam = "bern") },
    bern20big = { d <- mk(20, 200); u <- rnorm(20, 0, 0.8)
      d$y <- rbinom(nrow(d), 1, plogis(-2 + 3 * d$x + u[d$g]))
      list(d = d, fam = "bern") },
    bern2sd3 = { d <- mk(100, 2); u <- rnorm(100, 0, 3)
      d$y <- rbinom(nrow(d), 1, plogis(0.5 + 0.5 * d$x + u[d$g]))
      list(d = d, fam = "bern") },
    bin500 = { d <- mk(10, 5); u <- rnorm(10, 0, 0.4); d$N <- 500
      d$k <- rbinom(nrow(d), 500, plogis(0.3 + 0.5 * d$x + u[d$g]))
      d$nk <- d$N - d$k; list(d = d, fam = "bin") },
    bin50 = { d <- mk(100, 3); u <- rnorm(100, 0, 1); d$N <- 50
      d$k <- rbinom(nrow(d), 50, plogis(-1 + 0.5 * d$x + u[d$g]))
      d$nk <- d$N - d$k; list(d = d, fam = "bin") },
    binz = { # a group-level covariate that moves each mode a lot
      d <- mk(6, 4); z <- c(-3, -2, -1, 1, 2, 3); u <- rnorm(6, 0, 0.5)
      d$z <- z[d$g]; d$N <- 100
      d$k <- rbinom(nrow(d), 100, plogis(1.5 * d$z + 0.2 * d$x + u[d$g]))
      d$nk <- d$N - d$k; list(d = d, fam = "binz") },
    pois8 = { d <- mk(8, 10); u <- rnorm(8, 0, 0.3)
      d$y <- rpois(nrow(d), exp(5 + 0.5 * d$x + u[d$g]))
      list(d = d, fam = "pois") },
    pois200 = { d <- mk(200, 2); u <- rnorm(200, 0, 1)
      d$y <- rpois(nrow(d), exp(1 + 0.3 * d$x + u[d$g]))
      list(d = d, fam = "pois") },
    poissd2 = { d <- mk(50, 3); u <- rnorm(50, 0, 2)
      d$y <- rpois(nrow(d), exp(-1 + 0.3 * d$x + u[d$g]))
      list(d = d, fam = "pois") },
    gam10 = { d <- mk(10, 20); u <- rnorm(10, 0, 0.5)
      d$y <- rgamma(nrow(d), 10, rate = 10 / exp(1 + 0.3 * d$x + u[d$g]))
      list(d = d, fam = "gam") },
    gam100 = { d <- mk(100, 3); u <- rnorm(100, 0, 0.3)
      d$y <- rgamma(nrow(d), 2, rate = 2 / exp(1 + 0.3 * d$x + u[d$g]))
      list(d = d, fam = "gam") },
    stop("design ", des))
}
fo <- list(bern = y ~ x + (1 | g), bin = cbind(k, nk) ~ x + (1 | g),
           binz = cbind(k, nk) ~ z + x + (1 | g), pois = y ~ x + (1 | g),
           gam = y ~ x + (1 | g))
frm_fam <- list(bern = bernoulli(), bin = binomial(), binz = binomial(),
                pois = poisson(), gam = Gamma(link = "log"))
glm_fam <- list(bern = binomial(), bin = binomial(), binz = binomial(),
                pois = poisson())

# exact marginal log-likelihood: per group, the integrand centred at its
# mode and scaled by its curvature, so narrow groups are not missed
ref_ll <- function(fam, d, beta, sd, shape = NA) {
  X <- if (fam == "binz") cbind(1, d$z, d$x) else cbind(1, d$x)
  eta <- drop(X %*% beta)
  lf <- function(y, e, dd) switch(fam,
    bern = dbinom(y, 1, plogis(e), log = TRUE),
    bin = , binz = dbinom(dd$k, dd$N, plogis(e), log = TRUE),
    pois = dpois(y, exp(e), log = TRUE),
    gam = dgamma(y, shape, rate = shape / exp(e), log = TRUE))
  if (sd < 1e-10) return(sum(lf(d$y, eta, d)))
  sum(vapply(split(seq_len(nrow(d)), d$g), function(ix) {
    dd <- d[ix, ]
    lh <- function(u) {
      v <- sum(lf(dd$y, eta[ix] + u, dd)) + dnorm(u, 0, sd, log = TRUE)
      if (is.finite(v)) v else -Inf
    }
    lo <- -30 * sd - 30; hi <- 30 * sd + 30
    m <- optimize(function(u) -lh(u), c(lo, hi), tol = 1e-12)$minimum
    m <- optimize(function(u) -lh(u), c(m - 1, m + 1), tol = 1e-14)$minimum
    h <- 1e-4 * max(sd, 1e-3)
    c2 <- -(lh(m + h) - 2 * lh(m) + lh(m - h)) / h^2
    s <- 1 / sqrt(max(c2, 1e-12))
    l0 <- lh(m)
    f <- function(t) vapply(t, function(ti) exp(lh(m + s * ti) - l0), 0)
    I <- integrate(f, -Inf, Inf, rel.tol = 1e-12, subdivisions = 2000L)
    l0 + log(s) + log(I$value)
  }, 0))
}

for (des in designs) for (s in seeds) {
  g <- gen(des, s)
  fam <- g$fam
  if (arm == "glmer") {
    if (is.null(glm_fam[[fam]])) next
    t0 <- proc.time()[[3]]
    m <- suppressMessages(suppressWarnings(
      glmer(fo[[fam]], data = g$d, family = glm_fam[[fam]], nAGQ = 25)))
    tt <- proc.time()[[3]] - t0
    beta <- unname(fixef(m)); sd <- unname(getME(m, "theta"))
    rf <- ref_ll(fam, g$d, beta, sd)
    cat(sprintf(paste0("ROW\t%s\t%d\tglmer\tll\t%.12f\tref\t%.12f\t",
                       "dref\t%.3e\tsd\t%.10g\tbeta\t%s\tsing\t%s\t",
                       "secs\t%.2f\n"),
                des, s, as.numeric(logLik(m)), rf,
                as.numeric(logLik(m)) - rf, sd,
                paste(sprintf("%.10g", beta), collapse = ","),
                isSingular(m), tt))
    next
  }
  w <- character(); ms <- character()
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(
    frm(fo[[fam]], family = frm_fam[[fam]], data = g$d, quadrature = TRUE),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, class(x)[1]); invokeRestart("muffleMessage")
    }), error = function(e) e)
  tt <- proc.time()[[3]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("ROW\t%s\t%d\t%s\tERROR\t%s\n", des, s, arm,
                substr(gsub("\\s+", " ", conditionMessage(f)), 1, 200)))
    next
  }
  beta <- unname(f$estimates$beta); sd <- exp(unname(f$estimates$theta[1]))
  shape <- if (fam == "gam") exp(unname(f$estimates$betad[1])) else NA
  rf <- ref_ll(fam, g$d, beta, sd, shape)
  ll <- as.numeric(logLik(f))
  gr <- tryCatch(max(abs(f$obj$gr(f$opt$par))), error = function(e) NA)
  qs <- f$opt[["quad_settle"]]
  cat(sprintf(paste0("ROW\t%s\t%d\t%s\tll\t%.12f\tref\t%.12f\tdref\t%.3e\t",
                     "sd\t%.10g\tbeta\t%s\tshape\t%.10g\tconv\t%d\t",
                     "rounds\t%s\tsgain\t%s\tsgrad\t%s\tgr\t%.3e\t",
                     "warn\t%s\tmsg\t%s\tpar\t%s\tsecs\t%.2f\n"),
              des, s, arm, ll, rf, ll - rf, sd,
              paste(sprintf("%.10g", beta), collapse = ","), shape,
              f$opt$convergence,
              if (is.null(qs)) "-" else qs[["rounds"]],
              if (is.null(qs)) "-" else format(qs[["gain"]], digits = 3),
              if (is.null(qs)) "-" else format(qs[["grad"]], digits = 3),
              gr,
              if (length(w)) paste(w, collapse = ",") else "-",
              if (length(ms)) paste(ms, collapse = ",") else "-",
              paste(sprintf("%a", f$opt$par), collapse = ","), tt))
}
cat("END\n")
