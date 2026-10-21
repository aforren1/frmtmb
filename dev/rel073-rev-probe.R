# Release reviewer, 0.73.0: cross-lane probes on a build.
#   Rscript dev/rel073-rev-probe.R <lib> <probe>
# Logs: dev/rel073-rev-log/probe-<probe>-<libtag>.txt (by the caller).
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]; probe <- a[2]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "| probe", probe, "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
# run, report the error message or the warnings
try_fit <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(label, ": ERROR:", gsub("\\s+", " ", substr(r, 1, 400)), "\n")
  } else {
    cat(label, ": ok", if (length(w)) paste0("| ", length(w), " warnings: ",
      paste(unique(substr(gsub("\\s+", " ", w), 1, 120)), collapse = " || ")),
      "\n")
  }
  invisible(r)
}
ll <- function(f) as.numeric(logLik(f))
mu_of <- function(f, X) as.vector(exp(X %*% f$estimates$beta))

switch(probe,
# P1: negbinomial2 censored and truncated against stats:: at the fit's
# own estimates (an independent reference, not negbinomial's code)
nb2ref = {
  set.seed(20261021)
  n <- 300
  d <- data.frame(x = rnorm(n))
  d$y <- rnbinom(n, mu = exp(1.2 + 0.4 * d$x), size = 2.5)
  X <- cbind(1, d$x)
  d$cc <- ifelse(d$y >= 7, "right", ifelse(d$y <= 1, "left", "none"))
  d$y2 <- d$y + 3L
  d$cc[seq(3, n, by = 17)] <- "interval"
  f <- q(frm(bf(y | cens(cc, y2) ~ x), family = brmsfamily("negbinomial2"),
             data = d))
  mu <- mu_of(f, X); size <- exp(-f$estimates$betad[[1]])
  ref <- ifelse(d$cc == "none", dnbinom(d$y, size = size, mu = mu, log = TRUE),
    ifelse(d$cc == "right", pnbinom(d$y, size, mu = mu, lower.tail = FALSE,
                                    log.p = TRUE),
    ifelse(d$cc == "left", pnbinom(d$y, size, mu = mu, log.p = TRUE),
           log(pnbinom(d$y2, size, mu = mu) - pnbinom(d$y, size, mu = mu)))))
  cat(sprintf("cens nb2: logLik %.10f ref %.10f diff %.2e | sigma %.5f\n",
              ll(f), sum(ref), ll(f) - sum(ref), 1 / size))
  # brms's reading: right at k is P(Y > k); the 0.72.0 reading P(Y >= k)
  ref_old <- ifelse(d$cc == "right", pnbinom(d$y - 1, size, mu = mu,
                    lower.tail = FALSE, log.p = TRUE), ref)
  cat(sprintf("  the inclusive reading at these estimates: %.4f\n",
              sum(ref_old)))
  dt <- d[d$y >= 2 & d$y <= 30, ]
  Xt <- cbind(1, dt$x)
  for (bd in list(c(2, Inf), c(2, 30))) {
    fm <- if (is.finite(bd[2])) bf(y | trunc(lb = 2, ub = 30) ~ x) else
      bf(y | trunc(lb = 2) ~ x)
    f <- q(frm(fm, family = brmsfamily("negbinomial2"), data = dt))
    mu <- mu_of(f, Xt); size <- exp(-f$estimates$betad[[1]])
    lo <- pnbinom(bd[1] - 1, size, mu = mu)
    hi <- if (is.finite(bd[2])) pnbinom(bd[2], size, mu = mu) else 1
    ref <- dnbinom(dt$y, size, mu = mu, log = TRUE) - log(hi - lo)
    # the truncated mean, summed exactly over the window
    tm <- vapply(seq_along(mu), function(i) {
      k <- bd[1]:(if (is.finite(bd[2])) bd[2] else 5000)
      p <- dnbinom(k, size, mu = mu[i]); sum(k * p) / sum(p)
    }, 0)
    ft <- fitted(f)[, 1]
    cat(sprintf(paste0("trunc nb2 [%g, %g]: logLik %.10f ref %.10f diff ",
                       "%.2e | fitted truncated mean max rel err %.2e\n"),
                bd[1], bd[2], ll(f), sum(ref), ll(f) - sum(ref),
                max(abs(ft / tm - 1))))
  }
  # refusals and acceptances on negbinomial2, as on negbinomial
  d0 <- d; d0$y2 <- d0$y
  d0$cc <- "interval"
  try_fit("nb2 interval y2 == y", frm(bf(y | cens(cc, y2) ~ x),
          family = brmsfamily("negbinomial2"), data = d0))
  d1 <- d; d1$yb <- d1$y + 0.5; d1$cc <- "right"
  try_fit("nb2 fractional cens bound", frm(bf(yb | cens(cc) ~ x),
          family = brmsfamily("negbinomial2"), data = d1))
  try_fit("nb2 trunc(lb = 0)", frm(bf(y | trunc(lb = 0) ~ x),
          family = brmsfamily("negbinomial2"), data = d))
  try_fit("nb2 trunc(lb = 3) with y below", frm(bf(y | trunc(lb = 3) ~ x),
          family = brmsfamily("negbinomial2"), data = d))
  try_fit("nb2 trunc(lb = 5, ub = 5)", frm(bf(y | trunc(lb = 5, ub = 5) ~ x),
          family = brmsfamily("negbinomial2"), data = d[d$y == 5, ]))
  cmp <- frm_compat()
  r <- cmp[cmp$feature_b %in% c("negbinomial2", "dirichlet_multinomial") &
             cmp$feature_a %in% c("cens()", "trunc()") |
           cmp$feature_a %in% c("negbinomial2", "dirichlet_multinomial") &
             cmp$feature_b %in% c("cens()", "trunc()"), ]
  print(r[, intersect(c("feature_a", "feature_b", "status", "note"),
                      names(r))], right = FALSE)
},
# P1b: dirichlet_multinomial with cens() or trunc()
dmcens = {
  set.seed(2)
  n <- 120
  d <- data.frame(x = rnorm(n), N = 20L)
  P <- t(apply(cbind(0, 0.5 + 0.3 * d$x, -0.2 + 0.1 * d$x), 1,
               function(e) exp(e) / sum(exp(e))))
  Y <- t(vapply(seq_len(n), function(i) rmultinom(1, 20, P[i, ]), numeric(3)))
  d$Y <- Y; d$cc <- "none"
  try_fit("DM plain", f <- frm(bf(Y | trials(N) ~ x),
          family = dirichlet_multinomial(), data = d))
  try_fit("DM cens", frm(bf(Y | trials(N) + cens(cc) ~ x),
          family = dirichlet_multinomial(), data = d))
  try_fit("DM trunc", frm(bf(Y | trials(N) + trunc(lb = 0) ~ x),
          family = dirichlet_multinomial(), data = d))
  try_fit("DM weights + obs", {
    d$obs <- factor(seq_len(n)); d$w <- 1.5
    frm(bf(Y | trials(N) + weights(w) ~ x + (1 | obs)),
        family = dirichlet_multinomial(), data = d)
  })
  try_fit("DM quadrature (1 | g)", {
    d$g <- factor(rep(1:12, length.out = n))
    frm(bf(Y | trials(N) ~ x + (1 | g)), family = dirichlet_multinomial(),
        data = d, quadrature = TRUE)
  })
},
# P2: the weights refusal with the window code (cens, trunc)
wcens = {
  set.seed(20261021)
  n <- 80
  d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)), w = 1.5)
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  d$gy <- rgamma(n, 3, 3 / exp(0.3 + 0.2 * d$x))
  d$cnone <- "none"
  d$cmix <- ifelse(d$y > quantile(d$y, 0.8), "right", "none")
  d$yc <- pmin(d$y, quantile(d$y, 0.8))
  d$call <- "right"
  for (cs in list(
    list("gaussian, no cens (refused)", bf(y | weights(w) ~ x + (1 | obs)),
         gaussian()),
    list("gaussian, cens all none", bf(y | weights(w) + cens(cnone) ~ x +
                                         (1 | obs)), gaussian()),
    list("gaussian, 20% right censored", bf(yc | weights(w) + cens(cmix) ~ x +
                                             (1 | obs)), gaussian()),
    list("gaussian, all right censored", bf(yc | weights(w) + cens(call) ~ x +
                                             (1 | obs)), gaussian()),
    list("gaussian, trunc(lb = -10)", bf(y | weights(w) + trunc(lb = -10) ~ x +
                                           (1 | obs)), gaussian()),
    list("gaussian, trunc(lb = -10), w = 1", bf(y | trunc(lb = -10) ~ x +
                                           (1 | obs)), gaussian()),
    list("Gamma, cens all none", bf(gy | weights(w) + cens(cnone) ~ x +
                                      (1 | obs)), Gamma(link = "log")),
    list("Gamma, trunc(ub = 50)", bf(gy | weights(w) + trunc(ub = 50) ~ x +
                                       (1 | obs)), Gamma(link = "log")),
    list("student, cens all none", bf(y | weights(w) + cens(cnone) ~ x +
                                        (1 | obs)), student()))) {
    r <- try_fit(cs[[1]], f <- frm(cs[[2]], family = cs[[3]], data = d))
    if (!inherits(r, "err")) {
      dp <- f$estimates$betad
      cat("   logLik", format(ll(f), digits = 8), "| betad",
          paste(format(dp, digits = 4), collapse = " "), "| conv",
          f$opt$convergence, "\n")
    }
  }
},
# P3: quadrature on the censored families, one block per row, and the
# nested refusal before any tape, with time and the process's memory
quadcens = {
  set.seed(20261021)
  ng <- 30; n <- ng * 6
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = 6)),
                  h = factor(rep(1:6, ng)))
  d$gb <- factor(paste(d$g, rep(1:2, length.out = n)))
  u <- rnorm(ng, 0, 0.6)[d$g]
  d$yb <- plogis(0.2 + 0.4 * d$x + u); d$yb <- rbeta(n, d$yb * 8,
                                                     (1 - d$yb) * 8)
  d$cb <- ifelse(d$yb > 0.8, "right", "none"); d$ybc <- pmin(d$yb, 0.8)
  d$k <- rnbinom(n, mu = exp(1 + 0.3 * d$x + u), size = 3)
  d$ck <- ifelse(d$k >= 8, "right", "none"); d$kc <- pmin(d$k, 8)
  d$gy <- rgamma(n, 3, 3 / exp(0.3 + 0.2 * d$x + u))
  d$cg <- ifelse(d$gy > 3, "right", "none"); d$gyc <- pmin(d$gy, 3)
  cases <- list(
    list("Beta cens (1 | g)", bf(ybc | cens(cb) ~ x + (1 | g)), Beta()),
    list("negbinomial2 cens (1 | g)", bf(kc | cens(ck) ~ x + (1 | g)),
         brmsfamily("negbinomial2")),
    list("negbinomial2 trunc(ub = 40) (1 | g)", bf(k | trunc(ub = 40) ~ x +
                                                  (1 | g)),
         brmsfamily("negbinomial2")),
    list("Gamma cens (1 | g)", bf(gyc | cens(cg) ~ x + (1 | g)),
         Gamma(link = "log")))
  for (cs in cases) {
    t0 <- proc.time()[[3]]
    fl <- q(frm(cs[[2]], family = cs[[3]], data = d))
    t1 <- proc.time()[[3]]
    r <- try_fit(paste(cs[[1]], "quadrature"),
                 fq <- frm(cs[[2]], family = cs[[3]], data = d,
                           quadrature = TRUE))
    t2 <- proc.time()[[3]]
    if (!inherits(r, "err")) {
      cat(sprintf(paste0("   Laplace %.6f (%.1fs) quadrature %.6f (%.1fs) ",
                         "sd laplace %.4f quad %.4f\n"),
                  ll(fl), t1 - t0, ll(fq), t2 - t1,
                  exp(fl$estimates$theta[1]), exp(fq$estimates$theta[1])))
    }
  }
  for (cs in list(
    list("Beta cens (1 | g/gb) quadrature", bf(ybc | cens(cb) ~ x +
                                               (1 | g/gb)), Beta()),
    list("nb2 cens (1 | g) + (1 | h) quadrature", bf(kc | cens(ck) ~ x +
                                                 (1 | g) + (1 | h)),
         brmsfamily("negbinomial2")))) {
    t0 <- proc.time()[[3]]
    try_fit(cs[[1]], frm(cs[[2]], family = cs[[3]], data = d,
                         quadrature = TRUE))
    cat(sprintf("   refused after %.2fs\n", proc.time()[[3]] - t0))
  }
  cat("memory (Mb used by R, max):", paste(format(gc()[, 6]), collapse = " "),
      "\n")
},
# P6: simulate()'s inverse fallback for negbinomial2 under trunc()
simtrunc = {
  set.seed(20261021)
  n <- 200
  d <- data.frame(x = rnorm(n))
  sz <- 4
  d$y <- rnbinom(n, mu = exp(2.5 + 0.8 * d$x), size = sz)
  # an upper bound far below the mean on many rows: rejection keeps a
  # few percent there, so some rows reach the inverse fallback
  ub <- 6
  dt <- d[d$y <= ub, ]
  cat("rows kept:", nrow(dt), "\n")
  for (fam in list(negbinomial(), brmsfamily("negbinomial2"))) {
    f <- q(frm(bf(y | trunc(ub = 6) ~ x), family = fam, data = dt))
    X <- cbind(1, dt$x)
    mu <- mu_of(f, X)
    size <- if (identical(fam$family, "negbinomial2"))
      exp(-f$estimates$betad[[1]]) else exp(f$estimates$betad[[1]])
    mass <- pnbinom(ub, size, mu = mu)
    set.seed(7)
    S <- simulate(f, nsim = as.integer(Sys.getenv("REV_NSIM", "400")))
    S <- as.matrix(S)
    # the exact truncated mean and the simulated one, per row
    tm <- vapply(seq_along(mu), function(i) {
      k <- 0:ub; p <- dnbinom(k, size, mu = mu[i]); sum(k * p) / sum(p)
    }, 0)
    sm <- rowMeans(S)
    tv <- vapply(seq_along(mu), function(i) {
      k <- 0:ub; p <- dnbinom(k, size, mu = mu[i]); p <- p / sum(p)
      sum(k^2 * p) - sum(k * p)^2
    }, 0)
    z <- (sm - tm) / sqrt(tv / ncol(S))
    lowm <- mass < 0.05
    cat(sprintf(paste0("%-13s size %.4f | rows with mass < 5%%: %d | out of ",
                       "window: %d | z of the row means: max |z| %.2f, ",
                       "mean z %.3f; on the low-mass rows max |z| %.2f, ",
                       "mean %.3f\n"),
                fam$family, size, sum(lowm), sum(S > ub | S < 0),
                max(abs(z)), mean(z), if (any(lowm)) max(abs(z[lowm])) else NA,
                if (any(lowm)) mean(z[lowm]) else NA))
    # a whole-sample chi-square on the pooled low-mass rows
    if (any(lowm)) {
      pr <- rowMeans(vapply(which(lowm), function(i) {
        p <- dnbinom(0:ub, size, mu = mu[i]); p / sum(p)
      }, numeric(ub + 1)))
      ob <- tabulate(S[lowm, ] + 1L, ub + 1L)
      ex <- pr * sum(ob)
      cat("   low-mass rows pooled: chi-square", format(sum((ob - ex)^2 / ex),
          digits = 5), "on", ub, "df\n")
    }
  }
},
# P8: dirichlet and dirichlet_multinomial draws at small phi
dirdraws = {
  set.seed(3)
  n <- 200
  d <- data.frame(x = rnorm(n), N = 30L)
  P <- t(apply(cbind(0, 0.5 + 0.3 * d$x, -0.2), 1,
               function(e) exp(e) / sum(exp(e))))
  # overdispersed counts at a small phi
  Y <- t(vapply(seq_len(n), function(i) {
    p <- rgamma(3, P[i, ] * 0.05); p <- p / sum(p)
    if (any(!is.finite(p))) p <- c(1, 0, 0)
    rmultinom(1, 30, p)
  }, numeric(3)))
  d$Y <- Y
  r <- try_fit("DM fit phi small", f <- frm(bf(Y | trials(N) ~ x),
                        family = dirichlet_multinomial(), data = d))
  if (!inherits(r, "err")) {
    cat("   phi", format(exp(f$estimates$betad), digits = 4), "\n")
    s <- try_fit("DM simulate", S <- simulate(f, nsim = 20, seed = 1))
    if (!inherits(s, "err")) cat("   NA in simulate:", sum(is.na(unlist(S))),
                                 "\n")
    try_fit("DM predict", pr <- predict(f))
  }
},
stop("unknown probe ", probe))
cat("END\n")
