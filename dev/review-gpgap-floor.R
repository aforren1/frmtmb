# Review of lane gpgap, item 3: where do the two nuggets' fits begin to
# differ, on seeds and designs that are not the lane's, against the warning.
# Closed-form marginal maxima (sigma conditions the covariance, so either
# nugget is readable), frm() itself, and the warning's class.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-floor.R <log> <arm>
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
arm <- commandArgs(trailingOnly = TRUE)[1]

fit_closed <- function(x, y, s0, A, nugget, rep_pos = TRUE) {
  u <- sort(unique(x))
  D2 <- outer(u, u, "-")^2
  J <- match(x, u)
  n <- length(y)
  ybar <- mean(y)
  nll <- function(p) {
    s2 <- exp(2 * p[2])
    K <- s2 * exp(-D2 / (2 * exp(2 * p[3])))
    diag(K) <- diag(K) + nugget(s2)
    V <- K[J, J]
    diag(V) <- diag(V) + exp(2 * p[4])
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
  }
  w <- diff(range(x))
  best <- NULL
  for (st in list(c(ybar, log(A), log(w / 4), log(s0)),
                  c(ybar, log(A) + 0.7, log(w / 3), log(s0)),
                  c(ybar, log(A) - 0.5, log(w / 6), log(s0) - 0.5),
                  c(ybar, log(A) + 1.5, log(w / 2.5), log(s0)))) {
    o <- stats::nlminb(st, nll,
                       lower = c(ybar - 100 * A, log(A) - 6, log(w) - 7,
                                 log(s0) - 25),
                       upper = c(ybar + 100 * A, log(A) + 6, log(w) + 3,
                                 log(s0) + 8),
                       control = list(iter.max = 800, eval.max = 1200))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  list(par = best$par, obj = best$objective, nll = nll)
}

num_hess_se <- function(f, p, i) {
  # the profile-free curvature of log sigma, by central differences on the
  # full Hessian, on the marginal's own scale
  h <- 1e-3
  k <- length(p)
  H <- matrix(0, k, k)
  for (a in seq_len(k)) for (b in a:k) {
    ea <- eb <- numeric(k)
    ea[a] <- h; eb[b] <- h
    H[a, b] <- H[b, a] <- (f(p + ea + eb) - f(p + ea - eb) - f(p - ea + eb) +
                             f(p - ea - eb)) / (4 * h^2)
  }
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || V[i, i] <= 0) return(NA_real_)
  sqrt(V[i, i])
}

one <- function(seed, n, ratio, A = 1, xs = 6, tag = "", rep = 1L) {
  set.seed(seed)
  x <- round(runif(n %/% rep, 0, xs), 6)
  x <- rep(x, each = rep)
  x <- x[sample(length(x))]
  s <- ratio * A
  y <- A * sin(2 * x * 6 / xs) + s * rnorm(length(x))
  a <- fit_closed(x, y, s, A, function(s2) 1e-6 * s2)
  b <- fit_closed(x, y, s, A, function(s2) 1e-12)
  sb <- exp(b$par[4]); sf <- exp(a$par[4])
  se <- num_hess_se(b$nll, b$par, 4)
  warned <- character()
  fr <- withCallingHandlers(
    tryCatch(frm(bf(y ~ gp(x)), family = gaussian(), data = data.frame(x, y)),
             error = function(e) NULL),
    warning = function(w) {
      warned <<- c(warned, class(w)[1])
      invokeRestart("muffleWarning")
    })
  sg <- if (is.null(fr)) NA else as.numeric(sigma(fr))
  sdg <- if (is.null(fr)) NA else {
    bk <- Filter(function(k) identical(k$covstruct, "gp"),
                 fr$frame$re_blocks)[[1]]
    exp(fr$estimates$theta[bk$theta_idx][[1]])
  }
  cat(sprintf(paste0("FLOOR %s seed %d n %d A %g xs %g rep %d ratio %.1e | ",
                     "brms: sigma %.4e sdgp %.4g rho %.4g ll %.3f | ",
                     "frmnug: sigma %.4e sdgp %.4g rho %.4g ll %.3f | ",
                     "sigma/brms %.4f, |dlog|/SE %.2f, brms-at-frm loss %.3g, ",
                     "sigma_b/sdgp_b %.2e | frm(): sigma %.4e sdgp %.4g ",
                     "ll %.3f conv %s, sigma/sdgp %.2e warn %s nugwarn %s\n"),
              tag, seed, length(y), A, xs, rep, ratio,
              sb, exp(b$par[2]), exp(b$par[3]), -b$obj,
              sf, exp(a$par[2]), exp(a$par[3]), -a$obj,
              sf / sb, abs(log(sf / sb)) / se, b$nll(a$par) - b$obj,
              sb / exp(b$par[2]),
              sg, sdg, if (is.null(fr)) NA else -fr$opt$objective,
              if (is.null(fr)) "NA" else fr$opt$convergence, sg / sdg,
              paste(unique(warned), collapse = "/"),
              "frmtmb_gp_nugget" %in% warned))
}

ratios <- c(0.1, 0.03, 0.02, 0.012, 0.01, 0.008, 0.005, 0.003, 0.002, 0.0015,
            0.001, 3e-4)
switch(arm,
  seeds = for (sd in c(101L, 202L, 303L, 404L)) {
    for (r in ratios) one(sd, 60, r, tag = "seeds")
  },
  n = for (nn in c(30L, 120L, 300L)) {
    for (r in ratios) one(111L, nn, r, tag = "n")
  },
  scale = {
    for (r in ratios) one(121L, 60, r, A = 1000, tag = "Ascale")
    for (r in ratios) one(122L, 60, r, A = 1e-3, tag = "Ascale")
    for (r in ratios) one(123L, 60, r, xs = 6000, tag = "xscale")
    for (r in ratios) one(124L, 60, r, xs = 0.06, tag = "xscale")
  },
  rep = {
    # replicated positions: the nugget is one draw per POSITION, so it cannot
    # take its variance from sigma, which replication measures directly
    for (r in c(0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 3e-4)) {
      one(131L, 600, r, rep = 30L, tag = "rep30")
    }
    for (r in c(0.01, 0.005, 0.003, 0.001)) {
      one(132L, 240, r, rep = 4L, tag = "rep4")
    }
  })
