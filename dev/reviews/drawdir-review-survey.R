# Reviewer of lane drawdir: a false-alarm survey of the per-direction draw
# check on designs the lane did not use. Per fit: predict() (200 draws,
# seed 5) on the lane library, the cached verdict, what was said. Where
# the per-direction reading fires (the summed one silent), both arms at
# 2000 draws, seed 1: the lane's (reweighted) and the normal law's (the
# verdict's fire set FALSE in the cache, which is the 0.83.0 behavior on
# a fit the summed check does not fire on), against the truth: in-sample
# MAE of the predictive mean against the true conditional mean, and the
# probability the oracle predictive law of a new level puts inside the
# central 90% and 98% intervals, with their widths.
# Usage: Rscript dev/reviews/drawdir-review-survey.R <out.tsv> <group> [...]
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
out <- args[1]
groups <- args[-1]
fq <- function(expr) suppressMessages(suppressWarnings(expr))

# ---- designs: each returns list(fit, mu (true E[y] per row), nd, onew
# (a function n -> oracle draws of y at nd's one row, a new level))
ilogit <- stats::plogis
cum_mean <- function(eta, th) {
  # E[category] of a cumulative logit with thresholds th
  P <- sapply(th, function(t) ilogit(t - eta))
  P <- cbind(P, 1)
  pk <- cbind(P[, 1], P[, -1, drop = FALSE] - P[, -ncol(P), drop = FALSE])
  drop(pk %*% seq_len(ncol(pk)))
}
cum_draw <- function(eta, th) {
  u <- stats::rlogis(length(eta))
  1L + rowSums(outer(eta + u, th, ">"))
}

design <- function(kind, seed) {
  set.seed(seed)
  if (kind %in% c("bfew4", "bfew6")) {
    G <- if (kind == "bfew4") 4L else 6L
    m <- 20L
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u <- stats::rnorm(G, 0, 0.8)
    eta <- -0.3 + 0.7 * d$x + u[d$g]
    d$y <- stats::rbinom(nrow(d), 1, ilogit(eta))
    fit <- fq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d))
    return(list(fit = fit, mu = ilogit(eta),
                nd = data.frame(x = 1, g = factor("new")),
                onew = function(n) stats::rbinom(n, 1, ilogit(
                  -0.3 + 0.7 + stats::rnorm(n, 0, 0.8)))))
  }
  if (kind == "binfew") {
    G <- 5L
    m <- 6L
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m), N = 8L)
    u <- stats::rnorm(G, 0, 0.6)
    eta <- 0.2 + 0.5 * d$x + u[d$g]
    d$y <- stats::rbinom(nrow(d), 8, ilogit(eta))
    fit <- fq(frm(bf(y | trials(N) ~ x + (1 | g)), family = binomial(),
                  data = d))
    return(list(fit = fit, mu = 8 * ilogit(eta),
                nd = data.frame(x = 1, g = factor("new"), N = 8L),
                onew = function(n) stats::rbinom(n, 8, ilogit(
                  0.7 + stats::rnorm(n, 0, 0.6)))))
  }
  if (kind == "xg") {
    d <- expand.grid(a = factor(1:6), b = factor(1:8), r = 1:2)
    d$x <- stats::rnorm(nrow(d))
    ua <- stats::rnorm(6, 0, 0.6)
    ub <- stats::rnorm(8, 0, 0.3)
    mu <- 1 + 0.5 * d$x + ua[d$a] + ub[d$b]
    d$y <- mu + stats::rnorm(nrow(d))
    fit <- fq(frm(bf(y ~ x + (1 | a) + (1 | b)), data = d))
    return(list(fit = fit, mu = mu,
                nd = data.frame(x = 1, a = factor("new"), b = factor("new")),
                onew = function(n) stats::rnorm(n, 1.5,
                                                sqrt(0.36 + 0.09 + 1))))
  }
  if (kind == "xp") {
    d <- expand.grid(a = factor(1:5), b = factor(1:12))
    d$x <- stats::rnorm(nrow(d))
    ua <- stats::rnorm(5, 0, 0.4)
    ub <- stats::rnorm(12, 0, 0.3)
    mu <- exp(1 + 0.3 * d$x + ua[d$a] + ub[d$b])
    d$y <- stats::rpois(nrow(d), mu)
    fit <- fq(frm(bf(y ~ x + (1 | a) + (1 | b)), family = poisson(),
                  data = d))
    return(list(fit = fit, mu = mu,
                nd = data.frame(x = 0, a = factor("new"), b = factor("new")),
                onew = function(n) stats::rpois(n, exp(
                  1 + stats::rnorm(n, 0, 0.4) + stats::rnorm(n, 0, 0.3)))))
  }
  if (kind == "ordg") {
    G <- 8L
    m <- 12L
    th <- c(-1, 0.3, 1.5)
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u <- stats::rnorm(G, 0, 0.7)
    eta <- 0.6 * d$x + u[d$g]
    d$y <- cum_draw(eta, th)
    fit <- fq(frm(bf(y ~ x + (1 | g)), family = cumulative(), data = d))
    return(list(fit = fit, mu = cum_mean(eta, th),
                nd = data.frame(x = 1, g = factor("new")),
                onew = function(n) cum_draw(0.6 + stats::rnorm(n, 0, 0.7),
                                            th)))
  }
  if (kind == "ordf") {
    th <- c(-2.5, -1, 0, 1, 2.5)
    d <- data.frame(x = stats::rnorm(60))
    eta <- 0.8 * d$x
    d$y <- cum_draw(eta, th)
    fit <- fq(frm(bf(y ~ x), family = cumulative(), data = d))
    return(list(fit = fit, mu = cum_mean(eta, th), nd = NULL))
  }
  if (kind == "nbglm") {
    d <- data.frame(x = stats::rnorm(40))
    mu <- exp(1 + 0.4 * d$x)
    d$y <- stats::rnbinom(40, mu = mu, size = 1)
    fit <- fq(frm(bf(y ~ x), family = negbinomial(), data = d))
    return(list(fit = fit, mu = mu, nd = NULL))
  }
  if (kind == "nbglmm") {
    G <- 8L
    m <- 8L
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u <- stats::rnorm(G, 0, 0.5)
    mu <- exp(1 + 0.4 * d$x + u[d$g])
    d$y <- stats::rnbinom(nrow(d), mu = mu, size = 2)
    fit <- fq(frm(bf(y ~ x + (1 | g)), family = negbinomial(), data = d))
    return(list(fit = fit, mu = mu,
                nd = data.frame(x = 0, g = factor("new")),
                onew = function(n) stats::rnbinom(n, mu = exp(
                  1 + stats::rnorm(n, 0, 0.5)), size = 2)))
  }
  if (kind == "dsig") {
    G <- 10L
    m <- 6L
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u <- stats::rnorm(G, 0, 0.5)
    mu <- 1 + 0.5 * d$x + u[d$g]
    d$y <- mu + stats::rnorm(nrow(d), 0, exp(0.3 * d$x))
    fit <- fq(frm(bf(y ~ x + (1 | g), sigma ~ x), data = d))
    return(list(fit = fit, mu = mu,
                nd = data.frame(x = 1, g = factor("new")),
                onew = function(n) stats::rnorm(n, 1.5,
                                                sqrt(0.25 + exp(0.6)))))
  }
  if (kind == "dphi") {
    d <- data.frame(x = stats::rnorm(60))
    mu <- ilogit(0.2 + 0.5 * d$x)
    phi <- exp(2 + 0.5 * d$x)
    d$y <- stats::rbeta(60, mu * phi, (1 - mu) * phi)
    fit <- fq(frm(bf(y ~ x, phi ~ x), family = Beta(), data = d))
    return(list(fit = fit, mu = mu, nd = NULL))
  }
  if (kind == "prs") {
    G <- 10L
    m <- 8L
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u0 <- stats::rnorm(G, 0, 0.8)
    u1 <- stats::rnorm(G, 0, 0.3)
    mu <- 1 + 0.5 * d$x + u0[d$g] + u1[d$g] * d$x
    d$y <- mu + stats::rnorm(nrow(d))
    pr <- prior(normal(0, 2), class = b) + prior(normal(0, 1), class = sd) +
      prior(lkj(2), class = cor)
    fit <- fq(frm(bf(y ~ x + (1 + x | g)), data = d, prior = pr))
    return(list(fit = fit, mu = mu,
                nd = data.frame(x = 2, g = factor("new")),
                onew = function(n) stats::rnorm(n, 2,
                                                sqrt(0.64 + 4 * 0.09 + 1))))
  }
  if (kind == "prb") {
    d <- data.frame(x1 = stats::rnorm(40), x2 = stats::rnorm(40))
    eta <- -0.5 + 2.5 * d$x1 + 0.5 * d$x2
    d$y <- stats::rbinom(40, 1, ilogit(eta))
    fit <- fq(frm(bf(y ~ x1 + x2), family = bernoulli(), data = d,
                  prior = prior(normal(0, 1.5), class = b)))
    return(list(fit = fit, mu = ilogit(eta), nd = NULL))
  }
  if (kind == "sigid") {
    d <- data.frame(x = stats::runif(40))
    mu <- 1 + 0.5 * d$x
    d$y <- mu + stats::rnorm(40, 0, 0.6 + 0.5 * d$x)
    fit <- fq(frm(bf(y ~ x, sigma ~ x), data = d,
                  family = brmsfamily("gaussian", link_sigma = "identity")))
    return(list(fit = fit, mu = mu, nd = NULL))
  }
  if (kind %in% c("slp", "slpw")) {
    # sleepstudy-like: 18 subjects, days 0 to 9, lme4's estimates
    # (slope sd 5.9, or 2 for slpw)
    d <- expand.grid(Days = 0:9, Subject = factor(1:18))
    ss <- if (kind == "slp") 5.92 else 2
    S <- matrix(c(24.74^2, 0.066 * 24.74 * ss, 0.066 * 24.74 * ss, ss^2), 2)
    U <- matrix(stats::rnorm(36), 18) %*% chol(S)
    mu <- 251.4 + 10.47 * d$Days + U[d$Subject, 1] + U[d$Subject, 2] * d$Days
    d$Reaction <- mu + stats::rnorm(nrow(d), 0, 25.59)
    fit <- fq(frm(bf(Reaction ~ Days + (Days | Subject)), data = d))
    v <- c(1, 9)
    return(list(fit = fit, mu = mu,
                nd = data.frame(Days = 9, Subject = factor("new")),
                onew = function(n) stats::rnorm(n, 251.4 + 94.23, sqrt(
                  drop(v %*% S %*% v) + 25.59^2))))
  }
  if (kind %in% c("prsp", "brs")) {
    # random slopes outside the gaussian family
    G <- if (kind == "prsp") 10L else 12L
    m <- if (kind == "prsp") 10L else 15L
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u0 <- stats::rnorm(G, 0, 0.5)
    u1 <- stats::rnorm(G, 0, if (kind == "prsp") 0.2 else 0.3)
    eta <- (if (kind == "prsp") 1 else -0.2) + 0.4 * d$x + u0[d$g] +
      u1[d$g] * d$x
    if (kind == "prsp") {
      mu <- exp(eta)
      d$y <- stats::rpois(nrow(d), mu)
      fam <- poisson()
    } else {
      mu <- ilogit(eta)
      d$y <- stats::rbinom(nrow(d), 1, mu)
      fam <- bernoulli()
    }
    fit <- fq(frm(bf(y ~ x + (1 + x | g)), family = fam, data = d))
    return(list(fit = fit, mu = mu, nd = NULL))
  }
  if (kind == "nbpois") {
    # negbinomial on poisson data: the shape runs toward its end
    d <- data.frame(x = stats::rnorm(100))
    mu <- exp(1 + 0.4 * d$x)
    d$y <- stats::rpois(100, mu)
    fit <- fq(frm(bf(y ~ x), family = negbinomial(), data = d))
    return(list(fit = fit, mu = mu, nd = NULL))
  }
  if (kind == "xg3") {
    d <- expand.grid(a = factor(1:3), b = factor(1:4), r = 1:4)
    d$x <- stats::rnorm(nrow(d))
    ua <- stats::rnorm(3, 0, 0.5)
    ub <- stats::rnorm(4, 0, 0.5)
    mu <- 1 + 0.5 * d$x + ua[d$a] + ub[d$b]
    d$y <- mu + stats::rnorm(nrow(d))
    fit <- fq(frm(bf(y ~ x + (1 | a) + (1 | b)), data = d))
    return(list(fit = fit, mu = mu,
                nd = data.frame(x = 1, a = factor("new"), b = factor("new")),
                onew = function(n) stats::rnorm(n, 1.5, sqrt(0.5 + 1))))
  }
  if (kind == "ordrs") {
    G <- 10L
    m <- 15L
    th <- c(-1, 0.3, 1.5)
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    u0 <- stats::rnorm(G, 0, 0.6)
    u1 <- stats::rnorm(G, 0, 0.3)
    eta <- 0.6 * d$x + u0[d$g] + u1[d$g] * d$x
    d$y <- cum_draw(eta, th)
    fit <- fq(frm(bf(y ~ x + (1 + x | g)), family = cumulative(), data = d))
    return(list(fit = fit, mu = cum_mean(eta, th), nd = NULL))
  }
  if (kind == "slp10") {
    d <- expand.grid(Days = 0:9, Subject = factor(1:10))
    S <- matrix(c(24.74^2, 0.066 * 24.74 * 5.92, 0.066 * 24.74 * 5.92,
                  5.92^2), 2)
    U <- matrix(stats::rnorm(20), 10) %*% chol(S)
    mu <- 251.4 + 10.47 * d$Days + U[d$Subject, 1] + U[d$Subject, 2] * d$Days
    d$Reaction <- mu + stats::rnorm(nrow(d), 0, 25.59)
    fit <- fq(frm(bf(Reaction ~ Days + (Days | Subject)), data = d))
    v <- c(1, 9)
    return(list(fit = fit, mu = mu,
                nd = data.frame(Days = 9, Subject = factor("new")),
                onew = function(n) stats::rnorm(n, 251.4 + 94.23, sqrt(
                  drop(v %*% S %*% v) + 25.59^2))))
  }
  if (kind == "sleep") {
    d <- lme4::sleepstudy
    fit <- fq(frm(bf(Reaction ~ Days + (Days | Subject)), data = d))
    return(list(fit = fit, mu = NULL,
                nd = data.frame(Days = 9, Subject = factor("new")),
                onew = NULL))
  }
  stop("unknown design ", kind)
}

seeds_of <- list(bfew4 = 1:15, bfew6 = 1:15, binfew = 1:10, xg = 1:15,
                 xp = 1:10, ordg = 1:15, ordf = 1:10, nbglm = 1:10,
                 nbglmm = 1:10, dsig = 1:10, dphi = 1:5, prs = 1:15,
                 prb = 1:10, sigid = 1:5, slp = 1:30, slpw = 1:20,
                 sleep = 1L, prsp = 1:15, brs = 1:15, nbpois = 1:10,
                 xg3 = 1:10, ordrs = 1:10, slp10 = 1:20)

said_of <- function(expr) {
  said <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      k <- if (grepl("Pareto k is", msg)) {
        sub(".*Pareto k is ([-0-9.e]+).*", "k=\\1", msg)
      } else ""
      said <<- c(said, paste0(class(w)[1], k))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      msg <- conditionMessage(m)
      k <- if (grepl("Pareto k [0-9]", msg)) {
        sub(".*Pareto k ([-0-9.e]+).*", "k=\\1", msg)
      } else ""
      said <<- c(said, paste0(class(m)[1], k))
      invokeRestart("muffleMessage")
    })
  list(val = val, said = paste(unique(said), collapse = ","))
}

unfire <- function(fit) {
  chk <- fit$cache[["draw_check"]]
  chk$fire <- FALSE
  assign("draw_check", chk, envir = fit$cache)
}

arm <- function(r, nd) {
  set.seed(1)
  pin <- fq(predict(r$fit, ndraws = 2000, summary = FALSE))
  mae <- if (!is.null(r$mu)) mean(abs(colMeans(pin) - r$mu)) else NA
  cov <- c(NA, NA, NA, NA)
  if (!is.null(nd)) {
    set.seed(1)
    pn <- fq(predict(r$fit, newdata = nd, allow_new_levels = TRUE,
                     ndraws = 2000, summary = FALSE))[, 1]
    if (!is.null(r$onew)) {
      set.seed(99)
      o <- r$onew(2e5)
      q90 <- stats::quantile(pn, c(0.05, 0.95), names = FALSE)
      q98 <- stats::quantile(pn, c(0.01, 0.99), names = FALSE)
      cov <- c(mean(o >= q90[1] & o <= q90[2]), q90[2] - q90[1],
               mean(o >= q98[1] & o <= q98[2]), q98[2] - q98[1])
    }
  }
  c(mae = mae, cov90 = cov[1], w90 = cov[2], cov98 = cov[3], w98 = cov[4])
}

rows <- list()
for (g in groups) {
  for (s in seeds_of[[g]]) {
    id <- paste0(g, ":", s)
    r <- tryCatch(design(g, s), error = function(e) {
      cat(id, "design error:", conditionMessage(e), "\n")
      NULL
    })
    if (is.null(r) || is.null(r$fit)) next
    fit <- r$fit
    t0 <- proc.time()[[3]]
    pr <- tryCatch(said_of({
      set.seed(5)
      predict(fit, ndraws = 200, summary = FALSE)
    }), error = function(e) list(val = NULL, said = paste("ERROR",
                                                           conditionMessage(e))))
    t1 <- proc.time()[[3]]
    chk <- fit$cache[["draw_check"]]
    gv <- function(x) if (is.null(x)) NA else x
    fired_dir <- isTRUE(chk$fire) && !is.na(gv(chk$dir))
    m_lane <- m_norm <- rep(NA, 5)
    if (fired_dir) {
      m_lane <- tryCatch(arm(r, r$nd), error = function(e) rep(NA, 5))
      unfire(fit)
      m_norm <- tryCatch(arm(r, r$nd), error = function(e) rep(NA, 5))
    }
    row <- data.frame(id = id, np = length(fit$opt$par),
                      verdict = !is.null(chk), fire = isTRUE(chk$fire),
                      excess = gv(chk$excess), rel = gv(chk$rel),
                      dir = gv(chk$dir), dir_par = gv(chk$dir_par),
                      dir_n = gv(chk$dir_n), t_first = t1 - t0,
                      said = pr$said,
                      mae_lane = m_lane[1], mae_norm = m_norm[1],
                      cov90_lane = m_lane[2], cov90_norm = m_norm[2],
                      w90_lane = m_lane[3], w90_norm = m_norm[3],
                      cov98_lane = m_lane[4], cov98_norm = m_norm[4],
                      w98_lane = m_lane[5], w98_norm = m_norm[5])
    rows[[length(rows) + 1L]] <- row
    cat(sprintf(paste0("%-10s np %2d fire %-5s ex %7.3f rel %6.3f | dir %8.3f",
                       " (%s, %s) | t %.2f | %s"),
                id, row$np, row$fire, row$excess, row$rel, row$dir,
                row$dir_par, row$dir_n, row$t_first, row$said))
    if (fired_dir) {
      cat(sprintf(paste0(" || mae lane %.4f norm %.4f | cov90 %.3f %.3f",
                         " w90 %.3f %.3f | cov98 %.3f %.3f w98 %.3f %.3f"),
                  m_lane[1], m_norm[1], m_lane[2], m_norm[2], m_lane[3],
                  m_norm[3], m_lane[4], m_norm[4], m_lane[5], m_norm[5]))
    }
    cat("\n")
  }
}
utils::write.table(do.call(rbind, rows), out, sep = "\t", row.names = FALSE,
                   quote = FALSE)
cat("DONE", length(rows), "\n")
