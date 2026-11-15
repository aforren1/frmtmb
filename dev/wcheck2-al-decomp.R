# Lane wcheck2, backlog 4949: where the censored rows' cost goes in
# asym_laplace's Laplace objective. The construction of
# dev/reviews/wcheck-review-misc.R part cens (seed 20261111, n = 80):
# 36 rows left censored, 36 right, 8 uncensored, weights 1,
# yac | cens(ca) ~ 1 + (1 | obs).
#
# 1. The TMB objective's profile along log sigma (the rest re-optimized),
#    as the review measured it.
# 2. At each profile point, the objective decomposed row by row from the
#    inner modes TMB found (obj$env$last.par): each row's log
#    probability, its curvature, and where its mode sits relative to its
#    bound in units of sigma.
# 3. The same Laplace approximation with each row's inner mode found by
#    a careful 1-D maximization ("exact inner mode"), and the exact
#    marginal (adaptive quadrature over each row's effect), at the same
#    outer parameters, and both profiled along log sigma themselves.
#   Rscript dev/wcheck2-al-decomp.R <base|lane|lib> [both|right72]
# "right72" is instead the construction of dev/wcheck-cens-path.R (data
# of dev/wcheck-cases.R part cens, seed 20261021): 72 of 80 rows right
# censored, the mean held 2 below the bound and the quantile at 0.5.
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
construction <- if (length(a) > 1L) a[2] else "both"
n <- 80
if (construction == "both") {
  set.seed(20261111)
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                  g = factor(rep(1:16, each = 5)), t = rep(1:5, 16))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  o <- order(d$y)
  d$ca <- "none"
  d$ca[o[1:36]] <- "left"
  d$ca[o[45:80]] <- "right"
  lo <- d$y[o[36]]
  hi <- d$y[o[45]]
  d$yac <- pmin(pmax(d$y, lo), hi)
  hold <- NULL
} else {
  set.seed(20261021)
  d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  d$ya <- 1 + 0.5 * d$x + rnorm(n)
  o <- order(d$ya, decreasing = TRUE)
  d$ca <- "none"
  d$ca[o[seq_len(72)]] <- "right"
  hi <- d$ya[o[72]]
  lo <- hi
  d$yac <- pmin(d$ya, hi)
  # the mean 2 below the bound, the quantile at the median
  hold <- c(m = hi - 2, lq = 0)
}
cat(sprintf("left bound %.4f, right bound %.4f, %d left, %d right\n", lo, hi,
            sum(d$ca == "left"), sum(d$ca == "right")))
fa <- bf(yac | cens(ca) ~ 1 + (1 | obs))
ob <- q(frm(fa, family = asym_laplace(), data = d, dry_run = "objective"))
obj <- if (!is.null(ob$obj)) ob$obj else ob
cat("par names:", paste(unique(names(obj$par)), collapse = " "), "\n")
cen <- ifelse(d$ca == "left", -1, ifelse(d$ca == "right", 1, 0))
yv <- d$yac

# One row's log probability as a function of its mean (doubles, careful
# logs), and its second derivative in the mean.
row_lp <- function(mu, yi, ci, s, p) {
  u <- (yi - mu) / s
  if (ci == 0) {
    return(log(p * (1 - p) / s) - u * (p - (u < 0)))
  }
  if (ci == 1) {   # P(Y > yi)
    if (u > 0) return(log(1 - p) - p * u)
    return(log1p(-p * exp((1 - p) * u)))
  }
  # P(Y < yi)
  if (u < 0) return(log(p) + (1 - p) * u)
  log1p(-(1 - p) * exp(-p * u))
}
row_curv <- function(mu, yi, ci, s, p) {
  u <- (yi - mu) / s
  if (ci == 0) return(0)
  if (ci == 1) {
    if (u > 0) return(0)
    e <- p * exp((1 - p) * u)
    return((1 - p)^2 * e / ((1 - e)^2 * s^2))
  }
  if (u < 0) return(0)
  e <- (1 - p) * exp(-p * u)
  p^2 * e / ((1 - e)^2 * s^2)
}
# The Laplace term of one row at inner mode b: log p(y | b) + log N(b; 0,
# sd) + log(2 pi) / 2 - log(H) / 2, H = curvature + 1 / sd^2.
row_laplace <- function(b, m, yi, ci, s, p, sd) {
  row_lp(m + b, yi, ci, s, p) + dnorm(b, 0, sd, log = TRUE) +
    0.5 * log(2 * pi) - 0.5 * log(row_curv(m + b, yi, ci, s, p) + 1 / sd^2)
}
# The row's exact inner mode: the maximum of log p(y | b) + log N(b; 0,
# sd). The function is concave; optimize() over a bracket that holds it.
row_mode <- function(m, yi, ci, s, p, sd) {
  f <- function(b) row_lp(m + b, yi, ci, s, p) + dnorm(b, 0, sd, log = TRUE)
  if (ci == 0) return(yi - m)
  # the mode is within a few hundred sigma of the bound, or at 0
  ctr <- yi - m
  br <- sort(c(0, ctr - 2000 * s * sign(ci), ctr + 2000 * s * sign(ci)))
  br <- c(min(br) - 1, max(br) + 1)
  optimize(f, br, maximum = TRUE, tol = 1e-15 * max(1, abs(ctr)))$maximum
}
row_exact <- function(m, yi, ci, s, p, sd) {
  f <- function(b) {
    vapply(b, function(bb) row_lp(m + bb, yi, ci, s, p), 0) +
      dnorm(b, 0, sd, log = TRUE)
  }
  bm <- row_mode(m, yi, ci, s, p, sd)
  fm <- f(bm)
  g <- function(b) exp(f(b) - fm)
  # the integrand's scale: sigma near the kink, sd far from it
  pts <- sort(unique(c(-Inf, bm - 10 * sd, bm - 50 * s, bm, bm + 50 * s,
                       bm + 10 * sd, Inf)))
  tot <- 0
  for (k in seq_len(length(pts) - 1L)) {
    r <- tryCatch(integrate(g, pts[k], pts[k + 1L], rel.tol = 1e-10,
                            subdivisions = 2000L)$value,
                  error = function(e) NA_real_)
    tot <- tot + r
  }
  fm + log(tot)
}
unpack <- function(par) {
  nm <- names(par)
  list(m = par[nm == "beta"][1L], ls = par[nm == "betad"][1L],
       lq = par[nm == "betad"][2L], lsd = par[nm == "theta"][1L])
}
# -log of the three approximations at outer parameters `par`
three <- function(par) {
  u <- unpack(par)
  s <- exp(u$ls)
  p <- plogis(u$lq)
  sd <- exp(u$lsd)
  tmb <- obj$fn(par)
  lp <- obj$env$last.par
  bh <- lp[obj$env$random]
  rows <- vapply(seq_len(n), function(i) {
    bt <- bh[i]
    be <- row_mode(u$m, yv[i], cen[i], s, p, sd)
    c(tmb_term = row_laplace(bt, u$m, yv[i], cen[i], s, p, sd),
      tmb_curv_sd2 = row_curv(u$m + bt, yv[i], cen[i], s, p) * sd^2,
      tmb_u = (yv[i] - u$m - bt) / s,
      ex_term = row_laplace(be, u$m, yv[i], cen[i], s, p, sd),
      ex_curv_sd2 = row_curv(u$m + be, yv[i], cen[i], s, p) * sd^2,
      ex_u = (yv[i] - u$m - be) / s,
      exact = row_exact(u$m, yv[i], cen[i], s, p, sd))
  }, c(tmb_term = 0, tmb_curv_sd2 = 0, tmb_u = 0, ex_term = 0,
       ex_curv_sd2 = 0, ex_u = 0, exact = 0))
  list(tmb = -tmb, tmb_sum = sum(rows["tmb_term", ]),
       lap = sum(rows["ex_term", ]), exact = sum(rows["exact", ]),
       rows = rows)
}

gs <- -c(2, 4, 6, 8, 10, 12, 14)
pp <- obj$par
isd <- which(names(pp) == "betad")[1L]
held <- if (is.null(hold)) integer() else
  c(which(names(pp) == "beta")[1L], which(names(pp) == "betad")[2L])
if (length(held)) pp[held] <- unname(hold)
free <- setdiff(seq_along(pp), c(isd, held))
prof <- list()
for (t in gs) {
  pp[isd] <- t
  fr <- function(r) {
    q2 <- pp
    q2[free] <- r
    v <- tryCatch(obj$fn(q2), error = function(e) NA_real_)
    if (!is.finite(v)) 1e10 else v
  }
  op <- nlminb(pp[free], fr, control = list(iter.max = 500, eval.max = 800))
  if (op$objective < 1e9) pp[free] <- op$par
  th <- three(pp)
  u <- unpack(pp)
  rw <- th$rows
  cat(sprintf(paste0("log sigma %g: TMB %.4f (row sum %.4f) | exact inner ",
                     "mode %.4f | exact marginal %.4f | m %.4f, quantile ",
                     "%.4f, log sd %.3f\n"),
              t, th$tmb, th$tmb_sum, th$lap, th$exact, u$m, plogis(u$lq),
              u$lsd))
  for (ci in c(-1, 1)) {
    k <- cen == ci
    cat(sprintf(paste0("    %s rows: TMB mode u/sigma in [%.4g, %.4g], ",
                       "curv*sd^2 max %.3g, term sum %.4f | exact mode u ",
                       "in [%.4g, %.4g], curv*sd^2 in [%.3g, %.3g], term ",
                       "sum %.4f\n"),
                if (ci < 0) "left " else "right",
                min(rw["tmb_u", k]), max(rw["tmb_u", k]),
                max(rw["tmb_curv_sd2", k]), sum(rw["tmb_term", k]),
                min(rw["ex_u", k]), max(rw["ex_u", k]),
                min(rw["ex_curv_sd2", k]), max(rw["ex_curv_sd2", k]),
                sum(rw["ex_term", k])))
  }
  k <- cen == 0
  cat(sprintf("    uncensored: TMB term sum %.4f, exact-mode term sum %.4f\n",
              sum(rw["tmb_term", k]), sum(rw["ex_term", k])))
  prof[[length(prof) + 1L]] <- c(t = t, tmb = th$tmb, lap = th$lap,
                                 exact = th$exact)
}
pm <- do.call(rbind, prof)
cat("\nAlong the TMB profile path, per unit of -log sigma:\n")
for (col in c("tmb", "lap", "exact")) {
  cat(sprintf("  %-5s slope %s\n", col,
              paste(format(diff(pm[, col]) / diff(-pm[, "t"]), digits = 3),
                    collapse = " ")))
}

# The Laplace approximation with exact inner modes, and the exact
# marginal, each profiled along log sigma on its own (m, quantile and
# log sd re-optimized).
own_profile <- function(which) {
  st <- c(m = 0, lq = 0, lsd = 0)
  vapply(gs, function(t) {
    f <- function(v) {
      if (!is.null(hold)) v <- c(hold[["m"]], hold[["lq"]], v[length(v)])
      s <- exp(t); p <- plogis(v[2]); sd <- exp(v[3])
      tot <- sum(vapply(seq_len(n), function(i) {
        if (which == "lap") {
          b <- row_mode(v[1], yv[i], cen[i], s, p, sd)
          row_laplace(b, v[1], yv[i], cen[i], s, p, sd)
        } else row_exact(v[1], yv[i], cen[i], s, p, sd)
      }, 0))
      if (!is.finite(tot)) 1e10 else -tot
    }
    best <- NULL
    starts <- list(st, c(hi, 0, 0), c(lo, 0, 0), c((lo + hi) / 2, 0, 3),
                   c(hi, 0, 5), c((lo + hi) / 2, 0, -2))
    # with the mean and quantile held, only the effects' log sd is free
    if (!is.null(hold)) starts <- lapply(starts, function(s0) s0[3L])
    for (st0 in starts) {
      op <- tryCatch(nlminb(st0, f), error = function(e) NULL)
      if (!is.null(op) && (is.null(best) || op$objective < best$objective)) {
        best <- op
      }
    }
    st <<- best$par
    -best$objective
  }, 0)
}
for (which in c("lap", "exact")) {
  v <- own_profile(which)
  show_profile(sprintf("own profile, %s", if (which == "lap")
    "Laplace with exact inner modes" else "exact marginal"), -gs, v)
}
cat("DONE\n")
