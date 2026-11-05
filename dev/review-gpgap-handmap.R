# Review of lane gpgap, item 2: the order claim, by hand. Nothing here calls
# brms_gp_order() or stan_pars_from_fit(). brms's side is read off
# make_standata() (Xgp, Jgp, Igp, Cgp, dmax) and the marginal is written out
# in R; the joint maps frmtmb's b to brms's positions by matching raw
# coordinates. Seed 7301.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-handmap.R handmap.txt [case]
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
H <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-stan.R")) {
  sys.source(file.path("tests/testthat", h), envir = H)
}
only <- commandArgs(trailingOnly = TRUE)

set.seed(7301)
n <- 40
base <- data.frame(
  x = sample(seq(0, 8, by = 0.25), n, TRUE),
  z = sample(seq(0, 8, by = 0.5), n, TRUE),
  f = factor(sample(c("a", "b", "c"), n, TRUE)),
  w = sample(c(0, 0.5, 1.5, 2), n, TRUE))
base$xu <- round(base$x + runif(n, -0.1, 0.1), 3)
base$y <- 2 * sin(0.8 * base$x) + 0.5 * cos(base$z) + rnorm(n, 0, 0.3)
small <- data.frame(x = c(5.1, 0.3, 7.7, 2.2, 3.9, 0.9, 6.4, 1.6, 4.8, 8.8,
                          2.7, 6.0))
small$y <- sin(small$x) + c(0.1, -0.2, 0.05, 0.15, -0.3, 0.2, -0.05, 0.1,
                            -0.1, 0.25, 0, -0.15)
d_nax <- base; d_nax$x[c(1, 17)] <- NA
d_nay <- base; d_nay$y[c(2, 9)] <- NA
d_ulp <- base
d_ulp$x[2] <- d_ulp$x[1] * (1 + 4 * .Machine$double.eps)

cases <- list(
  small_unsorted = list(y ~ gp(x), small),
  dup_x = list(y ~ gp(x), base),
  gr_false_distinct = list(y ~ gp(xu, gr = FALSE), base),
  gr_false_dup = list(y ~ gp(x, gr = FALSE), base),
  by_factor = list(y ~ gp(x, by = f), base),
  by_factor_cmc_false = list(y ~ gp(x, by = f, cmc = FALSE), base),
  by_numeric = list(y ~ gp(x, by = w), base),
  twod_iso = list(y ~ gp(x, z), base),
  twod_aniso = list(y ~ gp(x, z, iso = FALSE), base),
  na_x = list(y ~ gp(x), d_nax),
  na_y = list(y ~ gp(x), d_nay),
  scale_false = list(y ~ gp(x, scale = FALSE), base),
  ulp_dup = list(y ~ gp(x), d_ulp),
  two_terms = list(y ~ gp(x) + gp(z), base)
)
if (length(only)) cases <- cases[only]

# All sub-GP data of term i off the standata, with no help from frmtmb.
term_parts <- function(sdat, i) {
  pi <- paste0("_", i)
  K <- sdat[[paste0("Kgp", pi)]]
  byf <- !is.null(sdat[[paste0("Igp", pi, "_1")]])
  lapply(seq_len(K), function(j) {
    sfx <- if (byf) paste0(pi, "_", j) else pi
    list(sfx = sfx,
         X = as.matrix(sdat[[paste0("Xgp", sfx)]]),
         J = sdat[[paste0("Jgp", sfx)]],
         I = sdat[[paste0("Igp", sfx)]],
         C = if (byf) sdat[[paste0("Cgp", sfx)]] else sdat[[paste0("Cgp", pi)]],
         dmax = sdat[[paste0("dmax", sfx)]] %||% 1)
  })
}

kern <- function(X, sd, rho, nugget) {
  D <- ncol(X)
  rho <- rep_len(rho, D)
  Q <- 0
  for (k in seq_len(D)) Q <- Q + outer(X[, k], X[, k], "-")^2 / rho[k]^2
  K <- sd^2 * exp(-Q / 2)
  K + diag(nugget(sd), nrow(K))
}

# The marginal covariance of y in brms's parameterization.
marg_cov <- function(sdat, N, blocks, theta, nugget) {
  V <- matrix(0, N, N)
  nt <- sum(grepl("^Kgp_[0-9]+$", names(sdat)))
  h <- 0L
  for (i in seq_len(nt)) {
    for (pt in term_parts(sdat, i)) {
      h <- h + 1L
      th <- theta[blocks[[h]]$theta_idx]
      sd <- exp(th[[1]])
      rho <- exp(th[-1]) / pt$dmax
      K <- kern(pt$X, sd, rho, nugget)
      rows <- if (is.null(pt$I)) seq_len(N) else as.integer(pt$I)
      J <- if (is.null(pt$J)) seq_along(rows) else as.integer(pt$J)
      Cv <- if (is.null(pt$C)) rep(1, length(rows)) else as.numeric(pt$C)
      V[rows, rows] <- V[rows, rows] + outer(Cv, Cv) * K[J, J, drop = FALSE]
    }
  }
  V
}

ll_gauss <- function(y, mu, V) {
  R <- chol(V)
  z <- backsolve(R, y - mu, transpose = TRUE)
  -0.5 * sum(z^2) - sum(log(diag(R))) - 0.5 * length(y) * log(2 * pi)
}

rel_nug <- function(sd) 1e-6 * sd^2
brms_nug <- function(sd) 1e-12

for (nm in names(cases)) {
  fo <- cases[[nm]][[1]]
  d <- cases[[nm]][[2]]
  cat("\n=== ", nm, ": ", deparse(fo), "\n", sep = "")
  r <- tryCatch({
    bfrm <- do.call(frmtmb::bf, list(fo))
    fit <- suppressWarnings(frm(bfrm, family = gaussian(), data = d))
    bb <- brms::bf(fo)
    pr <- H$brms_flat_prior(bb, data = d, family = gaussian())
    sdat <- suppressMessages(brms::make_standata(bb, data = d, family = gaussian(),
                                                 prior = pr))
    code <- brms::make_stancode(bb, data = d, family = gaussian(), prior = pr)
    bks <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)
    keep <- stats::complete.cases(d[, intersect(all.vars(fo), names(d)), drop = FALSE])
    y <- as.numeric(sdat$Y)
    stopifnot(length(y) == sum(keep))
    N <- length(y)
    b0 <- fit$estimates$beta[[1]]
    sg <- as.numeric(sigma(fit))
    th <- fit$estimates$theta
    Sy <- diag(sg^2, N)
    llf <- ll_gauss(y, b0, marg_cov(sdat, N, bks, th, rel_nug) + Sy)
    llb <- ll_gauss(y, b0, marg_cov(sdat, N, bks, th, brms_nug) + Sy)
    # two points off the optimum, through the tape
    par0 <- fit$opt$par
    pts <- lapply(c(0.05, -0.07), function(e) {
      p <- par0 + e
      list(tape = -fit$obj$fn(p),
           hand = {
             pn <- names(p)
             beta <- p[pn == "beta"]
             theta <- p[pn == "theta"]
             lsg <- p[pn == "sigma" | grepl("sigma", pn)]
             NA
           })
    })
    cat(sprintf(paste0("  N %d, positions frm %s, brms %s, sdgp %.4f, ",
                       "sigma %.4f, conv %d\n"),
                N, paste(vapply(bks, function(b) nrow(b$aux_D2[[1]]), 0L),
                         collapse = ","),
                paste(vapply(unlist(lapply(
                  seq_len(sum(grepl("^Kgp_[0-9]+$", names(sdat)))),
                  term_parts, sdat = sdat), recursive = FALSE),
                  function(p) nrow(p$X), 0L), collapse = ","),
                exp(th[[bks[[1]]$theta_idx[1]]]), sg, fit$opt$convergence))
    cat(sprintf(paste0("  marginal: frm logLik %.8f | hand at frm nugget %.8f ",
                       "(diff %.2e) | hand at brms nugget %.8f (diff %.2e)\n"),
                as.numeric(logLik(fit)), llf, llf - as.numeric(logLik(fit)),
                llb, llb - as.numeric(logLik(fit))))

    # fitted at the training rows against the closed-form conditional mean
    V <- marg_cov(sdat, N, bks, th, rel_nug) + Sy
    Knn <- V - Sy
    # the nugget's own share is not the smooth field: both readings
    Knn0 <- marg_cov(sdat, N, bks, th, function(sd) 0)
    alpha <- solve(V, y - b0)
    cond_with <- b0 + as.numeric(Knn %*% alpha)
    cond_wo <- b0 + as.numeric(Knn0 %*% alpha)
    fv <- as.numeric(fitted(fit)[, "Estimate"])
    cat(sprintf(paste0("  fitted at the rows: max |frm - E[f+nug|y]| %.2e, ",
                       "max |frm - E[f|y] no nugget| %.2e, range y %.2f; ",
                       "spearman with row order sanity %.4f\n"),
                max(abs(fv - cond_with)), max(abs(fv - cond_wo)),
                diff(range(y)), cor(fv, y, method = "spearman")))

    # the joint: frmtmb's b mapped by raw coordinates to brms's positions
    prg <- list(own = code,
                frm = gsub("cov[n, n] += 1e-12;",
                           "cov[n, n] += 1e-6 * square(sdgp);", code,
                           fixed = TRUE))
    nt <- sum(grepl("^Kgp_[0-9]+$", names(sdat)))
    pars <- list(Intercept = b0, sigma = sg)
    logdet <- 0
    h <- 0L
    for (i in seq_len(nt)) {
      pts_i <- term_parts(sdat, i)
      Ks <- length(pts_i)
      Dls <- if (length(bks[[h + 1L]]$theta_idx) == 2L) 1L else
        ncol(pts_i[[1]]$X)
      sdgp <- numeric(Ks)
      lsc <- matrix(0, Ks, Dls)
      for (j in seq_len(Ks)) {
        h <- h + 1L
        bk <- bks[[h]]
        pt <- pts_i[[j]]
        tb <- th[bk$theta_idx]
        sdgp[j] <- exp(tb[[1]])
        lsc[j, ] <- exp(tb[-1]) / pt$dmax
        # raw coordinates through the frame's own record of positions,
        # matched to brms's by value
        gi <- Filter(function(g) identical(g$block_id, bk$block_id %||% h),
                     fit$frame$linpreds[[1]]$gps)
        if (!length(gi)) gi <- list(fit$frame$linpreds[[1]]$gps[[h]])
        posf <- as.matrix(gi[[1]]$positions)
        raw <- pt$X * pt$dmax
        io <- vapply(seq_len(nrow(raw)), function(q) {
          which.min(rowSums(abs(sweep(posf, 2, raw[q, ]))))
        }, 1L)
        stopifnot(!anyDuplicated(io) || nm == "gr_false_dup")
        Lm <- t(chol(kern(pt$X, sdgp[j], lsc[j, ], rel_nug)))
        fb <- as.numeric(fit$estimates$b[bk$b_idx])[io]
        zv <- forwardsolve(Lm, fb)
        logdet <- logdet + sum(log(diag(Lm)))
        pars[[paste0("zgp", pt$sfx)]] <- array(zv, length(zv))
      }
      pars[[paste0("sdgp_", i)]] <- array(sdgp, Ks)
      pars[[paste0("lscale_", i)]] <- lsc
    }
    ours <- -fit$obj$env$f(fit$obj$env$last.par.best) + logdet
    vals <- vapply(prg, function(cd) {
      sf <- suppressMessages(rstan::sampling(H$brms_stan_model(cd), data = sdat,
                                             chains = 0))
      u <- rstan::unconstrain_pars(sf, pars)
      rstan::log_prob(sf, u, adjust_transform = FALSE) - ours
    }, 0)
    cat(sprintf(paste0("  joint by hand: brms's own line %.6f; brms's program ",
                       "at frmtmb's line %.3e (tol 1e-6*|ours| = %.1e)\n"),
                vals[["own"]], vals[["frm"]], 1e-6 * abs(ours)))
    # and the lane's translator against the hand map
    tp <- H$stan_pars_from_fit(fit, sdat, code)
    zn <- grep("^zgp_", names(tp), value = TRUE)
    cat(sprintf("  translator z vs hand z: max diff %.2e\n",
                max(vapply(zn, function(q) {
                  max(abs(as.numeric(tp[[q]]) - as.numeric(pars[[q]])))
                }, 0))))
    TRUE
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    FALSE
  })
}

# Readers: new positions, against the closed form, on frmtmb's own side, in a
# shuffled newdata.
cat("\n=== readers at new positions and in shuffled row order\n")
d <- base
fit <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
bk <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)[[1]]
tb <- fit$estimates$theta[bk$theta_idx]
sd <- exp(tb[[1]]); rho <- exp(tb[[2]])
b0 <- fit$estimates$beta[[1]]; sg <- as.numeric(sigma(fit))
u <- sort(unique(d$x))
Kuu <- sd^2 * exp(-outer(u, u, "-")^2 / (2 * rho^2)) + diag(1e-6 * sd^2, length(u))
J <- match(d$x, u)
V <- Kuu[J, J] + diag(sg^2, n)
alpha <- solve(V, d$y - b0)
newx <- c(1.13, 7.21, 0.5, 3.37, 8.0, -0.5, 9.0, d$x[5])
Ks <- sd^2 * exp(-outer(newx, d$x, "-")^2 / (2 * rho^2))
closed <- b0 + as.numeric(Ks %*% alpha)
fn <- as.numeric(fitted(fit, newdata = data.frame(x = newx))[, "Estimate"])
print(round(rbind(x = newx, frm = fn, closed = closed), 5))
cat(sprintf("  max |frm - closed form| at new positions: %.2e\n",
            max(abs(fn - closed))))
idx <- sample(n)
f1 <- as.numeric(fitted(fit)[, "Estimate"])
f2 <- as.numeric(fitted(fit, newdata = d[idx, ])[, "Estimate"])
cat(sprintf("  fitted(newdata = shuffled training rows) vs fitted()[shuffle]: %.2e\n",
            max(abs(f2 - f1[idx]))))
re <- tryCatch(ranef(fit), error = function(e) conditionMessage(e))
cat("  ranef names: ", if (is.list(re)) paste(names(re), collapse = ",") else re, "\n")
if (is.list(re) && length(re)) {
  r1 <- re[[1]]
  if (is.array(r1) || is.data.frame(r1)) {
    print(utils::head(rownames(r1), 4)); print(dim(r1))
    # the first ranef row is frmtmb's first SORTED position
    bv <- as.numeric(fit$estimates$b[bk$b_idx])
    print(round(c(first_b = bv[1], first_ranef = as.numeric(r1[1, 1]),
                  sorted_u1 = u[1], data_x1 = d$x[1]), 5))
  }
}
cat("  cnms: ", paste(utils::head(bk$cnms, 3), collapse = ", "), " ...\n")
cat("  coef names with gp: ", paste(grep("gp", names(fixef(fit)[, 1]),
                                         value = TRUE), collapse = ","), "\n")
