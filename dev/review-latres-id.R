# Attack 1: the identity, with no translator. Hand map (beta, ar, ma, sderr)
# to brms's Stan parameters at chosen points that are NOT frm()'s estimates,
# compute the Laplace value of brms's own program in zerr by Newton on
# rstan::grad_log_prob, and compare with frm()'s obj$fn at the same point.
# Also: Cov(err) of the block against brms's Stan functions by hand.
# Usage: Rscript review-latres-id.R <case> <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
case <- args[[1]]
outf <- args[[2]]
suppressPackageStartupMessages({library(rstan); library(brms)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
stan_model_of <- hel$brms_stan_model

# data shapes
mk <- function(shape) {
  set.seed(4242)
  d <- switch(shape,
    sorted = expand.grid(t = 1:4, g = factor(1:8)),
    unsorted = {
      x <- expand.grid(t = 1:4, g = factor(1:8))
      x[sample(nrow(x)), ]
    },
    ragged = {
      n <- c(1, 2, 3, 5, 1, 7, 2, 4)
      x <- do.call(rbind, lapply(seq_along(n), function(i) {
        data.frame(t = seq_len(n[i]), g = i)
      }))
      x$g <- factor(x$g)
      x[sample(nrow(x)), ]
    },
    gaps = {
      x <- do.call(rbind, lapply(1:6, function(i) {
        data.frame(t = sort(sample(1:12, 5)), g = i)
      }))
      x$g <- factor(x$g)
      x
    },
    single = data.frame(t = 1:20, g = factor(1:20)))
  rownames(d) <- NULL
  d$x <- rnorm(nrow(d))
  e <- rnorm(nrow(d), 0, 0.5)
  d$y <- rpois(nrow(d), exp(0.6 + 0.4 * d$x + e))
  d$bin <- rbinom(nrow(d), 1, plogis(0.3 + 0.4 * d$x + e))
  d$pos <- rgamma(nrow(d), 3, 3 / exp(0.6 + 0.4 * d$x + e))
  d
}

cases <- list(
  ar = list(term = "ar(t, g, cov = TRUE)", shape = "sorted", fam = poisson()),
  ma = list(term = "ma(t, g, cov = TRUE)", shape = "sorted", fam = poisson()),
  arma = list(term = "arma(t, g, cov = TRUE)", shape = "sorted", fam = poisson()),
  cosy = list(term = "cosy(t, g)", shape = "sorted", fam = poisson()),
  unstr = list(term = "unstr(t, g)", shape = "sorted", fam = poisson()),
  ar_uns = list(term = "ar(t, g, cov = TRUE)", shape = "unsorted", fam = poisson()),
  ma_uns = list(term = "ma(t, g, cov = TRUE)", shape = "unsorted", fam = poisson()),
  unstr_uns = list(term = "unstr(t, g)", shape = "unsorted", fam = poisson()),
  ar_rag = list(term = "ar(t, g, cov = TRUE)", shape = "ragged", fam = poisson()),
  arma_rag = list(term = "arma(t, g, cov = TRUE)", shape = "ragged", fam = poisson()),
  cosy_rag = list(term = "cosy(t, g)", shape = "ragged", fam = poisson()),
  ar_gap = list(term = "ar(t, g, cov = TRUE)", shape = "gaps", fam = poisson()),
  ar_single = list(term = "ar(t, g, cov = TRUE)", shape = "single", fam = poisson()),
  ar_notime = list(term = "ar(gr = g, cov = TRUE)", shape = "unsorted", fam = poisson()),
  ar_nogr = list(term = "ar(t, cov = TRUE)", shape = "single", fam = poisson()),
  cosy_nogr = list(term = "cosy(t)", shape = "single", fam = poisson()),
  ma_notime = list(term = "ma(gr = g, cov = TRUE)", shape = "ragged", fam = poisson()),
  ar_gamma = list(term = "ar(t, g, cov = TRUE)", shape = "sorted", fam = Gamma(link = "log"), y = "pos"),
  ar_bern = list(term = "ar(t, g, cov = TRUE)", shape = "ragged", fam = bernoulli(), y = "bin"),
  unstr_rag = list(term = "unstr(t, g)", shape = "ragged", fam = poisson()),
  ma_rag = list(term = "ma(t, g, cov = TRUE)", shape = "ragged", fam = poisson()),
  arma_uns = list(term = "arma(t, g, cov = TRUE)", shape = "unsorted", fam = poisson()),
  cosy_uns = list(term = "cosy(t, g)", shape = "unsorted", fam = poisson()),
  condar = list(term = "ar(t, g)", shape = "unsorted", fam = poisson()),
  condar2 = list(term = "ar(t, g, p = 2)", shape = "ragged", fam = poisson()),
  condar_gamma = list(term = "ar(t, g)", shape = "sorted", fam = Gamma(link = "log"), y = "pos"),
  condar_gamma_inv = list(term = "ar(t, g)", shape = "sorted", fam = Gamma(), y = "pos"),
  condar_probit = list(term = "ar(t, g)", shape = "ragged", fam = bernoulli("probit"), y = "bin"))
cs <- cases[[case]]
d <- mk(cs$shape)
y <- cs$y %||% "y"
form <- as.formula(paste(y, "~ x +", cs$term))
pri <- c(set_prior("", class = "sderr"), set_prior("", class = "Intercept"))
if (identical(cs$fam$family, "Gamma")) pri <- c(pri, set_prior("", class = "shape"))
if (grepl("^unstr", case)) pri <- c(pri, set_prior("", class = "cortime"))
if (grepl("^(ar|ma|cosy)", cs$term)) {
  # flat coefficients are the default on brms's side
}
code <- brms::stancode(form, d, family = cs$fam, prior = pri)
sdat <- brms::standata(form, d, family = cs$fam, prior = pri)
writeLines(code, sub("[.]tsv$", paste0("-", case, ".stan"), outf))
m <- stan_model_of(code)
sf <- rstan::sampling(m, data = sdat, chains = 0, refresh = 0)

# the chosen points: not fitted, different per index
pts <- list(
  list(b0 = 0.3, b = 0.2, s = 0.7, a = 0.5, m = -0.4, r = 0.3),
  list(b0 = 0.9, b = -0.3, s = 0.35, a = -0.6, m = 0.7, r = 0.6),
  list(b0 = -0.2, b = 0.5, s = 1.3, a = 0.85, m = 0.2, r = 0.15))
shape_par <- 2.5

pars_of <- function(p) {
  xm <- mean(d$x)
  fam <- cs$fam$family
  pl <- list(b = array(p$b, 1), Intercept = p$b0 + p$b * xm, sderr = p$s)
  N <- length(sdat$Y)
  pl$zerr <- rep(0, N)
  nm <- rstan::get_stancode(sf)
  if (grepl("vector<[^>]*>\\[Kar\\] ar", nm)) pl$ar <- array(if (case == "condar2" ) c(p$a, -0.2) else p$a, sdat$Kar)
  if (grepl("vector<[^>]*>\\[Kma\\] ma", nm)) pl$ma <- array(p$m, sdat$Kma)
  if (grepl("real<[^>]*> cosy", nm)) pl$cosy <- p$r
  if (grepl("Lcortime", nm) && grepl("cholesky_factor_corr", nm)) {
    C <- p$r^abs(outer(1:sdat$n_unique_t, 1:sdat$n_unique_t, "-"))
    pl$Lcortime <- t(chol(C))
  }
  if (fam == "Gamma") pl$shape <- shape_par
  pl
}

# which upars are zerr
pl0 <- pars_of(pts[[1]])
u0 <- rstan::unconstrain_pars(sf, pl0)
pl1 <- pl0
pl1$zerr <- rep(1, length(pl0$zerr))
iz <- which(rstan::unconstrain_pars(sf, pl1) != u0)
stopifnot(length(iz) == length(pl0$zerr))

laplace_brms <- function(p) {
  pl <- pars_of(p)
  u <- rstan::unconstrain_pars(sf, pl)
  n <- length(iz)
  z <- rep(0, n)
  lp <- function(z) {
    u[iz] <- z
    tryCatch(as.numeric(rstan::log_prob(sf, u, adjust_transform = FALSE)),
             error = function(e) -Inf)
  }
  gr <- function(z) {
    u[iz] <- z
    tryCatch(as.numeric(rstan::grad_log_prob(sf, u, adjust_transform = FALSE)[iz]),
             error = function(e) rep(NA_real_, n))
  }
  hs <- function(z, h = 1e-5) {
    H <- matrix(0, n, n)
    for (i in seq_len(n)) {
      e <- rep(0, n)
      e[i] <- h
      H[i, ] <- (gr(z + e) - gr(z - e)) / (2 * h)
    }
    (H + t(H)) / 2
  }
  for (it in 1:200) {
    g <- gr(z)
    if (max(abs(g)) < 1e-8) break
    H <- hs(z)
    stp <- solve(-H, g)
    f0 <- lp(z)
    a <- 1
    while (a > 1e-9 && (!is.finite(lp(z + a * stp)) || lp(z + a * stp) < f0)) a <- a / 2
    z <- z + a * stp
  }
  H <- hs(z)
  lp(z) + n / 2 * log(2 * pi) - 0.5 * as.numeric(determinant(-H)$modulus)
}

# frm() side: the same point in its charts
fit0 <- msg_of(frm(form, data = d, family = cs$fam))
if (!is.null(fit0$err)) {
  cat(case, "frm ERROR:", fit0$err, "\n")
  write.table(data.frame(case = case, pt = NA, brms = NA, frm = NA, diff = fit0$err),
              outf, sep = "\t", row.names = FALSE, quote = FALSE)
  quit("no")
}
fit <- fit0$val
cat(case, "frm warnings:", fit0$warn, "\n")
bk <- fi("frame_latres")(fit$frame)[[1]]
a <- bk$aux
th_of <- function(p) {
  inv_pacf <- function(r) r / sqrt(1 - r^2)
  q <- function(r) log((r + 1) / 2) - log(1 - (r + 1) / 2)
  switch(a$kind,
    ar1 = c(log(p$s), inv_pacf(p$a)),
    ma1 = c(log(p$s), inv_pacf(p$m)),
    arma1 = c(log(p$s), inv_pacf(p$a), inv_pacf(p$m)),
    condar = c(log(p$s), q(if (case == "condar2") c(p$a, -0.2) else p$a)),
    cosy = c(log(p$s), qlogis(p$r)),
    unstr = {
      dd <- a$d
      C <- p$r^abs(outer(1:dd, 1:dd, "-"))
      c(log(p$s), fi("us_theta_cor")(C))
    })
}
par_names <- names(fit$obj$par)
fn_at <- function(p, f = fit) {
  par <- f$obj$env$last.par.best[f$obj$env$lfixed()]
  par[] <- NA
  nm <- names(par)
  ib <- which(nm == "beta")
  par[ib] <- c(p$b0, p$b)
  it <- which(nm == "theta")
  par[it] <- th_of(p)
  rest <- setdiff(seq_along(par), c(ib, it))
  if (length(rest)) {
    # the family's shape: frm's internal log shape
    par[rest] <- log(shape_par)
  }
  f$obj$env$last.par.best[f$obj$env$random] <- 0
  -f$obj$fn(par)
}

rows <- list()
for (i in seq_along(pts)) {
  lb <- laplace_brms(pts[[i]])
  lf <- tryCatch(fn_at(pts[[i]]), error = function(e) NA_real_)
  rows[[i]] <- data.frame(case = case, pt = i, brms = lb, frm = lf, diff = lf - lb)
  cat(sprintf("%-14s pt %d brms %.8f frm %.8f diff %.3e\n", case, i, lb, lf, lf - lb))
}
out <- do.call(rbind, rows)

# deliberately wrong models, on brms's program side: wrong correlation, wrong
# scale, MA sign, row order. The comparison must break.
wrong <- function(tag, p_frm, p_brms) {
  lb <- laplace_brms(p_brms)
  lf <- fn_at(p_frm)
  data.frame(case = case, pt = tag, brms = lb, frm = lf, diff = lf - lb)
}
p1 <- pts[[1]]
w <- rbind(
  wrong("wrong_scale", within(p1, s <- s * 1.1), p1),
  if (grepl("ar", cs$term)) wrong("wrong_ar", within(p1, a <- a + 0.1), p1),
  if (grepl("ma", cs$term)) wrong("ma_flip", within(p1, m <- -m), p1),
  if (grepl("cosy", cs$term)) wrong("wrong_cosy", within(p1, r <- r + 0.1), p1))
print(w)
out <- rbind(out, w)
# wrong row order within a group: refit frm() on data whose time labels are
# reversed within group, same chart point
if (!grepl("single|nogr", case) && !is.null(d$t) && !grepl("notime", case)) {
  dr <- d
  dr$t <- ave(dr$t, dr$g, FUN = function(z) rev(sort(z))[rank(z)])
  # reversing a stationary symmetric law changes nothing for ar and ma: permute
  dr$t <- ave(d$t, d$g, FUN = function(z) if (length(z) > 2) sample(z) else z)
  fr <- msg_of(frm(form, data = dr, family = cs$fam))
  if (is.null(fr$err)) {
    fw <- fr$val
    lf <- tryCatch(fn_at(p1, fw), error = function(e) NA_real_)
    lb <- laplace_brms(p1)
    out <- rbind(out, data.frame(case = case, pt = "wrong_order", brms = lb, frm = lf, diff = lf - lb))
  } else cat("wrong order fit:", fr$err, "\n")
}
write.table(out, outf, sep = "\t", row.names = FALSE, quote = FALSE)
print(out)
