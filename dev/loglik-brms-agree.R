# Lane loglik: frmtmb.sample's log_lik() against brms 2.23.0's log_lik()
# at the SAME parameter vectors, three per case: frmtmb draws written by
# hand (no sampler), and a brms fit whose three chains sit at those
# vectors (Fixed_param). Every residual covariance term, gaussian and
# student, ragged groups (no interior gap for the ar() family, where
# brms reads lags by position; one for cosy() and unstr(), where both
# readings agree), rows shuffled, sigma ~ x. Seed 4.
# Output: dev/loglik-brms-agree.txt
source("C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/loglik-brms-lib.R")
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"),
    find.package("frmtmb.spatial"), "\n")
set.seed(4)
sizes <- c(7, 5, 6, 4, 7, 3, 6, 5)
d <- do.call(rbind, lapply(seq_along(sizes), function(g) {
  data.frame(g = g, t = seq_len(sizes[g]))
}))
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d)) + rnorm(8)[d$g]
d <- d[sample(nrow(d)), ]
rownames(d) <- NULL
dgap <- d[!(d$g == 2 & d$t == 3) & !(d$g == 5 & d$t == 1), ]
rownames(dgap) <- NULL
n <- nrow(d)
V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
dg <- 1 / sqrt(rowSums(P))
W <- P * outer(dg, dg)   # symmetric: real eigenvalues for brms's Stan
d2 <- list(V = V, W = W)

draws_at <- function(fit, vals) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  k <- max(lengths(vals))
  m <- matrix(0, k, length(lab), dimnames = list(NULL, lab))
  for (nm in names(vals)) {
    stopifnot(nm %in% lab)
    m[, nm] <- vals[[nm]]
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}
pinv <- function(r) r / sqrt(1 - r^2)   # frmtmb's pacf chart, inverted
b0 <- c(1.1, 0.9, 1.3)
b1 <- c(0.4, 0.6, 0.2)
sg <- c(0.8, 1.2, 1.0)
nus <- c(5, 12, 3.5)
three <- function(f) lapply(1:3, f)
cmp <- function(nm, llf, bfit) {
  llb <- log_lik(bfit)
  cat(sprintf("%-22s %d x %d  max|frmtmb - brms| / max|brms| = %.3e\n",
              nm, nrow(llf), ncol(llf),
              max(abs(llf - llb)) / max(abs(llb))))
}

# ar(), ma(), arma(): brms's sigma is the innovation sd, frmtmb's the
# marginal one; gamma0 is the ARMA(1,1) variance at unit innovations
g0 <- function(a, m) (1 + 2 * a * m + m^2) / (1 - a^2)
tcases <- list(
  ar = list(f = y ~ x + ar(t, g, cov = TRUE), a = c(0.5, -0.3, 0.8),
            m = c(0, 0, 0)),
  ma = list(f = y ~ x + ma(t, g, cov = TRUE), a = c(0, 0, 0),
            m = c(0.4, -0.6, 0.2)),
  arma = list(f = y ~ x + arma(t, g, cov = TRUE), a = c(0.5, -0.2, 0.3),
              m = c(0.3, 0.5, -0.4)))
for (nm in names(tcases)) {
  cs <- tcases[[nm]]
  for (fam in c("gaussian", "student")) {
    fit <- frm(cs$f, family = get(fam)(), data = d, dry_run = "objective")
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
    th <- list()
    if (nm != "ma") th$a <- pinv(cs$a)
    if (nm != "ar") th$m <- pinv(cs$m)
    for (k in seq_along(th)) vals[[paste0("thetaac_", k)]] <- th[[k]]
    if (fam == "student") vals$nu <- nus
    llf <- log_lik(draws_at(fit, vals))
    inits <- three(function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s] / sqrt(g0(cs$a[s], cs$m[s])))
      if (nm != "ma") p$ar <- array(cs$a[s], 1)
      if (nm != "ar") p$ma <- array(cs$m[s], 1)
      if (fam == "student") p$nu <- nus[s]
      p
    })
    bf_ <- brms::bf(stats::as.formula(deparse(cs$f)))
    cmp(paste(nm, fam), llf, ll_fixed_fit(bf_, brms::brmsfamily(fam), d, inits))
  }
}

# sigma ~ x under ar(): log sigma shifts by log sqrt(1 - ar^2)
for (fam in c("gaussian", "student")) {
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), family = get(fam)(),
             data = d, dry_run = "objective")
  a <- c(0.5, -0.4, 0.2)
  bs0 <- c(-0.2, 0.1, 0)
  bs1 <- c(0.3, -0.2, 0.1)
  vals <- list(b_Intercept = b0, b_x = b1, b_sigma_Intercept = bs0,
               b_sigma_x = bs1, thetaac_1 = pinv(a))
  if (fam == "student") vals$nu <- nus
  llf <- log_lik(draws_at(fit, vals))
  inits <- three(function(s) {
    p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
              ar = array(a[s], 1), b_sigma = array(bs1[s], 1),
              Intercept_sigma = bs0[s] + bs1[s] * mean(d$x) +
                0.5 * log(1 - a[s]^2))
    if (fam == "student") p$nu <- nus[s]
    p
  })
  cmp(paste("ar sigma~x", fam), llf,
      ll_fixed_fit(brms::bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x),
                   brms::brmsfamily(fam), d, inits))
}

# cosy() on ragged data with gaps: one correlation for every pair, so
# positions and levels agree. frmtmb's cosy interval is (-1/(d-1), 1)
for (fam in c("gaussian", "student")) {
  fit <- frm(y ~ x + cosy(t, g), family = get(fam)(), data = dgap,
             dry_run = "objective")
  dd <- fit$frame$autocor$y$d
  a <- 1 / (dd - 1)
  rho <- c(0.4, 0.1, 0.7)
  vals <- list(b_Intercept = b0, b_x = b1, sigma = sg,
               thetaac_1 = qlogis((rho + a) / (1 + a)))
  if (fam == "student") vals$nu <- nus
  llf <- log_lik(draws_at(fit, vals))
  inits <- three(function(s) {
    p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(dgap$x),
              sigma = sg[s], cosy = rho[s])
    if (fam == "student") p$nu <- nus[s]
    p
  })
  cmp(paste("cosy gaps", fam), llf,
      ll_fixed_fit(brms::bf(y ~ x + cosy(t, g)), brms::brmsfamily(fam), dgap, inits))
}

# unstr() on ragged data with gaps: brms carries the time index
for (fam in c("gaussian", "student")) {
  fit <- frm(y ~ x + unstr(t, g), family = get(fam)(), data = dgap,
             dry_run = "objective")
  ac <- fit$frame$autocor$y
  set.seed(40)
  ths <- three(function(s) rnorm(ac$npar, 0, 0.5))
  vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
  for (k in seq_len(ac$npar)) {
    vals[[paste0("thetaac_", k)]] <- vapply(ths, `[`, 0, k)
  }
  if (fam == "student") vals$nu <- nus
  llf <- log_lik(draws_at(fit, vals))
  inits <- three(function(s) {
    p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(dgap$x),
              sigma = sg[s],
              Lcortime = as.matrix(frmtmb:::us_chol_L(ths[[s]], ac$d)))
    if (fam == "student") p$nu <- nus[s]
    p
  })
  cmp(paste("unstr gaps", fam), llf,
      ll_fixed_fit(brms::bf(y ~ x + unstr(t, g)), brms::brmsfamily(fam), dgap, inits))
}

# fcor() and sar(): one group of every row
for (term in c("fcor(V)", "sar(W)", "sar(W, type = 'error')")) {
  for (fam in c("gaussian", "student")) {
    fo <- stats::as.formula(paste("y ~ x +", term))
    fit <- frm(fo, family = get(fam)(), data = d, data2 = d2,
               dry_run = "objective")
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
    rho <- c(0.4, -0.3, 0.6)
    if (term != "fcor(V)") {
      aux <- fit$frame$autocor$y$aux
      vals$thetaac_1 <- qlogis((rho - aux$lo) / (aux$hi - aux$lo))
    }
    if (fam == "student") vals$nu <- nus
    llf <- log_lik(draws_at(fit, vals))
    inits <- three(function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s])
      if (term == "sar(W)") p$lagsar <- rho[s]
      if (grepl("error", term)) p$errorsar <- rho[s]
      if (fam == "student") p$nu <- nus[s]
      p
    })
    cmp(paste(term, fam), llf,
        ll_fixed_fit(brms::bf(fo), brms::brmsfamily(fam), d, inits, data2 = d2))
  }
}
