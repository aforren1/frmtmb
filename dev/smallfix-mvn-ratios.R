# The quantities frmtmb.latent's test-mixture-mvn.R compares, on the
# right answer and on a constructed wrong one, so each absolute
# tolerance can become a ratio to something the run measures and still
# fail on a wrong answer. Same data, seeds and fits as the test file.
#   Rscript dev/smallfix-mvn-ratios.R <lib;lib>
# Log: dev/smallfix-out/mvn-ratios.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "| latent", format(packageVersion("frmtmb.latent")), "\n")
eps <- .Machine$double.eps
row <- function(test, what, right, wrong, bound) {
  cat(sprintf("%-4s %-34s right %-10.3g wrong %-10.3g bound %-8.3g %s\n",
              test, what, right, wrong, bound,
              if (right < bound && wrong >= bound) "SEPARATES" else "FAILS"))
}
rel <- function(x, y) max(abs(unname(x) - unname(y))) / max(abs(unname(y)))
se_of <- function(fit, nm, coef) {
  fx <- fixef(fit)
  fx[paste0(nm, "_", coef), "Est.Error"]
}
est_of <- function(fit, nm, coef) {
  fixef(fit)[paste0(nm, "_", coef), "Estimate"]
}

sim_mvnmix_data <- function(seed = 42, n = 400, p1 = 0.4,
                            m1 = c(0, 0), m2 = c(3, 4)) {
  set.seed(seed)
  cl <- rbinom(n, 1, p1)
  L1 <- t(chol(matrix(c(1, 0.5, 0.5, 1), 2)))
  L2 <- t(chol(matrix(c(0.5, -0.2, -0.2, 0.8), 2)))
  E <- matrix(rnorm(2 * n), 2)
  Y <- t(ifelse(matrix(cl == 1, 2, n, byrow = TRUE),
                m1 + L1 %*% E, m2 + L2 %*% E))
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  list(dd = dd, cl = cl)
}

## T1-T3: "mixture_mvn matches direct ML"
sim <- sim_mvnmix_data()
fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = sim$dd)
Y <- sim$dd$Y
ldmvn2 <- function(Y, m, S) {
  det_ <- S[1, 1] * S[2, 2] - S[1, 2]^2
  d1 <- Y[, 1] - m[1]
  d2 <- Y[, 2] - m[2]
  q <- (S[2, 2] * d1^2 - 2 * S[1, 2] * d1 * d2 + S[1, 1] * d2^2) / det_
  -log(2 * pi) - 0.5 * log(det_) - 0.5 * q
}
nll <- function(p) {
  L1 <- matrix(c(exp(p[5]), p[7], 0, exp(p[6])), 2)
  L2 <- matrix(c(exp(p[8]), p[10], 0, exp(p[9])), 2)
  pi1 <- stats::plogis(p[11])
  -sum(log(pi1 * exp(ldmvn2(Y, p[1:2], L1 %*% t(L1))) +
             (1 - pi1) * exp(ldmvn2(Y, p[3:4], L2 %*% t(L2)))))
}
op <- suppressWarnings(
  stats::optim(c(0, 0, 3, 4, 0, 0, 0.5, log(0.7), log(0.8), -0.3,
                 stats::qlogis(0.4)),
               nll, method = "BFGS",
               control = list(reltol = 1e-13, maxit = 5000)))
ll <- as.numeric(logLik(fit))
# wrong: the direct ML stopped short, its mixing logit 0.05 off
pw <- op$par
pw[11] <- pw[11] + 0.05
cat("T1 ll", format(ll, digits = 15), "direct", format(-op$value, digits = 15),
    "wrong", format(-nll(pw), digits = 15), "\n")
row("T1", "|ll - direct| / |ll|", abs(ll + op$value) / abs(ll),
    abs(ll + nll(pw)) / abs(ll), 1e-7)
P <- mixture_probs(fit)
# wrong: posterior class weights left unnormalized by one ulp per row
# is invisible; a missed normalization by the weights' own sum is not
Pw <- P * (1 + 1e-9)
row("T2", "max|rowSums(P) - 1| / eps", max(abs(rowSums(P) - 1)) / eps,
    max(abs(rowSums(Pw) - 1)) / eps, 16)
fv <- fitted(fit)[, "Estimate", ]
se_mean <- apply(Y, 2, stats::sd) / sqrt(nrow(Y))
# wrong: the mixture mean at mixing weights 1% off the estimate
pr <- colMeans(P)
m1 <- c(est_of(fit, "mu1d1", "Intercept"), est_of(fit, "mu1d2", "Intercept"))
m2 <- c(est_of(fit, "mu2d1", "Intercept"), est_of(fit, "mu2d2", "Intercept"))
fvw <- colMeans(fv) + 0.01 * (m1 - m2)
row("T3", "max|colMeans diff| / se(mean)",
    max(abs(colMeans(fv) - colMeans(Y)) / se_mean),
    max(abs(fvw - colMeans(Y)) / se_mean), 1e-3)

## T4: "mixture_mvn recovers the faithful clusters"
Yf <- data.matrix(datasets::faithful)
df <- data.frame(row = seq_len(nrow(Yf)))
df$Y <- Yf
ff <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = df)
ref <- rbind(c(2.04, 54.5), c(4.29, 80.0))
zf <- function(f) {
  m <- rbind(c(est_of(f, "mu1d1", "Intercept"), est_of(f, "mu1d2", "Intercept")),
             c(est_of(f, "mu2d1", "Intercept"), est_of(f, "mu2d2", "Intercept")))
  s <- rbind(c(se_of(f, "mu1d1", "Intercept"), se_of(f, "mu1d2", "Intercept")),
             c(se_of(f, "mu2d1", "Intercept"), se_of(f, "mu2d2", "Intercept")))
  o <- order(m[, 1])
  print(cbind(m[o, ], s[o, ]))
  max(abs(m[o, ] - ref) / s[o, ])
}
# wrong: the spherical model, the right family with the wrong covariance
fw <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = "EII"), data = df)
row("T4", "max|mean - ref| / se", zf(ff), zf(fw), 0.5)

## T5: slopes, "class means can depend on covariates" (seed 43)
set.seed(43)
n <- 300
x <- rnorm(n)
cl <- rbinom(n, 1, 0.5)
L <- t(chol(matrix(c(0.49, 0.15, 0.15, 0.49), 2)))
E <- L %*% matrix(rnorm(2 * n), 2)
Y5 <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + E[1, ],
            ifelse(cl == 1, 1, 5) - 0.8 * x + E[2, ])
d5 <- data.frame(x = x)
d5$Y <- Y5
zslope <- function(f) {
  z <- c(sapply(c("mu1d1", "mu2d1"), function(nm)
    abs(est_of(f, nm, "x") - 1.2) / se_of(f, nm, "x")),
    sapply(c("mu1d2", "mu2d2"), function(nm)
      abs(est_of(f, nm, "x") + 0.8) / se_of(f, nm, "x")))
  print(round(z, 3))
  max(z)
}
f_x <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = d5)
# wrong: the covariate's values shuffled, so each slope estimates 0
d5w <- d5
set.seed(1)
d5w$x <- sample(d5w$x)
f_w <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = d5w)
row("T5", "max|slope - truth| / se", zslope(f_x), zslope(f_w), 4)

## T5b: the same under every covariance model (seed 47)
mvn_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VVV")
set.seed(47)
x <- rnorm(n)
cl <- rbinom(n, 1, 0.5)
Y7 <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + rnorm(n, 0, 0.7),
            ifelse(cl == 1, 1, 5) - 0.8 * x + rnorm(n, 0, 0.7))
d7 <- data.frame(x = x)
d7$Y <- Y7
ctl <- frmtmb_control(restarts = 3)
zs <- numeric(0)
lr <- numeric(0)
for (m in mvn_models) {
  fx7 <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2, model = m), data = d7,
             control = ctl)
  zs[m] <- zslope(fx7)
}
d7w <- d7
set.seed(2)
d7w$x <- sample(d7w$x)
fw7 <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2, model = "VVV"), data = d7w,
           control = ctl)
row("T5b", "max over models |slope-truth|/se", max(zs), zslope(fw7), 4)

## T6: the gated model nests the flat one (seed 44)
sim <- sim_mvnmix_data(seed = 44, n = 150)
dd <- sim$dd
dd$x <- rnorm(nrow(dd))
f_gate <- frm(bf(Y ~ 1, theta1 ~ x) + mixture_mvn(K = 2, D = 2), data = dd)
f_flat <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
lg <- as.numeric(logLik(f_gate))
lf <- as.numeric(logLik(f_flat))
cat("T6 gate", format(lg, digits = 15), "flat", format(lf, digits = 15),
    "\n")
# wrong: a gated maximum 1e-3 below the flat one it contains
row("T6", "(ll_flat - ll_gate)/|ll_flat|", max(0, lf - lg) / abs(lf),
    1e-3 / abs(lf), 1e-8)

## T7, T8: the shared point and the nesting chains (seed 3)
set.seed(3)
n <- 200
g <- rbinom(n, 1, 0.5)
Y3 <- cbind(ifelse(g == 1, 0, 3), ifelse(g == 1, 0, 3)) +
  matrix(rnorm(2 * n, sd = 0.9), n, 2)
d3 <- data.frame(row = seq_len(n))
d3$Y <- Y3
ls <- log(0.9)
point <- list(
  EII = list(sigmaraw = ls),
  VII = list(sigmaraw1 = ls, sigmaraw2 = ls),
  EEI = list(sigmaraw = c(ls, ls)),
  VEI = list(sigmavol1 = ls, sigmavol2 = ls, sigmashape = 0),
  EVI = list(sigmavol = ls, sigmashape1 = 0, sigmashape2 = 0),
  VVI = list(sigmaraw1 = c(ls, ls), sigmaraw2 = c(ls, ls)),
  EEE = list(sigmaraw = c(ls, ls, 0)),
  VVV = list(sigmaraw1 = c(ls, ls, 0), sigmaraw2 = c(ls, ls, 0)))
mvn_start_vector <- function(fit, st) {
  tpl <- fit$frame$par_template
  v <- unlist(lapply(names(tpl), function(nm) {
    if (lengths(tpl[nm]) == 0L) numeric(0) else st[[nm]]
  }), use.names = FALSE)
  stopifnot(length(v) == length(fit$obj$par))
  v
}
llp <- numeric(0)
llf <- numeric(0)
llw <- NA
for (m in mvn_models) {
  st <- c(list(beta = c(0, 0, 3, 3), betad = 0.2), point[[m]])
  fm <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m), data = d3,
            start = st)
  llp[m] <- -fm$obj$fn(mvn_start_vector(fm, st))
  llf[m] <- as.numeric(logLik(fm))
  if (m == "VEI") {
    # wrong: a point the models do not share, the shape 1e-3 off
    sw <- st
    sw$sigmashape <- 1e-3
    llw <- -fm$obj$fn(mvn_start_vector(fm, sw))
  }
}
print(llp, digits = 17)
row("T7", "max|ll_m - ll_EII| / (|ll| eps)",
    max(abs(llp - llp[["EII"]])) / (abs(llp[["EII"]]) * eps),
    abs(llw - llp[["EII"]]) / (abs(llp[["EII"]]) * eps), 1024)
chains <- list(c("EII", "VII", "VEI", "VVI", "VVV"),
               c("EII", "EEI", "VEI", "VVI", "VVV"),
               c("EII", "EEI", "EVI", "VVI", "VVV"),
               c("EII", "EEI", "EEE", "VVV"))
gap <- 0
for (ch in chains) for (i in seq_len(length(ch) - 1L)) {
  gap <- max(gap, (llf[[ch[i]]] - llf[[ch[i + 1L]]]) / abs(llf[[ch[i]]]))
}
print(llf, digits = 17)
# wrong: VVV's maximum 1e-3 below EEE's
row("T8", "max (ll_sub - ll_super)/|ll_sub|", max(gap, 0),
    1e-3 / abs(llf[["EEE"]]), 1e-8)

## T9, T10: mclust agreement
mclust_extras <- function(model, Sig, K, D) {
  lsd <- vapply(Sig, function(S) 0.5 * log(diag(S)), numeric(D))
  us <- function(S) us_theta(S)
  per <- function(pre, v) stats::setNames(v, paste0(pre, seq_len(K)))
  free <- function(k) (lsd[, k] - mean(lsd[, k]))[seq_len(D - 1L)]
  switch(model,
    EII = list(sigmaraw = mean(lsd)),
    VII = per("sigmaraw", as.list(colMeans(lsd))),
    EEI = list(sigmaraw = rowMeans(lsd)),
    VEI = c(per("sigmavol", as.list(colMeans(lsd))),
            list(sigmashape = free(1))),
    EVI = c(list(sigmavol = mean(lsd)),
            per("sigmashape", lapply(seq_len(K), free))),
    VVI = per("sigmaraw", lapply(seq_len(K), function(k) lsd[, k])),
    EEE = list(sigmaraw = us(Sig[[1]])),
    VVV = per("sigmaraw", lapply(Sig, us)))
}
mc_rows <- function(Y, K, D, lab) {
  dd <- data.frame(row = seq_len(nrow(Y)))
  dd$Y <- Y
  mclustBIC <- mclust::mclustBIC
  r <- matrix(NA, 0, 4)
  for (m in mvn_models) {
    mc <- mclust::Mclust(Y, G = K, modelNames = m, verbose = FALSE,
                         control = mclust::emControl(tol = 1e-13,
                                                     itmax = 1e6))
    Sig <- lapply(seq_len(K), function(k) mc$parameters$variance$sigma[, , k])
    st <- c(list(beta = as.vector(mc$parameters$mean),
                 betad = log(mc$parameters$pro[-K] / mc$parameters$pro[K])),
            mclust_extras(m, Sig, K, D))
    st <- lapply(st, unname)
    fit <- frm(bf(Y ~ 1) + mixture_mvn(K = K, D = D, model = m),
               data = dd, start = st)
    fx <- fixef_by_dpar(fit)
    ourM <- t(vapply(seq_len(K), function(k) vapply(seq_len(D), function(j)
      unname(fx[[paste0("mu", k, "d", j)]]), numeric(1)), numeric(D)))
    sg <- fit$spec$responses[[1]]$family$mix$sigma
    ex <- fit$estimates[fit$frame$extra_names]
    ourS <- lapply(seq_len(K), function(k) sg(ex, k))
    mcM <- t(mc$parameters$mean)
    ord <- order(ourM[, 1])
    mo <- order(mcM[, 1])
    r <- rbind(r, c(
      abs(as.numeric(logLik(fit)) - mc$loglik) / abs(mc$loglik),
      rel(ourM[ord, ], mcM[mo, ]), rel(unlist(ourS[ord]), unlist(Sig[mo])),
      rel(mixture_probs(fit)[, ord], mc$z[, mo])))
    rownames(r)[nrow(r)] <- paste(lab, m)
  }
  colnames(r) <- c("ll", "means", "covs", "probs")
  r
}
Yf2 <- unname(data.matrix(datasets::faithful))
Yf2[, 2] <- Yf2[, 2] / 10
set.seed(11)
n <- 600
g <- sample(1:3, n, TRUE)
mus <- rbind(c(0, 0, 0), c(3, 1, -2), c(-2, 3, 2))
Y3d <- mus[g, ] + matrix(rnorm(3 * n, sd = 0.8), n, 3) %*%
  chol(matrix(c(1, .3, .1, .3, 1, -.2, .1, -.2, 1), 3))
R <- rbind(mc_rows(Yf2, 2L, 2L, "faithful"), mc_rows(Y3d, 3L, 3L, "sim3"))
print(signif(R, 3))
# wrong: each model's frmtmb fit against mclust's next model in the list
mclustBIC <- mclust::mclustBIC
mcv <- mclust::Mclust(Yf2, G = 2, modelNames = "EII", verbose = FALSE)
mvv <- mclust::Mclust(Yf2, G = 2, modelNames = "VVV", verbose = FALSE)
w_ll <- abs(mcv$loglik - mvv$loglik) / abs(mvv$loglik)
w_m <- rel(t(mcv$parameters$mean), t(mvv$parameters$mean))
w_s <- rel(mcv$parameters$variance$sigma, mvv$parameters$variance$sigma)
w_z <- rel(mcv$z, mvv$z)
row("T9", "|ll - mclust ll| / |mclust ll|", max(R[, "ll"]), w_ll, 1e-9)
row("T10", "means rel", max(R[, "means"]), w_m, 1e-5)
row("T10", "covariances rel", max(R[, "covs"]), w_s, 1e-5)
row("T10", "probabilities rel", max(R[, "probs"]), w_z, 1e-5)
