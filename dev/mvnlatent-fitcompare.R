# Fits must be unchanged by the move. This fits every mixture_mvn()
# model of the moved tests (frmtmb's test-mvn-mixture.R, the blocks of
# test-quadrature-defects.R, test-unpinned-seams.R and
# test-prior-location-dpar.R, and frmtmb.sample's test-simulators.R),
# on the data those tests build with the same seeds, and saves what a
# comparison needs: logLik, estimates, standard errors, predict(),
# fitted(), mixture_probs() and simulate().
#
#   Rscript dev/mvnlatent-fitcompare.R <lib or "base"> <out.rds>
#
# "base" reads rellib-r11 (frmtmb 0.73.0 with mixture_mvn() in core);
# a lib reads the lane's build, where frmtmb.latent supplies it.
a <- commandArgs(trailingOnly = TRUE)
base <- c("C:/Users/adf44/source/r/rellib-r11",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(a[1], "base")) base else c(a[1], base))
suppressMessages(library(frmtmb))
if (!"mixture_mvn" %in% getNamespaceExports("frmtmb")) {
  suppressMessages(library(frmtmb.latent))
}
where <- environmentName(environment(mixture_mvn))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "; mixture_mvn from", where, "\n")

res <- list()
grab <- function(id, fit) {
  se <- tryCatch(sqrt(diag(as.matrix(vcov(fit)))),
                 error = function(e) conditionMessage(e))
  fx <- tryCatch(fixef(fit), error = function(e) conditionMessage(e))
  pr <- tryCatch(predict(fit), error = function(e) conditionMessage(e))
  ft <- tryCatch(fitted(fit), error = function(e) conditionMessage(e))
  mp <- tryCatch(mixture_probs(fit), error = function(e) conditionMessage(e))
  sm <- tryCatch(simulate(fit, nsim = 2, seed = 11),
                 error = function(e) conditionMessage(e))
  res[[id]] <<- list(
    logLik = as.numeric(logLik(fit)), df = attr(logLik(fit), "df"),
    est = unlist(fit$estimates), opt = fit$opt$par, se = se, fixef = fx,
    predict = pr, fitted = ft, probs = mp,
    simulate = lapply(sm, unclass), conv = fit$opt$convergence)
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
mvn_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VVV")

# test-mvn-mixture.R, "matches direct ML"
sim <- sim_mvnmix_data()
grab("directML", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = sim$dd))

# "recovers the faithful clusters"
Y <- data.matrix(datasets::faithful)
dd <- data.frame(row = seq_len(nrow(Y)))
dd$Y <- Y
grab("faithful", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd))

# "class means can depend on covariates"
set.seed(43)
n <- 300
x <- rnorm(n)
cl <- rbinom(n, 1, 0.5)
L <- t(chol(matrix(c(0.49, 0.15, 0.15, 0.49), 2)))
E <- L %*% matrix(rnorm(2 * n), 2)
Y <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + E[1, ],
           ifelse(cl == 1, 1, 5) - 0.8 * x + E[2, ])
dd <- data.frame(x = x)
dd$Y <- Y
grab("cov_x", frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd))
grab("cov_0", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd))

# "dpar overrides and covariate gating work"
sim <- sim_mvnmix_data(seed = 44, n = 150)
dd <- sim$dd
dd$x <- rnorm(nrow(dd))
grab("override", frm(bf(Y ~ x, mu2d1 ~ 1) + mixture_mvn(K = 2, D = 2),
                     data = dd))
grab("gate", frm(bf(Y ~ 1, theta1 ~ x) + mixture_mvn(K = 2, D = 2),
                 data = dd))
grab("flat", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd))

# "validation and guards"
sim <- sim_mvnmix_data(seed = 45, n = 80)
grab("guards", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = sim$dd))

# "covariance models size the template correctly"
sim <- sim_mvnmix_data(seed = 46, n = 120)
for (m in mvn_models) {
  grab(paste0("size_", m),
       frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m), data = sim$dd))
}

# "covariance models agree at a shared parameter point"
set.seed(3)
n <- 200
g <- rbinom(n, 1, 0.5)
Y <- cbind(ifelse(g == 1, 0, 3), ifelse(g == 1, 0, 3)) +
  matrix(rnorm(2 * n, sd = 0.9), n, 2)
dd <- data.frame(row = seq_len(n))
dd$Y <- Y
ls <- log(0.9)
point <- list(
  EII = list(sigmaraw = ls),
  VII = list(sigmaraw1 = ls, sigmaraw2 = ls),
  EEI = list(sigmaraw = c(ls, ls)),
  VEI = list(sigmavol1 = ls, sigmavol2 = ls, sigmashape = 0),
  EVI = list(sigmavol = ls, sigmashape1 = 0, sigmashape2 = 0),
  VVI = list(sigmaraw1 = c(ls, ls), sigmaraw2 = c(ls, ls)),
  EEE = list(sigmaraw = c(ls, ls, 0)),
  VVV = list(sigmaraw1 = c(ls, ls, 0), sigmaraw2 = c(ls, ls, 0))
)
for (m in mvn_models) {
  st <- c(list(beta = c(0, 0, 3, 3), betad = 0.2), point[[m]])
  grab(paste0("point_", m),
       frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m), data = dd,
           start = st))
}

# "class means take covariates under every model"
set.seed(47)
n <- 300
x <- rnorm(n)
cl <- rbinom(n, 1, 0.5)
Y <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + rnorm(n, 0, 0.7),
           ifelse(cl == 1, 1, 5) - 0.8 * x + rnorm(n, 0, 0.7))
dd <- data.frame(x = x)
dd$Y <- Y
ctl <- frmtmb_control(restarts = 3)
for (m in mvn_models) {
  grab(paste0("every_x_", m),
       frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2, model = m), data = dd,
           control = ctl))
  grab(paste0("every_0_", m),
       frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m), data = dd,
           control = ctl))
}

# "matches mclust for every covariance model": the starts come from
# mclust, which is the same package in both arms
us_theta_cor <- function(C) {
  Lc <- t(chol(C))
  L <- Lc / diag(Lc)
  L[lower.tri(L)]
}
mclust_extras <- function(model, Sig, K, D) {
  lsd <- vapply(Sig, function(S) 0.5 * log(diag(S)), numeric(D))
  us <- function(S) {
    c(log(sqrt(diag(S))), us_theta_cor(stats::cov2cor(S)))
  }
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
    VVV = per("sigmaraw", lapply(Sig, us))
  )
}
mclustBIC <- mclust::mclustBIC
mclust_fits <- function(tag, Y, K, D) {
  dd <- data.frame(row = seq_len(nrow(Y)))
  dd$Y <- Y
  for (m in mvn_models) {
    mc <- mclust::Mclust(Y, G = K, modelNames = m, verbose = FALSE,
                         control = mclust::emControl(tol = 1e-13,
                                                     itmax = 1e6))
    Sig <- lapply(seq_len(K),
                  function(k) mc$parameters$variance$sigma[, , k])
    st <- c(list(beta = as.vector(mc$parameters$mean),
                 betad = log(mc$parameters$pro[-K] / mc$parameters$pro[K])),
            mclust_extras(m, Sig, K, D))
    st <- lapply(st, unname)
    grab(paste0(tag, "_", m),
         frm(bf(Y ~ 1) + mixture_mvn(K = K, D = D, model = m), data = dd,
             start = st))
  }
}
Y <- unname(data.matrix(datasets::faithful))
Y[, 2] <- Y[, 2] / 10
mclust_fits("mclust2", Y, 2L, 2L)
set.seed(11)
n <- 600
g <- sample(1:3, n, TRUE)
mus <- rbind(c(0, 0, 0), c(3, 1, -2), c(-2, 3, 2))
Y3 <- mus[g, ] + matrix(rnorm(3 * n, sd = 0.8), n, 3) %*%
  chol(matrix(c(1, .3, .1, .3, 1, -.2, .1, -.2, 1), 3))
mclust_fits("mclust3", Y3, 3L, 3L)

# test-quadrature-defects.R, "mixture_mvn() refuses REML and profile too"
set.seed(5)
n <- 200
z <- rbinom(n, 1, 0.4)
dz <- data.frame(x = rnorm(n))
dz$Y <- cbind(rnorm(n, ifelse(z == 1, 3, 0)),
              rnorm(n, ifelse(z == 1, -2, 0)))
grab("quadrature_defects", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2),
                               data = dz))

# test-unpinned-seams.R
set.seed(42)
n <- 200
cl <- rbinom(n, 1, 0.4)
x <- rnorm(n)
Y <- cbind(ifelse(cl == 1, 0, 3) + 0.5 * x + rnorm(n),
           ifelse(cl == 1, 0, 4) + rnorm(n))
dd <- data.frame(x = x)
dd$Y <- Y
grab("unpinned", frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd))

# frmtmb.sample's test-simulators.R
set.seed(46)
n <- 300L
cl <- stats::rbinom(n, 1L, 0.4)
L1 <- t(chol(matrix(c(1, 0.5, 0.5, 1), 2L)))
L2 <- t(chol(matrix(c(0.5, -0.2, -0.2, 0.8), 2L)))
E <- matrix(stats::rnorm(2L * n), 2L)
Y <- t(ifelse(matrix(cl == 1L, 2L, n, byrow = TRUE),
              c(0, 0) + L1 %*% E, c(3, 4) + L2 %*% E))
dd <- data.frame(row = seq_len(n))
dd$Y <- Y
grab("sample_sim", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd))

cat("fits:", length(res), "\n")
saveRDS(res, a[2])
