# Second surface probe: each claim that frmtmb does something drmTMB does
# not (and the reverse) is run in BOTH packages here, so the related-work
# text rests on a run, not on either package's documentation.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
probe <- function(label, expr) {
  cl <- substitute(expr)
  cat("\n=== ", label, "\n", sep = "")
  cat(paste(deparse(cl, width.cutoff = 100), collapse = "\n"), "\n")
  res <- tryCatch(suppressWarnings(eval(cl, parent.frame())),
                  error = function(e) e)
  if (inherits(res, "error")) {
    cat("  REFUSED/ERROR:", substr(gsub("\n", " ", conditionMessage(res)), 1,
                                   400), "\n")
    return(invisible(NULL))
  }
  if (inherits(res, c("drmTMB", "frmtmb_fit"))) {
    ll <- tryCatch(format(as.numeric(logLik(res)), digits = 12),
                   error = function(e) paste("unavailable:",
                                             conditionMessage(e)))
    cat("  FITTED. logLik =", substr(ll, 1, 120), " convergence =",
        res$opt$convergence, "(", res$opt$message, ")\n")
  } else if (inherits(res, "ggplot")) {
    cat("  RETURNED a ggplot with", length(res$layers), "layers\n")
  } else print(res)
  invisible(res)
}
d <- sim_grouped(101)
d$yc3 <- d$yc
set.seed(7)
d$yzi <- ifelse(runif(nrow(d)) < 0.2, 0L, d$yc)
d$y3 <- rnorm(nrow(d))
d$ybin <- rbinom(nrow(d), 1, plogis(0.3 * d$x))
D <- function(...) drmTMB::drmTMB(...)

cat("\n##### A. frmtmb fits these; drmTMB 0.7.0 refusals are in surface.log\n")
probe("frm: nonlinear formula",
  frm(bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
      family = gaussian(), data = d, start = list(beta = c(1, 0.3))))
probe("frm: three responses with rescor",
  frm(mvbf(bf(y ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE),
      family = gaussian(), data = d))
probe("frm: mixed-family multivariate, gaussian + negbinomial, shared |p|",
  frm(mvbf(bf(y2 ~ x + (1 | p | g), family = gaussian()),
           bf(yc ~ x + (1 | p | g), family = negbinomial())), data = d))
probe("frm: |p| across mu and shape, negbinomial",
  frm(bf(yc ~ x + (1 | p | g), shape ~ 1 + (1 | p | g)),
      family = negbinomial(), data = d))
probe("frm: |p| across mu and phi, Beta",
  frm(bf(yb ~ x + (1 | p | g), phi ~ 1 + (1 | p | g)), family = Beta(),
      data = d))
probe("frm: correlated slope block, Beta",
  frm(bf(yb ~ x + (1 + x | g), phi ~ 1), family = Beta(), data = d))
probe("frm: student sigma random effect",
  frm(bf(yt ~ x + (1 | g), sigma ~ 1 + (1 | g), nu ~ 1), family = student(),
      data = d))
probe("frm: zero-inflation random effect",
  frm(bf(yzi ~ x, shape ~ 1, zi ~ 1 + (1 | g)),
      family = zero_inflated_negbinomial(), data = d))
probe("frm: REML, Beta",
  frm(bf(yb ~ x + (1 | g), phi ~ 1), family = Beta(), data = d, REML = TRUE))
probe("frm: REML, negbinomial",
  frm(bf(yc ~ x + (1 | g), shape ~ 1), family = negbinomial(), data = d,
      REML = TRUE))
probe("frm: penalized smooth s(x)",
  frm(bf(y ~ s(x) + (1 | g)), family = gaussian(), data = d))
# A logistic log-density: smooth, so nlminb's convergence code means
# something (a Laplace density's kink gave "false convergence").
probe("frm: custom family (logistic log-density)",
  frm(bf(y ~ x), data = d,
      family = custom_family("mylogistic", dpars = c("mu", "sigma"),
        links = list(mu = "identity", sigma = "log"), type = "continuous",
        lpdf = function(y, dpars, aterms)
          -(y - dpars$mu) / dpars$sigma - log(dpars$sigma) -
            2 * log1p(exp(-(y - dpars$mu) / dpars$sigma)))))
probe("frm: nbinom2 mu RE and shape RE together",
  frm(bf(yc ~ x + (1 | g), shape ~ 1 + (1 | g)), family = negbinomial(),
      data = d))

cat("\n##### B. drmTMB extras, run in drmTMB and tried in frmtmb\n")
f_g <- D(dbf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d)
probe("drm: worm_plot()", drmTMB::worm_plot(f_g, seed = 1))
probe("drm: centile_chart()", drmTMB::centile_chart(f_g, covariate = "x"))
probe("drm: check_drm()", drmTMB::check_drm(f_g))
probe("drm: anova() LRT between two fits",
  anova(D(dbf(y ~ x, sigma ~ z), family = gaussian(), data = d), f_g))
probe("drm: mspl estimator, binomial random intercept",
  D(dbf(ybin ~ x + (1 | g)), family = binomial(), data = d,
    estimator = "mspl"))
probe("drm: binary missing predictor, impute_model(binomial)",
  { dm <- d; dm$bx <- rbinom(nrow(dm), 1, plogis(dm$z))
    dm$y <- dm$y + 0.5 * dm$bx; dm$bx[seq(1, nrow(dm), 7)] <- NA
    D(dbf(y ~ x + mi(bx), sigma ~ 1), family = gaussian(), data = dm,
      impute = list(bx = drmTMB::impute_model(bx ~ z, family = binomial())),
      missing = drmTMB::miss_control(predictor = "model")) })
probe("frm: binary missing predictor, bf(bx | mi() ~ z) + bernoulli",
  { dm <- d; dm$bx <- rbinom(nrow(dm), 1, plogis(dm$z))
    dm$y <- dm$y + 0.5 * dm$bx; dm$bx[seq(1, nrow(dm), 7)] <- NA
    frm(bf(y ~ x + mi(bx), family = gaussian()) +
          bf(bx | mi() ~ z, family = bernoulli()) + set_rescor(FALSE),
        data = dm) })
probe("drm: missing response included, gaussian",
  { dm <- d; dm$y[seq(1, nrow(dm), 9)] <- NA
    D(dbf(y ~ x + (1 | g), sigma ~ 1), family = gaussian(), data = dm,
      missing = drmTMB::miss_control(response = "include")) })
probe("drm: bivariate labelled all-four block mu1, mu2, sigma1, sigma2",
  D(dbf(mu1 = y ~ x + (1 | p | g), mu2 = y2 ~ x + (1 | p | g),
        sigma1 = ~ 1 + (1 | p | g), sigma2 = ~ 1 + (1 | p | g)),
    family = drmTMB::biv_gaussian(), data = d))
probe("frm: same all-four block",
  frm(mvbf(bf(y ~ x + (1 | p | g), sigma ~ 1 + (1 | p | g)),
           bf(y2 ~ x + (1 | p | g), sigma ~ 1 + (1 | p | g))) +
        set_rescor(TRUE), family = gaussian(), data = d))
probe("drm: rho12 regression rho12 ~ x",
  D(dbf(mu1 = y ~ x, mu2 = y2 ~ x, rho12 = ~ x),
    family = drmTMB::biv_gaussian(), data = d))

# Bipartite phylogenetic interaction: drmTMB builds the Kronecker
# precision from two trees; frmtmb needs the Kronecker covariance handed
# in through gr(cov =). Same model if both agree.
set.seed(8)
t1 <- ape::rcoal(8); t1$tip.label <- paste0("p", 1:8)
t2 <- ape::rcoal(6); t2$tip.label <- paste0("q", 1:6)
pairs <- expand.grid(plant = t1$tip.label, poll = t2$tip.label)
dn <- pairs[rep(seq_len(nrow(pairs)), 3), ]
dn$x <- rnorm(nrow(dn))
C1 <- ape::vcv(t1, corr = TRUE); C2 <- ape::vcv(t2, corr = TRUE)
K <- kronecker(C2, C1)
nm <- paste(pairs$plant, pairs$poll, sep = ":")
dimnames(K) <- list(nm, nm)
a <- as.vector(t(chol(K)) %*% rnorm(nrow(K))) * 0.7
dn$pair <- factor(paste(dn$plant, dn$poll, sep = ":"), levels = nm)
dn$y <- 0.3 + 0.4 * dn$x + a[as.integer(dn$pair)] + rnorm(nrow(dn), 0, 0.5)
dn$plant <- factor(dn$plant); dn$poll <- factor(dn$poll)
fpi <- probe("drm: phylo_interaction(1 | plant:poll)",
  D(dbf(y ~ x + phylo_interaction(1 | plant:poll, tree1 = t1, tree2 = t2),
        sigma ~ 1), family = gaussian(), data = dn))
fpk <- probe("frm: gr(pair, cov = kronecker(C2, C1))",
  frm(bf(y ~ x + (1 | gr(pair, cov = K))), family = gaussian(), data = dn,
      data2 = list(K = K)))
if (!is.null(fpi) && !is.null(fpk))
  cat("  logLik difference frm - drm:",
      format(as.numeric(logLik(fpk)) - as.numeric(logLik(fpi)), digits = 4),
      "\n")

# Known full sampling covariance: drmTMB meta_V(V = M) against frmtmb
# fcor(M), each checked against the hand-computed likelihood of
# y ~ N(Xb, M + s^2 I) and of y ~ N(Xb, s^2 M).
set.seed(9)
k <- 40
M <- diag(runif(k, 0.05, 0.2))
for (b in seq(1, k, 4)) M[b:(b + 3), b:(b + 3)][upper.tri(diag(4)) |
  lower.tri(diag(4))] <- 0.03
dmeta <- data.frame(x = rnorm(k))
dmeta$yi <- as.vector(0.2 + 0.3 * dmeta$x +
  t(chol(M + 0.04 * diag(k))) %*% rnorm(k))
fmv <- probe("drm: meta_V(V = M) with a full covariance",
  D(dbf(yi ~ x + meta_V(V = M), sigma ~ 1), family = gaussian(), data = dmeta))
ffc <- probe("frm: fcor(M)",
  frm(bf(yi ~ x + fcor(M)), family = gaussian(), data = dmeta,
      data2 = list(M = M)))
X <- cbind(1, dmeta$x)
ll_add <- function(p) mvtnorm::dmvnorm(dmeta$yi, X %*% p[1:2],
                                       M + exp(2 * p[3]) * diag(k), log = TRUE)
ll_scl <- function(p) mvtnorm::dmvnorm(dmeta$yi, X %*% p[1:2],
                                       exp(2 * p[3]) * M, log = TRUE)
o_add <- optim(c(0, 0, -1), function(p) -ll_add(p), method = "BFGS")
o_scl <- optim(c(0, 0, 0), function(p) -ll_scl(p), method = "BFGS")
cat("  hand ML, Cov = M + s^2 I:", format(-o_add$value, digits = 10),
    "\n  hand ML, Cov = s^2 M    :", format(-o_scl$value, digits = 10), "\n")
