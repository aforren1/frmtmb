# Reviewer probes: dirichlet_multinomial() through the post-fit paths,
# side by side with multinomial(). Data seed 2027, n = 300, 15 groups.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
dm_sim <- function(seed, n = 300, phi = 6, re = FALSE, Nlo = 5, Nhi = 30) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(15, n / 15),
                   f = factor(sample(c("u", "v"), n, TRUE)))
  u <- if (re) rnorm(15, 0, 0.5)[dd$g] else 0
  E <- cbind(0, 0.4 + 0.6 * dd$x + u, -0.3 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$N <- sample(Nlo:Nhi, n, TRUE)
  Y <- t(vapply(seq_len(n), function(i) {
    p <- rgamma(3, P[i, ] * phi)
    as.numeric(rmultinom(1, dd$N[i], p / sum(p)))
  }, numeric(3)))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}
dd <- dm_sim(2027)
DM <- dirichlet_multinomial()
fd <- frm(Y | trials(N) ~ x, family = DM, data = dd)
fm <- frm(Y | trials(N) ~ x, family = multinomial(), data = dd)
probe("summary dm", summary(fd))
probe("confint dm", confint(fd))
probe("fitted dims dm vs mn", list(dim(fitted(fd)), dim(fitted(fm)),
                                   dimnames(fitted(fd))[[3]],
                                   dimnames(fitted(fm))[[3]]))
probe("frm_linpred dm dims", dim(frm_linpred(fd)))
probe("frm_linpred mn dims", dim(frm_linpred(fm)))
probe("fitted dpar phi head", head(fitted(fd, dpar = "phi")))
nd <- data.frame(x = c(-1, 0, 2), N = c(2L, 50L, 1000L), f = "u")
probe("predict newdata trials", {
  set.seed(1)
  p <- predict(fd, newdata = nd, ndraws = 2000)
  list(dim = dim(p), est = p[, "Estimate", ],
       rowsum = rowSums(p[, "Estimate", ]))
})
probe("predict newdata trials, multinomial", {
  set.seed(1)
  p <- predict(fm, newdata = nd, ndraws = 2000)
  list(dim = dim(p), est = p[, "Estimate", ])
})
probe("predict newdata: dm spread exceeds mn at N = 1000 (sd of cat a)", {
  set.seed(2)
  pd <- predict(fd, newdata = nd, ndraws = 4000)
  pm <- predict(fm, newdata = nd, ndraws = 4000)
  c(dm = pd[3, "Est.Error", 1], mn = pm[3, "Est.Error", 1])
})
probe("simulate newdata trials", {
  s <- simulate(fd, nsim = 2, seed = 3, newdata = nd)
  lapply(s, rowSums)
})
probe("simulate newdata trials missing", simulate(fd, nsim = 1, seed = 3,
                                                  newdata = nd[, "x",
                                                               drop = FALSE]))
probe("conditional_effects dm", {
  ce <- conditional_effects(fd, effects = "x")
  list(names = names(ce), head = head(ce[[1]], 3))
})
probe("conditional_effects mn", {
  ce <- conditional_effects(fm, effects = "x")
  list(names = names(ce), head = head(ce[[1]], 3))
})
probe("emmeans dm", summary(emmeans::emmeans(fd, ~ 1)))
probe("emmeans mn", summary(emmeans::emmeans(fm, ~ 1)))
probe("marginaleffects dm avg_predictions",
      marginaleffects::avg_predictions(fd))
probe("marginaleffects mn avg_predictions",
      marginaleffects::avg_predictions(fm))
probe("anova mn vs dm (nested at phi -> Inf)", anova(fm, fd))
probe("drop1 dm", drop1(fd))
probe("refcat = 'b' logLik", c(logLik(fd), logLik(frm(Y | trials(N) ~ x,
  family = dirichlet_multinomial(refcat = "b"), data = dd))))
probe("refcat = NA (all categories predicted)", {
  m <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(refcat = NA),
           data = dd)
  list(ll = logLik(m), conv = m$opt$convergence,
       se = sqrt(diag(vcov(m))))
})
probe("refcat = NA, multinomial", {
  m <- frm(Y | trials(N) ~ x, family = multinomial(refcat = NA), data = dd)
  list(ll = logLik(m), conv = m$opt$convergence,
       se = sqrt(diag(vcov(m))))
})
probe("phi ~ z", {
  m <- frm(bf(Y | trials(N) ~ x, phi ~ z), family = DM, data = dd)
  list(ll = logLik(m), phi = fixef_by_dpar(m)$phi)
})
probe("trials(10) constant on equal-total data", {
  d2 <- dm_sim(5, Nlo = 10, Nhi = 10)
  logLik(frm(Y | trials(10) ~ x, family = DM, data = d2))
})
probe("weights", {
  dd$wt <- runif(nrow(dd), 0.5, 2)
  m <- frm(Y | trials(N) + weights(wt) ~ x, family = DM, data = dd)
  logLik(m)
})
probe("zero-trial rows", {
  d2 <- dd
  d2$Y[1:5, ] <- 0
  d2$N[1:5] <- 0L
  m <- frm(Y | trials(N) ~ x, family = DM, data = d2)
  c(ll = logLik(m), ll_drop = logLik(frm(Y | trials(N) ~ x, family = DM,
                                          data = d2[-(1:5), ])))
})
probe("a category with no counts", {
  d2 <- dd
  d2$N <- d2$N - d2$Y[, "c"]
  d2$Y[, "c"] <- 0
  m <- frm(Y | trials(N) ~ x, family = DM, data = d2)
  mm <- frm(Y | trials(N) ~ x, family = multinomial(), data = d2)
  list(dm_conv = m$opt$convergence, dm_ll = logLik(m),
       mn_conv = mm$opt$convergence, mn_ll = logLik(mm),
       dm_b = fixef_by_dpar(m))
})
probe("mvbf(dm, gaussian)", {
  dd$w <- 0.5 * dd$x + rnorm(nrow(dd))
  m <- frm(mvbf(bf(Y | trials(N) ~ x, family = DM),
                bf(w ~ x, family = gaussian())), data = dd)
  c(logLik(m), logLik(fd) + logLik(lm(w ~ x, dd)))
})
probe("mixture(dm, multinomial) refused",
      frm(Y | trials(N) ~ x, family = mixture(DM, multinomial()), data = dd))
probe("cens() refused", {
  dd$cc <- 0
  frm(Y | trials(N) + cens(cc) ~ x, family = DM, data = dd)
})
probe("group term VarCorr", {
  d3 <- dm_sim(54, re = TRUE)
  m <- frm(Y | trials(N) ~ x + (1 | g), family = DM, data = d3)
  VarCorr(m)
})
probe("K = 2 dm fit vs beta_binomial fit", {
  d4 <- dm_sim(9)
  d4$Y2 <- cbind(a = d4$N - d4$Y[, "b"], b = d4$Y[, "b"])
  d4$yb <- d4$Y[, "b"]
  m2 <- frm(Y2 | trials(N) ~ x, family = DM, data = d4)
  bb <- frm(yb | trials(N) ~ x, family = beta_binomial(), data = d4)
  c(dm = logLik(m2), bb = logLik(bb), phi_dm = exp(fixef_by_dpar(m2)$phi),
    phi_bb = exp(fixef_by_dpar(bb)$phi))
})
# phi runs off: what the user is told
set.seed(4)
n <- 300
d <- data.frame(x = rnorm(n))
E <- cbind(0, 0.3 + 0.5 * d$x, -0.2 - 0.4 * d$x)
P <- exp(E) / rowSums(exp(E))
d$N <- sample(5:40, n, TRUE)
d$Y <- t(vapply(seq_len(n), function(i) {
  as.numeric(rmultinom(1, d$N[i], P[i, ]))
}, numeric(3)))
colnames(d$Y) <- c("a", "b", "c")
probe("multinomial data: dm fit, messages and summary", {
  m <- frm(Y | trials(N) ~ x, family = DM, data = d)
  print(summary(m))
  list(conv = m$opt$convergence, msg = m$opt$message,
       logphi = fixef_by_dpar(m)$phi, se = sqrt(diag(vcov(m))))
})
probe("multinomial data: confint", confint(frm(Y | trials(N) ~ x, family = DM,
                                               data = d)))
probe("frm_compat dm rows", {
  tb <- frm_compat("dirichlet_multinomial")
  mb <- frm_compat("multinomial")
  key <- function(t) paste(t$feature_a, t$feature_b)
  k1 <- sub("^dirichlet_multinomial ", "", key(tb))
  k2 <- sub("^multinomial ", "", key(mb))
  m <- merge(data.frame(k = k1, dm = tb$status),
             data.frame(k = k2, mn = mb$status), all = TRUE)
  m[is.na(m$dm) | is.na(m$mn) | m$dm != m$mn, ]
})
probe("frm_compat nb2 vs nb differing rows", {
  tb <- frm_compat("negbinomial2")
  mb <- frm_compat("negbinomial")
  key <- function(t) paste(t$feature_a, t$feature_b)
  k1 <- sub("^negbinomial2 ", "", key(tb))
  k2 <- sub("^negbinomial ", "", key(mb))
  m <- merge(data.frame(k = k1, nb2 = tb$status),
             data.frame(k = k2, nb = mb$status), all = TRUE)
  m[is.na(m$nb2) | is.na(m$nb) | m$nb2 != m$nb, ]
})
