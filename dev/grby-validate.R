# gr(g, by = f): agreement with independent references.
#
# 1. brms's own density, written from its Stan code (make_stancode()
#    and make_standata() of brms 2.23.0, read in dev/grby-log/
#    brms-code.txt): r_j = diag(sd[, Jby[j]]) L[Jby[j]] z_j, so for a
#    gaussian response the marginal law is y ~ N(X beta, Z S Z' +
#    sigma^2 I) with S block diagonal, level j taking the covariance of
#    its by-level. Evaluated with mvtnorm at frmtmb's estimates, on
#    brms's own standata (J_1, Jby_1, Z_1_*), it must equal logLik():
#    the Laplace approximation is exact for a gaussian response.
# 2. lme4 and glmmTMB at their ML optimum, fitting one term per
#    by-level with indicator columns, (0 + fa + fa:x | g) + (0 + fb +
#    fb:x | g), which is the same model.
#
# Run: Rscript dev/grby-validate.R > dev/grby-log/validate.txt
.libPaths(c("/opt/rlib/lane-grby", "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages({
  library(lme4)
  library(glmmTMB)
  library(brms)
  library(frmtmb)   # last, so bf() and VarCorr() are frmtmb's
})
fmt <- function(x) formatC(x, digits = 10, format = "g")
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
cat("brms", format(packageVersion("brms")), " lme4",
    format(packageVersion("lme4")), " glmmTMB",
    format(packageVersion("glmmTMB")), "\n\n")

set.seed(11)
ng <- 45
d <- data.frame(g = factor(rep(1:ng, each = 8)))
d$f <- factor(c("a", "b", "c")[(as.integer(d$g) - 1L) %% 3L + 1L])
d$x <- rnorm(nrow(d))
sd0 <- rbind(a = c(1.0, 0.2), b = c(0.3, 0.8), c = c(0.6, 0.5))
rho <- c(a = 0.4, b = -0.5, c = 0)
u <- t(vapply(seq_len(ng), function(j) {
  k <- as.character(d$f[match(j, as.integer(d$g))])
  S <- diag(sd0[k, ]) %*% matrix(c(1, rho[k], rho[k], 1), 2) %*%
    diag(sd0[k, ])
  as.vector(t(chol(S)) %*% rnorm(2))
}, numeric(2)))
d$y <- 1 + 0.5 * d$x + u[d$g, 1] + u[d$g, 2] * d$x + rnorm(nrow(d), 0, 0.5)
for (k in c("a", "b", "c")) {
  d[[paste0("f", k)]] <- as.numeric(d$f == k)
  d[[paste0("x", k)]] <- d[[paste0("f", k)]] * d$x
}
cat("data: seed 11,", nrow(d), "rows,", ng, "levels of g in 3 levels of f\n\n")

# brms's density, independent of frmtmb, at frmtmb's estimates
brms_marginal_ll <- function(fit, form, d) {
  sdat <- brms::make_standata(form, data = d)
  vc <- VarCorr(fit)
  rn <- rownames(vc$g$sd)
  cvm <- if (is.null(vc$g$cov)) {
    diag(vc$g$sd[, "Estimate"]^2, length(rn))
  } else {
    vc$g$cov[, "Estimate", ]
  }
  byl <- paste0("f", c("a", "b", "c"))
  M <- sdat$M_1
  Zc <- sapply(seq_len(M), function(m) sdat[[paste0("Z_1_", m)]])
  Zc <- matrix(Zc, ncol = M)
  n <- sdat$N
  V <- diag(vc$residual__$sd[, "Estimate"]^2, n)
  lev <- sdat$J_1
  by <- sdat$Jby_1[lev]
  for (k in seq_len(sdat$Nby_1)) {
    cols <- grep(paste0(":", byl[k], "$"), rn)
    S <- cvm[cols, cols, drop = FALSE]
    for (j in unique(lev[by == k])) {
      i <- which(lev == j)
      V[i, i] <- V[i, i] + Zc[i, , drop = FALSE] %*% S %*%
        t(Zc[i, , drop = FALSE])
    }
  }
  X <- sdat$X
  b <- fixef(fit)[, "Estimate"]
  mvtnorm::dmvnorm(as.vector(sdat$Y), drop(X %*% b), V, log = TRUE)
}

## 1. random intercept by f -------------------------------------------
cat("== (1 | gr(g, by = f))\n")
f1 <- frm(bf(y ~ x + (1 | gr(g, by = f))) + gaussian(), data = d)
l1 <- lmer(y ~ x + (0 + fa | g) + (0 + fb | g) + (0 + fc | g), data = d,
           REML = FALSE)
t1 <- glmmTMB(y ~ x + (0 + fa | g) + (0 + fb | g) + (0 + fc | g), data = d)
bl1 <- brms_marginal_ll(f1, y ~ x + (1 | gr(g, by = f)), d)
cat("logLik frmtmb ", fmt(as.numeric(logLik(f1))), "\n")
cat("logLik brms density at frmtmb's estimates ", fmt(bl1), " diff ",
    fmt(as.numeric(logLik(f1)) - bl1), "\n")
cat("logLik lme4   ", fmt(as.numeric(logLik(l1))), "\n")
cat("logLik glmmTMB", fmt(as.numeric(logLik(t1))), "\n")
sd_f <- VarCorr(f1)$g$sd[, "Estimate"]
sd_l <- sapply(VarCorr(l1)[c("g", "g.1", "g.2")], function(v) attr(v, "stddev"))
cat("sd frmtmb ", fmt(sd_f), "\n")
cat("sd lme4   ", fmt(unname(sd_l)), "\n")
cat("max |rel diff| sd frmtmb vs lme4", fmt(max(abs(sd_f / sd_l - 1))), "\n\n")

## 2. intercept + slope by f ------------------------------------------
cat("== (1 + x | gr(g, by = f))\n")
f2 <- frm(bf(y ~ x + (1 + x | gr(g, by = f))) + gaussian(), data = d)
l2 <- lmer(y ~ x + (0 + fa + xa | g) + (0 + fb + xb | g) + (0 + fc + xc | g),
           data = d, REML = FALSE,
           control = lmerControl(check.conv.singular = "ignore"))
t2 <- glmmTMB(y ~ x + (0 + fa + xa | g) + (0 + fb + xb | g) +
                (0 + fc + xc | g), data = d)
bl2 <- brms_marginal_ll(f2, y ~ x + (1 + x | gr(g, by = f)), d)
cat("logLik frmtmb ", fmt(as.numeric(logLik(f2))), "\n")
cat("logLik brms density at frmtmb's estimates ", fmt(bl2), " diff ",
    fmt(as.numeric(logLik(f2)) - bl2), "\n")
cat("logLik lme4   ", fmt(as.numeric(logLik(l2))), "\n")
cat("logLik glmmTMB", fmt(as.numeric(logLik(t2))), "\n")
vc2 <- VarCorr(f2)$g
sd_f <- vc2$sd[, "Estimate"]
vl <- VarCorr(l2)[c("g", "g.1", "g.2")]
sd_l <- unlist(lapply(vl, function(v) attr(v, "stddev")))
cor_f <- c(vc2$cor["Intercept:fa", "Estimate", "x:fa"],
           vc2$cor["Intercept:fb", "Estimate", "x:fb"],
           vc2$cor["Intercept:fc", "Estimate", "x:fc"])
cor_l <- vapply(vl, function(v) attr(v, "correlation")[1, 2], 1)
cat("sd frmtmb ", fmt(sd_f), "\n")
cat("sd lme4   ", fmt(unname(sd_l)), "\n")
cat("cor frmtmb", fmt(cor_f), "\n")
cat("cor lme4  ", fmt(unname(cor_l)), "\n")
cat("max |rel diff| sd", fmt(max(abs(sd_f / sd_l - 1))),
    " max |diff| cor", fmt(max(abs(cor_f - cor_l))), "\n")
cat("fixef frmtmb", fmt(fixef(f2)[, "Estimate"]), "\n")
cat("fixef lme4  ", fmt(unname(lme4::fixef(l2))), "\n")
cat("variables:", paste(variables(f2), collapse = " "), "\n\n")

## 3. poisson, Laplace against glmmTMB's Laplace ----------------------
cat("== poisson (1 + x | gr(g, by = f))\n")
set.seed(12)
d$cnt <- rpois(nrow(d), exp(0.3 + 0.2 * d$x + u[d$g, 1] * 0.6 +
                             u[d$g, 2] * 0.4 * d$x))
f3 <- frm(bf(cnt ~ x + (1 + x | gr(g, by = f))) + poisson(), data = d)
t3 <- glmmTMB(cnt ~ x + (0 + fa + xa | g) + (0 + fb + xb | g) +
                (0 + fc + xc | g), data = d, family = poisson())
cat("logLik frmtmb ", fmt(as.numeric(logLik(f3))), "\n")
cat("logLik glmmTMB", fmt(as.numeric(logLik(t3))), "\n")
sd_f <- VarCorr(f3)$g$sd[, "Estimate"]
sd_t <- unlist(lapply(glmmTMB::VarCorr(t3)$cond, function(v) attr(v, "stddev")))
cat("sd frmtmb ", fmt(sd_f), "\n")
cat("sd glmmTMB", fmt(unname(sd_t)), "\n")
cat("max |rel diff| sd", fmt(max(abs(sd_f / sd_t - 1))), "\n\n")

## 4. a structured block per by-level: diag and ar1 --------------------
cat("== diag(1 + x | gr(g, by = f))\n")
f4 <- frm(bf(y ~ x + diag(1 + x | gr(g, by = f))) + gaussian(), data = d)
t4 <- glmmTMB(y ~ x + diag(0 + fa + xa | g) + diag(0 + fb + xb | g) +
                diag(0 + fc + xc | g), data = d)
cat("logLik frmtmb ", fmt(as.numeric(logLik(f4))), "\n")
cat("logLik glmmTMB", fmt(as.numeric(logLik(t4))), "\n\n")

cat("== ar1(0 + t | gr(g, by = f))\n")
set.seed(13)
da <- expand.grid(t = factor(1:4), g = factor(1:30))
da$f <- factor(ifelse(as.integer(da$g) <= 15, "a", "b"))
ua <- matrix(0, 30, 4)
for (j in 1:30) {
  r <- if (j <= 15) 0.7 else -0.2
  s <- if (j <= 15) 1 else 0.5
  ua[j, ] <- s * as.vector(t(chol(r^abs(outer(1:4, 1:4, "-")))) %*% rnorm(4))
}
da$y <- ua[cbind(as.integer(da$g), as.integer(da$t))] + rnorm(nrow(da), 0, 0.4)
da$fa <- as.numeric(da$f == "a")
da$fb <- as.numeric(da$f == "b")
f5 <- frm(bf(y ~ 1 + ar1(0 + t | gr(g, by = f))) + gaussian(), data = da)
t5 <- glmmTMB(y ~ 1 + ar1(0 + t:fa | g) + ar1(0 + t:fb | g), data = da)
cat("logLik frmtmb ", fmt(as.numeric(logLik(f5))), "\n")
cat("logLik glmmTMB", fmt(as.numeric(logLik(t5))), "\n")
th_f <- f5$estimates$theta
th_t <- t5$fit$par[names(t5$fit$par) == "theta"]
cat("theta frmtmb ", fmt(th_f), "\n")
cat("theta glmmTMB", fmt(unname(th_t)), "\n")
cat("max |diff| theta", fmt(max(abs(th_f - th_t))), "\n")

## 5. the other positional structures against glmmTMB -----------------
cat("\n== other structures, one glmmTMB term per by-level\n")
# glmmTMB's homtoep does not converge on this design (logLik NA), so
# homtoep is checked by factorization in section 6
for (cs in c("cs", "homcs", "toep", "hetar1")) {
  ff <- as.formula(sprintf("y ~ 1 + %s(0 + t | gr(g, by = f))", cs))
  ft <- as.formula(sprintf("y ~ 1 + %s(0 + t:fa | g) + %s(0 + t:fb | g)",
                           cs, cs))
  fr <- suppressWarnings(frm(bf(ff) + gaussian(), data = da))
  tm <- suppressWarnings(glmmTMB(ft, data = da))
  cat(sprintf("%-8s logLik frmtmb %s glmmTMB %s diff %s\n", cs,
              fmt(as.numeric(logLik(fr))), fmt(as.numeric(logLik(tm))),
              fmt(as.numeric(logLik(fr)) - as.numeric(logLik(tm)))))
}

## 6. structures glmmTMB cannot split: the likelihood factorizes -------
# With a fixed effect and a residual sd per by-level, nothing is shared
# across by-levels, so the by-split fit's logLik is the sum of the two
# subset fits' (an identity of the model; the residual is numerical).
cat("\n== factorization over by-levels: by-split fit vs two subset fits\n")
set.seed(14)
times <- c(0, 0.3, 0.4, 1.1)
dt <- expand.grid(k = 1:4, g = factor(1:40))
dt$f <- factor(ifelse(as.integer(dt$g) <= 20, "a", "b"))
dt$tim <- num_factor(times[dt$k])
dt$y <- rnorm(40, 0, ifelse(1:40 <= 20, 1, 0.4))[dt$g] +
  rnorm(nrow(dt), 0, 0.5)
dt$y <- dt$y + 0.6 * ave(rnorm(nrow(dt)), dt$g, FUN = cumsum)
specs <- list(
  homdiag = "homdiag(0 + factor(k) | %s)",
  homtoep = "homtoep(0 + factor(k) | %s)",
  ou = "ou(tim + 0 | %s)",
  exp = "exp(tim + 0 | %s)",
  gau = "gau(tim + 0 | %s)",
  mat = "mat(tim + 0 | %s)",
  us_student = "(1 | %s)"
)
for (nm in names(specs)) {
  grp_by <- if (nm == "us_student") {
    "gr(g, by = f, dist = \"student\")"
  } else "gr(g, by = f)"
  grp <- if (nm == "us_student") "gr(g, dist = \"student\")" else "g"
  fby <- as.formula(paste("y ~ 0 + f +", sprintf(specs[[nm]], grp_by)))
  fsub <- as.formula(paste("y ~ 1 +", sprintf(specs[[nm]], grp)))
  whole <- suppressWarnings(frm(bf(fby, sigma ~ 0 + f) + gaussian(),
                                data = dt))
  subs <- lapply(c("a", "b"), function(k) {
    suppressWarnings(frm(bf(fsub) + gaussian(),
                         data = droplevels(dt[dt$f == k, ])))
  })
  parts <- vapply(subs, function(s) as.numeric(logLik(s)), 1)
  # the by-split objective AT the subsets' estimates: the identity
  # itself, free of where either optimizer stopped
  est <- lapply(subs, `[[`, "estimates")
  p0 <- whole$estimates
  p0$beta[] <- vapply(est, function(e) e$beta[[1L]], 1)
  p0$betad[] <- vapply(est, function(e) e$betad[[1L]], 1)
  p0$theta[] <- unlist(lapply(est, `[[`, "theta"))
  at <- -whole$obj$fn(unlist(p0[names(whole$obj$par)[
    !duplicated(names(whole$obj$par))]]))
  cat(sprintf(paste0("%-10s by-split %s  sum of subsets %s  diff %s\n",
                     "           by-split objective at the subsets' ",
                     "estimates %s  diff %s\n"), nm,
              fmt(as.numeric(logLik(whole))), fmt(sum(parts)),
              fmt(as.numeric(logLik(whole)) - sum(parts)), fmt(at),
              fmt(at - sum(parts))))
}

## 7. mm(g1, g2, by = cbind(f1, f2)) against brms's density ------------
cat("\n== (1 | mm(g1, g2, by = cbind(f1, f2)))\n")
set.seed(3)
nm_ <- 400
dm <- data.frame(g1 = factor(sample(1:30, nm_, TRUE), levels = 1:30),
                 g2 = factor(sample(1:30, nm_, TRUE), levels = 1:30))
fl <- ifelse(1:30 <= 15, "a", "b")
dm$f1 <- factor(fl[dm$g1])
dm$f2 <- factor(fl[dm$g2])
dm$x <- rnorm(nm_)
um <- rnorm(30, 0, ifelse(1:30 <= 15, 1, 0.3))
dm$y <- 0.5 + 0.3 * dm$x + 0.5 * (um[dm$g1] + um[dm$g2]) +
  rnorm(nm_, 0, 0.5)
fm <- y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))
f7 <- frm(bf(fm) + gaussian(), data = dm)
sdat <- brms::make_standata(fm, data = dm)
vc <- VarCorr(f7)
s_by <- vc$mmg1g2$sd[, "Estimate"]
lev_sd <- s_by[sdat$Jby_1]
A <- matrix(0, nm_, sdat$N_1)
for (k in 1:2) {
  J <- sdat[[paste0("J_1_", k)]]
  w <- sdat[[paste0("W_1_", k)]] * sdat[[paste0("Z_1_1_", k)]]
  A[cbind(seq_len(nm_), J)] <- A[cbind(seq_len(nm_), J)] + w
}
V <- A %*% diag(lev_sd^2) %*% t(A) +
  diag(vc$residual__$sd[, "Estimate"]^2, nm_)
bl7 <- mvtnorm::dmvnorm(as.vector(sdat$Y),
                        drop(sdat$X %*% fixef(f7)[, "Estimate"]), V,
                        log = TRUE)
cat("brms Nby_1", sdat$Nby_1, " Jby_1 by-level of each pooled level:",
    sdat$Jby_1, "\n")
cat("logLik frmtmb ", fmt(as.numeric(logLik(f7))), "\n")
cat("logLik brms density at frmtmb's estimates ", fmt(bl7), "\n")
cat("diff", fmt(as.numeric(logLik(f7)) - bl7), "\n")
cat("variables:", paste(variables(f7), collapse = " | "), "\n")

## 8. quadrature: each level integrated alone --------------------------
cat("\n== bernoulli (1 | gr(g, by = f)), quadrature = TRUE\n")
set.seed(7)
db <- data.frame(g = factor(rep(1:40, each = 6)))
db$f <- factor(ifelse(as.integer(db$g) <= 20, "a", "b"))
db$x <- rnorm(nrow(db))
ub <- rnorm(40, 0, ifelse(1:40 <= 20, 1, 0.4))
db$yb <- rbinom(nrow(db), 1, plogis(0.2 + 0.5 * db$x + ub[db$g]))
f8 <- frm(bf(yb ~ x + (1 | gr(g, by = f))) + bernoulli(), data = db,
          quadrature = TRUE)
b8 <- fixef(f8)[, "Estimate"]
s8 <- VarCorr(f8)$g$sd[, "Estimate"]
eta <- drop(cbind(1, db$x) %*% b8)
ll <- 0
for (j in levels(db$g)) {
  i <- which(db$g == j)
  s <- s8[as.character(db$f[i[1]]) == c("a", "b")]
  ll <- ll + log(integrate(function(u) {
    vapply(u, function(v) {
      prod(dbinom(db$yb[i], 1, plogis(eta[i] + v)))
    }, 1) * dnorm(u, 0, s)
  }, -Inf, Inf, rel.tol = 1e-12)$value)
}
cat("logLik frmtmb quadrature", fmt(as.numeric(logLik(f8))), "\n")
cat("logLik integrate() at frmtmb's estimates", fmt(ll), "\n")
cat("diff", fmt(as.numeric(logLik(f8)) - ll), "\n")

## 9. an |ID| key across two responses, per by-level -------------------
cat("\n== mvbf(bf(y1 ~ 1 + (1 | p | gr(g, by = f))), bf(y2 ~ ...))\n")
set.seed(21)
d9 <- data.frame(g = factor(rep(1:40, each = 6)))
d9$f <- factor(ifelse(as.integer(d9$g) <= 20, "a", "b"))
u1 <- rnorm(40, 0, ifelse(1:40 <= 20, 1, 0.3))
u2 <- 0.6 * u1 + rnorm(40, 0, 0.4)
d9$y1 <- u1[d9$g] + rnorm(nrow(d9), 0, 0.5)
d9$y2 <- u2[d9$g] + rnorm(nrow(d9), 0, 0.5)
f9 <- frm(mvbf(bf(y1 ~ 1 + (1 | p | gr(g, by = f))),
               bf(y2 ~ 1 + (1 | p | gr(g, by = f)))) + gaussian(),
          data = d9)
dl <- rbind(data.frame(v = d9$y1, trait = "y1", g = d9$g, f = d9$f),
            data.frame(v = d9$y2, trait = "y2", g = d9$g, f = d9$f))
dl$trait <- factor(dl$trait)
for (tr in c("1", "2")) {
  for (k in c("a", "b")) {
    dl[[paste0("t", tr, k)]] <- as.numeric(dl$trait == paste0("y", tr) &
                                             dl$f == k)
  }
}
t9 <- glmmTMB(v ~ 0 + trait + (0 + t1a + t2a | g) + (0 + t1b + t2b | g),
              dispformula = ~ 0 + trait, data = dl)
cat("logLik frmtmb ", fmt(as.numeric(logLik(f9))), "\n")
cat("logLik glmmTMB, long format with trait x by-level indicators",
    fmt(as.numeric(logLik(t9))), "\n")
cat("diff", fmt(as.numeric(logLik(f9)) - as.numeric(logLik(t9))), "\n")
cat("variables:", paste(variables(f9), collapse = " | "), "\n")

## 10. REML against lme4's REML ------------------------------------------
cat("\n== REML (1 | gr(g, by = f))\n")
f10 <- frm(bf(y ~ x + (1 | gr(g, by = f))) + gaussian(), data = d,
           REML = TRUE)
l10 <- lmer(y ~ x + (0 + fa | g) + (0 + fb | g) + (0 + fc | g), data = d,
            REML = TRUE)
cat("REML logLik frmtmb", fmt(as.numeric(logLik(f10))), " lme4",
    fmt(as.numeric(logLik(l10))), " diff",
    fmt(as.numeric(logLik(f10)) - as.numeric(logLik(l10))), "\n")
sd_f <- VarCorr(f10)$g$sd[, "Estimate"]
sd_l <- sapply(VarCorr(l10)[c("g", "g.1", "g.2")],
               function(v) attr(v, "stddev"))
cat("sd frmtmb", fmt(sd_f), "\nsd lme4  ", fmt(unname(sd_l)), "\n")
cat("max |rel diff| sd", fmt(max(abs(sd_f / sd_l - 1))), "\n")
