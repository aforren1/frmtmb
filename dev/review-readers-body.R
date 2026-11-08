# Review of lane readers: items 5 and 8 through a nonlinear body.
# A: fitted() with allow_new_levels = TRUE on a body, both arms (seed 8201).
# B: item 5's first-order variance through a curved body, exp(a) with
#    a ~ me(xo, sx), against Monte Carlo over the draw (seed 41, 1e6 draws),
#    and against the same mean written as a log link (seed 8202).
# C: item 8's delta method with random effects in an nlpar, against draws
#    from the joint normal of fixed and random effects (seed 43, 2e5 draws).
source("dev/review-readers-setup.R")
tc <- function(expr) tryCatch(expr, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
  NULL
})

cat("== A. allow_new_levels = TRUE on a body\n")
set.seed(8201)
n <- 300; J <- 12
d <- data.frame(x = runif(n, 0, 3), g = factor(sample(J, n, TRUE)))
d$y <- (3 + rnorm(J, 0, 0.4)[d$g]) * exp(-0.8 * d$x) + rnorm(n, 0, 0.1)
fA0 <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE),
           family = gaussian(), data = d)
dA <- d
fA <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
          family = gaussian(), data = d)
nd <- data.frame(x = c(0.5, 2), g = d$g[1])
cat("no group term, allow_new_levels = TRUE:\n")
print(tc(fitted(fA0, newdata = nd, allow_new_levels = TRUE)))
cat("group term, seen level, allow_new_levels = TRUE:\n")
print(tc(fitted(fA, newdata = nd, allow_new_levels = TRUE)))
cat("group term, new level, allow_new_levels = TRUE:\n")
print(tc(fitted(fA, newdata = transform(nd, g = "new"),
                allow_new_levels = TRUE)))
cat("group term, seen level, default:\n")
print(tc(fitted(fA, newdata = nd)))
cat("conditional_effects, re_formula = NULL:\n")
print(tc(conditional_effects(fA, effects = "x", re_formula = NULL,
                             int_conditions = list(x = c(0.5, 2)))[[1]][,
  c("x", "estimate__", "se__", "lower__", "upper__")]))

if (Sys.getenv("REV_ARM") == "before") quit(save = "no")

cat("\n== B. exp(a), a ~ me(xo, sx), gaussian; seed 8202\n")
set.seed(8202)
n <- 600
d <- data.frame(l = rnorm(n), sx = 0.3)
d$xo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.2 + 0.8 * d$l) + rnorm(n, 0, 0.3)
fB <- frm(bf(y ~ exp(a), a ~ me(xo, sx), nl = TRUE), family = gaussian(),
          data = d)
fL <- frm(bf(y ~ me(xo, sx)), family = gaussian(link = "log"), data = d)
ba <- function(f, dp) {
  as.numeric(frm_linpred(f, newdata = data.frame(xo = 1, sx = 1), dpar = dp)) -
    as.numeric(frm_linpred(f, newdata = data.frame(xo = 0, sx = 1), dpar = dp))
}
b1 <- ba(fB, "a")
cat("slope of a on the latent:", b1, "\n")
set.seed(41)
z <- rnorm(1e6)
for (s in c(0.25, 0.5, 1, 1.5)) {
  nd <- data.frame(xo = 0.5, sx = s / b1)
  a0 <- as.numeric(frm_linpred(fB, newdata = nd, dpar = "a"))
  vB <- tc(fitted(fB, newdata = nd))
  ce <- tc(conditional_effects(fB, effects = "xo",
                               int_conditions = list(xo = 0.5),
                               conditions = data.frame(sx = s / b1))[[1]])
  vL <- fitted(fL, newdata = data.frame(xo = 0.5, sx = s / ba(fL, "mu")))
  mc <- exp(a0 + s * z)
  cat(sprintf(paste("link sd %.2f: body Est %.4f Err %.4f ce se__ %.4f |",
                    "MC mean %.4f sd %.4f | Est/MC %.3f Err/MC %.3f |",
                    "log-link fit Est %.4f Err %.4f\n"),
              s, vB[, "Estimate"], vB[, "Est.Error"], ce$se__, mean(mc),
              sd(mc), vB[, "Estimate"] / mean(mc), vB[, "Est.Error"] / sd(mc),
              vL[, "Estimate"], vL[, "Est.Error"]))
  cat(sprintf("   plug-in exp(a0) %.4f, ce estimate__ %.4f lower__ %.4f upper__ %.4f, MC 2.5/97.5%% %.4f %.4f, fitted Q %.4f %.4f
", exp(a0), ce$estimate__, ce$lower__, ce$upper__, quantile(mc, 0.025), quantile(mc, 0.975), vB[, "Q2.5"], vB[, "Q97.5"]))
}

cat("\n== C. a * exp(-b x), a ~ 1 + (1 | g): Jacobian over random effects\n")
jc <- frm_joint_cov(fA)
lab <- jc$labels
print(lab)
est <- c(fA$estimates$beta, fA$estimates$b)
cat("length est", length(est), " dim V", dim(jc$V), "\n")
d <- dA
nd <- data.frame(x = c(0.5, 2), g = d$g[c(1, 1)])
for (rf in list(NULL, NA)) {
  ff <- fitted(fA, newdata = nd, re_formula = rf)
  cat("re_formula", if (is.null(rf)) "NULL" else "NA", ": Est.Error",
      format(ff[, "Est.Error"], digits = 6), "\n")
}
ia <- grep("^a_Intercept$|b_a_Intercept|^a\\b", lab)
ib <- grep("^b_Intercept$|b_b_Intercept|^b\\b", lab)
ig <- grep(paste0("\\[", d$g[1], ",|\\[", d$g[1], "\\]"), lab)
cat("index a", ia, "b", ib, "g level", ig, ":", lab[c(ia, ib, ig)], "\n")
# Monte Carlo from the joint normal of the coefficients the rows load, read
# off frm_lp_basis(): A, V and coef_names. The body is evaluated by hand.
lb <- frm_lp_basis(fA, newdata = nd)
cat("coef_names:", lb$coef_names, "\n")
cn <- lb$coef_names
ca <- which(cn == "beta.a_(Intercept)")
cb <- which(cn == "beta.b_(Intercept)")
cg <- which(cn == paste0("b.a: 1 | g.", d$g[1], ".(Intercept)"))
cat("picked:", cn[c(ca, cb, cg)], "\n")
if (length(ca) == 1L && length(cb) == 1L && length(cg) == 1L) {
  jcf <- frm_joint_cov(fA)
  pos <- lb$coef_pos
  mu0 <- c(fA$estimates$beta, fA$estimates$b)[pos]
  if (length(mu0) != length(pos) || anyNA(mu0)) {
    mu0 <- NULL
  }
  # the centre is read from the body itself, so the estimate vector's
  # layout does not matter: a + u and b at the estimate
  fe <- fixef(fA)[, "Estimate"]
  re <- ranef(fA)$g
  ug <- if (length(dim(re)) == 3L) re[as.character(d$g[1]), "Estimate", 1L]
        else re[as.character(d$g[1]), 1L]
  ctr <- numeric(length(pos))
  ctr[ca] <- fe[grep("^a_", names(fe))]
  ctr[cb] <- fe[grep("^b_", names(fe))]
  ctr[cg] <- ug
  set.seed(43)
  Z <- MASS::mvrnorm(2e5, ctr, lb$V)
  sims <- vapply(nd$x, function(x) (Z[, ca] + Z[, cg]) * exp(-Z[, cb] * x),
                 numeric(2e5))
  ff <- fitted(fA, newdata = nd)
  cat("body at centre vs fitted:", ctr[ca] + ctr[cg], "*exp(-", ctr[cb],
      "x):", (ctr[ca] + ctr[cg]) * exp(-ctr[cb] * nd$x), "|",
      ff[, "Estimate"], "\n")
  cat("Est.Error", format(ff[, "Est.Error"], digits = 6), " MC sd",
      format(apply(sims, 2, sd), digits = 6), " ratio",
      format(ff[, "Est.Error"] / apply(sims, 2, sd), digits = 5),
      " (MC se of an sd ~", format(1 / sqrt(2 * 2e5), digits = 2), ")\n")
  # without the random effect's column: what a Jacobian over fixed effects
  # alone would read
  sims0 <- vapply(nd$x, function(x) {
    (Z[, ca] + ctr[cg]) * exp(-Z[, cb] * x)
  }, numeric(2e5))
  cat("MC sd with u_g held:", format(apply(sims0, 2, sd), digits = 6), "\n")
}
