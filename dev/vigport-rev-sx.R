# Reviewer: what frmtmb's `sx_1` is, against brms's.
#
#   Rscript dev/vigport-rev-sx.R
#
# 1. The factor between the two columns over several seeds: constant or
#    data-dependent, and its sign.
# 2. Why: brms builds Xs with smoothCon(diagonal.penalty = TRUE) and no
#    modCon; frmtmb with smoothCon(modCon = 3), default diagonal.penalty.
#    Rebuild both with mgcv alone and compare to standata()$Xs.
# 3. Whether the random part (and so sds()) is on the same scale.
# 4. The brms_distreg fit_smooth1 data under the lane's port_seed.
# 5. A prior on coef sx1_1: what it targets; frmtmb's variable names.
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), " brms",
    as.character(packageVersion("brms")), " mgcv",
    as.character(packageVersion("mgcv")), "\n")
cols <- function(dat, form = y ~ s(x1) + s(x2)) {
  f <- frm(bf(form), data = dat)
  sd <- brms::standata(brms::bf(form), data = dat)
  list(f = f, X = f$frame$linpreds[[1]]$X, Xs = sd$Xs, sd = sd)
}
cat("\n## 1. factor frmtmb column / brms column, by seed\n")
for (seed in 1:6) {
  set.seed(seed)
  dat <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
  cc <- cols(dat)
  out <- character()
  for (v in c("x1", "x2")) {
    a <- cc$X[, paste0("s(", v, ").fx1")]
    b <- cc$Xs[, paste0("s", v, "_1")]
    k <- stats::coef(stats::lm(a ~ 0 + b))
    res <- max(abs(a - k * b))
    out <- c(out, sprintf("%s: a = %+.4f * b (max resid %.1e)", v, k, res))
  }
  cat("seed", seed, ":", paste(out, collapse = "; "), "\n")
}

cat("\n## 2. rebuilt with mgcv alone, seed 1\n")
set.seed(1)
dat <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
cc <- cols(dat)
for (v in c("x1", "x2")) {
  spec <- eval(parse(text = sprintf("mgcv::s(%s)", v)))
  smA <- mgcv::smoothCon(spec, data = dat, absorb.cons = TRUE,
                         modCon = 3)[[1]]
  smB <- mgcv::smoothCon(spec, data = dat, absorb.cons = TRUE,
                         diagonal.penalty = TRUE)[[1]]
  rA <- mgcv::smooth2random(smA, names(dat), type = 2)
  rB <- mgcv::smooth2random(smB, names(dat), type = 2)
  cat(sprintf(paste0("%s: frmtmb route == frmtmb column: %s; ",
                     "brms route == standata Xs: %s\n"), v,
              isTRUE(all.equal(unname(rA$Xf[, 1]),
                               unname(cc$X[, paste0("s(", v, ").fx1")]))),
              isTRUE(all.equal(unname(rB$Xf[, 1]),
                               unname(cc$Xs[, paste0("s", v, "_1")])))))
  # the random part: compare Z Z' (rotation-free) between the routes
  ZA <- do.call(cbind, rA$rand); ZB <- do.call(cbind, rB$rand)
  zb <- cc$sd[[paste0("Zs_", match(v, c("x1", "x2")), "_1")]]
  cat(sprintf(paste0("   random part: ncol %d vs %d; ",
                     "max|ZA ZA' - ZB ZB'| / max|ZB ZB'| = %.2e; ",
                     "brms Zs == ZB: %s\n"), ncol(ZA), ncol(ZB),
              max(abs(tcrossprod(ZA) - tcrossprod(ZB))) /
                max(abs(tcrossprod(ZB))),
              isTRUE(all.equal(unname(zb), unname(ZB)))))
}

cat("\n## 3. fixed-part coefficients: frmtmb and lm on each basis\n")
f <- cc$f
fx <- fixef(f)
print(round(fx[, 1:2], 4))
# Same fitted values, so a coefficient on frmtmb's column times the
# factor k is the coefficient on brms's column.
for (v in c("x1", "x2")) {
  a <- cc$X[, paste0("s(", v, ").fx1")]
  b <- cc$Xs[, paste0("s", v, "_1")]
  k <- stats::coef(stats::lm(a ~ 0 + b))
  e <- fx[paste0("s", v, "_1"), "Estimate"]
  cat(sprintf("%s: frmtmb sx_1 %.4f; on brms's column that is %.4f\n",
              v, e, e * k))
}

cat("\n## 4. a prior on coef sx1_1, normal(0, 0.1), on frmtmb\n")
pr <- set_prior("normal(0, 0.1)", class = "b", coef = "sx1_1")
fp <- frm(bf(y ~ s(x1) + s(x2)), data = dat, prior = pr)
a <- cc$X[, "s(x1).fx1"]; b <- cc$Xs[, "sx1_1"]
k <- stats::coef(stats::lm(a ~ 0 + b))
e <- fixef(fp)["sx1_1", "Estimate"]
cat(sprintf(paste0("penalized frmtmb sx1_1 = %.4f, i.e. %.4f on brms's ",
                   "column; brms's prior sd 0.1 on its column is sd %.3f ",
                   "on frmtmb's\n"), e, e * k, 0.1 / abs(k)))
cat("brms prior table coef names:",
    paste(unique(brms::default_prior(brms::bf(y ~ s(x1) + s(x2)),
                                     data = dat)$coef), collapse = " "),
    "\n")

cat("\n## 5. frmtmb's variable names for the smooth's fixed part\n")
cat(grep("sx", variables(f), value = TRUE), "\n")

cat("\n## 6. brms_distreg fit_smooth1 data, port_seed of its expression\n")
source("C:/Users/adf44/source/r/frmtmb-wt-vigport/dev/brms-port/shim.R")
set.seed(port_seed("brms_distreg.11.1"))
ds <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
form <- y ~ s(x1) + s(x2) + (1 | fac)
fd <- frm(bf(form, sigma ~ s(x0) + (1 | fac)), data = ds)
sdd <- brms::standata(brms::bf(form, sigma ~ s(x0) + (1 | fac)), data = ds)
Xd <- fd$frame$linpreds[[1]]$X
for (v in c("x1", "x2")) {
  a <- Xd[, paste0("s(", v, ").fx1")]
  b <- sdd$Xs[, paste0("s", v, "_1")]
  k <- stats::coef(stats::lm(a ~ 0 + b))
  e <- fixef(fd)[paste0("s", v, "_1"), "Estimate"]
  cat(sprintf("%s: a = %+.4f * b; frmtmb %.3f -> brms column %.3f\n",
              v, k, e, e * k))
}
