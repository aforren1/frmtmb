# Lane postfit2: every new or changed function against brms 2.23.0, on
# compiled brms fits of the same models AT THE SAME PARAMETERS.
#
#   Rscript dev/postfit2-brms-compare.R > dev/postfit2-log/brms-compare.txt
#
# brms is sampled with algorithm = "fixed_param" and one chain per
# parameter vector, each chain initialized at a vector mapped from
# frmtmb (the ML estimate, or one frm_sample() draw), so each brms draw
# IS that vector. A smooth's coefficients are mapped by projecting
# frmtmb's term, evaluated at the data, on brms's columns for the same
# mgcv basis; the residual of that projection is printed and must be at
# rounding level, or the two bases differ and nothing below would mean
# anything. Seeds: data 1, frm_sample 20260929, posterior_average 5.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/pf2rev-stan-cache")
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
say <- function(...) cat(sprintf(...), "\n", sep = "")
maxdiff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
# rstan reads a length-1 R vector as a scalar, so a Stan vector needs a
# 1-d array
arr <- function(v) array(v, dim = length(v))
strip <- function(d) {
  d <- as.data.frame(d)
  rownames(d) <- NULL
  d
}

## ---- A. make_conditions() and update_adterms(), no fit needed ------
say("== A. make_conditions(), update_adterms(): identical() to brms")
ep <- brms::epilepsy
ep$fz <- factor(ifelse(ep$zAge > 0, "old", "young"))
mc_cases <- list(
  list(ep, c("zBase", "zAge")),
  list(ep, c("zBase", "Trt")),
  list(ep, c("fz", "zAge"), digits = 1),
  list(ep, "Trt", incl_vars = FALSE),
  list(ep, c("zAge", "zBase"), sep = " | ")
)
for (cs in mc_cases) {
  a <- do.call(brms::make_conditions, cs)
  b <- do.call(frmtmb::make_conditions, cs)
  say("make_conditions(%s): identical %s, dim %s",
      paste(deparse(cs[-1]), collapse = ""), identical(a, b),
      paste(dim(b), collapse = "x"))
}
ua_cases <- list(
  list(y | trials(size) ~ x, ~ trials(10)),
  list(y | trials(size) ~ x, ~ weights(w)),
  list(y | trials(size) ~ x, ~ weights(w), action = "replace"),
  list(y ~ x, ~ trials(10)),
  list(y | se(s, sigma = TRUE) + weights(w) ~ x + (1 | g), ~ se(s2)),
  list(y | resp_se(s) ~ x, ~ se(s3) + cens(c)),
  list(y | cens(c, y2) ~ x, ~ weights(w), action = "replace"),
  list(y | weights(w) ~ x, ~ 1, action = "replace"),
  list(log(y) | trunc(lb = 0) ~ s(x), ~ trunc(ub = 5))
)
for (cs in ua_cases) {
  a <- do.call(brms::update_adterms, cs)
  b <- do.call(frmtmb::update_adterms, cs)
  say("update_adterms(%s, %s%s) = %s; identical %s",
      deparse1(cs[[1]]), deparse1(cs[[2]]),
      if (length(cs) > 2) paste0(", ", cs[[3]]) else "",
      deparse1(b), identical(a, b))
}

## ---- B. the models ---------------------------------------------------
set.seed(1)
n <- 150
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(2 * pi * d$x) + d$z * (d$f == "b") + rnorm(n, 0, 0.3)
d$y2 <- 2 * exp(0.8 * d$z) * (1 + 0.3 * (d$f == "b")) + rnorm(n, 0, 0.3)

# frmtmb's contribution of every smooth basis at the data rows, in
# sm_info order: list of (value, parts)
frm_sm_values <- function(fit, est = fit$estimates) {
  lp <- fit$frame$linpreds[[1]]
  cvec <- frmtmb:::coef_b(fit, est[["b"]])
  lapply(lp$smooths, function(si) {
    v <- 0
    for (k in si$block_ids) {
      bk <- fit$frame$re_blocks[[k]]
      v <- v + as.vector(lp$Z[, bk$c_idx, drop = FALSE] %*% cvec[bk$c_idx])
    }
    v + as.vector(as.matrix(lp$X[, si$xf_idx, drop = FALSE]) %*%
                    est[["beta"]][si$xf_idx])
  })
}

# brms init for one frmtmb parameter vector of a linear gaussian model
# with smooths: projection of each smooth basis on brms's columns
brms_init <- function(fit, sdata, est, sigma) {
  vals <- frm_sm_values(fit, est)
  lp <- fit$frame$linpreds[[1]]
  init <- list()
  ksx <- 0L
  bs <- numeric(0)
  resid <- 0
  for (k in seq_along(vals)) {
    si <- lp$smooths[[k]]
    zn <- grep(paste0("^Zs_", k, "_"), names(sdata), value = TRUE)
    Z <- do.call(cbind, sdata[zn])
    xs <- ksx + seq_len(si$nf)
    ksx <- ksx + si$nf
    M <- cbind(Z, sdata$Xs[, xs, drop = FALSE])
    cf <- qr.solve(M, vals[[k]])
    resid <- max(resid, maxdiff(M %*% cf, vals[[k]]) /
                   max(abs(vals[[k]])))
    pos <- 0L
    for (j in seq_along(zn)) {
      nj <- ncol(sdata[[zn[j]]])
      init[[paste0("zs_", k, "_", j)]] <- arr(cf[pos + seq_len(nj)])
      pos <- pos + nj
    }
    init[[paste0("sds_", k)]] <- arr(rep(1, length(zn)))
    bs <- c(bs, cf[pos + seq_len(si$nf)])
  }
  init$bs <- arr(bs)
  X <- sdata$X
  bcols <- setdiff(colnames(X), "Intercept")
  beta <- est[["beta"]]
  names(beta) <- colnames(lp$X)
  init$b <- arr(unname(beta[bcols]))
  init$Intercept <- unname(beta["(Intercept)"] +
                             sum(colMeans(X[, bcols, drop = FALSE]) *
                                   beta[bcols]))
  init$sigma <- sigma
  attr(init, "proj_resid") <- resid
  init
}

brm_fixed <- function(form, data, inits, family = gaussian()) {
  suppressMessages(suppressWarnings(
    brm(form, data = data, family = family, algorithm = "fixed_param",
        chains = length(inits), iter = 1, warmup = 0, init = inits,
        refresh = 0, seed = 1, backend = "rstan", silent = 2)
  ))
}

cmp_frames <- function(label, a, b, cols = c("estimate__", "se__",
                                              "lower__", "upper__")) {
  say("%s: keys brms [%s] frmtmb [%s]", label,
      paste(names(a), collapse = ", "), paste(names(b), collapse = ", "))
  for (k in names(a)) {
    fa <- strip(a[[k]])
    fb <- strip(b[[k]])
    gd <- vapply(setdiff(names(fa), c(cols, "cond__")), function(v) {
      if (is.numeric(fa[[v]])) maxdiff(fa[[v]], fb[[v]]) else
        as.numeric(!identical(as.character(fa[[v]]),
                              as.character(fb[[v]])))
    }, 1)
    say("  %s: rows %d / %d; columns identical %s; grid max diff %.3g%s",
        k, nrow(fa), nrow(fb), identical(names(fa), names(fb)), max(gd),
        if (any(gd > 1e-12)) {
          paste0(" (in ", paste(names(gd)[gd > 1e-12], collapse = ", "),
                 ": brms ", paste(signif(unique(fa[[names(gd)[gd > 1e-12][1]]]),
                                         4)[1:2], collapse = " "),
                 " frmtmb ", paste(signif(unique(fb[[names(gd)[gd > 1e-12][1]]]),
                                          4)[1:2], collapse = " "), ")")
        } else "")
    for (cc in intersect(cols, names(fa))) {
      say("    %-10s max |brms - frmtmb| = %.3g (scale %.3g)", cc,
          maxdiff(fa[[cc]], fb[[cc]]), max(abs(fa[[cc]])))
    }
    say("    cond__ equal %s; attr names brms [%s] frmtmb [%s]",
        identical(as.character(fa$cond__), as.character(fb$cond__)),
        paste(sort(setdiff(names(attributes(a[[k]])),
                           c("row.names", "names", "class"))),
              collapse = ","),
        paste(sort(setdiff(names(attributes(b[[k]])),
                           c("row.names", "names", "class"))),
              collapse = ","))
    for (at in c("effects", "response", "surface")) {
      say("    attr %-8s brms %s | frmtmb %s", at,
          paste(attr(a[[k]], at), collapse = " "),
          paste(attr(b[[k]], at), collapse = " "))
    }
  }
}

## ---- B1. conditional_smooths() at the ML estimate --------------------
say("\n== B1. conditional_smooths(): ML fit vs brms at the ML estimate")
f1 <- frm(bf(y ~ f + s(x) + s(z, by = f)), family = gaussian(), data = d)
sig1 <- frm_linpred(f1, dpar = "sigma", type = "response")[1]
form1 <- brms::bf(y ~ f + s(x) + s(z, by = f))
sd1 <- standata(form1, data = d)
i1 <- brms_init(f1, sd1, f1$estimates, sig1)
say("projection residual (relative): %.3g", attr(i1, "proj_resid"))
b1 <- brm_fixed(form1, d, list(i1))
cmp_frames("default", conditional_smooths(b1), conditional_smooths(f1))
cmp_frames("surface = FALSE, int_conditions",
           conditional_smooths(b1, surface = FALSE, smooths = "s(z, by = f)",
                               int_conditions = list(z = c(0.2, 0.5))),
           conditional_smooths(f1, surface = FALSE, smooths = "s(z, by = f)",
                               int_conditions = list(z = c(0.2, 0.5))))
say("brms s3: %s", tryCatch(conditional_smooths(b1, smooths = "s3"),
                           error = conditionMessage))
say("frmtmb s3: %s", tryCatch(conditional_smooths(f1, smooths = "s3"),
                             error = conditionMessage))

## ---- B2. t2() surface and too_far --------------------------------------
say("\n== B2. t2(x, z): surface and too_far")
f2 <- frm(bf(y ~ t2(x, z)), family = gaussian(), data = d)
sig2 <- frm_linpred(f2, dpar = "sigma", type = "response")[1]
form2 <- brms::bf(y ~ t2(x, z))
sd2 <- standata(form2, data = d)
i2 <- brms_init(f2, sd2, f2$estimates, sig2)
say("projection residual (relative): %.3g", attr(i2, "proj_resid"))
b2 <- brm_fixed(form2, d, list(i2))
for (tf in c(0, 0.1)) {
  cmp_frames(sprintf("t2 surface, resolution 20, too_far %g", tf),
             conditional_smooths(b2, resolution = 20, too_far = tf),
             conditional_smooths(f2, resolution = 20, too_far = tf))
}
cmp_frames("t2 surface = FALSE",
           conditional_smooths(b2, surface = FALSE, resolution = 20),
           conditional_smooths(f2, surface = FALSE, resolution = 20))
cmp_frames("conditional_effects x:z surface, too_far 0.2",
           conditional_effects(b2, "x:z", surface = TRUE, resolution = 15,
                               too_far = 0.2),
           conditional_effects(f2, "x:z", surface = TRUE, resolution = 15,
                               too_far = 0.2),
           cols = "estimate__")

## ---- B3. select_points -------------------------------------------------
say("\n== B3. select_points: the observations each keeps")
f3 <- frm(bf(y ~ f + x + z), family = gaussian(), data = d)
form3 <- brms::bf(y ~ f + x + z)
be3 <- fixef(f3)[, 1]
X3 <- standata(form3, data = d)$X
i3 <- list(b = arr(unname(be3[-1])),
           Intercept = unname(be3[1] + sum(colMeans(X3[, -1]) * be3[-1])),
           sigma = frm_linpred(f3, dpar = "sigma", type = "response")[1])
b3 <- brm_fixed(form3, d, list(i3))
for (sp in c(0, 0.05, 0.1, 0.3)) {
  pa <- attr(conditional_effects(b3, "f", select_points = sp)[[1]],
             "points")
  pb <- attr(conditional_effects(f3, "f", select_points = sp)[[1]],
             "points")
  key <- function(p) sort(paste(p$f, signif(p$resp__, 12)))
  # brms also measures the RESPONSE against the value its conditions
  # hold it at (its mean); frmtmb measures predictors only. Applying
  # brms's extra filter to frmtmb's rows must give brms's rows exactly.
  yh <- mean(d$y)
  uy <- abs((pb$resp__ - min(d$y)) / diff(range(d$y)) -
              (yh - min(d$y)) / diff(range(d$y)))
  pby <- if (sp > 0) pb[uy <= sp, , drop = FALSE] else pb
  say("select_points %.2f: brms %d rows, frmtmb %d rows, same rows %s; frmtmb plus brms's response filter %d rows, same as brms %s",
      sp, nrow(pa), nrow(pb), identical(key(pa), key(pb)), nrow(pby),
      identical(key(pa), key(pby)))
}
cea <- conditional_effects(b3, "x")
ceb <- conditional_effects(f3, "x")
say("estimate__ of x at the ML estimate: max diff %.3g",
    maxdiff(cea[[1]]$estimate__, ceb[[1]]$estimate__))

## ---- B4. the nonlinear Wald band ---------------------------------------
say("\n== B4. nonlinear predictor: estimate vs brms, band vs Monte Carlo")
f4 <- frm(bf(y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1, nl = TRUE),
          family = gaussian(), data = d)
form4 <- brms::bf(y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1, nl = TRUE)
lpa <- f4$frame$linpreds[["y2.a"]]
lpb <- f4$frame$linpreds[["y2.b"]]
i4 <- list(b_a = arr(unname(f4$estimates$beta[lpa$idx])),
           b_b = arr(unname(f4$estimates$beta[lpb$idx])),
           sigma = frm_linpred(f4, dpar = "sigma", type = "response")[1])
b4 <- brm_fixed(form4, d, list(i4),
                family = gaussian())
ce4b <- conditional_effects(b4, "z")
ce4f <- conditional_effects(f4, "z")
say("brms estimate__ vs frmtmb estimate__, z: max diff %.3g (scale %.3g)",
    maxdiff(ce4b[[1]]$estimate__, ce4f[[1]]$estimate__),
    max(abs(ce4f[[1]]$estimate__)))
ce4b2 <- conditional_effects(b4, "z:f")
ce4f2 <- conditional_effects(f4, "z:f")
say("brms estimate__ vs frmtmb estimate__, z:f: max diff %.3g",
    maxdiff(ce4b2[[1]]$estimate__, ce4f2[[1]]$estimate__))
# Monte Carlo over the fit's covariance: coefficient vectors from
# N(chat, V), the curve recomputed at each, its SD per grid point
nd <- ce4f[[1]][, c("z", "f")]
lb <- frm_lp_basis(f4, newdata = nd, re_formula = NA)
pm <- frmtmb:::joint_pos_map(frmtmb:::get_joint_cov(f4))
chat <- vapply(lb$coef_pos, function(k) {
  f4$estimates[[pm$comp[k]]][pm$idx[k]]
}, 0)
set.seed(20260929)
R <- 4000
L <- t(chol(lb$V))
sims <- matrix(NA_real_, R, nrow(nd))
fz <- f4
for (r in seq_len(R)) {
  cr <- chat + as.vector(L %*% rnorm(length(chat)))
  est <- f4$estimates
  for (j in seq_along(lb$coef_pos)) {
    k <- lb$coef_pos[j]
    est[[pm$comp[k]]][pm$idx[k]] <- cr[j]
  }
  fz$estimates <- est
  sims[r, ] <- frm_linpred(fz, newdata = nd, type = "link", dpar = "mu",
                           re_formula = NA)
}
mc_sd <- apply(sims, 2, sd)
ratio <- ce4f[[1]]$se__ / mc_sd
# the standard error of a sample SD is about sd / sqrt(2 (R - 1))
zsc <- (ce4f[[1]]$se__ - mc_sd) / (mc_sd / sqrt(2 * (R - 1)))
say("R = %d, seed 20260929: se__/MC sd over %d grid points: min %.4f, median %.4f, max %.4f",
    R, length(ratio), min(ratio), median(ratio), max(ratio))
say("  |z| of (se__ - MC sd) against the MC error of an SD: max %.2f, points beyond 2: %d, beyond 3: %d",
    max(abs(zsc)), sum(abs(zsc) > 2), sum(abs(zsc) > 3))
mc_lo <- apply(sims, 2, quantile, 0.025)
mc_hi <- apply(sims, 2, quantile, 0.975)
say("  band ends vs MC 2.5%%/97.5%% quantiles: max |diff| / MC width %.4f",
    max(abs(c(ce4f[[1]]$lower__ - mc_lo, ce4f[[1]]$upper__ - mc_hi)) /
          (mc_hi - mc_lo)))

## ---- C. draws: spaghetti, conditional_smooths, posterior_average ------
say("\n== C. draws at the same parameter vectors")
ds1 <- suppressWarnings(suppressMessages(
  frm_sample(f1, chains = 1, iter = 300, warmup = 150, seed = 20260929,
             refresh = 0)))
K <- 10
rows <- round(seq(1, ndraws(ds1), length.out = K))
dsK <- ds1
dsK$draws <- ds1$draws[rows, , drop = FALSE]
idx <- frmtmb.sample:::draws_par_index(ds1$fit)
initsK <- lapply(rows, function(i) {
  fi <- frmtmb.sample:::draws_fit_at(ds1, i, idx)
  frmtmb.sample:::draws_fit_at(ds1, i, idx)
  sg <- frm_linpred(fi, dpar = "sigma", type = "response")[1]
  brms_init(f1, sd1, fi$estimates, sg)
})
say("projection residual over the %d draws (relative): max %.3g", K,
    max(vapply(initsK, attr, 1, "proj_resid")))
bK <- brm_fixed(form1, d, initsK)
say("brms draws: %d", ndraws(bK))
cmp_frames("conditional_smooths on draws",
           conditional_smooths(bK), conditional_smooths(dsK))
sa <- conditional_smooths(bK, spaghetti = TRUE)
sb <- conditional_smooths(dsK, spaghetti = TRUE)
for (k in names(sa)) {
  pa <- attr(sa[[k]], "spaghetti")
  pb <- attr(sb[[k]], "spaghetti")
  if (is.null(pa) && is.null(pb)) {
    say("  %s: no spaghetti in either", k)
    next
  }
  say("  %s spaghetti: rows %d / %d, names identical %s, estimate__ max diff %.3g, sample__ identical %s",
      k, nrow(pa), nrow(pb), identical(names(pa), names(pb)),
      maxdiff(pa$estimate__, pb$estimate__),
      identical(as.character(pa$sample__), as.character(pb$sample__)))
}
for (eff in c("x", "x:f")) {
  ca <- conditional_effects(bK, eff, spaghetti = TRUE)
  cb <- conditional_effects(dsK, eff, spaghetti = TRUE)
  pa <- strip(attr(ca[[1]], "spaghetti"))
  pb <- strip(attr(cb[[1]], "spaghetti"))
  say("conditional_effects(%s, spaghetti): rows %d / %d, names brms [%s] frmtmb [%s]",
      eff, nrow(pa), nrow(pb), paste(names(pa), collapse = ","),
      paste(names(pb), collapse = ","))
  say("  estimate__ max diff %.3g; sample__ identical %s; band estimate__ max diff %.3g",
      maxdiff(pa$estimate__, pb$estimate__),
      identical(as.character(pa$sample__), as.character(pb$sample__)),
      maxdiff(ca[[1]]$estimate__, cb[[1]]$estimate__))
}
say("brms spaghetti + surface: %s",
    tryCatch(conditional_effects(bK, "x", spaghetti = TRUE, surface = TRUE),
             error = conditionMessage))
say("frmtmb spaghetti + surface: %s",
    tryCatch(conditional_effects(dsK, "x", spaghetti = TRUE,
                                 surface = TRUE),
             error = conditionMessage))

say("\n== C2. posterior_average(): same draws, weights and seed")
vars <- c("b_fb", "sigma", "b_Intercept")
say("variables present: brms %s, frmtmb %s",
    all(vars %in% variables(bK)), all(vars %in% variables(dsK)))
for (w in list(c(0.3, 0.7), c(1, 0), c(0.5, 0.5))) {
  pa <- posterior_average(bK, bK, variable = vars, weights = w, seed = 5)
  pb <- posterior_average(dsK, dsK, variable = vars, weights = w, seed = 5)
  say("weights %s: dim brms %s frmtmb %s; names identical %s; max diff %.3g; attr ndraws brms %s frmtmb %s; weights identical %s",
      paste(w, collapse = "/"), paste(dim(pa), collapse = "x"),
      paste(dim(pb), collapse = "x"), identical(names(pa), names(pb)),
      maxdiff(as.matrix(pa), as.matrix(pb)),
      paste(attr(pa, "ndraws"), collapse = "/"),
      paste(attr(pb, "ndraws"), collapse = "/"),
      identical(unname(attr(pa, "weights")), unname(attr(pb, "weights"))))
}
pa <- posterior_average(bK, bK, variable = vars, weights = c(1, 3),
                        ndraws = 7, seed = 5)
pb <- posterior_average(dsK, dsK, variable = vars, weights = c(1, 3),
                        ndraws = 7, seed = 5)
say("ndraws = 7: brms %s, frmtmb %s, max diff %.3g",
    paste(dim(pa), collapse = "x"), paste(dim(pb), collapse = "x"),
    maxdiff(as.matrix(pa), as.matrix(pb)))
say("names(attr(, 'weights')): brms %s | frmtmb %s",
    paste(names(attr(pa, "weights")), collapse = ","),
    paste(names(attr(pb, "weights")), collapse = ","))
pa <- tryCatch(posterior_average(bK, bK, variable = "nope", weights = c(1, 1)),
               error = conditionMessage)
pb <- tryCatch(posterior_average(dsK, dsK, variable = "nope",
                                 weights = c(1, 1)),
               error = conditionMessage)
say("unknown variable: brms '%s' | frmtmb '%s'", pa, pb)
pa <- tryCatch(posterior_average(bK, bK, variable = c("sigma", "nope"),
                                 weights = c(1, 1), missing = 0),
               error = conditionMessage)
pb <- tryCatch(posterior_average(dsK, dsK, variable = c("sigma", "nope"),
                                 weights = c(1, 1), missing = 0),
               error = conditionMessage)
say("missing = 0, a variable no model has: brms '%s' | frmtmb '%s'", pa, pb)
# a variable one model lacks: brms's example fit b3 has b_x, which the
# smooth model does not; the frmtmb side takes the same draws of f3
d3 <- suppressWarnings(suppressMessages(
  frm_sample(f3, chains = 1, iter = 60, warmup = 50, seed = 7,
             refresh = 0)))
i3K <- lapply(seq_len(ndraws(d3)), function(i) {
  fi <- frmtmb.sample:::draws_fit_at(d3, i,
                                     frmtmb.sample:::draws_par_index(d3$fit))
  be <- fi$estimates$beta
  list(b = arr(be[-1]),
       Intercept = be[1] + sum(colMeans(X3[, -1]) * be[-1]),
       sigma = frm_linpred(fi, dpar = "sigma", type = "response")[1])
})
b3K <- brm_fixed(form3, d, i3K)
pa <- posterior_average(bK, b3K, variable = c("sigma", "b_x"),
                        weights = c(2, 3), missing = 0, seed = 5)
pb <- posterior_average(dsK, d3, variable = c("sigma", "b_x"),
                        weights = c(2, 3), missing = 0, seed = 5)
say("missing = 0 across two models: dims %s / %s, max diff %.3g, ndraws %s / %s",
    paste(dim(pa), collapse = "x"), paste(dim(pb), collapse = "x"),
    maxdiff(as.matrix(pa), as.matrix(pb)),
    paste(attr(pa, "ndraws"), collapse = "/"),
    paste(attr(pb, "ndraws"), collapse = "/"))
say("done")
