# Agreement between frmtmb and drmTMB on every model both can fit.
# Run: Rscript dev/drmtmb-agree.R > dev/drmtmb-log/agree.log
# Seeds: sim_grouped(101), sim_pedigree(202), sim_meta(303), tree 404.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
cat("frmtmb", format(packageVersion("frmtmb")), " drmTMB",
    format(packageVersion("drmTMB")), " RTMB", format(packageVersion("RTMB")),
    " TMB", format(packageVersion("TMB")), "\n")

d <- sim_grouped(101)
ids <- function(k) lapply(seq_len(k), function(i) lin(i, i))
res <- list()
add <- function(r) {
  r <- tryCatch(r, error = function(e) {
    cat("\n########## COMPARISON FAILED:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(r)) return(invisible())
  fmt_block(r); res[[r$summary$model]] <<- r
}
dfit <- function(...) drmTMB::drmTMB(...)

# (a) Gaussian location-scale, random intercept in mu.
add(compare_fits("a gaussian mu RE, ML",
  dfit(dbf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d),
  frm(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d),
  ids(5)))
add(compare_fits("a gaussian mu RE, REML",
  dfit(dbf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d,
       REML = TRUE),
  frm(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d,
      REML = TRUE),
  ids(3)))

# (b) Random intercepts in mu and sigma, independent and correlated.
add(compare_fits("b gaussian mu+sigma RE, ML",
  dfit(dbf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
       data = d),
  frm(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
      data = d),
  ids(6)))
# REML with a sigma random effect is NOT compared here: drmTMB then also
# integrates out beta_sigma, a different criterion from frmtmb's. See
# dev/drmtmb-reml-nesting.R.
add(compare_fits("b gaussian mu+sigma RE correlated |p|, ML",
  dfit(dbf(y ~ x + (1 | p | g), sigma ~ z + (1 | p | g)),
       family = gaussian(), data = d),
  frm(bf(y ~ x + (1 | p | g), sigma ~ z + (1 | p | g)), family = gaussian(),
      data = d),
  c(ids(4), list(lin(5, 5), lin(6, 7), rho_re_map(7, 6)))))

# (c) Beta: drmTMB sigma = 1/sqrt(phi), so log phi = -2 log sigma.
add(compare_fits("c beta mu RE, ML",
  dfit(dbf(yb ~ x + (1 | g), sigma ~ 1), family = drmTMB::beta(), data = d),
  frm(bf(yb ~ x + (1 | g), phi ~ 1), family = Beta(), data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), lin(4, 4))))

# (d) NB2: drmTMB sigma = 1/sqrt(shape), so log shape = -2 log sigma.
add(compare_fits("d nbinom2 mu RE, ML",
  dfit(dbf(yc ~ x + (1 | g), sigma ~ 1), family = drmTMB::nbinom2(),
       data = d),
  frm(bf(yc ~ x + (1 | g), shape ~ 1), family = negbinomial(), data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), lin(4, 4))))
add(compare_fits("d nbinom2 sigma RE, ML",
  dfit(dbf(yc ~ x, sigma ~ 1 + (1 | g)), family = drmTMB::nbinom2(),
       data = d),
  frm(bf(yc ~ x, shape ~ 1 + (1 | g)), family = negbinomial(), data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), log2x(4, 4))))

# (e) Bivariate Gaussian with residual correlation.
add(compare_fits("e bivariate gaussian rho12, ML",
  dfit(dbf(mu1 = y ~ x, mu2 = y2 ~ x, rho12 = ~ 1),
       family = drmTMB::biv_gaussian(), data = d),
  frm(mvbf(bf(y ~ x), bf(y2 ~ x)) + set_rescor(TRUE), family = gaussian(),
      data = d),
  c(ids(6), list(rho_re_map(7, 7)))))
add(compare_fits("e bivariate gaussian rho12 + correlated RE |p|, ML",
  dfit(dbf(mu1 = y ~ x + (1 | p | g), mu2 = y2 ~ x + (1 | p | g),
           rho12 = ~ 1), family = drmTMB::biv_gaussian(), data = d),
  frm(mvbf(bf(y ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
        set_rescor(TRUE), family = gaussian(), data = d),
  c(ids(6), list(rho_re_map(10, 7), lin(7, 8), lin(8, 9), rho_re_map(9, 10)))))
add(compare_fits("e bivariate gaussian rho12 + correlated RE |p|, REML",
  dfit(dbf(mu1 = y ~ x + (1 | p | g), mu2 = y2 ~ x + (1 | p | g),
           rho12 = ~ 1), family = drmTMB::biv_gaussian(), data = d,
       REML = TRUE),
  frm(mvbf(bf(y ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
        set_rescor(TRUE), family = gaussian(), data = d, REML = TRUE),
  list(lin(1, 1), lin(2, 2), rho_re_map(6, 3), lin(3, 4), lin(4, 5),
       rho_re_map(5, 6))))
# Boundary case: y2b shares y's group effect, so the group correlation
# is 1 in the generating model.
add(compare_fits("e boundary: bivariate RE correlation at 1, ML",
  dfit(dbf(mu1 = y ~ x + (1 | p | g), mu2 = y2b ~ x + (1 | p | g),
           rho12 = ~ 1), family = drmTMB::biv_gaussian(), data = d),
  frm(mvbf(bf(y ~ x + (1 | p | g)), bf(y2b ~ x + (1 | p | g))) +
        set_rescor(TRUE), family = gaussian(), data = d),
  c(ids(6), list(rho_re_map(10, 7), lin(7, 8), lin(8, 9), rho_re_map(9, 10)))))

# (f) Animal model: drmTMB animal(A = A) against frmtmb gr(cov = A).
pd <- sim_pedigree(202)
dd <- pd$data; A <- pd$A
add(compare_fits("f animal model A, ML",
  dfit(dbf(y ~ x + animal(1 | id, A = A), sigma ~ 1), family = gaussian(),
       data = dd),
  frm(bf(y ~ x + (1 | gr(id, cov = A))), family = gaussian(), data = dd,
      data2 = list(A = A)),
  ids(4)))
add(compare_fits("f animal model A, REML",
  dfit(dbf(y ~ x + animal(1 | id, A = A), sigma ~ 1), family = gaussian(),
       data = dd, REML = TRUE),
  frm(bf(y ~ x + (1 | gr(id, cov = A))), family = gaussian(), data = dd,
      data2 = list(A = A), REML = TRUE),
  ids(2)))
# drmTMB's own pedigree builder against the A built here.
ped <- pd$ped
f_ped <- dfit(dbf(y ~ x + animal(1 | id, pedigree = ped), sigma ~ 1),
              family = gaussian(), data = dd)
f_A <- dfit(dbf(y ~ x + animal(1 | id, A = A), sigma ~ 1),
            family = gaussian(), data = dd)
cat("\nanimal: drmTMB pedigree route minus A route, logLik:",
    format(as.numeric(logLik(f_ped)) - as.numeric(logLik(f_A)),
           digits = 6), "\n")

# (f2) Phylogenetic intercept: drmTMB phylo(tree) against frmtmb
# gr(cov = ape::vcv(tree)), with both the covariance and the
# correlation form, because the tree is not unit height.
set.seed(404)
tree <- ape::rcoal(50)
tree$tip.label <- paste0("sp", seq_len(50))
Vt <- ape::vcv(tree)
cat("tree height (diag of vcv):", format(range(diag(Vt)), digits = 8), "\n")
dp <- data.frame(sp = factor(rep(tree$tip.label, each = 4),
                             levels = tree$tip.label), x = rnorm(200))
a <- as.vector(t(chol(Vt)) %*% rnorm(50)) * 0.8
dp$y <- 0.5 + 0.3 * dp$x + a[as.integer(dp$sp)] + rnorm(200, 0, 0.6)
fd_ph <- dfit(dbf(y ~ x + phylo(1 | sp, tree = tree), sigma ~ 1),
              family = gaussian(), data = dp)
Ccor <- ape::vcv(tree, corr = TRUE)
# drmTMB scales the tree to unit height, so its SD is sqrt(h) times
# the SD that multiplies the raw vcv(tree): theta_f = d - log(h) / 2.
h <- diag(Vt)[[1]]
add(compare_fits("f2 phylo vs gr(cov = vcv(tree)), ML", fd_ph,
  frm(bf(y ~ x + (1 | gr(sp, cov = Vt))), family = gaussian(), data = dp,
      data2 = list(Vt = Vt)),
  c(ids(3), list(lin(4, 4, 1, -log(h) / 2)))))
add(compare_fits("f2 phylo vs gr(cov = vcv(tree, corr = TRUE)), ML", fd_ph,
  frm(bf(y ~ x + (1 | gr(sp, cov = Ccor))), family = gaussian(), data = dp,
      data2 = list(Ccor = Ccor)), ids(4)))

# (g) Student-t. nu links differ; see nu_map().
add(compare_fits("g student mu RE, ML",
  dfit(dbf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1), family = drmTMB::student(),
       data = d),
  frm(bf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1), family = student(), data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3), nu_map(4, 4), lin(5, 5))))

# (h) Cumulative logit. Both: first cutpoint raw, then log increments.
add(compare_fits("h cumulative logit mu RE, ML",
  dfit(dbf(yo ~ x + (1 | g)), family = drmTMB::cumulative_logit(), data = d),
  frm(bf(yo ~ x + (1 | g)), family = cumulative("logit"), data = d),
  list(lin(1, 1), lin(3, 2), lin(4, 3), lin(5, 4), lin(2, 5))))

# (i) Skew-normal, fixed effects: drmTMB nu against brms alpha.
add(compare_fits("i skew_normal, ML",
  dfit(dbf(ysn ~ x, sigma ~ 1, nu ~ 1), family = drmTMB::skew_normal(),
       data = d),
  frm(bf(ysn ~ x, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = d),
  ids(4)))

# (j) Meta-analysis with known sampling variances. frmtmb's se() with
# sigma = TRUE adds the residual SD that drmTMB's sigma models; metafor
# fits the same model and is the third opinion.
m <- sim_meta(303)
fd_m <- dfit(dbf(yi ~ x + meta_V(V = vi), sigma ~ 1), family = gaussian(),
             data = m)
ff_m <- frm(bf(yi | se(sei, sigma = TRUE) ~ x), family = gaussian(), data = m)
rma <- metafor::rma(yi, vi, mods = ~ x, data = m, method = "ML")
add(compare_fits("j meta-analysis meta_V vs se(sigma = TRUE), ML", fd_m,
  ff_m, ids(3),
  extra = c(metafor_ll = as.numeric(logLik(rma)),
            metafor_tau = sqrt(rma$tau2),
            drm_sigma = exp(fd_m$opt$par[[3]]),
            frm_sigma = exp(ff_m$opt$par[[3]]))))

summ <- do.call(rbind, lapply(res, `[[`, "summary"))
cat("\n\n########## SUMMARY\n")
op <- options(width = 250)
print(summ[, c("model", "ll_frm", "ll_drm", "ll_diff", "obj_gap_at_drm_opt",
               "obj_gap_at_frm_opt", "obj_gap_off_opt",
               "max_abs_diff_over_se",
               "max_abs_log_se_ratio", "grad_frm", "grad_drm")],
      digits = 4, row.names = FALSE)
options(op)
saveRDS(res, file.path("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev",
                      "drmtmb-log/agree.rds"))
