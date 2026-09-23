# Second batch of shared models. Run:
# Rscript dev/drmtmb-agree2.R > dev/drmtmb-log/agree2.log
# Seed: sim_batch2(505); mi() drops 60 xm values after that draw.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree2-data.R")
d <- sim_batch2(505)
res <- list()
add <- function(r) {
  r <- tryCatch(r, error = function(e) {
    cat("\n########## COMPARISON FAILED:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(r)) return(invisible())
  fmt_block(r); res[[r$summary$model]] <<- r
}
ids <- function(k) lapply(seq_len(k), function(i) lin(i, i))
dfit <- function(...) drmTMB::drmTMB(...)

# Beta-binomial: drmTMB phi = 1 / sigma^2.
add(compare_fits("k beta_binomial mu RE, ML",
  dfit(dbf(cbind(ys, nt - ys) ~ x + (1 | g), sigma ~ 1),
       family = drmTMB::beta_binomial(), data = d),
  frm(bf(ys | trials(nt) ~ x + (1 | g), phi ~ 1), family = beta_binomial(),
      data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), lin(4, 4))))
# Tweedie: drmTMB phi = sigma^2; power 1 + plogis(eta) in both.
add(compare_fits("l tweedie mu RE, ML",
  dfit(dbf(ytw ~ x + (1 | g), sigma ~ 1, nu ~ 1), family = drmTMB::tweedie(),
       data = d),
  frm(bf(ytw ~ x + (1 | g), phi ~ 1, power ~ 1), family = tweedie(), data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, 2), lin(4, 4), lin(5, 5))))
add(compare_fits("m lognormal mu RE, ML",
  dfit(dbf(yln ~ x + (1 | g), sigma ~ 1), family = drmTMB::lognormal(),
       data = d),
  frm(bf(yln ~ x + (1 | g), sigma ~ 1), family = lognormal(), data = d),
  ids(4)))
add(compare_fits("m lognormal sigma RE, ML",
  dfit(dbf(yln ~ x, sigma ~ 1 + (1 | g)), family = drmTMB::lognormal(),
       data = d),
  frm(bf(yln ~ x, sigma ~ 1 + (1 | g)), family = lognormal(), data = d),
  ids(4)))
# Gamma: drmTMB sigma is the CV, shape = 1 / sigma^2.
add(compare_fits("n Gamma(log) mu RE, ML",
  dfit(dbf(ygam ~ x + (1 | g), sigma ~ 1), family = Gamma(link = "log"),
       data = d),
  frm(bf(ygam ~ x + (1 | g), shape ~ 1), family = Gamma(link = "log"),
      data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), lin(4, 4))))
add(compare_fits("n Gamma(log) sigma RE, ML",
  dfit(dbf(ygam ~ x, sigma ~ 1 + (1 | g)), family = Gamma(link = "log"),
       data = d),
  frm(bf(ygam ~ x, shape ~ 1 + (1 | g)), family = Gamma(link = "log"),
      data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), log2x(4, 4))))
add(compare_fits("o zero-inflated nbinom2, ML",
  dfit(dbf(yzi ~ x, sigma ~ 1, zi ~ 1), family = drmTMB::nbinom2(), data = d),
  frm(bf(yzi ~ x, shape ~ 1, zi ~ 1), family = zero_inflated_negbinomial(),
      data = d),
  list(lin(1, 1), lin(2, 2), lin(3, 3, -2), lin(4, 4))))
# drmTMB's sd(g) ~ w is log sd_g = c0 + c1 w_g. frmtmb has no such
# grammar; the nl formula y = b0 + exp(c1 w) * zz, zz ~ N(0, s^2) is the
# same model with c0 = log s.
add(compare_fits("p sd(g) ~ w vs nl exp(lsd) * zz, ML",
  dfit(dbf(y ~ x + (1 | g), sigma ~ 1, sd(g) ~ w), family = gaussian(),
       data = d),
  frm(bf(y ~ b0 + exp(lsd) * zz, b0 ~ x, lsd ~ 0 + w, zz ~ 0 + (1 | g),
         nl = TRUE), family = gaussian(), data = d),
  list(lin(1, 1), lin(2, 2), lin(4, 3), lin(5, 4), lin(3, 5))))
dm <- d
set.seed(506)
dm$xm[sample(nrow(dm), 60)] <- NA
add(compare_fits("q mi() gaussian missing predictor, ML",
  dfit(dbf(y ~ z + mi(xm), sigma ~ 1), family = gaussian(), data = dm,
       impute = list(xm = xm ~ w2),
       missing = drmTMB::miss_control(predictor = "model")),
  frm(bf(y ~ z + mi(xm)) + bf(xm | mi() ~ w2) + set_rescor(FALSE),
      family = gaussian(), data = dm),
  list(lin(1, 1), lin(2, 2), lin(3, 3), lin(6, 4), lin(4, 5), lin(5, 6),
       lin(7, 7))))

s <- do.call(rbind, lapply(res, `[[`, "summary"))
cat("\n\n########## SUMMARY\n")
op <- options(width = 250)
print(s[, c("model", "ll_frm", "ll_drm", "ll_diff", "obj_gap_at_drm_opt",
            "obj_gap_at_frm_opt", "obj_gap_off_opt", "max_abs_diff_over_se",
            "max_abs_log_se_ratio")], digits = 4, row.names = FALSE)
options(op)
saveRDS(res, file.path("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev",
                      "drmtmb-log/agree2.rds"))
