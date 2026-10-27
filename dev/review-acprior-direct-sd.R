# Reviewer: the sd and tail quantiles of the ar posterior, frm_sample
# after against the direct brms fit (dev/review-acprior-direct.R), as
# ratios to posterior's MC standard errors.
# Log: dev/review-acprior-log/direct-sd.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
D <- readRDS(file.path(wt, "dev", "review-acprior-log", "direct-brms.rds"))$par
A <- readRDS(file.path(wt, "dev", "acprior-out", "agree-after.rds"))$ar$par
C <- readRDS(file.path(wt, "dev", "acprior-out", "agree-after.rds"))$cosy$par
Rb <- readRDS(file.path(wt, "dev", "acprior-out", "agree-brms.rds"))$cosy$par
st <- function(x) {
  m <- matrix(x, ncol = 4)
  c(sd = sd(x), mcse_sd = unname(posterior::mcse_sd(m)),
    q05 = unname(quantile(x, 0.05)), mcse_q05 = unname(posterior::mcse_quantile(m, 0.05)),
    q95 = unname(quantile(x, 0.95)), mcse_q95 = unname(posterior::mcse_quantile(m, 0.95)))
}
cmp <- function(lab, a, b) {
  sa <- st(a); sb <- st(b)
  for (s in c("sd", "q05", "q95")) {
    z <- (sa[[s]] - sb[[s]]) / sqrt(sa[[paste0("mcse_", s)]]^2 +
                                      sb[[paste0("mcse_", s)]]^2)
    cat(sprintf("%s %s: frm %.4f brms %.4f z %.2f\n", lab, s, sa[[s]],
                sb[[s]], z))
  }
}
cmp("ar (direct)", A, D)
cmp("cosy", C, Rb)
