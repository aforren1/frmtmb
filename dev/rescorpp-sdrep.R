## Lane rescorpp: the sd-ratio z of dev/rescorpp-compare.R over 30
## predictive seeds at FIXED parameter draws, so E[z | draws] is read off
## directly. Data g3 (dev/rescorpp-data.R), frm_sample seed 21, 2000
## draws; predictive seeds 1..30.
##   Rscript dev/rescorpp-sdrep.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-data.R")
d <- rescorpp_data("g3")
f <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE),
         family = gaussian(), data = d)
ds <- suppressWarnings(suppressMessages(frm_sample(
  f, chains = 2, iter = 2000, seed = 21, refresh = 0)))
E <- posterior_epred(ds)
tg <- vapply(c("y1", "y2", "y3"), function(r) {
  s <- posterior_epred(ds, dpar = "sigma", resp = r)
  sqrt(mean(s[, 1]^2))
}, 1)
zst <- function(e, t) {
  dd <- apply(e, 2, sd) - t
  mean(dd) / (sd(dd) / sqrt(length(dd)))
}
z <- t(vapply(1:30, function(k) {
  set.seed(k)
  e <- posterior_predict(ds) - E
  c(joint = vapply(1:3, function(j) zst(e[, , j], tg[j]), 1),
    single1 = zst(posterior_predict(ds, resp = "y1") - E[, , 1], tg[1]))
}, numeric(4)))
print(round(colMeans(z), 3))
print(round(apply(z, 2, sd), 3))
# a pure-R draw at the same parameters, as the control
dr <- as.data.frame(posterior::as_draws_df(ds))
zc <- vapply(1:30, function(k) {
  set.seed(100 + k)
  e <- dr$sigma_y1 * matrix(rnorm(nrow(dr) * nrow(d)), nrow(dr))
  zst(e, tg[1])
}, 1)
cat("control (rnorm at the draws' sigma_y1): mean z", round(mean(zc), 3),
    " sd", round(sd(zc), 3), "\n")
