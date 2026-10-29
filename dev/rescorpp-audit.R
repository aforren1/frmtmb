## Lane rescorpp: every other reader of set_rescor(TRUE) draws and ML
## fits, measured on one arm.
##   Rscript dev/rescorpp-audit.R <lane|base>
## Data g3 and st of dev/rescorpp-data.R; frm_sample(chains = 2,
## iter = 1000, seed = 21).
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-data.R")
d3 <- rescorpp_data("g3")
dst <- rescorpp_data("st")
f3 <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE),
          family = gaussian(), data = d3)
fst <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE),
           family = student(), data = dst)
smp <- function(f) suppressWarnings(suppressMessages(frm_sample(
  f, chains = 2, iter = 1000, seed = 21, refresh = 0)))
ds3 <- smp(f3)
dsst <- smp(fst)

## ---- the loo family: weights from log_lik(resp =) ----------------------
show("loo(ds3) elpd_loo", loo(ds3)$estimates["elpd_loo", ])
show("loo(ds3, resp = 'y1') elpd_loo", loo(ds3, resp = "y1")$estimates[
  "elpd_loo", ])
show("sum of colMeans(log_lik(resp = 'y1')) | of log_lik()",
     c(sum(colMeans(log_lik(ds3, resp = "y1"))),
       sum(colMeans(log_lik(ds3)))))
show("loo_predict(ds3) dim", dim(loo_predict(ds3)))
show("loo_predict(ds3, resp = 'y1') length", length(loo_predict(ds3,
                                                                resp = "y1")))
show("loo_epred(ds3) dimnames3", dimnames(loo_epred(ds3))[[3]])
show("loo_R2(ds3)", loo_R2(ds3))
show("waic(ds3) elpd", waic(ds3)$estimates["elpd_waic", ])

## ---- pp_check: one response, its margin ---------------------------------
show("pp_check(ds3) without resp", class(pp_check(ds3, ndraws = 5)))
show("pp_check(ds3, resp = 'y2')", class(pp_check(ds3, resp = "y2",
                                                  ndraws = 5)))

## ---- conditional_effects(method = "posterior_predict") -----------------
## the band at x = 0 against posterior_predict()'s own quantiles there;
## for student() the margin must carry the shared nu
for (nm in c("g3", "st")) {
  ds <- if (nm == "g3") ds3 else dsst
  set.seed(5)
  ce <- conditional_effects(ds, effects = "x", method = "posterior_predict",
                            resp = "y1", int_conditions = list(x = 0))
  ce1 <- ce[[1]]
  set.seed(6)
  pp0 <- posterior_predict(ds, newdata = data.frame(x = 0), resp = "y1")
  q <- quantile(pp0[, 1], c(0.025, 0.975))
  show(paste0(nm, ": ce band at x = 0 | posterior_predict quantiles"),
       rbind(ce = c(ce1$lower__[1], ce1$upper__[1]), pp = q))
  show(paste0(nm, ": conditional_effects() default names"),
       names(conditional_effects(ds, method = "posterior_predict")))
}

## ---- core on the ML fit ---------------------------------------------------
set.seed(7)
pc <- predict(f3, summary = FALSE, ndraws = 3000, resp = c("y1", "y3"))
show("core predict(resp = c('y1','y3')) dim", dim(pc))
ec <- fitted(f3, resp = c("y1", "y3"))[, "Estimate", ]
rc <- vapply(seq_len(nrow(d3)), function(i) {
  cor(pc[, i, 1] - ec[i, 1], pc[, i, 2] - ec[i, 2])
}, 1)
show("core predict(resp = c('y1','y3')): mean residual cor | rescor",
     c(mean(rc), rescor_matrix(f3)["y1", "y3"]))
set.seed(8)
pst <- predict(fst, summary = FALSE, ndraws = 3000)
est <- fitted(fst)[, "Estimate", ]
show("core predict() student: mean spearman residual cor | rescor",
     c(mean(vapply(seq_len(nrow(dst)), function(i) {
       cor(pst[, i, 1] - est[i, 1], pst[, i, 2] - est[i, 2],
           method = "spearman")
     }, 1)), rescor_matrix(fst)[1, 2]))
show("core simulate(f3)", dim(simulate(f3)))
show("core frm_simulate(mvbf rescor)", {
  s <- frm_simulate(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE),
                    data = d3, family = gaussian(), seed = 1)
  str(s)
  s
})

## ---- missing responses: mi() and subset() with rescor -------------------
dmi <- d3
dmi$y1[1:5] <- NA
show("mi() with rescor", frm(mvbf(bf(y1 | mi() ~ x), bf(y2 ~ x)) +
                               set_rescor(TRUE), data = dmi))
dsub <- d3
dsub$s <- rep(c(TRUE, FALSE), length.out = nrow(dsub))
show("subset() with rescor", frm(mvbf(bf(y1 | subset(s) ~ x), bf(y2 ~ x)) +
                                   set_rescor(TRUE), data = dsub))
