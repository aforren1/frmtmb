## Lane rescorpp: frmtmb.sample's predictive draws and log_lik() of a
## set_rescor(TRUE) fit against brms 2.23.0 on the same data.
##   Rscript dev/rescorpp-compare.R <lane|base> <g3|gs|st|g2f>
## Draws: frm_sample(chains = 2, iter = 2000, seed = 21); brms's own fit
## is dev/rescorpp-brmsfit.R's, seed 21. Predictive seeds are set below.
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-data.R")
w <- commandArgs(TRUE)[2]
suppressMessages(requireNamespace("brms"))
fd <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-brmsfits"
b <- readRDS(file.path(fd, paste0(w, ".rds")))
d <- rescorpp_data(rescorpp_dname(w))
fam <- if (rescorpp_family(w) == "student") student() else gaussian()
f <- frm(eval(rescorpp_frm_forms[[w]]), family = fam, data = d)
ds <- suppressWarnings(suppressMessages(frm_sample(
  f, chains = 2, iter = 2000, seed = 21, refresh = 0)))
S <- ndraws(ds)
cat("model", w, " draws", S, " rows", nrow(d), "\n")
resps <- names(f$spec$responses)
K <- length(resps)
rescor <- isTRUE(f$spec$rescor)

## ---- brms's own code on frmtmb's draws (the transplant) -------------
fr <- as.matrix(posterior::as_draws_matrix(ds))
bm <- as.matrix(posterior::as_draws_matrix(b))[seq_len(S), , drop = FALSE]
read <- grep("^(b_|sigma|rescor__|nu$)", colnames(bm), value = TRUE)
show("brms variables read, absent from frmtmb's draws",
     setdiff(read, colnames(fr)))
for (nm in intersect(read, colnames(fr))) bm[, nm] <- fr[, nm]
dm_frm <- posterior::as_draws_matrix(bm)
dm_brm <- posterior::subset_draws(posterior::as_draws_matrix(b),
                                  draw = seq_len(S))
prep_from <- function(x, draws, resp = NULL) {
  # brms:::prepare_predictions.brmsfit() with the draws handed in
  x <- brms:::restructure(x)
  x <- brms:::exclude_terms(x, incl_autocor = TRUE, offset = TRUE,
                            smooths_only = FALSE)
  resp <- brms:::validate_resp(resp, x)
  sdata <- brms::standata(x, resp = resp, internal = TRUE)
  bframe <- brms:::brmsframe(brms:::update_re_terms(x$formula, NULL),
                             data = x$data)
  prep_re <- brms:::prepare_predictions_re_global(
    bframe = bframe, draws = draws, sdata = sdata, resp = resp,
    old_reframe = x$ranef, sample_new_levels = "uncertainty")
  brms:::prepare_predictions(bframe, draws = draws, sdata = sdata,
                             prep_re = prep_re, resp = resp,
                             sample_new_levels = "uncertainty", nug = NULL,
                             new = FALSE, oos = NULL, stanvars = x$stanvars)
}
ll_brms <- function(draws, resp = NULL) {
  out <- rstantools::log_lik(prep_from(b, draws, resp))
  unname(out)
}
show("control: transplant route with brms's own draws vs log_lik(b)",
     max(abs(ll_brms(dm_brm) - brms::log_lik(b)[seq_len(S), ])))
rel <- function(a, r) max(abs(a - r)) / max(abs(r))

## ---- log_lik ----------------------------------------------------------
L <- show("log_lik(ds) vs brms at frmtmb's draws: max abs / max |ref|",
          rel(log_lik(ds), ll_brms(dm_frm)))
show("log_lik(ds, resp = 'y1') vs brms at frmtmb's draws",
     rel(log_lik(ds, resp = "y1"), ll_brms(dm_frm, "y1")))
if (w == "g3") {
  show("log_lik(ds, resp = c('y1','y3')) vs brms at frmtmb's draws",
       rel(log_lik(ds, resp = c("y1", "y3")),
           ll_brms(dm_frm, c("y1", "y3"))))
  show("log_lik(ds, resp = c('y3','y1')) vs brms, order reversed",
       rel(log_lik(ds, resp = c("y3", "y1")),
           ll_brms(dm_frm, c("y3", "y1"))))
}
if (w == "st") {
  # brms's one-response margin of a Student-t rescor model reads nu = 1;
  # the margin of the multivariate t is a t with the shared nu
  dr <- as.data.frame(posterior::as_draws_df(ds))
  m1 <- outer(dr$b_y1_Intercept, rep(1, nrow(d))) + outer(dr$b_y1_x, d$x)
  z <- (matrix(d$y1, S, nrow(d), byrow = TRUE) - m1) / dr$sigma_y1
  show("st: log_lik(ds, resp = 'y1') vs t(shared nu) margin",
       rel(log_lik(ds, resp = "y1"),
           dt(z, dr$nu, log = TRUE) - log(dr$sigma_y1)))
}
show("log_lik(ds, combine = FALSE) dim", dim(log_lik(ds, combine = FALSE)))

## ---- shapes -------------------------------------------------------------
shape <- function(a) list(dim = dim(a), names3 = dimnames(a)[3][[1]])
show("posterior_predict(ds) shape", shape(posterior_predict(ds, ndraws = 5)))
show("posterior_epred(ds) shape", shape(posterior_epred(ds, ndraws = 5)))
show("posterior_linpred(ds) shape", shape(posterior_linpred(ds, ndraws = 5)))
show("predictive_error(ds) shape", shape(predictive_error(ds, ndraws = 5)))
show("fitted(ds) shape", shape(fitted(ds, ndraws = 5)))
show("predict(ds) shape", shape(predict(ds, ndraws = 5)))
show("residuals(ds) shape", shape(residuals(ds, ndraws = 5)))
show("residuals(ds, type = 'pearson') shape",
     shape(residuals(ds, type = "pearson", ndraws = 5)))
show("predictive_interval(ds)", predictive_interval(ds, ndraws = 5))
show("bayes_R2(ds)", bayes_R2(ds))

## ---- the predictive residual correlation --------------------------------
## Per row i, the correlation over draws of the predictive residuals
## e = prediction - epred of responses j and k. Given the parameter
## draws, its expectation is target_i = mean(w r s_j s_k) /
## sqrt(mean(w s_j^2) mean(w s_k^2)), w = nu / (nu - 2) for student()
## and 1 otherwise. z uses the spread over rows as the MCSE.
dr <- as.data.frame(posterior::as_draws_df(ds))
sig <- function(r) {
  s <- posterior_epred(ds, dpar = "sigma", resp = r)
  if (is.null(dim(s))) matrix(s, S, nrow(d)) else s
}
wgt <- if (w == "st") dr$nu / (dr$nu - 2) else rep(1, S)
target <- function(j, k) {
  rc <- if (rescor) dr[[paste0("rescor__", resps[j], "__", resps[k])]] else 0
  sj <- sig(resps[j]); sk <- sig(resps[k])
  colMeans(wgt * rc * sj * sk) /
    sqrt(colMeans(wgt * sj^2) * colMeans(wgt * sk^2))
}
rowcor <- function(e, j, k) {
  vapply(seq_len(dim(e)[2]), function(i) cor(e[, i, j], e[, i, k]), 1)
}
zs <- function(a, b) {
  dlt <- a - b
  c(mean_a = mean(a), mean_b = mean(b), ratio = mean(a) / mean(b),
    z = mean(dlt) / (sd(dlt) / sqrt(length(dlt))))
}
set.seed(101)
t0 <- proc.time()[["elapsed"]]
P <- posterior_predict(ds)
cat("posterior_predict(ds) seconds", proc.time()[["elapsed"]] - t0, "\n")
E <- posterior_epred(ds)
show("posterior_predict(ds) dim", dim(P))
if (length(dim(P)) == 3L) {
  e <- P - E
  set.seed(102)
  pb <- rstantools::posterior_predict(prep_from(b, dm_frm))
  eb <- pb - E
  for (j in 1:(K - 1)) for (k in (j + 1):K) {
    tg <- target(j, k)
    cf <- rowcor(e, j, k)
    cb <- rowcor(eb, j, k)
    show(paste0("pair ", j, k, ": frmtmb residual cor vs its target"),
         zs(cf, tg))
    show(paste0("pair ", j, k, ": brms code on frmtmb's draws vs target"),
         zs(cb, tg))
    show(paste0("pair ", j, k, ": frmtmb vs brms code, same draws"),
         zs(cf, cb))
    if (w == "st") {
      # one chi-square per row shared by the responses makes |e_j| and
      # |e_k| depend even through the scale; brms's draw has the same
      show(paste0("pair ", j, k, ": cor of squared residuals, frmtmb vs ",
                  "brms code"),
           zs(vapply(seq_len(dim(e)[2]), function(i) {
             cor(e[, i, j]^2, e[, i, k]^2, method = "spearman")
           }, 1), vapply(seq_len(dim(e)[2]), function(i) {
             cor(eb[, i, j]^2, eb[, i, k]^2, method = "spearman")
           }, 1)))
    }
  }
  # each response's spread against its expected predictive sd
  for (j in seq_len(K)) {
    sj <- sig(resps[j])
    show(paste0("resp ", j, ": sd(residual) / sqrt(mean(w sigma^2))"),
         zs(apply(e[, , j], 2, sd), sqrt(colMeans(wgt * sj^2))))
  }
  # a linear combination: Var(e_1 + e_2) / its expectation
  vs <- apply(e[, , 1] + e[, , 2], 2, var)
  s1 <- sig(resps[1]); s2 <- sig(resps[2])
  rc <- if (rescor) dr[[paste0("rescor__", resps[1], "__", resps[2])]] else 0
  ev <- colMeans(wgt * (s1^2 + s2^2 + 2 * rc * s1 * s2))
  show("Var(e1 + e2) / its expectation", zs(vs, ev))
}
## core predict() on the ML fit, for the record
if (rescor) {
  set.seed(103)
  pc <- predict(f, summary = FALSE, ndraws = 2000)
  ec <- fitted(f)[, "Estimate", ]
  rc <- vapply(seq_len(nrow(d)), function(i) {
    cor(pc[, i, 1] - ec[i, 1], pc[, i, 2] - ec[i, 2])
  }, 1)
  show("core predict(): mean residual cor rows | rescor_matrix",
       c(mean(rc), rescor_matrix(f)[1, 2]))
}
