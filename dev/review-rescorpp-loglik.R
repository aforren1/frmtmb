## Review rescorpp: log_lik(resp =) margins and subsets on 3-response
## gaussian and student rescor models against mvtnorm at the draws;
## loo / waic / loo_compare / moment_match on resp.
##   Rscript dev/review-rescorpp-loglik.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
source(file.path(wt, "dev/rescorpp-data.R"))
library(mvtnorm)
dg <- rescorpp_data("g3")
set.seed(41); n <- 200
x <- rnorm(n)
R0 <- matrix(c(1, .7, .4, .7, 1, .2, .4, .2, 1), 3)
E <- MASS::mvrnorm(n, rep(0, 3), R0) / sqrt(rchisq(n, 6) / 6)
dt3 <- data.frame(x = x, y1 = 1 + x + E[, 1], y2 = 2 * E[, 2],
                  y3 = 0.3 * x + 0.5 * E[, 3])
f3 <- mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE)
out <- list()
for (w in c("gauss", "student")) {
  d <- if (w == "gauss") dg else dt3
  fit <- quiet(frm(f3, family = if (w == "gauss") gaussian() else student(),
                   data = d))
  ds <- quiet(frm_sample(fit, chains = 2, iter = 1000, refresh = 0, seed = 21))
  dr <- as.data.frame(posterior::as_draws_df(ds))
  S <- nrow(dr)
  rs <- c("y1", "y2", "y3")
  mu <- lapply(setNames(rs, rs), function(r) posterior_epred(ds, resp = r))
  if (w == "student") {
    # epred of a student is mu where nu > 1; read mu as linpred instead
    mu <- lapply(setNames(rs, rs), function(r) posterior_linpred(ds, resp = r))
  }
  sg <- lapply(setNames(rs, rs), function(r)
    posterior_epred(ds, resp = r, dpar = "sigma"))
  Y <- as.matrix(d[, rs])
  use <- round(seq(1, S, length.out = 60))
  hand <- function(sub) {
    t(vapply(use, function(s) {
      Rm <- diag(length(sub))
      for (a in seq_along(sub)) for (bb in seq_along(sub)) if (a != bb) {
        nm <- paste0("rescor__", min(sub[a], sub[bb]), "__", max(sub[a], sub[bb]))
        Rm[a, bb] <- dr[[nm]][s]
      }
      vapply(seq_len(n), function(i) {
        sgi <- vapply(sub, function(r) sg[[r]][s, i], 1)
        Sg <- diag(sgi, length(sub)) %*% Rm %*% diag(sgi, length(sub))
        m <- vapply(sub, function(r) mu[[r]][s, i], 1)
        if (w == "gauss") {
          dmvnorm(Y[i, sub], m, Sg, log = TRUE)
        } else {
          dmvt(Y[i, sub], m, Sg, df = dr$nu[s], log = TRUE)
        }
      }, 1)
    }, numeric(n)))
  }
  subs <- list("y1", "y2", "y3", c("y1", "y3"), c("y2", "y3"), c("y3", "y1"),
               c("y1", "y2", "y3"), c("y3", "y2", "y1"))
  for (sub in subs) {
    ll <- tryCatch(log_lik(ds, resp = sub), error = function(e) e)
    if (inherits(ll, "error")) {
      cat(w, paste(sub, collapse = ","), "ERROR", conditionMessage(ll), "\n")
      next
    }
    ref <- hand(sub)
    cat(sprintf("%-8s resp=%-9s dim=%s  max|diff|/max|ref| = %.3g\n", w,
                paste(sub, collapse = ","), paste(dim(ll), collapse = "x"),
                max(abs(ll[use, ] - ref)) / max(abs(ref))))
  }
  ll0 <- log_lik(ds)
  cat(w, "no resp == resp all three:", identical(ll0, log_lik(ds, resp = rs)),
      "; combine = FALSE identical:", identical(ll0, log_lik(ds, combine = FALSE)),
      "; combine = FALSE resp y1 dim:",
      paste(dim(log_lik(ds, combine = FALSE, resp = "y1")), collapse = "x"), "\n")
  ## the brms spelling of a name, and a duplicated name
  show(paste(w, "log_lik duplicated resp"),
       dim(log_lik(ds, resp = c("y1", "y1"))))
  ## loo / waic by resp against the loo package on the same matrix
  for (r in list(NULL, "y1", c("y1", "y3"))) {
    lw <- waic(ds, resp = r)
    lr <- suppressWarnings(loo::waic(log_lik(ds, resp = r)))
    lo <- suppressWarnings(loo(ds, resp = r))
    llr <- log_lik(ds, resp = r)
    chain <- rep(1:2, each = S / 2)
    lo_ref <- suppressWarnings(loo::loo(llr, r_eff = loo::relative_eff(
      exp(llr), chain_id = chain)))
    cat(sprintf("%-8s resp=%-6s waic %.4f (loo pkg %.4f)  loo %.4f (loo pkg %.4f)\n",
                w, paste(r, collapse = ","), lw$estimates["elpd_waic", 1],
                lr$estimates["elpd_waic", 1], lo$estimates["elpd_loo", 1],
                lo_ref$estimates["elpd_loo", 1]))
  }
  ## loo_compare across resp choices: brms warns "Not all models have
  ## the same y variable"
  l1 <- suppressWarnings(loo(ds, resp = "y1"))
  l2 <- suppressWarnings(loo(ds, resp = "y2"))
  l12 <- suppressWarnings(loo(ds, resp = c("y1", "y2")))
  show(paste(w, "loo_compare(y1, y2)"), loo::loo_compare(l1, l2))
  show(paste(w, "loo_compare(y1, joint y1y2)"), loo::loo_compare(l1, l12))
  show(paste(w, "loo(moment_match = TRUE, resp = y1)"),
       loo(ds, resp = "y1", moment_match = TRUE))
  show(paste(w, "reloo"), loo(ds, resp = "y1", reloo = TRUE))
  show(paste(w, "loo_R2 by resp"), loo_R2(ds))
  show(paste(w, "loo_predict resp y1 dim"), dim(loo_predict(ds, resp = "y1")))
  show(paste(w, "loo_predict no resp dim"), dim(loo_predict(ds)))
  show(paste(w, "kfold"), kfold(ds, resp = "y1", K = 2))
  show(paste(w, "loo_subsample"), loo_subsample(ds, resp = "y1"))
  out[[w]] <- ds
}
