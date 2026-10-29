## Review rescorpp: layout, dimnames and refusals of every draws method
## against brms 2.23.0 fits on the same data.
##   Rscript dev/review-rescorpp-layout.R lane <g3|g2f|st|gs|mix|re>
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
w <- commandArgs(TRUE)[2]
source(file.path(wt, "dev/rescorpp-data.R"))
source(file.path(wt, "dev/review-rescorpp-data.R"))
suppressMessages(loadNamespace("brms"))
own <- w %in% c("mix", "re", "ord")
d <- if (own) rv_data(w) else rescorpp_data(rescorpp_dname(w))
b <- readRDS(if (own) file.path(rv_out, paste0("brms-", w, ".rds")) else
  file.path(wt, "dev/rescorpp-brmsfits", paste0(w, ".rds")))
fit <- quiet(if (own) frm(rv_frm_form(w), data = d) else
  frm(eval(rescorpp_frm_forms[[w]]), family = rescorpp_family(w), data = d))
ds <- quiet(frm_sample(fit, chains = 2, iter = 1000, refresh = 0, seed = 21))
saveRDS(ds, file.path(rv_out, paste0("ds-", w, "-", arm, ".rds")))
cat("frmtmb draws:", posterior::ndraws(ds), " brms draws:", posterior::ndraws(b), "\n")
resps <- names(fit$spec$responses)
cat("responses:", resps, "\n")

desc <- function(x) {
  if (inherits(x, "error")) return(paste("ERROR:", conditionMessage(x)))
  if (inherits(x, "ggplot")) return("ggplot")
  if (inherits(x, "brms_conditional_effects")) {
    return(paste("ce:", paste(names(x), collapse = ",")))
  }
  if (inherits(x, "loo")) {
    return(paste("loo elpd", paste(names(x$estimates[, 1]), collapse = ","),
                 "pointwise", paste(dim(x$pointwise), collapse = "x")))
  }
  dn <- dimnames(x)
  dns <- if (is.null(dn)) "NULL" else paste(vapply(dn, function(e) {
    if (is.null(e)) "-" else if (length(e) > 6) paste0("[", length(e), "]")
    else paste(e, collapse = ",")
  }, ""), collapse = " | ")
  paste0(class(x)[1], " dim=", paste(dim(x) %||% length(x), collapse = "x"),
         " dn=", dns)
}
try_q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
}
nd <- d[1:5, , drop = FALSE]
cases <- list(
  pp = function(o, ...) posterior_predict(o, ndraws = 30, ...),
  pp_nd = function(o, ...) posterior_predict(o, newdata = nd, ndraws = 30, ...),
  epred = function(o, ...) posterior_epred(o, ndraws = 30, ...),
  epred_sigma = function(o, ...) posterior_epred(o, ndraws = 30, dpar = "sigma", ...),
  linpred = function(o, ...) posterior_linpred(o, ndraws = 30, ...),
  linpred_tr = function(o, ...) posterior_linpred(o, ndraws = 30, transform = TRUE, ...),
  perr = function(o, ...) predictive_error(o, ndraws = 30, ...),
  fitted = function(o, ...) fitted(o, ndraws = 30, ...),
  fitted_lin = function(o, ...) fitted(o, ndraws = 30, scale = "linear", ...),
  fitted_nosum = function(o, ...) fitted(o, ndraws = 30, summary = FALSE, ...),
  predict = function(o, ...) predict(o, ndraws = 30, ...),
  predict_nosum = function(o, ...) predict(o, ndraws = 30, summary = FALSE, ...),
  resid = function(o, ...) residuals(o, ndraws = 30, ...),
  resid_pearson = function(o, ...) residuals(o, ndraws = 30, type = "pearson", ...),
  pint = function(o, ...) predictive_interval(o, ndraws = 30, ...),
  loglik = function(o, ...) log_lik(o, ndraws = 30, ...),
  loglik_nocomb = function(o, ...) log_lik(o, ndraws = 30, combine = FALSE, ...),
  ppcheck = function(o, ...) pp_check(o, ndraws = 5, ...),
  bR2 = function(o, ...) bayes_R2(o, ndraws = 30, ...),
  ce_pred = function(o, ...) conditional_effects(o, method = "posterior_predict", ndraws = 30, ...),
  loo_pred = function(o, ...) loo_predict(o, ...),
  loo_ep = function(o, ...) loo_epred(o, ...),
  loo_pint = function(o, ...) loo_predictive_interval(o, ...)
)
rsets <- list(none = NULL, first = resps[1], two = resps[1:2],
              rev = rev(resps[1:2]))
if (length(resps) == 3) rsets$skip <- resps[c(1, 3)]
same <- 0; diff <- 0
for (cn in names(cases)) {
  for (rn in names(rsets)) {
    rr <- rsets[[rn]]
    set.seed(1)
    fb <- try_q(if (is.null(rr)) cases[[cn]](b) else cases[[cn]](b, resp = rr))
    set.seed(1)
    ff <- try_q(if (is.null(rr)) cases[[cn]](ds) else cases[[cn]](ds, resp = rr))
    db <- desc(fb); df <- desc(ff)
    tag <- if (identical(db, df)) "SAME" else {
      if (grepl("^ERROR", db) && grepl("^ERROR", df)) "BOTH-ERR" else "DIFF"
    }
    if (tag == "SAME") same <- same + 1 else diff <- diff + 1
    cat(sprintf("%-14s resp=%-5s %-8s\n   brms : %s\n   frmtmb: %s\n", cn, rn,
                tag, substr(db, 1, 220), substr(df, 1, 220)))
  }
}
cat("\nSUMMARY", w, "same", same, "not same", diff, "\n")

## draw-order and seed semantics on frmtmb's joint draw
if (isTRUE(fit$spec$rescor)) {
  set.seed(5); a1 <- posterior_predict(ds, ndraws = 40)
  set.seed(5); a2 <- posterior_predict(ds, ndraws = 40)
  cat("\nseeded repeat identical:", identical(a1, a2), "\n")
  e95 <- posterior_epred(ds, draw_ids = c(9, 5))
  e59 <- posterior_epred(ds, draw_ids = c(5, 9))
  cat("epred draw_ids c(9,5) == c(5,9) reversed:",
      identical(e95, e59[2:1, , , drop = FALSE]), "\n")
  set.seed(3); p95 <- posterior_predict(ds, draw_ids = c(9, 5))
  set.seed(3)
  ok <- rep(TRUE, nrow(d))
  man <- lapply(c(9, 5), function(id) {
    sh <- frmtmb.sample:::draws_fit_at(ds, id)
    predict_rescor_draw(sh, sh$spec$responses, eval_dpars(sh), ok)
  })
  cat("pp draw_ids c(9,5) equals manual draws in that order:",
      all(vapply(1:2, function(k) identical(unname(p95[k, , 1]), man[[k]][[1]]), NA)),
      "\n")
  rows <- frmtmb.sample:::draws_subsample(ds, 40, NULL)
  e40 <- posterior_epred(ds, ndraws = 40)
  cat("ndraws = 40 reads draws_subsample's rows (epred check):",
      identical(e40, posterior_epred(ds, draw_ids = rows)), "\n")
  r1 <- try_q(posterior_predict(ds, newdata = nd[1, , drop = FALSE], draw_ids = 3))
  rb <- try_q(posterior_predict(b, newdata = nd[1, , drop = FALSE], draw_ids = 3))
  cat("one row, one draw: frmtmb", desc(r1), " brms", desc(rb), "\n")
  pe1 <- try_q(predictive_error(ds, draw_ids = 3))
  pb1 <- try_q(predictive_error(b, draw_ids = 3))
  cat("predictive_error one draw: frmtmb", desc(pe1), " brms", desc(pb1), "\n")
}
