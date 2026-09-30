# Lane sampfix: exact-text replacements for the nits round (S2): every
# laplace refusal names the function the user called.
WT <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/R/"
edit <- function(file, old, new) {
  f <- paste0(WT, file)
  s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, " times in ", file, ":\n", old)
  s <- sub(old, new, s, fixed = TRUE)
  writeLines(s, f)
}

edit("methods-draws.R",
"  out <- if (identical(scale, \"response\")) {
    posterior_epred(object, newdata = newdata, re_formula = re_formula,
                    resp = resp, dpar = dpar, nlpar = nlpar,
                    ndraws = ndraws, draw_ids = draw_ids, ...)
  } else {
    posterior_linpred(object, newdata = newdata, re_formula = re_formula,
                      resp = resp, dpar = dpar, nlpar = nlpar,
                      ndraws = ndraws, draw_ids = draw_ids, ...)
  }",
"  out <- draws_as_caller(\"fitted()\", if (identical(scale, \"response\")) {
    posterior_epred(object, newdata = newdata, re_formula = re_formula,
                    resp = resp, dpar = dpar, nlpar = nlpar,
                    ndraws = ndraws, draw_ids = draw_ids, ...)
  } else {
    posterior_linpred(object, newdata = newdata, re_formula = re_formula,
                      resp = resp, dpar = dpar, nlpar = nlpar,
                      ndraws = ndraws, draw_ids = draw_ids, ...)
  })")

edit("methods-draws.R",
"  out <- posterior_predict(object, newdata = newdata,
                           re_formula = re_formula, transform = transform,
                           resp = resp, negative_rt = negative_rt,
                           ndraws = ndraws, draw_ids = draw_ids, ...)",
"  out <- draws_as_caller(\"predict()\", posterior_predict(
    object, newdata = newdata, re_formula = re_formula,
    transform = transform, resp = resp, negative_rt = negative_rt,
    ndraws = ndraws, draw_ids = draw_ids, ...))")

edit("methods-draws.R",
"  out <- predictive_error(object, newdata = newdata,
                          re_formula = re_formula, method = method,
                          resp = resp, ndraws = ndraws,
                          draw_ids = draw_ids, ...)",
"  out <- draws_as_caller(\"residuals()\", predictive_error(
    object, newdata = newdata, re_formula = re_formula, method = method,
    resp = resp, ndraws = ndraws, draw_ids = draw_ids, ...))")

edit("methods-draws.R",
"    pp <- posterior_predict(object, newdata = newdata,
                            re_formula = re_formula, resp = resp,
                            ndraws = ndraws, draw_ids = draw_ids, ...)",
"    pp <- draws_as_caller(\"residuals()\", posterior_predict(
      object, newdata = newdata, re_formula = re_formula, resp = resp,
      ndraws = ndraws, draw_ids = draw_ids, ...))")

edit("methods-draws.R",
"  yrep <- if (identical(method, \"posterior_epred\")) {
    posterior_epred(object, newdata = newdata, re_formula = re_form,
                    resp = resp, ndraws = ndraws, draw_ids = draw_ids)
  } else {
    posterior_predict(object, newdata = newdata, re_formula = re_form,
                      resp = resp, ndraws = ndraws,
                      draw_ids = draw_ids)
  }",
"  yrep <- draws_as_caller(\"predictive_error()\",
                          if (identical(method, \"posterior_epred\")) {
    posterior_epred(object, newdata = newdata, re_formula = re_form,
                    resp = resp, ndraws = ndraws, draw_ids = draw_ids)
  } else {
    posterior_predict(object, newdata = newdata, re_formula = re_form,
                      resp = resp, ndraws = ndraws,
                      draw_ids = draw_ids)
  })")

edit("methods-draws.R",
"  yrep <- posterior_predict(object, newdata = newdata, resp = resp,
                            re_formula = re_form, ndraws = ndraws)",
"  yrep <- draws_as_caller(\"predictive_interval()\", posterior_predict(
    object, newdata = newdata, resp = resp, re_formula = re_form,
    ndraws = ndraws))")

edit("methods-draws.R",
"  yrep <- pred(object, newdata = newdata, resp = resp, draw_ids = rows,
               re_formula = re_formula)",
"  yrep <- draws_as_caller(\"pp_check()\", pred(
    object, newdata = newdata, resp = resp, draw_ids = rows,
    re_formula = re_formula))")

# pp_check's weights quote log_lik() by design, so log_lik() names itself
edit("methods-draws.R",
"  ll <- tryCatch(do.call(log_lik, lla), error = function(e) {",
"  ll <- tryCatch(draws_as_caller(\"log_lik()\", do.call(log_lik, lla),
                                 re_formula = FALSE, override = TRUE),
                 error = function(e) {")

edit("loo.R",
"    ep <- posterior_epred(object, resp = resps[r], ndraws = ndraws)",
"    ep <- draws_as_caller(\"bayes_R2()\", re_formula = FALSE,
                          posterior_epred(object, resp = resps[r],
                                          ndraws = ndraws))")

edit("loo.R",
"  ll <- log_lik(x, ndraws = ndraws, resp = resp)",
"  ll <- draws_as_caller(what, re_formula = FALSE,
                        log_lik(x, ndraws = ndraws, resp = resp))")

edit("loo.R",
"  crit <- lapply(models, function(m) {
    if (identical(criterion, \"loo\")) loo(m) else waic(m)
  })",
"  crit <- draws_as_caller(\"loo_compare()\", re_formula = FALSE,
                          lapply(models, function(m) {
    if (identical(criterion, \"loo\")) loo(m) else waic(m)
  }))")

edit("loo.R",
"  frm_stop(what, \" needs draws of the random effects, and these draws come \",
           \"from frm_sample(laplace = TRUE), which integrates them out \",
           \"instead of sampling them. The pointwise log-density is the one \",
           \"CONDITIONAL on each draw's own group-level values, so there is \",
           \"nothing left to condition on. Resample without laplace = TRUE\",
           call. = FALSE)",
"  frm_stop(draws_caller(what), \" needs draws of the random effects, and \",
           \"these draws come from frm_sample(laplace = TRUE), which \",
           \"integrates them out instead of sampling them. The pointwise \",
           \"log-density (log_lik()) is the one CONDITIONAL on each draw's \",
           \"own group-level values, so there is nothing left to condition \",
           \"on. Sample without laplace = TRUE\", call. = FALSE)")

edit("draws-brms.R",
"    frm_stop(\"ranef() and coef() have no draws of the group-level \",
             \"coefficients of '\", g, \"': these draws come from \",
             \"frm_sample(laplace = TRUE), which integrates the random \",
             \"effects out. Sample with laplace = FALSE for them\",
             call. = FALSE)",
"    frm_stop(draws_caller(\"ranef()\"), \" has no draws of the group-level \",
             \"coefficients of '\", g, \"': these draws come from \",
             \"frm_sample(laplace = TRUE), which integrates the random \",
             \"effects out. Sample without laplace = TRUE for them\",
             call. = FALSE)")

edit("methods-draws.R",
"  fe <- fixef(object, summary = FALSE)
  co <- ranef(object, summary = FALSE)",
"  fe <- fixef(object, summary = FALSE)
  co <- draws_as_caller(\"coef()\", re_formula = FALSE,
                        ranef(object, summary = FALSE))")
cat("edits done\n")
