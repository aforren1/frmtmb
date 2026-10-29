# One pair's measurement, from spec to verdict. Sourced after
# dev/equiv-lib.R and dev/equiv-specs.R.

eq_has_re <- function(fit) length(fit$obj$env$random) > 0L

# NA is a spec whose families live inside its bf() calls, so neither
# package is given a family argument
eq_fam_text <- function(spec, side) {
  own <- spec[[paste0("family_", side)]]
  if (!is.null(own)) return(own)
  if (identical(spec$family, NA)) return(NA_character_)
  spec$family %||% "gaussian()"
}

eq_frm_fit <- function(spec, ctrl = NULL) {
  bform <- eq_eval(spec$bf_frm %||% spec$bf, "frmtmb")
  fam <- eq_fam_text(spec, "frm")
  fam <- if (is.na(fam)) NULL else eq_eval(fam, "frmtmb")
  prior <- eq_eval(spec$prior_frm %||% spec$prior, "frmtmb")
  a <- c(list(bform, data = spec$data), spec$frm_args)
  if (!is.null(fam)) a$family <- fam
  if (!is.null(prior)) a$prior <- prior
  if (length(spec$data2)) a$data2 <- spec$data2
  cc <- c(spec$ctrl, ctrl)
  if (length(cc)) a$control <- do.call(frmtmb::frmtmb_control, cc)
  do.call(frmtmb::frm, a)
}

# A method compared on both sides at the same parameter values. Each
# side runs on its own, so a method both refuse is a recorded refusal
# and one that only frmtmb refuses is a recorded mismatch.
eq_both <- function(fa, fb, cmp) {
  a <- eq_try(fa())
  b <- eq_try(fb())
  if (!is.null(a$error) && !is.null(b$error)) {
    return(list(ok = TRUE, refused_both = TRUE, metric = NA_real_,
                what = "both refuse the method",
                detail = eq_clip(c(a$error, b$error))))
  }
  if (!is.null(a$error)) {
    return(list(ok = FALSE, frm_refused = TRUE, metric = NA_real_,
                what = "frmtmb refuses the method; brms answers",
                detail = eq_clip(a$error)))
  }
  if (!is.null(b$error)) {
    return(list(ok = NA, brms_refused = TRUE, metric = NA_real_,
                what = "brms refuses the method; frmtmb answers",
                detail = eq_clip(b$error)))
  }
  out <- tryCatch(cmp(a$value, b$value), error = function(e) {
    list(ok = NA, metric = NA_real_,
         what = paste("the comparison itself failed:",
                      conditionMessage(e)))
  })
  out$warn_frm <- eq_clip(a$warnings, 200)
  out
}

eq_est_of <- function(x) {
  if (length(dim(x)) == 3L) x[, "Estimate", ] else x[, "Estimate"]
}
eq_first <- function(x) if (length(dim(x)) == 3L) x[1, , ] else x[1, ]
eq_exact <- function(what) {
  function(a, b) {
    d <- eq_rel(a, b)
    list(ok = is.finite(d) && d < eq_tol_exact, metric = d, what = what)
  }
}

# The comparison of one method at the frmtmb estimate. `bfit` is the
# Fixed_param brmsfit at that estimate; `bfit_n` the same with 4000
# draws, all equal, for a method that simulates.
eq_method1 <- function(m, a, fit, bfit, bfit_n) {
  switch(m,
    fitted = eq_both(
      function() eq_est_of(do.call(fitted, c(list(fit), a))),
      function() eq_first(do.call(brms::posterior_epred,
                                  c(list(bfit), a))),
      eq_exact("fitted Estimate vs posterior_epred at the estimate")),
    residuals = eq_both(
      function() eq_est_of(do.call(residuals,
                                   c(list(fit, type = "response"), a))),
      function() eq_first(suppressWarnings(do.call(residuals, c(
        list(bfit, method = "posterior_epred", type = "ordinary",
             summary = FALSE), a)))),
      eq_exact("response residual vs brms ordinary residual")),
    predict = , simulate = eq_both(
      function() {
        if (m == "predict") {
          do.call(predict, c(list(fit, summary = FALSE, ndraws = 4000,
                                  propagate_error = FALSE), a))
        } else {
          s <- do.call(simulate, c(list(fit, nsim = 4000, seed = 7), a))
          # an ordinal response simulates as a factor; its codes are brms's
          t(vapply(s, function(v) as.numeric(v), numeric(nrow(s))))
        }
      },
      function() do.call(brms::posterior_predict, c(list(bfit_n), a)),
      function(sa, sb) {
        if (length(dim(sa)) == 3L) sa <- matrix(sa, nrow = dim(sa)[1])
        if (length(dim(sb)) == 3L) sb <- matrix(sb, nrow = dim(sb)[1])
        cmp <- eq_sim_compare(sa, sb)
        if (is.null(cmp$crit)) {
          return(list(ok = FALSE, metric = Inf,
                      what = paste("draw shapes differ:", cmp$note)))
        }
        list(ok = isTRUE(cmp$ok) && isTRUE(cmp$const_agree),
             metric = max(cmp$maxz_mean, cmp$maxz_var) / cmp$crit,
             what = sprintf(paste0("per-row mean and variance z, 4000 ",
                                   "draws a side at the estimate; max |z| ",
                                   "%.2f (mean) %.2f (var), crit %.2f"),
                            cmp$maxz_mean, cmp$maxz_var, cmp$crit))
      }),
    emmeans = eq_both(
      function() as.data.frame(do.call(emmeans::emmeans, c(list(fit), a))),
      function() as.data.frame(do.call(emmeans::emmeans, c(list(bfit), a))),
      function(ea, eb) {
        col <- function(e) intersect(c("emmean", "rate", "prob",
                                       "response"), names(e))[1]
        d <- eq_rel(ea[[col(ea)]], eb[[col(eb)]])
        list(ok = is.finite(d) && d < eq_tol_exact, metric = d,
             what = paste("emmeans estimate", col(ea), "vs", col(eb)))
      }),
    hypothesis_profile = eq_both(
      function() hypothesis(fit, a$h, method = "profile")$hypothesis,
      function() brms::hypothesis(bfit, a$h)$hypothesis,
      function(ha, hb) {
        d <- eq_rel(ha$Estimate, hb$Estimate)
        list(ok = is.finite(d) && d < eq_tol_exact, metric = d,
             what = "hypothesis(method = 'profile') Estimate vs brms",
             ci = rbind(c(ha$CI.Lower, ha$CI.Upper)))
      }),
    confint_profile = eq_both(
      function() confint(fit, parm = a$parm, method = "profile"),
      function() TRUE,
      function(ci, b) {
        list(ok = all(is.finite(ci)), metric = NA_real_,
             what = "profile interval; read against the posterior (L4)",
             ci = ci, check_only = TRUE)
      }),
    frm_curve = , frm_curve_contrast = eq_both(
      function() {
        if (m == "frm_curve") {
          frm_curve(fit, newdata = a$newdata, simultaneous = FALSE)$.estimate
        } else {
          frm_curve(fit, newdata = a$newdata, contrast = a$contrast,
                    simultaneous = FALSE)$.estimate
        }
      },
      function() {
        lb <- function(nd) eq_first(brms::posterior_linpred(
          bfit, newdata = nd, re_formula = NA))
        if (m == "frm_curve") lb(a$newdata) else
          lb(a$newdata) - lb(a$contrast)
      },
      eq_exact("frm_curve .estimate vs posterior_linpred(re_formula = NA)")),
    list(ok = NA, metric = NA_real_,
         what = paste("no comparison written for method", m)))
}

eq_method <- function(spec, fit, bfit, bfit_n) {
  out <- list()
  for (m in spec$method) {
    a <- if (!is.null(spec$margs[[m]]) && is.list(spec$margs[[m]])) {
      spec$margs[[m]]
    } else if (length(spec$method) == 1L) {
      spec$margs %||% list()
    } else {
      list()
    }
    out[[m]] <- eq_method1(m, a, fit, bfit, bfit_n)
  }
  out
}

eq_measure <- function(key, row) {
  rec <- list(key = key, status = row$status)
  mk <- eq_specs[[key]]
  if (is.null(mk)) {
    rec$verdict <- "cannot measure"
    rec$reason <- "no spec written for this pair"
    return(rec)
  }
  spec <- mk()
  if (length(spec$data2)) spec$brms_args$data2 <- spec$data2
  # a prior spelled differently on the two sides (an alias test): the
  # brms text is the one the prior machinery reads as "a prior is set"
  if (is.null(spec$prior) && !is.null(spec$prior_brms)) {
    spec$prior <- spec$prior_brms
  }
  rec$spec <- spec[setdiff(names(spec), c("data", "data2"))]
  rec$nobs <- NROW(spec$data)
  if (!is.null(spec$skip)) {
    rec$verdict <- "cannot measure"
    rec$reason <- spec$skip
    return(rec)
  }

  ff <- eq_try(eq_frm_fit(spec))
  rec$frm_error <- ff$error
  rec$frm_warnings <- ff$warnings
  fit <- if (is.null(ff$error)) ff$value
  if (!is.null(fit) && !inherits(fit, "frmtmb_fit")) {
    rec$frm_error <- paste("frm() returned class", class(fit)[1])
    fit <- NULL
  }

  bform <- family <- prior <- NULL
  if (is.null(spec$nobrms)) {
    bb <- eq_try({
      bform <- eq_eval(spec$bf_brms %||% spec$bf, "brms")
      ft <- eq_fam_text(spec, "brms")
      family <- if (is.na(ft)) NULL else eq_eval(ft, "brms")
      prior <- eq_brms_prior(spec, bform, family, spec$data)
      do.call(brms::make_stancode,
              c(list(bform, family = family),
                eq_brms_args(spec, spec$data, prior)))
    })
    rec$brms_error <- bb$error
    rec$brms_warnings <- bb$warnings
  } else {
    rec$brms_error <- paste("no brms analogue:", spec$nobrms)
  }
  frm_ok <- is.null(rec$frm_error)
  brms_ok <- is.null(rec$brms_error)
  rec$frm_ok <- frm_ok
  rec$brms_ok <- brms_ok

  if (!frm_ok && !brms_ok) {
    if (!is.null(spec$nobrms)) {
      rec$verdict <- "cannot measure"
      rec$reason <- paste0("frmtmb refuses; no brms analogue (",
                           spec$nobrms, ")")
    } else if (row$status == "refused") {
      rec$verdict <- "refusal correct"
      rec$reason <- "both refuse"
      rec$agree <- TRUE
    } else {
      # the table claims the pair is reachable, and neither package
      # accepts any spelling of it that the spec could write
      rec$verdict <- "compat row wrong"
      rec$reason <- "both refuse; the table does not say refused"
      rec$agree <- TRUE
      rec$proposal <- "refused"
    }
    return(rec)
  }
  if (!frm_ok) {
    # make_stancode() succeeding is not brms fitting: a nonlinear body
    # that calls an unknown function is written into the program and
    # fails only at the compile, so brms "works" is confirmed by one
    # (cached) compile before frmtmb's refusal is called a divergence
    cc <- eq_try({
      args <- eq_brms_args(spec, spec$data, prior)
      code <- do.call(brms::make_stancode,
                      c(list(bform, family = family), args))
      brms_stan_model(code)
      TRUE
    })
    if (!is.null(cc$error)) {
      rec$brms_error <- paste("brms's program does not compile:",
                              eq_clip(cc$error, 200))
      rec$brms_ok <- FALSE
      if (row$status == "refused") {
        rec$verdict <- "refusal correct"
        rec$reason <- "both refuse (brms at the compile)"
        rec$agree <- TRUE
      } else {
        rec$verdict <- "compat row wrong"
        rec$reason <- "both refuse (brms at the compile); not refused here"
        rec$agree <- TRUE
        rec$proposal <- "refused"
      }
      return(rec)
    }
    if (!is.null(spec$decision)) {
      rec$verdict <- "divergence by decision"
      rec$reason <- paste("frmtmb refuses, brms fits;", spec$decision)
    } else {
      rec$verdict <- "defect"
      rec$reason <- "frmtmb refuses or errors where brms builds the model"
    }
    rec$agree <- FALSE
    if (row$status != "refused") rec$proposal <- "refused"
    return(rec)
  }
  if (!brms_ok) {
    if (row$status == "refused") {
      rec$verdict <- "compat row wrong"
      rec$proposal <- "works"
    } else {
      rec$verdict <- "cannot measure"
    }
    rec$reason <- if (!is.null(spec$nobrms)) {
      paste("frmtmb fits; no brms analogue:", spec$nobrms)
    } else {
      "frmtmb fits; brms refuses the same spelling"
    }
    return(rec)
  }

  # both build the model: measure
  level <- character()
  ok_all <- TRUE
  # An estimation mode that changes the estimator but not the model
  # (importance =) keeps its own objective, which is not the joint
  # density, so the model identity is read on the same spec without the
  # mode and the mode's estimate goes to L4 against the posterior.
  spec_l1 <- spec
  fit_l1 <- fit
  if (!is.null(spec$l1_without)) {
    spec_l1$frm_args[spec$l1_without] <- NULL
    f0 <- eq_try(eq_frm_fit(spec_l1))
    fit_l1 <- if (is.null(f0$error)) f0$value
    rec$l1_without <- spec$l1_without
  }
  joint <- if (is.null(fit_l1)) eq_has_re(fit) else eq_has_re(fit_l1)
  rec$joint <- joint
  refit <- function(k) {
    eq_frm_fit(spec_l1, ctrl = list(
      optimizer = eq_step_optimizer(c(1, -1)[k]), restarts = 0,
      check_se = "ignore"))
  }
  if (!isTRUE(spec$no_l1) && !is.null(fit_l1)) {
    l1 <- eq_try(eq_l1(spec, bform, family, spec$data, fit_l1, prior, joint,
                       refit = if (isTRUE(spec$no_points)) NULL else refit))
    if (is.null(l1$error)) {
      v <- l1$value
      # under a prior the constant is free: brms normalizes a density
      # truncated at a bound (minus log of its mass, log 2 for a half
      # normal) and frm()'s penalty need not, which moves no estimate.
      # What must hold is that it is the SAME at all three points.
      cst <- if (!is.null(spec$prior)) v$at$const else spec$const %||% 0
      rec$l1_const_free <- !is.null(spec$prior)
      tol <- 1e-6 * max(1, abs(v$at$ours))
      rec$l1_const <- v$at$const
      rec$l1_grad <- v$at$grad
      rec$l1_ours <- v$at$ours
      pts <- vapply(v$extra, function(e) {
        if (!is.null(e$error)) NA_real_ else e$const
      }, 0)
      rec$l1_extra_const <- pts
      rec$l1_extra_err <- unlist(lapply(v$extra, `[[`, "error"))
      val_ok <- abs(v$at$const - cst) < tol
      grad_ok <- is.finite(v$at$grad) && v$at$grad < 1e-3
      pts_ok <- all(is.na(pts) | abs(pts - cst) < tol)
      rec$l1_points <- 1L + sum(!is.na(pts))
      rec$l1_val_ok <- val_ok
      rec$l1_grad_ok <- grad_ok
      rec$l1_pts_ok <- pts_ok
      rec$l1_ok <- val_ok && grad_ok && pts_ok
      if (!rec$l1_ok && !is.null(spec$l1_decision)) {
        # the two programs differ by a recorded design choice (the exact
        # gp() nugget), so L1 is kept as a number and the verdict is
        # read from the level the choice does not touch
        rec$l1_by_decision <- spec$l1_decision
      } else {
        level <- c(level, "L1")
        ok_all <- ok_all && rec$l1_ok
      }
    } else {
      rec$l1_error <- l1$error
    }
  }
  if (!is.null(spec$method)) {
    bf1 <- eq_try(eq_fixed(spec, bform, family, spec$data, fit))
    needs_n <- any(spec$method %in% c("predict", "simulate"))
    bfn <- if (needs_n && is.null(bf1$error)) {
      eq_try(eq_fixed(spec, bform, family, spec$data, fit, ndraws = 4000))
    } else bf1
    if (is.null(bf1$error) && is.null(bfn$error)) {
      mm <- eq_method(spec, fit, bf1$value, bfn$value)
      rec$l2 <- lapply(mm, function(x) x[setdiff(names(x), "ci")])
      oks <- vapply(mm, function(x) as.logical(x$ok)[1], NA)
      rec$l2_ok <- if (anyNA(oks)) NA else all(oks)
      rec$l2_metric <- suppressWarnings(max(vapply(
        mm, function(x) as.numeric(x$metric)[1], 0), na.rm = TRUE))
      rec$l2_what <- paste(vapply(names(mm), function(n) {
        paste0(n, ": ", mm[[n]]$what,
               if (!is.null(mm[[n]]$detail)) paste0(" [", mm[[n]]$detail,
                                                    "]"))
      }, ""), collapse = " ; ")
      rec$l2_refused_both <- any(vapply(mm, function(x)
        isTRUE(x$refused_both), NA))
      # a finiteness check (confint_profile) compares nothing with brms;
      # its comparison is the posterior, at L4
      compared <- !vapply(mm, function(x) isTRUE(x$check_only), NA)
      if (!anyNA(oks)) {
        ok_all <- ok_all && rec$l2_ok
        if (any(compared)) {
          level <- c(level, "L2")
          rec$method_reached <- TRUE
        }
      }
      cis <- Filter(Negate(is.null), lapply(mm, `[[`, "ci"))
      if (length(cis)) rec$ci <- cis[[1]]
    } else {
      rec$l2_error <- bf1$error %||% bfn$error
    }
  }
  # L4: a profile interval against the flat-prior posterior, and any
  # pair the translator could not map
  want_l4 <- isTRUE(spec$l4) || !is.null(rec$ci) ||
    !is.null(spec$l1_without) ||
    (!length(level) && !is.null(rec$l1_error))
  if (want_l4 && !joint && !length(level) && !is.null(rec$l1_error) &&
        is.null(rec$ci) && is.null(spec$l1_without) && !isTRUE(spec$l4)) {
    # no map, no random effects: compare the two maxima instead
    want_l4 <- FALSE
    o4 <- eq_try(eq_l4_opt(spec, bform, family, spec$data, prior, fit))
    if (is.null(o4$error)) {
      rec$l4 <- o4$value
      rec$l4_route <- "optimizing"
      rec$l4_ok <- isTRUE(o4$value$ok)
      level <- c(level, "L4")
      ok_all <- ok_all && rec$l4_ok
    } else {
      rec$l4_error <- o4$error
    }
  }
  if (want_l4) {
    s4 <- eq_try(eq_brms_sample(spec, bform, family, spec$data, prior))
    if (is.null(s4$error)) {
      post <- s4$value
      rec$l4_rhat <- suppressWarnings(max(brms::rhat(post), na.rm = TRUE))
      l4 <- eq_try(eq_l4(spec, fit, post, rec))
      # the same reading for the fit without the mode is the control:
      # how far an ML point lands from this posterior when nothing is
      # wrong
      if (!is.null(spec$l1_without) && !is.null(fit_l1)) {
        l4c <- eq_try(eq_l4(spec, fit_l1, post, list()))
        if (is.null(l4c$error)) rec$l4_control <- l4c$value$z
      }
      if (is.null(l4$error) && is.finite(rec$l4_rhat) &&
            rec$l4_rhat > 1.05) {
        # an unconverged posterior is no reference, so the level is not
        # reached, whatever the comparison says
        rec$l4 <- l4$value
        rec$l4_error <- sprintf("brms posterior not converged (max Rhat %.2f)",
                                rec$l4_rhat)
      } else if (is.null(l4$error)) {
        rec$l4 <- l4$value
        rec$l4_ok <- isTRUE(l4$value$ok)
        level <- c(level, "L4")
        ok_all <- ok_all && rec$l4_ok
        # a profile interval read against the posterior is the method's
        # own comparison
        if (!is.null(rec$ci)) rec$method_reached <- TRUE
      } else {
        rec$l4_error <- l4$error
      }
    } else {
      rec$l4_error <- s4$error
    }
  }
  rec$level <- paste(level, collapse = "+")
  if (!length(level)) {
    rec$verdict <- "cannot measure"
    rec$reason <- paste("no level completed:",
                        eq_clip(c(rec$l1_error, rec$l2_error, rec$l2_what,
                                  rec$l4_error), 250))
    return(rec)
  }
  # a method pair whose method half reached no level is not measured,
  # whatever the model half says: the model alone would score the pair
  # on one of its two features (p029 did, before this rule)
  if (!is.null(spec$method) && !isTRUE(rec$method_reached)) {
    rec$verdict <- "cannot measure"
    rec$reason <- paste("the method half reached no level:",
                        eq_clip(c(rec$l2_error, rec$l2_what), 250))
    rec$model_ok <- ok_all
    return(rec)
  }
  rec$equivalent <- ok_all
  rec$agree <- ok_all
  if (row$status == "refused") {
    rec$verdict <- "compat row wrong"
    rec$reason <- paste("frmtmb fits a pair the table refuses,",
                        if (ok_all) "and agrees with brms" else
                          "and DISAGREES with brms")
    rec$proposal <- if (ok_all) "works" else "broken"
  } else {
    rec$verdict <- if (ok_all) "equivalent" else "defect"
    rec$reason <- if (ok_all) "every level agrees" else
      "a level disagrees; see the metrics"
  }
  rec
}

# L4: the frm() estimate against brms's flat-prior posterior. For a
# profile interval, its ends against the posterior quantiles; otherwise
# the fitted means against the posterior mean of posterior_epred. The
# tolerance is a ratio to the posterior sd: 0.25 for an interval end
# and 0.5 for a fitted mean. A wrong model misses by many sds, and an
# ML point and a flat-prior posterior at n = 200 meet well inside both.
eq_l4 <- function(spec, fit, post, rec) {
  a <- spec$margs[[spec$method[1]]] %||% spec$margs %||% list()
  if (!is.null(rec$ci) && "hypothesis_profile" %in% spec$method) {
    a <- spec$margs[["hypothesis_profile"]] %||% a
    hd <- brms::hypothesis(post, a$h)$samples[[1]]
    q <- stats::quantile(hd, c(0.025, 0.975))
    z <- (as.numeric(rec$ci) - q) / stats::sd(hd)
    return(list(ok = all(abs(z) < 0.25), z = z,
                what = "profile hypothesis CI ends vs posterior quantiles"))
  }
  if (!is.null(rec$ci)) {
    a <- spec$margs[["confint_profile"]] %||% a
    dr <- as.matrix(post, variable = a$bpar)
    ci <- as.matrix(rec$ci)
    pc <- which(colnames(ci) %in% c("lwr", "upr") |
                  grepl("%", colnames(ci), fixed = TRUE))
    ci <- ci[, if (length(pc) == 2L) pc else seq_len(2L), drop = FALSE]
    z <- numeric()
    for (j in seq_along(a$bpar)) {
      q <- stats::quantile(dr[, j], c(0.025, 0.975))
      z <- c(z, (ci[j, ] - q) / stats::sd(dr[, j]))
    }
    return(list(ok = all(abs(z) < 0.25), z = z,
                what = "profile CI ends vs posterior quantiles"))
  }
  ep <- brms::posterior_epred(post)
  if (length(dim(ep)) == 3L) ep <- matrix(ep, nrow = dim(ep)[1])
  fa <- fitted(fit)
  fa <- if (length(dim(fa)) == 3L) as.numeric(fa[, "Estimate", ]) else
    fa[, "Estimate"]
  z <- (fa - colMeans(ep)) / apply(ep, 2, stats::sd)
  z <- z[is.finite(z)]
  list(ok = max(abs(z)) < 0.5, z = stats::quantile(abs(z), c(0.5, 1)),
       what = "fitted at the ML point vs posterior mean of epred")
}

eq_flatten <- function(rec, row) {
  num <- function(x) if (is.null(x) || !length(x)) NA_real_ else
    as.numeric(x)[1]
  ov <- c(if ("rescor" %in% c(row$feature_a, row$feature_b)) "rescorpp",
          if ("me()" %in% c(row$feature_a, row$feature_b)) "meme")
  data.frame(
    id = row$id, feature_a = row$feature_a, kind_a = row$kind_a,
    feature_b = row$feature_b, kind_b = row$kind_b, status = row$status,
    draws = row$draws, weight = row$weight,
    verdict = rec$verdict %||% NA, level = rec$level %||% "",
    agree = rec$agree %||% NA,
    reason = eq_clip(rec$reason %||% ""),
    proposal = rec$proposal %||% "",
    frm_ok = rec$frm_ok %||% NA, brms_ok = rec$brms_ok %||% NA,
    frm_error = eq_clip(rec$frm_error %||% ""),
    brms_error = eq_clip(rec$brms_error %||% ""),
    l1_const = num(rec$l1_const), l1_grad = num(rec$l1_grad),
    l1_points = num(rec$l1_points), l1_ok = rec$l1_ok %||% NA,
    l1_extra = paste(signif(rec$l1_extra_const %||% NA, 3),
                     collapse = ","),
    l1_error = eq_clip(rec$l1_error %||% ""),
    l2_metric = num(rec$l2_metric), l2_ok = rec$l2_ok %||% NA,
    l2_what = eq_clip(rec$l2_what %||% "", 400),
    l2_error = eq_clip(rec$l2_error %||% ""),
    l4_ok = rec$l4_ok %||% NA,
    l4_z = eq_clip(if (is.null(rec$l4$z)) "" else
      paste(signif(rec$l4$z, 3), collapse = ",")),
    l4_rhat = num(rec$l4_rhat),
    l4_route = rec$l4_route %||% (if (!is.null(rec$l4)) "posterior" else ""),
    l4_control = paste(signif(rec$l4_control %||% NA, 3), collapse = ","),
    l1_without = rec$l1_without %||% "",
    l1_by_decision = eq_clip(rec$l1_by_decision %||% "", 120),
    method_reached = if (is.null(rec$spec$method)) NA else
      isTRUE(rec$method_reached),
    l4_error = eq_clip(rec$l4_error %||% ""),
    overlap = paste(ov, collapse = ","),
    frm_warn = eq_clip(rec$frm_warnings %||% "", 200),
    seconds = num(rec$seconds),
    frmtmb = rec$frmtmb_version %||% "",
    stringsAsFactors = FALSE)
}
