# The equivalence survey's measuring instrument.
#
# Sourced by dev/equiv-measure.R. It reuses the brms identity harness of
# the test suite (tests/testthat/helper-brms.R, helper-stan.R and
# helper-brms-methods.R) rather than restating it: the compile cache,
# the flat prior set, the parameter translator stan_pars_from_fit(), the
# ranef table, the inner-block index and the Fixed_param brmsfit are all
# theirs. What is added here is the part a survey needs and a test does
# not: every outcome is RECORDED instead of asserted, so one pair that
# fails a check yields a number and a verdict instead of an abort.

eq_root <- normalizePath(".", winslash = "/")
eq_helpers <- file.path(eq_root, "tests", "testthat",
                        c("helper-brms.R", "helper-stan.R",
                          "helper-brms-methods.R"))
for (h in eq_helpers) sys.source(h, envir = environment())

# Relative difference with a floor of 1 on the scale, the rule of the
# method tier's expect_exact_num(): a value near 0 is compared
# absolutely and a large one relatively.
eq_rel <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) return(Inf)
  if (!length(x)) return(NA_real_)
  if (any(is.na(x) != is.na(y))) return(Inf)
  ok <- !is.na(x)
  max(abs(x[ok] - y[ok]) / pmax(1, pmax(abs(x[ok]), abs(y[ok]))))
}
eq_tol_exact <- 1e-8

# A spec's code is text so that ONE spelling is evaluated by each
# package: frmtmb's bf() and set_prior() for frm(), brms's for brms.
eq_eval <- function(txt, pkg, data = NULL) {
  if (is.null(txt)) return(NULL)
  e <- new.env(parent = asNamespace(pkg))
  if (!is.null(data)) assign(".data", data, envir = e)
  # a spec may name one of the survey's own constructors (the custom
  # family); those live beside this function, not in either package
  home <- environment(sys.function())
  for (nm in ls(home, pattern = "^eq_bb2_")) {
    assign(nm, get(nm, envir = home), envir = e)
  }
  eval(parse(text = txt), envir = e)
}

# Conditions are collected rather than printed, so a pair's record says
# what each package said.
eq_try <- function(expr) {
  warns <- character()
  msgs <- character()
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(
      list(message = conditionMessage(e), class = class(e)),
      class = "eq_error")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(value = val, warnings = unique(warns), messages = unique(msgs),
       error = if (inherits(val, "eq_error")) val$message else NULL)
}
eq_clip <- function(x, n = 300) {
  x <- paste(x, collapse = " | ")
  x <- gsub("[\t\r\n]+", " ", x)
  if (nchar(x) > n) paste0(substr(x, 1, n), "...") else x
}

# The brms prior: the spec's rows where it has them and brms's flat ""
# everywhere else, so a prior the spec does not name is flat on both
# sides, as frm() leaves it.
eq_brms_prior <- function(spec, bform, family, data) {
  flat <- do.call(brms_flat_prior,
                  c(list(bform, data = data, family = family),
                    spec$brms_args))
  if (is.null(spec$prior)) return(flat)
  p <- eq_eval(spec$prior_brms %||% spec$prior, "brms")
  key <- function(t) paste(t$class, t$coef, t$group, t$resp, t$dpar,
                           t$nlpar, sep = "\r")
  keep <- !(key(flat) %in% key(p))
  # a class-level row overrides the flat class row and every coef row
  # under it, which brms would otherwise read as the more specific prior
  cls <- paste(p$class, p$resp, p$dpar, p$nlpar)[p$coef == "" &
                                                   p$group == ""]
  keep <- keep & !(paste(flat$class, flat$resp, flat$dpar, flat$nlpar) %in%
                     cls)
  rbind(p, flat[keep, , drop = FALSE])
}

eq_brms_args <- function(spec, data, prior) {
  c(list(data = data, prior = prior), spec$brms_args)
}

# ------------------------------------------------------------------ L1

# brms_lp_check() of helper-brms.R, with the prior made a parameter and
# the outcome returned instead of asserted. `points` adds evaluations
# away from the optimum: a fit whose optimizer stopped at the estimate
# plus a fixed step, so the translator reads that point as it reads the
# estimate. The value must keep the constant measured at the estimate;
# the gradient is read at the estimate only, where it must vanish.
eq_l1 <- function(spec, bform, family, data, fit, prior, joint,
                  refit = NULL) {
  args <- eq_brms_args(spec, data, prior)
  code <- do.call(brms::make_stancode,
                  c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  rtab <- do.call(brms_ranef_table,
                  c(list(bform, data, family, prior), spec$brms_args))
  mod <- brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  one <- function(f, away = FALSE) {
    pars <- stan_pars_from_fit(f, sdat, code, rtab)
    upars <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, upars, adjust_transform = FALSE,
                          gradient = FALSE)
    # away from the optimum the objective's best point is not the point
    # the fit reports, so it is evaluated there and its last point read
    if (away && joint) f$obj$fn(f$opt$par)
    ours <- if (joint) {
      env <- f$obj$env
      -env$f(if (away) env$last.par else env$last.par.best) +
        attr(pars, "logJ")
    } else if (!is.null(spec$prior)) {
      # the penalized objective, which carries the prior frm() applied
      -f$opt$objective
    } else {
      as.numeric(stats::logLik(f))
    }
    grad <- rstan::grad_log_prob(sf, upars, adjust_transform = FALSE)
    # frm()'s MAP carries the log-Jacobian of the map from each
    # prior-carrying parameter to its internal one (NEWS 0.36.0 on
    # class "sigma", 0.33.0 on "sdgp"/"lscale"; the ordered map on
    # thresholds); brms's program without adjust_transform does not.
    # Stan's own transform for those blocks is the same chart (log for
    # a positive parameter, ordered for thresholds), so its Jacobian is
    # added for the prior blocks alone: Stan's Jacobian is a sum over
    # blocks, so a block's share is the change in the total when that
    # block alone moves from the reference point.
    jadj <- 0
    if (length(pblocks)) {
      jt <- function(p) {
        u <- rstan::unconstrain_pars(sf, p)
        rstan::log_prob(sf, u, adjust_transform = TRUE) -
          rstan::log_prob(sf, u, adjust_transform = FALSE)
      }
      if (is.null(ref)) ref <<- pars
      mixed <- ref
      mixed[pblocks] <- pars[pblocks]
      jadj <- jt(mixed) - jt(ref)
      gt <- rstan::grad_log_prob(sf, upars, adjust_transform = TRUE)
      ix <- stan_upar_index(sf, pblocks)
      grad[ix] <- gt[ix]
    }
    if (joint) {
      idx <- brms_inner_index(sf, pars)
      grad <- if (length(idx)) grad[idx] else NA_real_
    }
    list(lp = lp, ours = ours, const = lp + jadj - ours, jadj = jadj,
         grad = max(abs(as.numeric(grad))))
  }
  ref <- NULL
  pblocks <- eq_prior_blocks(spec, brms_stan_par_names(code))
  at <- one(fit)
  extra <- list()
  if (is.function(refit)) {
    for (k in 1:2) {
      fk <- eq_try(refit(k))
      if (is.null(fk$error) && inherits(fk$value, "frmtmb_fit")) {
        ek <- eq_try(one(fk$value, away = TRUE))
        extra[[k]] <- if (is.null(ek$error)) ek$value else
          list(error = ek$error)
      } else {
        extra[[k]] <- list(error = fk$error %||% "refit failed")
      }
    }
  }
  list(at = at, extra = extra, code = code, sdat = sdat, sf = sf,
       rtab = rtab)
}

# The two away points: frm() run with an optimizer that converges and
# then steps by +/- 0.05 in every coordinate (the step a test of the
# value can see and a well-posed model survives). newton_polish() and
# the restarts would walk it back, so both are off for this refit.
eq_step_optimizer <- function(sign) {
  force(sign)
  function(par, fn, gr, lower, upper, control) {
    r <- stats::nlminb(par, fn, gr, lower = lower, upper = upper,
                       control = control)
    p <- r$par + sign * 0.05
    p <- pmin(pmax(p, lower), upper)
    list(par = p, objective = fn(p), convergence = 0L,
         message = "equiv: stepped off the optimum")
  }
}

# ------------------------------------------------------------------ L2

# The Fixed_param brmsfit of helper-brms-methods.R, at the frmtmb
# estimate, with the spec's extra brms arguments (data2) passed through.
eq_fixed <- function(spec, bform, family, data, fit, ndraws = 10) {
  do.call(brms_fixed_fit, c(list(bform, family, data, fit,
                                 ndraws = ndraws), spec$brms_args))
}

# Per-row Monte Carlo agreement of two simulators run at the same
# parameter values: the mean and the variance of each row, as z-scores
# against their Monte Carlo errors. With n rows, a max |z| under the
# two-sided Bonferroni point at 0.001 is agreement.
eq_sim_compare <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (ncol(a) != ncol(b)) return(list(maxz = Inf, note = "row count"))
  na <- nrow(a)
  nb <- nrow(b)
  ma <- colMeans(a)
  mb <- colMeans(b)
  va <- apply(a, 2, stats::var)
  vb <- apply(b, 2, stats::var)
  keep <- (va + vb) > 0
  zm <- (ma - mb)[keep] / sqrt(va[keep] / na + vb[keep] / nb)
  # the variance of a sample variance is (m4 - v^2) / n, with m4 the
  # row's fourth central moment, measured: the gaussian 2 v^2 / (n - 1)
  # is too small for a lognormal or a count row, which made agreeing
  # simulators look apart in the first run
  m4 <- function(s, m) colMeans(sweep(s, 2, m)^4)
  va4 <- pmax(m4(a, ma) - va^2, 0)
  vb4 <- pmax(m4(b, mb) - vb^2, 0)
  # a fair coin's sample variance has no spread to first order, so the
  # denominator is floored rather than left at zero
  den <- va4[keep] / na + vb4[keep] / nb
  den <- pmax(den, ((va + vb)[keep] / na)^2)
  zv <- (va - vb)[keep] / sqrt(den)
  crit <- stats::qnorm(1 - 0.001 / (4 * max(1, sum(keep))))
  list(maxz_mean = max(abs(zm)), maxz_var = max(abs(zv)), crit = crit,
       ok = max(abs(zm)) < crit && max(abs(zv)) < crit,
       const_rows = sum(!keep),
       const_agree = all(abs(ma[!keep] - mb[!keep]) < 1e-12))
}

# ------------------------------------------------------------------ L4

# brms sampled with flat priors, short chains, against the frm()
# estimate; read only where L1 and L2 cannot apply or for a method
# whose output is an interval (profile), which only a posterior
# answers. Compiles through the same cache.
eq_brms_sample <- function(spec, bform, family, data, prior,
                           iter = 2000, chains = 2, seed = 1) {
  args <- eq_brms_args(spec, data, prior)
  code <- do.call(brms::make_stancode,
                  c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  mod <- brms_stan_model(code)
  scaffold <- suppressMessages(do.call(brms::brm, c(
    list(bform, family = family, empty = TRUE), args)))
  sf <- suppressMessages(rstan::sampling(
    mod, data = sdat, chains = chains, iter = iter, seed = seed,
    refresh = 0, cores = 1))
  scaffold$fit <- sf
  brms::rename_pars(scaffold)
}

eq_quantile_z <- function(est, draws) {
  (est - mean(draws)) / stats::sd(draws)
}

# ------------------------------------------------------------ the flow

# ------------------------------------------- translator extensions
#
# Two block shapes the suite's translator was never given, added here
# rather than in tests/ because this lane changes no tracked test code.
# Each is a change of variables brms's Stan code states, so a wrong
# extension shows up as an L1 mismatch rather than as a quiet pass.

# mm(g1, g2) is brms's group "mmg1g2"; frmtmb labels the block by the
# call. brms names a grouping by pasting the call's text with its
# punctuation removed, so the same is done to frmtmb's label.
eq_orig_block_group <- brms_block_group
brms_block_group <- function(nm) {
  g <- eq_orig_block_group(nm)
  if (grepl("^(mm|gr)[(]", g)) {
    inner <- sub("^(mm|gr)[(](.*)[)]$", "\\2", g)
    args <- trimws(strsplit(inner, ",", fixed = TRUE)[[1]])
    args <- args[!grepl("=", args, fixed = TRUE)]
    g <- if (startsWith(g, "mm(")) {
      paste0("mm", paste(args, collapse = ""))
    } else {
      args[1]
    }
  }
  g
}

# gr(g, cov = A): brms draws r = (Lcov %x% diag(sd) L) z, so the
# standardized z the iid translation returns is multiplied by
# t(Lcov)^-1, and the change of variables adds M log|Lcov|.
#
# A mixture whose proportions have no predictor: brms declares
# simplex[K] theta; frmtmb carries theta1..theta(K-1) as intercepts on
# the logit (softmax, last component the reference) scale. The
# translator reads the parameter list off the code only, so the
# simplex is taken out of the list it reads and filled in here.
eq_orig_stan_pars <- stan_pars_from_fit
stan_pars_from_fit <- function(fit, sdat, code, rtab = NULL) {
  simplex_theta <- grepl("\n  simplex\\[[0-9]+\\] theta;", code)
  if (simplex_theta) {
    code <- sub("\n  simplex\\[[0-9]+\\] theta;[^\n]*", "", code)
  }
  pars <- eq_orig_stan_pars(fit, sdat, code, rtab)
  if (simplex_theta) {
    fe <- fixef_by_dpar(fit)
    a <- vapply(grep("^theta[0-9]+$", names(fe), value = TRUE),
                function(nm) unname(fe[[nm]][["(Intercept)"]]), 0)
    a <- c(a, 0)
    pars[["theta"]] <- exp(a) / sum(exp(a))
  }
  lj <- attr(pars, "logJ") %||% 0
  for (nm in grep("^Lcov_[0-9]+$", names(sdat), value = TRUE)) {
    i <- sub("^Lcov_", "", nm)
    zn <- paste0("z_", i)
    if (is.null(pars[[zn]])) next
    lc <- as.matrix(sdat[[nm]])
    z <- as.matrix(pars[[zn]])
    pars[[zn]] <- z %*% t(solve(lc))
    lj <- lj + nrow(z) * sum(log(diag(lc)))
  }
  attr(pars, "logJ") <- lj
  pars
}

# A VarCorr column keeps the call's spelling, mmc(x1, x2), where brms's
# coefficient is the same text with its punctuation removed, mmcx1x2.
eq_orig_frm_coef <- brms_frm_coef
brms_frm_coef <- function(cols, coef, dpar, resp = "", nlpar = "") {
  out <- tryCatch(eq_orig_frm_coef(cols, coef, dpar, resp, nlpar),
                  error = function(e) e)
  if (!inherits(out, "error")) return(out)
  bare <- gsub("[^A-Za-z0-9_.]", "", cols)
  hit <- cols[bare == coef]
  if (length(hit) == 1L) return(hit)
  stop(conditionMessage(out), call. = FALSE)
}

# A categorical or multinomial fit carries one (1 | g) block per
# non-reference category ("mub: 1 | g", "muc: 1 | g"), each with a
# column named (Intercept), so selecting blocks by grouping factor alone
# hands every category the first block's sd. brms keeps them as
# separate ids, one dpar each, so the selection also reads the dpar
# prefix of the label. The rest is brms_group_pars() unchanged.
eq_orig_group_pars <- brms_group_pars
brms_group_pars <- function(fit, sdat, rtab, i) {
  info <- brms_group_info(rtab, i)
  dp <- unique(info$dpars)
  blocks <- names(unclass(varcorr_matrices(fit)))
  pre <- paste0(dp, ": ")
  if (length(dp) == 1L && nzchar(dp) && dp != "mu" &&
        sum(startsWith(blocks, pre)) >= 1L &&
        sum(vapply(blocks, function(b) identical(brms_block_group(b),
                                                 info$group), NA)) > 1L) {
    vc <- varcorr_matrices
    on.exit(assign("varcorr_matrices", vc, envir = environment(
      brms_group_pars)), add = TRUE)
    assign("varcorr_matrices", function(f) {
      v <- unclass(vc(f))
      v[startsWith(names(v), pre)]
    }, envir = environment(brms_group_pars))
  }
  eq_orig_group_pars(fit, sdat, rtab, i)
}

# The Stan parameter blocks a spec's prior rows put a density on, for
# the Jacobian adjustment in eq_l1(). Classes on an unconstrained block
# (b, an unordered Intercept) are listed too: their Jacobian is zero.
eq_prior_blocks <- function(spec, stan_names) {
  if (is.null(spec$prior_brms %||% spec$prior)) return(character())
  p <- eq_eval(spec$prior_brms %||% spec$prior, "brms")
  # sds: a smooth's sd, on the log chart frmtmb gives every sd
  pat <- c(sd = "^sd_[0-9]+$", sigma = "^sigma$", Intercept = "^Intercept$",
           rescor = "^Lrescor$", sdgp = "^sdgp_", lscale = "^lscale_",
           cor = "^L_[0-9]+$", b = "^b(_.*)?$", nu = "^nu$",
           phi = "^phi$", shape = "^shape$", sds = "^sds_")
  out <- character()
  # a row with no density (a bound alone, set_prior("", lb = 0)) is a
  # constraint in frm() and carries no Jacobian (drop_flat_slots())
  p <- p[nzchar(p$prior), , drop = FALSE]
  for (cl in unique(p$class)) {
    if (!is.na(pat[cl])) {
      out <- c(out, grep(pat[[cl]], stan_names, value = TRUE))
    }
  }
  unique(out)
}

# L4 by optimization, for a model without random effects whose blocks
# the translator cannot map: brms's program under flat priors is
# maximized from several random starts (rstan::optimizing, no Jacobian,
# so its mode is the maximum likelihood point), and its best log
# density is compared with frmtmb's maximized log likelihood. A
# different model has a different maximum, so this needs no map
# between the two parameterizations. Mixtures and an improper flat nu
# leave a sampled posterior unconverged; the maximum does not care.
eq_l4_opt <- function(spec, bform, family, data, prior, fit, starts = 6) {
  args <- eq_brms_args(spec, data, prior)
  code <- do.call(brms::make_stancode, c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  mod <- brms_stan_model(code)
  vals <- vapply(seq_len(starts), function(s) {
    o <- tryCatch(suppressWarnings(rstan::optimizing(
      mod, data = sdat, seed = s, iter = 20000, as_vector = TRUE,
      init = if (s == 1L) "0" else "random", tol_rel_obj = 1e-12,
      tol_rel_grad = 1e-12)), error = function(e) NULL)
    if (is.null(o) || o$return_code != 0) NA_real_ else o$value
  }, 0)
  best <- suppressWarnings(max(vals, na.rm = TRUE))
  ours <- if (!is.null(spec$prior)) -fit$opt$objective else
    as.numeric(stats::logLik(fit))
  d <- best - ours - (spec$const %||% 0)
  list(ok = is.finite(d) && abs(d) < 1e-4 * max(1, abs(ours)),
       z = c(brms_max = best, frm_max = ours, diff = d),
       starts_ok = sum(is.finite(vals)),
       what = paste("brms's maximized log density under flat priors vs",
                    "frmtmb's maximized log likelihood"))
}
