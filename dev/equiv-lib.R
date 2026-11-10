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
                  refit = NULL, code_fix = NULL) {
  # a value no density reads has left the fit: brms's program on the
  # observed rows is the same model (brms_mi_observed())
  data <- brms_mi_observed(fit, data)
  args <- eq_brms_args(spec, data, prior)
  code <- do.call(brms::make_stancode,
                  c(list(bform, family = family), args))
  # brms's program with one recorded design choice made frmtmb's way (the
  # exact gp() nugget), so the identity reads every other line of it
  if (is.function(code_fix)) code <- code_fix(code)
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  # a change of brms's data that is a change of variables of its program
  # (dev/translator-bslogdet.R scales Xs so that eq_bs_logdet() is not 0);
  # no survey spec sets it
  if (is.function(spec$l1_sdat_fix)) sdat <- spec$l1_sdat_fix(sdat)
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
      as.numeric(stats::logLik(f)) + (attr(pars, "mi_ld") %||% 0)
    }
    grad <- brms_lp_grad(sf, upars, f, pars)
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
    out <- list(lp = lp, ours = ours, const = lp + jadj - ours, jadj = jadj,
                grad = max(abs(as.numeric(grad))))
    if (isTRUE(spec$l1_marginal)) {
      # last, because obj$fn() moves the objective's stored inner mode
      mg <- eq_laplace_stan(sf, pars)
      ours_m <- -f$obj$fn(f$opt$par)
      ld <- eq_bs_logdet(f, sdat)
      out$marg <- c(stan = mg$value, ours = ours_m, k = mg$k, bs_logdet = ld,
                    const = mg$value - ours_m - ld)
    }
    out
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

# Points that are not an optimum, for a model whose optimum the identity
# cannot read: an ordinal mixture that runs off to a degenerate boundary
# (Stan's gradient is NaN there and the densities saturate), or one whose
# optimum lies outside brms's parameter space. frm() is given an
# optimizer that returns its own starting values moved by a fixed
# offset, different for each k and in sign across coordinates, so every
# parameter, cs() offsets included, sits away from zero and from the
# others. The offset is deterministic, so a pair replays bitwise.
eq_jitter_optimizer <- function(scale, k) {
  force(scale)
  force(k)
  function(par, fn, gr, lower, upper, control) {
    p <- par + scale * sin(seq_along(par) * 1.7 + k)
    p <- pmin(pmax(p, lower), upper)
    list(par = p, objective = fn(p), convergence = 0L,
         message = "equiv: a fixed point off the start")
  }
}

# ------------------------------------------- REML: the marginal identity
#
# Under REML frm() maximizes the density integrated over the random
# effects AND mu's fixed effects. The joint identity above holds at
# whatever point the fit reports, so it says nothing about that integral,
# which is why a REML pair at L1 counted as weak evidence. For a gaussian
# response the joint is exactly quadratic in those blocks, so the Laplace
# value of brms's own program over the same blocks IS the integral, and it
# must equal minus frm()'s objective at every outer point. A spec asks
# for this with `l1_marginal = TRUE`.
eq_marg_pat <- "^(Intercept|b|bs|zs_[0-9]+_[0-9]+|z_[0-9]+)$"

# The Hessian comes from central differences of Stan's gradient: on an
# exact quadratic the gradient is linear, so the difference is exact for
# any step and a wide one keeps rounding out of it. The translated inner
# mode is a mode only to the inner tolerance, so one Newton step in
# Stan's own chart puts the value at the mode a quadratic has exactly.
eq_laplace_stan <- function(sf, pars, pattern = eq_marg_pat, h = 0.01) {
  u0 <- rstan::unconstrain_pars(sf, pars)
  ix <- brms_inner_index(sf, pars, pattern)
  g <- function(u) {
    rstan::grad_log_prob(sf, u, adjust_transform = FALSE)[ix]
  }
  H <- vapply(seq_along(ix), function(j) {
    up <- u0
    dn <- u0
    up[ix[j]] <- up[ix[j]] + h
    dn[ix[j]] <- dn[ix[j]] - h
    (g(up) - g(dn)) / (2 * h)
  }, numeric(length(ix)))
  H <- (H + t(H)) / 2
  u <- u0
  u[ix] <- u[ix] - solve(H, g(u0))
  lp <- rstan::log_prob(sf, u, adjust_transform = FALSE, gradient = FALSE)
  ld <- determinant(-H, logarithm = TRUE)
  if (ld$sign <= 0) stop("brms's inner Hessian is not negative definite")
  list(value = lp + length(ix) / 2 * log(2 * pi) - 0.5 * as.numeric(ld$modulus),
       k = length(ix))
}

# brms's unpenalized smooth columns (Xs) and frmtmb's (<term>.fx<j>) span
# one space in two bases, v_brms = A v_frm, so the integral over brms's
# coefficients is |det A| times frmtmb's: a constant of the design, which
# the marginal identity states rather than absorbs. REML integrates mu's
# coefficients only, so only mu's map enters.
eq_bs_logdet <- function(fit, sdat) {
  xb <- sdat[["Xs"]]
  if (is.null(xb) || !NCOL(xb)) return(0)
  fe <- brms_fe_of(fit, "mu")
  v <- fe[grepl("[.]fx[0-9]+$", names(fe))]
  xf <- as.matrix(brms_lp_of(fit, "mu")$X)[, names(v), drop = FALSE]
  a <- brms_basis_map(as.matrix(xb), xf, orthogonal = FALSE)
  as.numeric(determinant(a, logarithm = TRUE)$modulus)
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
  # the same in a block merged by |ID|, whose columns carry their linear
  # predictor: y.mu:mmc(x1, x2)
  dp <- if (nzchar(dpar)) dpar else "mu"
  pre <- sub(":.*$", "", cols)
  mine <- if (nzchar(resp)) {
    pre == paste0(resp, ".", dp)
  } else {
    pre == dp | endsWith(pre, paste0(".", dp))
  }
  hit <- cols[mine][gsub("[^A-Za-z0-9_.]", "",
                         sub("^[^:]*:", "", cols[mine])) == coef]
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

# ------------------------------------- translator extensions, round r081
#
# Rules for Stan parameters the suite's translator stops on, added here
# because this lane changes no tracked test code. The suite's translator
# reads the parameter list off the program's parameters block, so a
# parameter one of these rules owns is taken out of the program text it
# reads and filled in afterwards, as the helper already does for a
# mixture's simplex theta. A rule returns the value and its share of the
# map's log-Jacobian. Each is a change of variables brms's Stan code
# states, so a wrong rule shows up as an L1 mismatch, not a quiet pass.
eq_rules <- list()

# The program with the named parameters' declarations removed from its
# parameters block.
eq_strip_pars <- function(code, nms) {
  if (!length(nms)) return(code)
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  start <- grep("^parameters[[:space:]]*[{]", lines)[1]
  depth <- 0L
  drop <- integer()
  for (i in seq(start, length(lines))) {
    ln <- sub("//.*$", "", lines[[i]])
    depth <- depth + lengths(regmatches(ln, gregexpr("[{]", ln))) -
      lengths(regmatches(ln, gregexpr("[}]", ln)))
    for (nm in nms) {
      if (grepl(paste0("(^|[^A-Za-z0-9_])", nm, "[[:space:]]*;[[:space:]]*$"),
                ln)) {
        drop <- c(drop, i)
      }
    }
    if (i > start && depth <= 0L) break
  }
  if (length(drop) != length(nms)) {
    stop("could not find the declarations of ", paste(nms, collapse = ", "))
  }
  paste(lines[-drop], collapse = "\n")
}

# The inner blocks check C reads its gradient on: the suite's list, with
# a smooth's zs under a prefix (zs_b1_1_1, zs_y1_1_1) as zgp already is.
# Without it a smooth in an nl body or a multivariate response had no
# inner index, its gradient read NA, and an identity scored as a defect.
brms_inner_pat <- sub("zs_\\d+_\\d+", "zs_([A-Za-z0-9]+_)*\\d+_\\d+",
                      brms_inner_pat, fixed = TRUE)
stopifnot(grepl(brms_inner_pat, "zs_b1_1_1"),
          grepl(brms_inner_pat, "zs_y1_1_3"),
          grepl(brms_inner_pat, "zs_1_2"),
          !grepl(brms_inner_pat, "sds_b1_1"),
          !grepl(brms_inner_pat, "bs_b1"))

eq_r081_prev_stan_pars <- stan_pars_from_fit
stan_pars_from_fit <- function(fit, sdat, code, rtab = NULL) {
  need <- brms_stan_par_names(code)
  own <- character()
  rule_of <- list()
  for (nm in need) {
    for (r in eq_rules) {
      if (isTRUE(r$match(nm, fit, sdat, code))) {
        own <- c(own, nm)
        rule_of[[nm]] <- r
        break
      }
    }
  }
  pars <- eq_r081_prev_stan_pars(fit, sdat, eq_strip_pars(code, own), rtab)
  lj <- attr(pars, "logJ") %||% 0
  keep <- attributes(pars)[c("mi_ld", "mi_gone")]
  for (nm in own) {
    v <- rule_of[[nm]]$value(nm, fit, sdat, code, rtab, pars)
    pars[[nm]] <- v$value
    lj <- lj + (v$logJ %||% 0)
  }
  attr(pars, "logJ") <- lj
  for (a in names(keep)) attr(pars, a) <- keep[[a]]
  pars
}

# A smooth outside mu: an nl parameter's (zs_b1_1_1) or a response's in a
# multivariate model (zs_y1_1_1). The suite's rule reads its basis from
# Zs_<i>_<j>, which is mu's of a univariate model; brms names this one's
# Zs_<prefix>_<i>_<j>. The rest is the suite's rule.
eq_rules$zs_prefixed <- list(
  match = function(nm, fit, sdat, code) {
    grepl("^zs_", nm) && !identical(brms_idx_parts(nm, "zs")$dpar, "mu")
  },
  value = function(nm, fit, sdat, code, rtab, pars) {
    ip <- brms_idx_parts(nm, "zs")
    bk <- brms_smooth_term(fit, ip$dpar, ip$idx[[1]])[[ip$idx[[2]]]]
    sdv <- sqrt(brms_block_cov(fit, bk)[1, 1])
    zb <- sdat[[paste0("Zs_", ip$dpar, "_", ip$idx[[1]], "_", ip$idx[[2]])]]
    if (is.null(zb)) stop("brms's standata has no basis for ", nm)
    amat <- brms_basis_map(as.matrix(zb), brms_smooth_z(fit, ip$dpar, bk))
    v <- as.numeric(amat %*% brms_block_b(fit, bk)) / sdv
    list(value = array(v, length(v)), logJ = 0)
  })

# The linear-predictor frame of one dpar of a univariate fit.
eq_lp_frame <- function(fit, dp) {
  lps <- fit$frame[["linpreds"]]
  hit <- names(lps)[endsWith(names(lps), paste0(".", dp))]
  if (length(hit) != 1L) stop("no unique linear predictor for ", dp)
  lps[[hit]]
}

# An ordinal mixture under order = "none" (brms's default for ordinal
# components, which it announces): component k's thresholds are
# Intercept_mu<k>. frmtmb holds them in tau_raw<k>, read through that
# component's own family (its threshold structure and whether it stores
# increments). brms centers mu<k>'s design as for any ordinal model, so
# the shift enters with the ordinal sign (see the suite's Intercept rule).
eq_rules$ord_mix_intercept <- list(
  match = function(nm, fit, sdat, code) {
    grepl("^Intercept_mu[0-9]+$", nm) && !is.null(sdat[["nthres"]]) &&
      !is.null(family(fit)[["mix"]][["ord"]])
  },
  value = function(nm, fit, sdat, code, rtab, pars) {
    k <- as.integer(sub("^Intercept_mu", "", nm))
    ord <- family(fit)[["mix"]][["ord"]]
    if (isTRUE(ord$shared)) {
      stop("an ordinal mixture with shared thresholds (order = 'mu') has ",
           "no rule for ", nm)
    }
    v <- ord$views[[k]]
    raw <- fit$estimates[[ord$tau_names[k]]]
    tau <- if (!identical(v[["threshold"]] %||% "flexible", "flexible")) {
      ord_threshold_values(v, raw)
    } else {
      ord_tau_from_raw(raw, ordered = v[["family"]] %in%
                         c("cumulative", "hurdle_cumulative"))
    }
    dp <- paste0("mu", k)
    x <- sdat[[paste0("X_", dp)]]
    shift <- 0
    if (!is.null(x) && NCOL(x) > 0L) {
      fe <- brms_fe_of(fit, dp)
      cn <- colnames(x)
      shift <- sum(colMeans(x)[cn] * unname(fe[brms_coef_to_frm(cn)]))
    }
    list(value = as.numeric(tau - shift), logJ = 0)
  })

# Category-specific effects of one mixture component, bcs_mu<k>: brms's
# matrix[Kcs_mu<k>, nthres], one row per cs() column; frmtmb holds one
# bcs<j> vector per column, named in the component's linear-predictor
# frame. The suite's rule reads the unsuffixed Kcs and frmtmb's bcs<j>
# in number order, which would mix the components.
eq_rules$bcs_prefixed <- list(
  match = function(nm, fit, sdat, code) {
    grepl("^bcs_mu[0-9]+$", nm) && !is.null(family(fit)[["mix"]])
  },
  value = function(nm, fit, sdat, code, rtab, pars) {
    dp <- sub("^bcs_", "", nm)
    kcs <- sdat[[paste0("Kcs_", dp)]]
    nth <- sdat[["nthres"]]
    keys <- vapply(eq_lp_frame(fit, dp)[["cs"]], `[[`, "", "par")
    v <- unlist(lapply(keys, function(k) fit$estimates[[k]]))
    if (length(v) != kcs * nth) {
      stop("brms wants ", kcs * nth, " category-specific coefficients for ",
           dp, " but frmtmb holds ", length(v))
    }
    list(value = matrix(v, nrow = kcs, ncol = nth, byrow = TRUE), logJ = 0)
  })

# The generalized extreme value shape. brms declares an unconstrained
# tmp_xi and maps it into the interval that keeps every response inside
# the support at the current mu and sigma (fun_scale_xi.stan):
#   x = (Y - mu) / sigma, lb = min(-1/min(x), -1/max(x)),
#   ub = max(-1/min(x), -1/max(x)), xi = inv_logit(tmp_xi) (ub - lb) + lb.
# The inverse is read at frmtmb's mu and sigma. In a mixture brms does
# this per component over EVERY response, so a component's xi must keep
# all of y inside that component's support; a frmtmb estimate outside
# the interval has no brms value, which the rule says rather than clips.
# The value identity needs no Jacobian: log_prob without
# adjust_transform reads the density, and tmp_xi carries no prior.
eq_rules$tmp_xi <- list(
  match = function(nm, fit, sdat, code) grepl("^tmp_xi([0-9]+)?$", nm),
  value = function(nm, fit, sdat, code, rtab, pars) {
    k <- sub("^tmp_xi", "", nm)
    dp <- function(p) eval_dpars(fit)[[1L]][[paste0(p, k)]]
    mu <- as.numeric(dp("mu"))
    sigma <- as.numeric(dp("sigma"))
    xi <- unique(as.numeric(dp("xi")))
    if (length(xi) != 1L) {
      stop(nm, ": xi varies by row, which brms's tmp_xi cannot hold")
    }
    y <- as.numeric(sdat[["Y"]])
    x <- (y - mu) / sigma
    b <- c(-1 / min(x), -1 / max(x))
    lb <- min(b)
    ub <- max(b)
    p <- (xi - lb) / (ub - lb)
    if (!(p > 0 && p < 1)) {
      stop(sprintf(paste0("%s: frmtmb's xi%s = %.6g lies outside brms's ",
                          "interval (%.6g, %.6g), which keeps every ",
                          "response inside the component's support"),
                   nm, k, xi, lb, ub))
    }
    list(value = stats::qlogis(p), logJ = 0)
  })
