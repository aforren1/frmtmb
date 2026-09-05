# Helpers for the brms PRIOR-placement tier (test-brms-priors.R).
#
# The flat-prior tier next door proves that frmtmb's objective is the
# same function of the parameters as brms's Stan program when neither
# side carries a prior. This tier asks the next question: what happens
# when both sides carry the priors brms's own get_prior() supplies.
#
# It reuses the translator and the compile cache from helper-brms.R
# (stan_pars_from_fit(), brms_stan_model(), brms_standata(),
# brms_ranef_table(), brms_inner_index()) and adds only what a prior needs.
# testthat sources helper files in alphabetical order, and
# "helper-brms-priors.R" sorts BEFORE "helper-brms.R", so nothing here
# may call those at source time. Everything below is a function body.
#
# See dev/brms-priors-findings.md for the measurements these pin.

# ---------------------------------------------------------------------
# What frm(prior =) does with one row of a brms get_prior() table.
# ---------------------------------------------------------------------

bp_dist_kind <- function(s) trimws(sub("[(].*$", "", s))

# The two gates are ASKED, never restated here. A local list of the
# accepted classes or of the parsed densities would freeze this file
# against the package: a decision that widens either gate has to move
# the statuses below, and it can only do that if the statuses are an
# observation of frmtmb rather than a parallel copy of its rules.
bp_class_ok <- function(cls, dist) {
  !inherits(try(frmtmb:::check_brms_prior_class(cls, dist),
                silent = TRUE), "try-error")
}

bp_dist_ok <- function(dist) {
  !inherits(try(frmtmb:::parse_prior_dist(dist), silent = TRUE),
            "try-error")
}

# One row of a brms table read by frmtmb's own translator, with
# `source` flipped to "user" so that the default-drop is out of the way
# and the row's own fate is what shows. Going through as_priorlist()
# rather than re-spelling the row means the target question follows
# whatever the class gate decides, including a routing decision that
# does not exist today.
bp_row_priorlist <- function(g, i) {
  row <- as.data.frame(g[i, , drop = FALSE], stringsAsFactors = FALSE)
  row$source <- "user"
  frmtmb:::as_priorlist(structure(row,
                                  class = c("brmsprior", "data.frame")))
}

# Five outcomes, and the order matters because a row can fail more than
# one way and only the first is reported:
#   "flat slot"             the prior cell is empty; not a prior at all
#   "refused: class"        check_brms_prior_class() refuses the class
#   "refused: distribution" parse_prior_dist() does not know the density
#   "refused: no target"    the class is accepted but addresses nothing
#                           in THIS model
#   "honored"               carried through to a resolved prior entry
#
# The order is as_priorlist()'s own: it runs the class gate before it
# spells the row, so a row that fails both reports the class.
#
# The `source` column is not consulted, here or in the package. It
# records who BUILT a row rather than who wrote the density in it, and
# reading it used to drop every one of these rows before its own fate
# could be asked about.
bp_classify_rows <- function(gp, fit = NULL) {
  g <- as.data.frame(gp)
  g$kind <- ifelse(nzchar(g$prior), bp_dist_kind(g$prior), "")
  g$class_ok <- FALSE
  g$dist_ok <- FALSE
  g$status <- "flat slot"
  live <- nzchar(g$prior)
  for (i in which(live)) {
    g$class_ok[[i]] <- bp_class_ok(g$class[[i]], g$prior[[i]])
    g$dist_ok[[i]] <- bp_dist_ok(g$prior[[i]])
    g$status[[i]] <- if (!g$class_ok[[i]]) {
      "refused: class"
    } else if (!g$dist_ok[[i]]) {
      "refused: distribution"
    } else {
      "honored"
    }
  }
  if (!is.null(fit)) {
    for (i in which(g$status == "honored")) {
      bad <- inherits(try(resolve_prior_input(fit, bp_row_priorlist(g, i)),
                          silent = TRUE), "try-error")
      if (bad) g$status[[i]] <- "refused: no target"
    }
  }
  g
}

# One frmtmb set_prior() call per named row. The rows are re-spelled
# rather than handed over as a frame so that each row's own fate is
# separable from the table's.
#
# A brms class that the gate accepts but set_prior() has no name for
# cannot be guessed at here: the translator is asked for the spelling
# frm(prior = ) would use instead. That is how a distributional class
# reaches its routed spelling, which is the one thing about this file
# that must not be restated locally.
bp_frm_prior <- function(g, idx) {
  out <- NULL
  for (i in idx) {
    one <- try(set_prior(g$prior[[i]], class = g$class[[i]],
                         coef = g$coef[[i]], group = g$group[[i]],
                         resp = g$resp[[i]], dpar = g$dpar[[i]],
                         nlpar = g$nlpar[[i]]), silent = TRUE)
    if (inherits(one, "try-error")) one <- bp_row_priorlist(g, i)
    out <- if (is.null(out)) one else out + one
  }
  out
}

# Is this row one the translator routes to the NATURAL placement, that
# is, a distributional parameter's own class? Asked of the package, so
# that a change to the routing set moves the row sets below with it.
bp_is_natural_row <- function(g, i) {
  pl <- try(bp_row_priorlist(g, i), silent = TRUE)
  if (inherits(pl, "try-error") || is.null(pl)) return(FALSE)
  any(vapply(unclass(pl), function(s) isTRUE(s$natural), TRUE))
}

# frmtmb's OWN spelling for the same parameter: class = "Intercept" with
# dpar = the brms class, which sits on the LINK scale. This is the
# spelling ?set_prior documents as diverging from brms, and the reason
# it still exists is that flipping it would change what an existing
# frmtmb script means (see dev/priors-findings.md). Measuring it beside
# the routed one is what shows the size of that divergence.
bp_dpar_link_prior <- function(g, idx) {
  out <- NULL
  for (i in idx) {
    one <- set_prior(g$prior[[i]], class = "Intercept",
                     dpar = g$class[[i]], resp = g$resp[[i]])
    out <- if (is.null(out)) one else out + one
  }
  out
}

bp_prior_c <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  a + b
}

# A brms prior table with every row outside `keep` blanked, so the Stan
# program carries exactly the prior set named and nothing else.
bp_stan_prior <- function(gp, keep = integer()) {
  gp$prior[!seq_len(nrow(gp)) %in% keep] <- ""
  gp
}

# ---------------------------------------------------------------------
# The densities, written out, so a residual can be attributed rather
# than only measured.
# ---------------------------------------------------------------------

# The numeric arguments of a brms prior string. They are data dependent
# (brms scales its default student_t by the spread of the response), so
# a test that hard-coded them would be pinning this data set rather than
# the placement rule.
bp_hyper <- function(p) {
  as.numeric(strsplit(gsub("^[a-z_]+[(]|[)]$", "", p), ",")[[1]])
}

# Stan's student_t_lpdf, summed over a vector.
bp_st <- function(x, df, mu, s) {
  sum(stats::dt((x - mu) / s, df = df, log = TRUE) - log(s))
}

# What brms writes for a lower-bounded parameter: the same density minus
# one lccdf renormalizer per element. For mu = 0 that renormalizer is
# log(1/2) whatever the degrees of freedom and whatever the scale, so
# the half-t costs exactly log(2) per element against the untruncated
# density frmtmb evaluates.
bp_half_st <- function(x, df, mu, s) {
  bp_st(x, df, mu, s) - length(x) *
    stats::pt(-mu / s, df = df, lower.tail = FALSE, log.p = TRUE)
}

bp_half_t_const <- function(k) k * log(2)

# The number of half-t elements in one Stan program's lprior block.
# brms writes a lower-bounded prior as the density minus one lccdf per
# element, with the count as a literal or as a standata name.
bp_half_t_count <- function(code, sdat) {
  ln <- strsplit(code, "\n", fixed = TRUE)[[1]]
  ln <- ln[grepl("_lccdf", ln, fixed = TRUE)]
  if (!length(ln)) {
    return(0L)
  }
  mult <- sub("^.*- ([A-Za-z0-9_]+) [*].*$", "\\1", ln)
  sum(vapply(mult, function(m) {
    v <- suppressWarnings(as.integer(m))
    if (is.na(v)) as.integer(sdat[[m]]) else v
  }, integer(1)))
}

# ---------------------------------------------------------------------
# The measurement.
# ---------------------------------------------------------------------

# frmtmb tapes the prior INTO fit$obj, and logLik() reports the
# penalized objective, so the only way to separate likelihood from
# penalty is a second tape of the same model with no prior. `fit0` is
# that tape, evaluated at the PENALIZED fit's parameters.
bp_objective <- function(fit, fit0, joint = FALSE, logJ = 0) {
  par <- fit$obj$env$last.par.best
  pen <- -fit$obj$env$f(par)
  ll <- -fit0$obj$env$f(par)
  if (joint) {
    pen <- pen + logJ
    ll <- ll + logJ
  }
  list(pen = pen, ll = ll, logprior = pen - ll)
}

# frmtmb's log prior, entry by entry, exactly as neg_log_prior_fn()
# assembles it. The labels carry the component, the index and the
# SCALE, which is the field this whole tier is about.
bp_prior_entries <- function(fit, prior) {
  ents <- resolve_prior_input(fit, prior)$entries
  # parList()'s FIRST argument is the fixed subvector, not the whole
  # parameter vector. Passing last.par.best positionally works on a
  # model with no random effects and silently mis-fills every component
  # on a model with them, which showed up as an sd entry of -106 nats
  # where -1.7 belonged.
  plist <- fit$obj$env$parList(par = fit$obj$env$last.par.best)
  data.frame(
    comp = vapply(ents, function(e) e$comp, ""),
    idx = vapply(ents, function(e) paste(e$idx, collapse = ","), ""),
    scale = vapply(ents, function(e) e$scale, ""),
    kind = vapply(ents, function(e) e$dist$kind, ""),
    # the entry may carry a centering offset over OTHER parameters, so
    # the value is read the way neg_log_prior_fn() reads it rather than
    # from the entry's own slot alone
    centered = vapply(ents, function(e) !is.null(e$offset), TRUE),
    value = vapply(ents, function(e) {
      as.numeric(sum(prior_logdens(plist[[e$comp]][e$idx], e$dist,
                                   e$scale, e$link,
                                   frmtmb:::entry_offset(e, plist))))
    }, numeric(1)),
    stringsAsFactors = FALSE)
}

# An uncompiled-data stanfit for one program, through the shared cache.
bp_stanfit <- function(code, sdat) {
  suppressMessages(suppressWarnings(
    rstan::sampling(brms_stan_model(code), data = sdat, chains = 0)))
}

# Checks A and B with priors on both sides, at frmtmb's estimates.
#
#   A  log_prob under BOTH adjust_transform settings against frmtmb's
#      PENALIZED objective. With flat priors only one of the two could
#      match; with priors on transformed parameters the question is
#      which, and that is what this tier reports.
#   B  grad_log_prob under both settings. The setting whose gradient
#      vanishes names the density frmtmb actually maximized.
#
# `sf_flat` is the same program with every prior blanked. Subtracting
# its log_prob at the SAME point isolates the prior contribution from
# the likelihood without assuming the flat-prior tier's result.
bp_check <- function(fit, fit0, sf, sf_flat, sdat, code, rtab = NULL,
                     joint = FALSE) {
  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  logJ <- attr(pars, "logJ")
  if (is.null(logJ)) logJ <- 0
  upars <- rstan::unconstrain_pars(sf, pars)
  # brms_z_index() became brms_inner_index() with a pattern when the
  # tier grew smooth, GP, mi() and CAR blocks; this tier has only z
  zidx <- brms_inner_index(sf, pars, pattern = "^z_[0-9]+$")
  gF <- rstan::grad_log_prob(sf, upars, adjust_transform = FALSE)
  gT <- rstan::grad_log_prob(sf, upars, adjust_transform = TRUE)
  lpF <- rstan::log_prob(sf, upars, adjust_transform = FALSE)
  lpT <- rstan::log_prob(sf, upars, adjust_transform = TRUE)
  lpF0 <- rstan::log_prob(sf_flat, upars, adjust_transform = FALSE)
  lpT0 <- rstan::log_prob(sf_flat, upars, adjust_transform = TRUE)
  ob <- bp_objective(fit, fit0, joint, logJ)
  list(pars = pars, upars = upars, zidx = zidx, logJ = logJ,
       lpF = lpF, lpT = lpT, lpF_flat = lpF0, lpT_flat = lpT0,
       # what each side calls the prior, isolated from the likelihood
       stan_prior_F = lpF - lpF0, stan_prior_T = lpT - lpT0,
       frm_prior = ob$logprior,
       # the Jacobian sum Stan adds over every constrained parameter
       jac = lpT - lpF,
       pen = ob$pen, ll = ob$ll,
       dF = lpF - ob$pen, dT = lpT - ob$pen,
       gF = max(abs(gF)), gT = max(abs(gT)),
       gFz = if (length(zidx)) max(abs(gF[zidx])) else NA_real_,
       gTz = if (length(zidx)) max(abs(gT[zidx])) else NA_real_)
}

# Everything one shape needs, built once: the row classification, the
# three Stan programs, and two frmtmb fits.
#
#   hon   every row frm(prior = ) accepts, in the meaning the
#         translator gives it. This is what a ported brms script gets.
#   link  the same set with each distributional-class row respelled in
#         frmtmb's OWN class = "Intercept" + dpar = spelling, which sits
#         on the LINK scale. Built only where such a row exists, and
#         only to measure the documented divergence between the two
#         spellings; it is not what a translated table does.
#
# Nothing is asserted here. The assertions belong in the test file,
# where the number being pinned is next to the reason it holds.
bp_shape <- function(bform, family, data, frm_model, joint = FALSE,
                     ...) {
  gp <- brms::get_prior(bform, data = data, family = family, ...)
  fit0 <- frm(frm_model, data = data)
  g <- bp_classify_rows(gp, fit0)
  hon <- which(g$status == "honored")
  # the honored rows the translator puts on the parameter ITSELF. Read
  # off the translator rather than from a list of class names here, so
  # that the row set follows the routing decision instead of restating
  # it.
  dpi <- hon[vapply(hon, function(i) bp_is_natural_row(g, i), TRUE)]
  mk <- function(pr) {
    brms::make_stancode(bform, data = data, family = family,
                        prior = pr, ...)
  }
  sdat <- brms_standata(bform, data = data, family = family,
                        prior = gp, ...)
  rtab <- brms_ranef_table(bform, data, family, gp, ...)
  code <- list(full = mk(gp), hon = mk(bp_stan_prior(gp, hon)),
               flat = mk(bp_stan_prior(gp)))
  sf <- lapply(code, bp_stanfit, sdat = sdat)
  pl <- list(hon = bp_frm_prior(g, hon),
             link = bp_prior_c(bp_frm_prior(g, setdiff(hon, dpi)),
                               bp_dpar_link_prior(g, dpi)))
  fitof <- function(p) {
    if (is.null(p)) fit0 else frm(frm_model, data = data, prior = p)
  }
  fit <- lapply(pl, fitof)
  out <- list(rows = g, honored = hon, dpar = dpi, prior = pl,
              fit = fit, fit0 = fit0, code = code, sdat = sdat,
              rtab = rtab, sf = sf,
              hon = bp_check(fit$hon, fit0, sf$hon, sf$flat, sdat,
                             code$hon, rtab, joint))
  if (length(dpi)) {
    out$link <- bp_check(fit$link, fit0, sf$hon, sf$flat, sdat,
                         code$hon, rtab, joint)
  }
  out
}
