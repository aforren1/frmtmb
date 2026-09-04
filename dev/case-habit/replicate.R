## Replication of the response-preparation ("habit") model of
## Hardwick, Forrence, Krakauer & Haith (2019, Nat Hum Behav) in frmtmb,
## validated against the authors' own MATLAB implementation.
##
## Runnable end to end. It clones the reference repository into a scratch
## directory, translates the reference likelihood into R by hand, checks that
## translation against the archived MATLAB fits, then ports the model to frm().
##
## LICENSING: the reference repository carries no license file. Nothing it
## contains, code or data, may be copied into the frmtmb tree. This script
## derives everything from a scratch clone at run time and writes only
## numeric summaries into dev/case-habit/. See model-notes.md.
##
## Usage:
##   Rscript dev/case-habit/replicate.R [stage ...]
## with stages: validate, persubject, compare, corrected, hier, all.
##
## Needs R.matlab to read the archived .mat fits. That is a development-side
## tool only and must NOT be added to DESCRIPTION. Install it into a private
## library and put that library on the path, for example:
##   R_LIBS=/path/to/devlib Rscript dev/case-habit/replicate.R validate
## Use R_LIBS rather than R_LIBS_USER, which would displace the user library
## that frmtmb itself lives in.

suppressWarnings(suppressMessages({
  library(frmtmb)
  library(R.matlab)
}))

# ---------------------------------------------------------------- paths ----

SCRATCH <- Sys.getenv("HABIT_SCRATCH", unset = file.path(tempdir(), "habit"))
## Default the output directory to this script's own location so the tables
## land beside model-notes.md however the script is invoked.
OUT <- Sys.getenv("HABIT_OUT", unset = "")
if (!nzchar(OUT)) {
  args <- commandArgs(FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  OUT <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
}
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

UPSTREAM <- file.path(SCRATCH, "HabitTR")
CLEAN    <- file.path(SCRATCH, "HabitTR-clean")

ensure_clone <- function(dest, url, branch = NULL) {
  if (dir.exists(dest)) return(invisible(dest))
  args <- c("clone", "--depth", "1")
  if (!is.null(branch)) args <- c(args, "-b", branch)
  system2("git", c(args, url, shQuote(dest)), stdout = FALSE, stderr = FALSE)
  invisible(dest)
}
ensure_clone(UPSTREAM, "https://github.com/adrianhaith/HabitTR.git")
ensure_clone(CLEAN, "https://github.com/aforren1/HabitTR.git", branch = "clean_data")

CACHE <- function(name) file.path(SCRATCH, paste0("cache-", name, ".rds"))
cached <- function(name, expr) {
  f <- CACHE(name)
  if (file.exists(f)) return(readRDS(f))
  v <- force(expr); saveRDS(v, f); v
}

# ============================================================================
# 1. THE REFERENCE LIKELIHOOD, TRANSLATED BY HAND FROM MATLAB
#
#    habit_lik.m and getResponseProbs.m, with no frmtmb involvement. This is
#    the independent yardstick everything else is measured against.
#
#    params = [muA sigmaA qA muB sigmaB qB qInit rho]
#    (the header comment at habit_lik.m:10 gives a different order; it is
#    stale. Section 3 of model-notes.md shows this order is the right one.)
# ============================================================================

## Rows are response categories (1 correct, 2 habit, 3 other-single-key);
## columns are the four preparation states (neither, A only, B only, both).
hb_alpha <- function(p, model) {
  qA <- p[3]; qB <- p[6]; qI <- p[7]
  switch(model,
    "habit" = rbind(
      c(qI,       (1 - qA) / 3, qB,           qB),
      c(qI,       qA,           (1 - qB) / 3, (1 - qB) / 3),
      c(.5 - qI,  (1 - qA) / 3, (1 - qB) / 3, (1 - qB) / 3)),
    ## Replacing the "A only" column with the "neither" column makes PhiA drop
    ## out algebraically, which is what reduces this to a single process.
    "no-habit" = rbind(
      c(qI,      qI,      qB,           qB),
      c(qI,      qI,      (1 - qB) / 3, (1 - qB) / 3),
      c(.5 - qI, .5 - qI, (1 - qB) / 3, (1 - qB) / 3)),
    "flex-habit" = {
      rho <- p[8]
      rbind(
        c(qI,      rho * (1 - qA) / 3 + (1 - rho) * qI,        qB,           qB),
        c(qI,      rho * qA + (1 - rho) * qI,                  (1 - qB) / 3, (1 - qB) / 3),
        c(.5 - qI, rho * (1 - qA) / 3 + (1 - rho) * (.5 - qI), (1 - qB) / 3, (1 - qB) / 3))
    },
    stop("no such model: ", model))
}

hb_presponse <- function(RT, p, model) {
  PhiA <- pnorm(RT, p[1], p[2])
  PhiB <- pnorm(RT, p[4], p[5])
  hb_alpha(p, model) %*% rbind((1 - PhiA) * (1 - PhiB), PhiA * (1 - PhiB),
                               (1 - PhiA) * PhiB,       PhiA * PhiB)
}

## Their convention, reproduced deliberately: category 3 covers two keys but
## is scored with the probability of one. See model-notes.md section 1.
hb_LL <- function(RT, response, p, model) {
  pr <- hb_presponse(RT, p, model)
  sum(log(pr[cbind(response, seq_along(response))]))
}

## The published objective adds a ridge on both preparation-time SDs.
HB_AA <- 1000; HB_SLOPE0 <- 0.07
hb_nLL <- function(RT, response, p, model) {
  -hb_LL(RT, response, p, model) +
    HB_AA * (p[2] - HB_SLOPE0)^2 + HB_AA * (p[5] - HB_SLOPE0)^2
}

# ---------------------------------------------------------------- data -----

MODELS <- c("no-habit", "habit", "flex-habit")
COND   <- c("minimal", "4day", "20day")
LEV    <- c("correct", "habit", "other")

load_mat <- function(f) readMat(file.path(UPSTREAM, f))

## data(subject, cond) as archived by preprocess_data.m.
mat_trials <- function(dat, s, cc, which = c("remapped", "unchanged")) {
  which <- match.arg(which)
  key <- if (which == "remapped") c("RT", "response") else c("RT.unchanged", "response.unchanged")
  RT <- as.vector(dat[[key[1], s, cc]])
  if (!length(RT)) return(NULL)
  data.frame(t = RT, response = as.vector(dat[[key[2], s, cc]]))
}

as_frm_data <- function(tr) {
  data.frame(t = tr$t, resp = factor(LEV[tr$response], levels = LEV))
}

# ============================================================================
# 2. VALIDATE THE R REFERENCE AGAINST THE ARCHIVED MATLAB FITS
# ============================================================================

stage_validate <- function() {
  rows <- list()
  for (f in c("HabitModelFits.mat", "HabitModelFits_bads.mat")) {
    m <- load_mat(f); mo <- m$model; dat <- m$data
    for (mi in seq_along(MODELS)) for (cc in 1:3) for (s in 1:24) {
      tr <- mat_trials(dat, s, cc); if (is.null(tr)) next
      p <- mo[["paramsOpt", 1, mi]][s, , cc]
      rows[[length(rows) + 1]] <- data.frame(
        fitfile = f, model = MODELS[mi], cond = COND[cc], subj = s,
        LL_R = hb_LL(tr$t, tr$response, p, MODELS[mi]),
        LL_matlab = mo[["LLactual", 1, mi]][cc, s],
        nLL_R = hb_nLL(tr$t, tr$response, p, MODELS[mi]),
        nLL_matlab = mo[["LLopt", 1, mi]][cc, s])
    }
  }
  r <- do.call(rbind, rows)
  r$dLL <- r$LL_R - r$LL_matlab
  r$dnLL <- r$nLL_R - r$nLL_matlab
  agg <- do.call(rbind, lapply(split(r, list(r$fitfile, r$model)), function(x)
    data.frame(fitfile = x$fitfile[1], model = x$model[1], n = nrow(x),
               max_abs_dLL = max(abs(x$dLL)), max_abs_dnLL = max(abs(x$dnLL)))))
  list(per_fit = r, summary = agg)
}

# ============================================================================
# 3. THE frm() PORT
#
#    Two blockers at v0.47.0 shape this spelling:
#      * nl = TRUE demands a family with a single mu (R/parse.R:1225), and
#        categorical() has primary_dpars mu<level>, so bf(nl = TRUE) +
#        categorical() is refused outright (registered at R/compat.R:609).
#      * categorical() accepts the logit link only (R/families.R:4297) and
#        builds the softmax internally, so there is no identity-link hook to
#        hand it probabilities directly.
#
#    The way through is nlf(), which names a dpar and so bypasses the single-mu
#    gate, combined with a log-ratio parameterization: feeding
#    log(c2/c1) and log(2*c3/c1) to the two non-reference categories makes the
#    softmax return exactly (c1, c2, 2*c3), which is the properly normalized
#    version of their category probabilities.
# ============================================================================

## Shared preparation-state pieces. PhiA is omitted for no-habit because it
## cancels there, and leaving muA and sgA in the model would make them
## unidentified and the Hessian singular.
nlf_states <- function(model) {
  base <- list(nlf(PhiB ~ RTMB::pnorm((t - muB) / sgB)))
  if (model == "no-habit") return(base)
  c(list(nlf(PhiA ~ RTMB::pnorm((t - muA) / sgA))), base)
}

nlf_probs <- function(model) {
  switch(model,
    "no-habit" = list(
      nlf(c1 ~ qI * (1 - PhiB) + qB * PhiB),
      nlf(c2 ~ qI * (1 - PhiB) + ((1 - qB) / 3) * PhiB),
      nlf(c3 ~ (0.5 - qI) * (1 - PhiB) + ((1 - qB) / 3) * PhiB)),
    "habit" = list(
      nlf(c1 ~ qI * (1 - PhiA) * (1 - PhiB) + ((1 - qA) / 3) * PhiA * (1 - PhiB) + qB * PhiB),
      nlf(c2 ~ qI * (1 - PhiA) * (1 - PhiB) + qA * PhiA * (1 - PhiB) + ((1 - qB) / 3) * PhiB),
      nlf(c3 ~ (0.5 - qI) * (1 - PhiA) * (1 - PhiB) + ((1 - qA) / 3) * PhiA * (1 - PhiB) +
                ((1 - qB) / 3) * PhiB)),
    "flex-habit" = list(
      nlf(c1 ~ qI * (1 - PhiA) * (1 - PhiB) +
               (rho * (1 - qA) / 3 + (1 - rho) * qI) * PhiA * (1 - PhiB) + qB * PhiB),
      nlf(c2 ~ qI * (1 - PhiA) * (1 - PhiB) +
               (rho * qA + (1 - rho) * qI) * PhiA * (1 - PhiB) + ((1 - qB) / 3) * PhiB),
      nlf(c3 ~ (0.5 - qI) * (1 - PhiA) * (1 - PhiB) +
               (rho * (1 - qA) / 3 + (1 - rho) * (0.5 - qI)) * PhiA * (1 - PhiB) +
               ((1 - qB) / 3) * PhiB)))
}

## The two log-ratios the softmax consumes. Doubling c3 restores the second
## "other" key that their likelihood drops, so these probabilities normalize.
nlf_link <- function() list(nlf(muhabit ~ log(c2 / c1)), nlf(muother ~ log(2 * c3 / c1)))

MODEL_PARS <- list(
  "no-habit"   = c("muB", "sgB", "qB", "qI"),
  "habit"      = c("muA", "sgA", "qA", "muB", "sgB", "qB", "qI"),
  "flex-habit" = c("muA", "sgA", "qA", "muB", "sgB", "qB", "qI", "rho"))

## Box constraints of fit_habit_model.m:13-21.
HB_LOWER <- c(muA = 0,   sgA = .01, qA = .99,  muB = 0,  sgB = .01, qB = .5,    qI = .0001, rho = .0001)
HB_UPPER <- c(muA = .75, sgA = 100, qA = .999, muB = 10, sgB = 100, qB = .9999, qI = .499,  rho = .9999)
HB_INIT  <- c(muA = .4,  sgA = .05, qA = .99,  muB = .5, sgB = .05, qB = .95,   qI = .25,   rho = .95)

build_formula <- function(model) {
  pars <- MODEL_PARS[[model]]
  lin <- eval(parse(text = paste0("lf(", paste(pars, "~ 1", collapse = ", "), ")")))
  f <- bf(resp ~ 1)
  for (x in c(nlf_states(model), nlf_probs(model), nlf_link())) f <- f + x
  f + lin
}

## The ridge of habit_lik.m:47 is a Gaussian MAP penalty: 1/(2 sd^2) = 1000.
PEN_SD <- 1 / sqrt(2 * HB_AA)

## House style is to carry bounds in the prior vocabulary, and since v0.49
## that is available here: set_prior("", nlpar = "muA", lb = 0, ub = 0.75)
## addresses the nonlinear parameter and produces the same hard box the
## optimizer gets from it. It is a hard box, not a transform.
##
## Until v0.49 it did not. The bound was keyed by the design-matrix column
## name of the parameter's own sub-formula, which is "(Intercept)" for every
## one of them, so the bound named no outer parameter (the fit refused with
## "Unknown parameter(s) in bounds: (Intercept)") and, with several nonlinear
## parameters, all of them collided on that single key. The distribution path
## was never affected, which is why the two normal() penalties below always
## worked through nlpar.
##
## Note that lb/ub is a hard box and NOT a transform. That distinction
## matters for this stage: it reproduces a box-constrained fmincon run whose
## optima mostly sit exactly ON a bound (qA is at its limit for the large
## majority of participants, being unidentified from rho), and a logit
## transform would turn an attainable boundary into an asymptote. The
## hierarchical model below transforms inside each body because it wants a
## smooth interior; here the box is the point, and lb/ub gives one.
##
## Every parameter is declared intercept-only, so the bare names HB_LOWER
## carries and the prior's nlpar= spelling resolve to the same positions. The
## switch stays so the two can be checked against each other.
BOUNDS_VIA_PRIOR <- TRUE

build_prior <- function(model, penalize = TRUE) {
  pars <- MODEL_PARS[[model]]
  specs <- list()
  if (penalize) {
    for (p in intersect(c("sgA", "sgB"), pars)) {
      specs[[length(specs) + 1]] <-
        set_prior(paste0("normal(", HB_SLOPE0, ",", PEN_SD, ")"), nlpar = p)
    }
  }
  if (BOUNDS_VIA_PRIOR) {
    for (p in pars) {
      specs[[length(specs) + 1]] <-
        set_prior("", nlpar = p, lb = HB_LOWER[[p]], ub = HB_UPPER[[p]])
    }
  }
  if (!length(specs)) return(NULL)
  Reduce(`+`, specs)
}

fit_one <- function(d, model, penalize = TRUE, start = NULL) {
  pars <- MODEL_PARS[[model]]
  st <- list(beta = unname(if (is.null(start)) HB_INIT[pars] else start[pars]))
  args <- list(build_formula(model), family = categorical(levels = LEV),
               data = d, start = st, prior = build_prior(model, penalize))
  if (!BOUNDS_VIA_PRIOR) {
    args$lower <- HB_LOWER[pars]; args$upper <- HB_UPPER[pars]
  }
  do.call(frm, args)
}

## frmtmb's logLik() on a penalized fit returns LL + log prior density, so it
## is not their LLactual. Recover the comparable number by evaluating the
## hand-written reference at frmtmb's own estimates.
frm_params <- function(fit, model) {
  fx <- fixef(fit)
  p <- HB_INIT
  for (nm in MODEL_PARS[[model]]) p[nm] <- fx[[nm]][[1]]
  if (model == "no-habit") { p["muA"] <- 0; p["sgA"] <- HB_SLOPE0 }
  p
}

frm_LL <- function(fit, d, model) {
  p <- frm_params(fit, model)
  hb_LL(d$t, as.integer(d$resp), unname(p), model)
}

at_bound <- function(fit, model, tol = 1e-6) {
  p <- frm_params(fit, model)[MODEL_PARS[[model]]]
  nm <- names(p)
  paste(nm[abs(p - HB_LOWER[nm]) < tol | abs(p - HB_UPPER[nm]) < tol], collapse = ",")
}

# ---------------------------------------------- per-subject replication ----

stage_persubject <- function(fitfile = "HabitModelFits.mat", models = MODELS,
                             conds = 1:3, max_subj = 24) {
  m <- load_mat(fitfile); mo <- m$model; dat <- m$data
  rows <- list()
  for (mi in seq_along(MODELS)) {
    model <- MODELS[mi]; if (!(model %in% models)) next
    for (cc in conds) for (s in seq_len(max_subj)) {
      tr <- mat_trials(dat, s, cc); if (is.null(tr)) next
      d <- as_frm_data(tr)
      pm <- mo[["paramsOpt", 1, mi]][s, , cc]
      fit <- tryCatch(fit_one(d, model), error = function(e) e)
      if (inherits(fit, "error")) {
        rows[[length(rows) + 1]] <- data.frame(model = model, cond = COND[cc], subj = s,
          n = nrow(d), ok = FALSE, err = conditionMessage(fit))
        next
      }
      pf <- frm_params(fit, model)
      ll_frm <- frm_LL(fit, d, model)
      ll_mat <- mo[["LLactual", 1, mi]][cc, s]
      row <- data.frame(model = model, cond = COND[cc], subj = s, n = nrow(d), ok = TRUE,
                        err = NA_character_,
                        LL_frmtmb = ll_frm, LL_matlab = ll_mat, dLL = ll_frm - ll_mat,
                        bounds = at_bound(fit, model))
      for (nm in MODEL_PARS[[model]]) {
        row[[paste0(nm, "_frmtmb")]] <- unname(pf[nm])
        row[[paste0(nm, "_matlab")]] <- unname(pm[match(nm, names(HB_INIT))])
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  merge_rows(rows)
}

merge_rows <- function(lst) {
  nms <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) {
    for (n in setdiff(nms, names(x))) x[[n]] <- NA
    x[nms]
  }))
}

# ============================================================================
# 4. MODEL COMPARISON: habit vs no-habit, per subject, per group
#
#    Their AIC uses the unpenalized likelihood evaluated at the penalized
#    optimum (fit_habit_model.m:102-103). For a likelihood-ratio test the
#    penalty must either be absent or identical across the two models, and it
#    is not (habit carries a ridge on sgA that no-habit has no room for), so
#    the comparison fits are unpenalized.
# ============================================================================

## Nesting, which decides which comparisons are legitimate:
##   flex-habit contains habit    at rho = 1  (a boundary of the parameter space)
##   flex-habit contains no-habit at rho = 0  (also a boundary, and muA and sgA
##                                             become unidentified there)
##   habit and no-habit are NOT nested in each other: collapsing the "A only"
##   column of habit onto the "neither" column would need qA = qI = 0.25, and
##   qA is bounded at or above 0.99.
## So habit vs no-habit is an AIC comparison only, which is what the paper
## reports, and the likelihood-ratio test is habit vs flex-habit, which is what
## LR_test.m runs. Note that frmtmb's anova() does not itself verify nesting;
## it requires only equal n_obs, so choosing the legitimate pair is on the user.
stage_compare <- function(fitfile = "HabitModelFits.mat", conds = 1:3) {
  m <- load_mat(fitfile); dat <- m$data
  rows <- list()
  for (cc in conds) for (s in 1:24) {
    tr <- mat_trials(dat, s, cc); if (is.null(tr)) next
    d <- as_frm_data(tr)
    f0 <- tryCatch(fit_one(d, "no-habit",   penalize = FALSE), error = function(e) NULL)
    f1 <- tryCatch(fit_one(d, "habit",      penalize = FALSE), error = function(e) NULL)
    f2 <- tryCatch(fit_one(d, "flex-habit", penalize = FALSE), error = function(e) NULL)
    if (is.null(f0) || is.null(f1)) next
    a <- if (is.null(f2)) NULL else tryCatch(anova(f1, f2), error = function(e) NULL)
    rows[[length(rows) + 1]] <- data.frame(
      cond = COND[cc], subj = s, n = nrow(d),
      LL_nohabit = frm_LL(f0, d, "no-habit"),
      LL_habit   = frm_LL(f1, d, "habit"),
      LL_flex    = if (is.null(f2)) NA else frm_LL(f2, d, "flex-habit"),
      AIC_nohabit = AIC(f0), AIC_habit = AIC(f1),
      AIC_flex = if (is.null(f2)) NA else AIC(f2),
      ## positive favors habit, matching the sign of the paper's figures
      dAIC = AIC(f0) - AIC(f1),
      ## habit vs flex-habit, the one nested pair, as in LR_test.m
      LR_chisq = if (is.null(a)) NA else a[["Chisq"]][2],
      LR_df    = if (is.null(a)) NA else a[["Chi Df"]][2],
      LR_p     = if (is.null(a)) NA else a[["Pr(>Chisq)"]][2])
  }
  do.call(rbind, rows)
}

# ============================================================================
# 5. THE CORRECTED TRIAL SUBSETTING
#
#    preprocess_data.m:59-66 splits remapped from unchanged with
#    ismember(recodedX, revisedX), which matches on millisecond-rounded prep
#    time VALUES rather than trial indices. Trials whose prep time collides
#    across the two stimulus classes are admitted to both datasets. This stage
#    rebuilds the remapped dataset by trial identity from the cleaned CSVs and
#    refits, to measure what the defect cost.
# ============================================================================

CLEAN_FILES <- c(minimal = "data_1day.csv", `4day` = "data_5day.csv", `20day` = "data_20day.csv")
CLEAN_SESS  <- c(minimal = 4L, `4day` = 8L, `20day` = 26L)

load_clean <- function(cond) {
  f <- file.path(CLEAN, "clean_data", CLEAN_FILES[[cond]])
  d <- read.csv(f)
  ## Issue #1: three trials for subject 4 of experiment 2 were duds.
  if (cond == "20day") {
    for (st in list(c(10, 1073), c(15, 790), c(26, 758))) {
      i <- which(d$Subject == 104 & d$Session == st[1] & d$Trial == st[2])
      if (length(i)) d$PrepTime[i] <- NA
    }
  }
  d <- d[d$Session == CLEAN_SESS[[cond]] & d$IsFreeResp == 0 &
         d$IsToCriterion == 0 & d$IsHand == 0, ]
  ## PR #3: mishits, user errors, and the -999 coding-error sentinel.
  d <- d[!is.na(d$PrepTime) & !is.na(d$ResponseButton) &
         d$ResponseButton <= 4 & d$PrepTime < 5 & d$PrepTime > 0, ]
  d$response <- ifelse(d$ResponseButton == d$NewButton, 1L,
                ifelse(d$ResponseButton == d$OldButton, 2L, 3L))
  d$t <- d$PrepTime
  d
}

stage_corrected <- function(conds = names(CLEAN_FILES)) {
  rows <- list()
  for (cond in conds) {
    dd <- load_clean(cond)
    rm_ <- dd[dd$WillRemap == 1, ]
    for (sb in sort(unique(rm_$Subject))) {
      x <- rm_[rm_$Subject == sb, ]
      d <- data.frame(t = x$t, resp = factor(LEV[x$response], levels = LEV))
      f0 <- tryCatch(fit_one(d, "no-habit", penalize = FALSE), error = function(e) NULL)
      f1 <- tryCatch(fit_one(d, "habit",    penalize = FALSE), error = function(e) NULL)
      if (is.null(f0) || is.null(f1)) next
      rows[[length(rows) + 1]] <- data.frame(
        cond = cond, subject = sb, n = nrow(d),
        dAIC = AIC(f0) - AIC(f1),
        logLik_nohabit = as.numeric(logLik(f0)), logLik_habit = as.numeric(logLik(f1)))
    }
  }
  do.call(rbind, rows)
}

# ============================================================================
# 6. THE VALUE-ADD: ONE HIERARCHICAL FIT PER GROUP
#
#    Their procedure fits each subject in isolation, so it cannot borrow
#    strength or estimate how the preparation-time parameters covary across
#    people. This is the model asked for in HabitTR issue #2 (July 2017).
#
#    Random effects need an unconstrained scale, so the bounded parameters are
#    reparameterized by hand inside the bodies rather than box-constrained.
# ============================================================================

## Which working parameters each model needs, in the order their fixed effects
## enter `start`.
HIER_PARS <- list(
  "no-habit" = c("lmuB", "lsgB", "lqB", "lqI"),
  "habit"    = c("lmuA", "lmuB", "lsgA", "lsgB", "lqA", "lqB", "lqI"))

## Working scale for each bounded parameter, so random effects live somewhere
## unbounded. The box constraints of the per-subject fits become these
## transforms; the bounds themselves are the ones in HB_LOWER/HB_UPPER.
##
## An earlier version wrote muA as muB - exp(ldmu), to carry over the
## mu_A <= mu_B inequality that fmincon imposed as a linear constraint. That
## parameterization does not survive a random effect: the offset is only weakly
## identified, its variance component ran to a boundary, and two of three
## groups failed to converge. The means are therefore modeled directly, which
## is also what HabitTR issue #2 sketches. Ordering is checked after the fact
## rather than imposed.
HIER_TRANSFORMS <- list(
  muA = quote(lmuA),
  muB = quote(lmuB),
  sgA = quote(exp(lsgA)),
  sgB = quote(exp(lsgB)),
  qA  = quote(0.99 + 0.009 / (1 + exp(-lqA))),
  qB  = quote(0.5 + 0.4999 / (1 + exp(-lqB))),
  qI  = quote(0.499 / (1 + exp(-lqI))))

HIER_START <- c(lmuA = 0.4, lmuB = 0.5, lsgA = log(0.05), lsgB = log(0.05),
                lqA = 0, lqB = 2.2, lqI = 0)

## Preparation-time parameters carry correlated subject effects. The asymptote
## and guess parameters stay pooled: they are weakly identified per subject and
## giving them their own effects makes the block singular.
HIER_RE <- c("lmuA", "lmuB", "lsgA", "lsgB")

hier_formula <- function(model, re = "(1 | p | subject)", pooled = FALSE) {
  need <- switch(model,
    "no-habit" = c("muB", "sgB", "qB", "qI"),
    "habit"    = c("muB", "muA", "sgA", "sgB", "qA", "qB", "qI"))
  f <- bf(resp ~ 1)
  for (nm in need) {
    f <- f + eval(call("nlf", call("~", as.name(nm), HIER_TRANSFORMS[[nm]])))
  }
  for (x in c(nlf_states(model), nlf_probs(model), nlf_link())) f <- f + x
  terms <- vapply(HIER_PARS[[model]], function(p)
    if (!pooled && p %in% HIER_RE) paste0(p, " ~ 1 + ", re) else paste0(p, " ~ 1"),
    character(1))
  f + eval(parse(text = paste0("lf(", paste(terms, collapse = ", "), ")")))
}

hier_data <- function(dat, cc) {
  parts <- list()
  for (s in 1:24) {
    tr <- mat_trials(dat, s, cc); if (is.null(tr)) next
    parts[[length(parts) + 1]] <- data.frame(
      t = tr$t, resp = factor(LEV[tr$response], levels = LEV), subject = factor(s))
  }
  do.call(rbind, parts)
}

fit_hier <- function(d, model, pooled = FALSE, start = NULL) {
  st <- if (is.null(start)) HIER_START[HIER_PARS[[model]]] else start[HIER_PARS[[model]]]
  suppressWarnings(frm(hier_formula(model, pooled = pooled),
      family = categorical(levels = LEV), data = d,
      start = list(beta = unname(st))))
}

## Place the fixed effects at the median of the per-subject estimates. A
## non-linear model of this size does not find its own way from a generic
## guess; the first attempt at the correlated block stalled from one.
hier_start_from <- function(persubject, cond) {
  x <- persubject[persubject$ok & persubject$model == "habit" &
                  persubject$cond == cond, ]
  md <- function(n) stats::median(x[[paste0(n, "_frmtmb")]])
  sq <- function(v, lo, hi) stats::qlogis(min(max((v - lo) / (hi - lo), 1e-4), 1 - 1e-4))
  c(lmuA = md("muA"), lmuB = md("muB"), lsgA = log(md("sgA")), lsgB = log(md("sgB")),
    lqA = sq(md("qA"), 0.99, 0.999), lqB = sq(md("qB"), 0.5, 0.9999),
    lqI = sq(md("qI"), 0, 0.499))
}

hier_report <- function(fit) {
  if (inherits(fit, "error")) return(list(ok = FALSE, msg = conditionMessage(fit)))
  list(ok = TRUE, logLik = as.numeric(logLik(fit)), AIC = AIC(fit),
       df = attr(logLik(fit), "df"), convergence = fit$opt$convergence)
}

stage_hier <- function(conds = 1:3, fitfile = "HabitModelFits.mat") {
  m <- load_mat(fitfile); dat <- m$data
  ps <- read.csv(file.path(OUT, "table-persubject.csv"))
  out <- list()
  for (cc in conds) {
    d <- hier_data(dat, cc)
    st <- hier_start_from(ps, COND[cc])
    ## Complete pooling first, then the correlated block started from it.
    pooled <- tryCatch(fit_hier(d, "habit", pooled = TRUE, start = st),
                       error = function(e) e)
    st2 <- if (inherits(pooled, "error")) st else
      vapply(HIER_PARS[["habit"]], function(n) fixef(pooled)[[n]][[1]], 0)
    hb <- tryCatch(fit_hier(d, "habit", start = st2), error = function(e) e)
    nh <- tryCatch(fit_hier(d, "no-habit", start = st2[HIER_PARS[["no-habit"]]]),
                   error = function(e) e)
    out[[COND[cc]]] <- list(
      n = nrow(d), nsubj = nlevels(droplevels(d$subject)),
      start = st, pooled = pooled, habit = hb, nohabit = nh,
      report = list(pooled = hier_report(pooled), habit = hier_report(hb),
                    nohabit = hier_report(nh)))
  }
  out
}

# ============================================================================
# 7. BOOTSTRAP BANDS ON THE FITTED CURVES
#
#    "Still TODO: bootstrap CIs on curves" is the open item on HabitTR PR #3
#    (2017-09-01). frm_bootstrap() answers it directly: it simulates from the
#    fit, refits each draw, and applies FUN. Here FUN returns the three
#    category probabilities on a grid of preparation times, so the bootstrap
#    distribution is a band around the whole curve rather than around a
#    parameter.
# ============================================================================

BOOT_GRID <- seq(0.001, 1.2, by = 0.01)

## The pooled (complete-pooling) habit fit per group. Bootstrapping this is
## cheap because there are no random effects to re-integrate on every draw.
stage_boot <- function(conds = 1:3, nsim = 100, seed = 1,
                       fitfile = "HabitModelFits.mat") {
  m <- load_mat(fitfile); dat <- m$data
  ps <- read.csv(file.path(OUT, "table-persubject.csv"))
  out <- list()
  for (cc in conds) {
    d <- hier_data(dat, cc)
    fit <- tryCatch(fit_hier(d, "habit", pooled = TRUE,
                             start = hier_start_from(ps, COND[cc])),
                    error = function(e) e)
    if (inherits(fit, "error")) { out[[COND[cc]]] <- fit; next }
    curve_of <- function(f) {
      b <- fixef(f)
      g <- function(n) b[[n]][[1]]
      muB <- g("lmuB"); muA <- g("lmuA")
      p <- c(muA, exp(g("lsgA")), 0.99 + 0.009 / (1 + exp(-g("lqA"))),
             muB, exp(g("lsgB")), 0.5 + 0.4999 / (1 + exp(-g("lqB"))),
             0.499 / (1 + exp(-g("lqI"))), 1)
      pr <- hb_presponse(BOOT_GRID, p, "habit")
      ## Restore the second "other" key so the three curves sum to one.
      ## rbind gives 3 x k and as.vector reads column-major, so the result is
      ## interleaved by time: (cat1_t1, cat2_t1, cat3_t1, cat1_t2, ...), which
      ## is the layout boot_band() indexes. Do not transpose.
      as.vector(rbind(pr[1, ], pr[2, ], 2 * pr[3, ]))
    }
    bo <- tryCatch(frm_bootstrap(fit, FUN = curve_of, nsim = nsim, seed = seed),
                   error = function(e) e)
    out[[COND[cc]]] <- list(fit = fit, boot = bo, grid = BOOT_GRID)
  }
  out
}

boot_band <- function(bo, grid, probs = c(0.025, 0.975)) {
  k <- length(grid)
  keep <- bo$t[stats::complete.cases(bo$t), , drop = FALSE]
  cats <- c("correct", "habit", "other")
  do.call(rbind, lapply(seq_along(cats), function(j) {
    idx <- (seq_len(k) - 1) * 3 + j
    q <- apply(keep[, idx, drop = FALSE], 2, quantile, probs = probs, na.rm = TRUE)
    data.frame(category = cats[j], t = grid, est = bo$t0[idx],
               lo = q[1, ], hi = q[2, ])
  }))
}

# ---------------------------------------------------------------- driver ---

main <- function(stages) {
  if ("all" %in% stages)
    stages <- c("validate", "persubject", "pooled", "compare", "corrected",
                "hier", "boot")
  res <- list()
  if ("validate" %in% stages) {
    v <- cached("validate", stage_validate())
    write.csv(v$summary, file.path(OUT, "table-reference-validation.csv"), row.names = FALSE)
    print(v$summary); res$validate <- v
  }
  if ("persubject" %in% stages) {
    p <- cached("persubject", stage_persubject())
    write.csv(p, file.path(OUT, "table-persubject.csv"), row.names = FALSE)
    res$persubject <- p
  }
  if ("pooled" %in% stages) {
    p <- cached("persubject", stage_persubject())
    q <- p[p$ok, ]
    rows <- list()
    for (mo in MODELS) for (g in COND) {
      x <- q[q$model == mo & q$cond == g, ]
      if (!nrow(x)) next
      rows[[length(rows) + 1]] <- data.frame(model = mo, cond = g, n_subj = nrow(x),
        LL_matlab = sum(x$LL_matlab), LL_frmtmb = sum(x$LL_frmtmb),
        diff = sum(x$LL_frmtmb) - sum(x$LL_matlab),
        max_abs_per_subj = max(abs(x$dLL)))
    }
    ps <- do.call(rbind, rows)
    write.csv(ps, file.path(OUT, "table-pooled-summary.csv"), row.names = FALSE)
    print(ps); res$pooled <- ps
  }
  if ("compare" %in% stages) {
    cmp <- cached("compare", stage_compare())
    write.csv(cmp, file.path(OUT, "table-model-comparison.csv"), row.names = FALSE)
    res$compare <- cmp
  }
  if ("corrected" %in% stages) {
    cr <- cached("corrected", stage_corrected())
    write.csv(cr, file.path(OUT, "table-corrected-subsetting.csv"), row.names = FALSE)
    res$corrected <- cr
  }
  if ("hier" %in% stages) res$hier <- cached("hier", stage_hier())
  if ("boot" %in% stages) {
    b <- cached("boot", stage_boot())
    bands <- do.call(rbind, lapply(names(b), function(g) {
      z <- b[[g]]
      if (inherits(z, "error") || inherits(z$boot, "error")) return(NULL)
      cbind(cond = g, boot_band(z$boot, z$grid))
    }))
    if (!is.null(bands)) write.csv(bands, file.path(OUT, "table-bootstrap-bands.csv"),
                                   row.names = FALSE)
    res$boot <- b
  }
  invisible(res)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a)) main(a)
}
