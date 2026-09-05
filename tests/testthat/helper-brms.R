# Helpers for the brms cross-validation suite (test-brms-agreement.R).
#
# Two tiers. Tier 1 goes through brms's data-generating functions
# (make_standata, brmsterms, get_prior), which never touch Stan, so it
# runs whenever brms is installed and NOT_CRAN is set. Tier 2 fits Stan
# models and is opt-in only.

# Structural tier: brms must be installed, but nothing compiles.
skip_unless_brms <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("brms")
}

# Numeric tier: Stan compilation, minutes per model. Opt in with
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
skip_unless_brms_fit <- function() {
  skip_unless_brms()
  testthat::skip_if_not_installed("rstan")
  if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
    testthat::skip("set FRMTMB_BRMS_FIT_TESTS=true to run brms fit tests")
  }
}

# brms chatters about mixture ordering and dpar defaults; the design
# objects are what the tests read.
brms_standata <- function(...) {
  suppressMessages(brms::make_standata(...))
}

# Design matrices are compared by VALUE: brms names the intercept
# "Intercept" (no parentheses) so that it survives Stan's identifier
# rules, and drops matrix dimnames in places.
expect_design_equal <- function(x, y, tol = 1e-10) {
  x <- unname(as.matrix(x))
  y <- unname(as.matrix(y))
  testthat::expect_identical(dim(x), dim(y))
  testthat::expect_lt(max(abs(x - y)), tol)
}

# Column-space equality, for bases that span the same model but use a
# different parameterization (mgcv's diagonal.penalty reparameterization).
col_span_proj <- function(M) {
  M <- as.matrix(M)
  q <- qr(M)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  Q %*% t(Q)
}

expect_span_equal <- function(x, y, tol = 1e-8) {
  testthat::expect_lt(max(abs(col_span_proj(x) - col_span_proj(y))), tol)
}

# ---------------------------------------------------------------------
# Log-density identity harness (test-brms-likelihood.R).
#
# The claim is that frmtmb's objective is the same function of the
# parameters as the Stan program brms generates with flat priors, so the
# two are compared at a point rather than through their optimizers. See
# dev/brms-likelihood-tests.md.
# ---------------------------------------------------------------------

# Compiling one brms program costs about a minute and the program is a
# function of the formula and the family only, never of the data, so the
# cache is content addressed on the Stan code. The rstan version joins
# the key because a DSO built by one rstan is not loadable by another.
brms_stan_cache_dir <- function() {
  dir <- Sys.getenv("FRMTMB_STAN_CACHE", "")
  if (!nzchar(dir)) {
    dir <- testthat::test_path("..", "..", "dev", "stan-cache")
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  normalizePath(dir, winslash = "/", mustWork = TRUE)
}

# md5 over a file rather than a string keeps this on base R: the test
# suite must not gain a dependency to hash one character vector.
brms_stan_cache_key <- function(code) {
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  unname(tools::md5sum(f))
}

# Programs already in hand in THIS session, keyed by cache key.
#
# This is not an optimization, it is a correctness fix. A stanfit model
# object that was compiled in this session and is then re-read from its
# own RDS comes back with a DSO that will not initialize:
#   Failed to initialize module pointer:
#   Error in FUN(X[[i]], ...): NULL value passed for DllInfo
# So the first brms_stan_model() call for a program works and every
# later call in the same session fails. Reading from an RDS written by
# an EARLIER session is fine, which is why a fully warm run never saw
# it and a run from an empty cache errored on the two programs this
# file uses twice.
.brms_stan_models <- new.env(parent = emptyenv())

# stan_model(save_dso = TRUE) is the default and is what makes the saved
# object usable from a later session, so the on-disk cache is a plain
# RDS. The session cache above sits in front of it so that a program is
# read from disk at most once per process.
brms_stan_model <- function(code) {
  key <- brms_stan_cache_key(code)
  hit <- .brms_stan_models[[key]]
  if (!is.null(hit)) {
    return(hit)
  }
  path <- file.path(brms_stan_cache_dir(), paste0(key, ".rds"))
  if (file.exists(path)) {
    mod <- try(readRDS(path), silent = TRUE)
    if (!inherits(mod, "try-error")) {
      .brms_stan_models[[key]] <- mod
      return(mod)
    }
  }
  mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(mod, path)
  .brms_stan_models[[key]] <- mod
  mod
}

# Empty strings are brms's spelling for an improper flat prior, so this
# strips every default and leaves the log posterior equal to the log
# likelihood up to constants brms cannot be told to drop.
brms_flat_prior <- function(bform, data, family, ...) {
  p <- brms::get_prior(bform, data = data, family = family, ...)
  p$prior <- ""
  p
}

# The names declared in the Stan parameters block, in declaration order.
# The translator walks these so that an untranslated parameter is named
# in the failure instead of surfacing as an rstan error about lengths.
brms_stan_par_names <- function(code) {
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  lines <- sub("//.*$", "", lines)
  start <- grep("^parameters\\s*\\{", lines)
  if (!length(start)) {
    return(character())
  }
  depth <- 0L
  out <- character()
  for (i in seq(start[[1]], length(lines))) {
    ln <- lines[[i]]
    depth <- depth + lengths(regmatches(ln, gregexpr("\\{", ln))) -
      lengths(regmatches(ln, gregexpr("\\}", ln)))
    m <- regmatches(ln, regexpr("([A-Za-z_][A-Za-z0-9_]*)\\s*;\\s*$", ln))
    if (length(m)) {
      out <- c(out, trimws(sub(";.*$", "", m)))
    }
    if (i > start[[1]] && depth <= 0L) {
      break
    }
  }
  out
}

# brms strips the parentheses from the intercept so that the name is a
# legal Stan identifier; frmtmb keeps model.matrix's spelling.
brms_coef_to_frm <- function(x) {
  ifelse(x == "Intercept", "(Intercept)", x)
}

# Population-level effects for one linear predictor. The suffix is the
# dpar or nlpar brms appended to the parameter name; a bare "b" is mu.
brms_dpar_of <- function(sfx) if (is.na(sfx) || !nzchar(sfx)) "mu" else sfx

brms_X_of <- function(sdat, sfx) {
  nm <- if (identical(sfx, "mu")) "X" else paste0("X_", sfx)
  if (is.null(sdat[[nm]])) {
    stop("no brms design matrix ", nm, " for linear predictor ", sfx)
  }
  sdat[[nm]]
}

# The columns of X that brms's `b` covers. Kc counts them, so comparing
# it with ncol(X) decides whether an intercept column was dropped: an
# ordinal model has no intercept column at all and Kc == K, every other
# model has one and Kc == K - 1.
brms_Xc_cols <- function(sdat, sfx) {
  x <- brms_X_of(sdat, sfx)
  kc <- sdat[[if (identical(sfx, "mu")) "Kc" else paste0("Kc_", sfx)]]
  cn <- colnames(x)
  if (!is.null(kc) && length(cn) == kc + 1L) {
    cn <- cn[-1]
  }
  cn
}

# frmtmb's coefficient vector for one linear predictor, indexed by
# brms's column names so that the order is brms's, not frmtmb's.
brms_fe_of <- function(fit, sfx) {
  fe <- fixef(fit)
  if (is.null(fe[[sfx]])) {
    stop("frmtmb fit has no linear predictor named ", sfx)
  }
  fe[[sfx]]
}

# brms's own group-level table, from a brmsfit skeleton. brm(empty =
# TRUE) builds every design object and stops before Stan, so it costs
# about a third of a second and is the only public source of two things
# standata does not carry: the coefficient order inside a group and the
# level labels J_<id> indexes.
#
# The first attempt read the order off get_prior()'s "sd" rows instead.
# Those rows are sorted alphabetically, so (Days | Subject) came back as
# Days then Intercept while Stan wants Intercept then Days. Check B
# caught it: the value was wrong by 1149 nats and the gradient on the z
# block was 295 where it should have been zero.
brms_ranef_table <- function(bform, data, family, prior, ...) {
  suppressMessages(brms::brm(bform, data = data, family = family,
                             prior = prior, empty = TRUE, ...))$ranef
}

# One group-level block, keyed by the id that names it in standata
# (J_<id>, sd_<id>, z_<id>, L_<id>).
brms_group_info <- function(rtab, id) {
  rows <- rtab[rtab$id == id, , drop = FALSE]
  if (!nrow(rows)) {
    stop("brms has no group-level block with id ", id)
  }
  rows <- rows[order(rows$cn), , drop = FALSE]
  group <- rows$group[[1]]
  list(id = id, group = group, coefs = rows$coef, dpars = rows$dpar,
       nlpars = rows$nlpar,
       labels = as.character(attr(rtab, "levels")[[group]]))
}

# Every frmtmb ranef column that belongs to grouping factor `gvar`, as
# one levels-by-coefficients matrix. Separate blocks on the same factor
# ((1 | g) + (0 + x | g)) are distinct brms IDs but arrive here as
# distinct frmtmb blocks too, so they are pooled and then selected.
brms_ranef_block <- function(fit, gvar) {
  re <- ranef(fit)
  keep <- vapply(names(re), function(nm) {
    identical(brms_block_group(nm), gvar)
  }, logical(1))
  if (!any(keep)) {
    stop("frmtmb fit has no random effect on grouping factor ", gvar)
  }
  blocks <- re[keep]
  out <- do.call(cbind, blocks)
  colnames(out) <- unlist(lapply(blocks, colnames))
  out
}

# The grouping factor named by a frmtmb ranef or VarCorr block label.
# A plain block is "Days | Subject"; a block merged across dpars with
# |ID| is "1 | g + sigma: 1 | g [ID]", so the trailing marker comes off
# before the last bar-separated segment is read.
brms_block_group <- function(nm) {
  nm <- sub("[[][^]]*[]][[:space:]]*$", "", nm)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  trimws(parts[[length(parts)]])
}

# frmtmb's ordinal thresholds on the scale brms declares them.
#
# The storage convention is not the same for all four ordinal families,
# and the family object does not carry the flag: R/families.R passes
# ordered = TRUE for cumulative and sratio, which store
# (tau_1, log increments), and ordered = FALSE for cratio and acat,
# which store the thresholds themselves. Reading the family name here
# duplicates that one fact deliberately. If frmtmb ever changes a
# convention, checks A and B fail loudly, which is the point.
#
# Getting this wrong is not a crash. cratio came back 0.83 nats out
# with a gradient of 9, and acat 9.6 nats out, because the ordered
# transform was applied to thresholds that were already thresholds.
brms_ord_thresholds <- function(fit) {
  fam <- family(fit)[["family"]]
  ord_tau_from_raw(fit$estimates[["tau_raw"]],
                   ordered = fam %in% c("cumulative", "sratio"))
}

# frmtmb's column name for one brms group-level coefficient. An
# unmerged block names its columns exactly as model.matrix does; a
# block merged across dpars or responses prefixes them, as in
# "y.sigma:(Intercept)", because the same coefficient name then appears
# once per linear predictor.
brms_frm_coef <- function(cols, coef, dpar) {
  want <- brms_coef_to_frm(coef)
  if (want %in% cols) {
    return(want)
  }
  suffix <- paste0(if (nzchar(dpar)) dpar else "mu", ":", want)
  hit <- cols[cols == suffix | endsWith(cols, paste0(".", suffix))]
  if (length(hit) != 1L) {
    stop("frmtmb has no unique group-level column for ", coef,
         " under linear predictor ", if (nzchar(dpar)) dpar else "mu")
  }
  hit
}

# The plan's z rule. brms builds r = (diag(sd) L z)^T per group, so the
# inverse is z = solve(diag(sd) %*% L, t(r)), with r in brms's level
# order. The Jacobian of that map is returned with it because check C
# needs it and it is a property of the map, not of the check.
brms_group_pars <- function(fit, sdat, rtab, i) {
  info <- brms_group_info(rtab, i)
  n_lev <- sdat[[paste0("N_", i)]]
  m <- sdat[[paste0("M_", i)]]
  if (length(info$coefs) != m) {
    stop("group ", i, " declares ", m, " coefficients in standata but ",
         length(info$coefs), " in the brms ranef table")
  }
  if (length(info$labels) != n_lev) {
    stop("group ", i, " declares ", n_lev, " levels in standata but ",
         length(info$labels), " in the brms ranef table")
  }
  blocks <- unclass(VarCorr(fit))
  keep <- vapply(names(blocks), function(nm) {
    identical(brms_block_group(nm), info$group)
  }, logical(1))
  if (!any(keep)) {
    stop("frmtmb VarCorr has no block on grouping factor ", info$group)
  }
  sigma <- as.matrix(Reduce(brms_blockdiag, blocks[keep]))
  want <- mapply(brms_frm_coef, coef = info$coefs, dpar = info$dpars,
                 MoreArgs = list(cols = colnames(sigma)),
                 USE.NAMES = FALSE)
  sigma <- sigma[want, want, drop = FALSE]
  sd <- sqrt(diag(sigma))
  corr <- sigma / tcrossprod(sd)
  lmat <- t(chol(corr))
  lmat[upper.tri(lmat)] <- 0
  r <- brms_ranef_block(fit, info$group)
  rwant <- mapply(brms_frm_coef, coef = info$coefs, dpar = info$dpars,
                  MoreArgs = list(cols = colnames(r)), USE.NAMES = FALSE)
  r <- r[info$labels, rwant, drop = FALSE]
  z <- solve(diag(sd, nrow = length(sd)) %*% lmat, t(r))
  list(sd = sd, L = lmat, z = z, n_lev = n_lev,
       logJ = n_lev * (sum(log(sd)) + sum(log(diag(lmat)))))
}

# Reduce() over one block must not drop it to a bare matrix operation.
brms_blockdiag <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  out <- matrix(0, nrow(a) + nrow(b), ncol(a) + ncol(b))
  dimnames(out) <- list(c(rownames(a), rownames(b)),
                        c(colnames(a), colnames(b)))
  out[seq_len(nrow(a)), seq_len(ncol(a))] <- a
  out[nrow(a) + seq_len(nrow(b)), ncol(a) + seq_len(ncol(b))] <- b
  out
}

# The translator. Returns the constrained parameter list that
# rstan::unconstrain_pars() accepts, built from frmtmb's estimates alone.
#
# It is driven by the Stan parameters block rather than by frmtmb's
# parameter vector, so a shape brms declares and this has no rule for is
# an error naming the parameter, never a silently omitted term.
#
# Deviation from the plan's signature: the Stan code and brms's ranef
# table join `fit` and `sdat` as arguments. The declared parameter set is
# only in the code, and standata carries neither the level labels nor the
# coefficient order inside a group.
stan_pars_from_fit <- function(fit, sdat, code, rtab = NULL) {
  need <- brms_stan_par_names(code)
  links <- family(fit)$links
  out <- list()
  jac <- 0
  for (nm in need) {
    if (grepl("^b(_(.+))?$", nm)) {
      sfx <- brms_dpar_of(sub("^b_?", "", nm))
      fe <- brms_fe_of(fit, sfx)
      cn <- brms_Xc_cols(sdat, sfx)
      out[[nm]] <- array(unname(fe[brms_coef_to_frm(cn)]), length(cn))
    } else if (grepl("^Intercept(_(.+))?$", nm)) {
      sfx <- brms_dpar_of(sub("^Intercept_?", "", nm))
      x <- brms_X_of(sdat, sfx)
      fe <- brms_fe_of(fit, sfx)
      cn <- brms_Xc_cols(sdat, sfx)
      shift <- sum(colMeans(x)[cn] * unname(fe[brms_coef_to_frm(cn)]))
      if (identical(sfx, "mu") && !is.null(sdat[["nthres"]])) {
        # ordinal: Intercept is the threshold vector, X carries no
        # intercept column, and the centering enters with the opposite
        # sign because the thresholds are compared to mu, not added to
        # it. brms's generated quantity is
        #   b_Intercept = Intercept + dot_product(means_X, b)
        out[[nm]] <- as.numeric(
          brms_ord_thresholds(fit) - shift)
      } else {
        # brms centers X inside the Stan program, so its Intercept is
        # the intercept of the centered fit, not the one frmtmb reports:
        #   b_Intercept = Intercept - dot_product(means_X, b)
        out[[nm]] <- as.numeric(fe[["(Intercept)"]] + shift)
      }
    } else if (grepl("^bsp(_(.+))?$", nm)) {
      sfx <- brms_dpar_of(sub("^bsp_?", "", nm))
      fe <- brms_fe_of(fit, sfx)
      known <- brms_coef_to_frm(colnames(brms_X_of(sdat, sfx)))
      sp <- fe[setdiff(names(fe), known)]
      out[[nm]] <- array(unname(sp), length(sp))
    } else if (identical(nm, "ordered_Intercept")) {
      # brms identifies a mixture by ordering the component intercepts:
      # ordered_Intercept[k] IS Intercept_mu<k>, declared ordered so
      # that label switching cannot happen. The entries are the centered
      # intercepts, exactly as for any other dpar.
      ks <- grep("^mu[0-9]+$", names(links), value = TRUE)
      ks <- ks[order(as.integer(sub("^mu", "", ks)))]
      out[[nm]] <- vapply(ks, function(k) {
        fe <- brms_fe_of(fit, k)
        x <- if (is.null(sdat[[paste0("X_", k)]])) NULL else {
          brms_X_of(sdat, k)
        }
        shift <- if (is.null(x)) 0 else {
          cn <- brms_Xc_cols(sdat, k)
          sum(colMeans(x)[cn] * unname(fe[brms_coef_to_frm(cn)]))
        }
        as.numeric(fe[["(Intercept)"]] + shift)
      }, numeric(1), USE.NAMES = FALSE)
      if (is.unsorted(out[[nm]])) {
        stop("frmtmb's mixture components are not in brms's increasing ",
             "intercept order, so the component labels differ")
      }
    } else if (grepl("^bcs(_(.+))?$", nm)) {
      # brms declares bcs as matrix[Kcs, nthres]: one row per
      # category-specific covariate, one column per threshold. frmtmb
      # stores it as one bcs<j> vector of length nthres per covariate.
      kcs <- sdat[["Kcs"]]
      nth <- sdat[["nthres"]]
      keys <- grep("^bcs[0-9]+$", names(fit$estimates), value = TRUE)
      keys <- keys[order(as.integer(sub("^bcs", "", keys)))]
      v <- unlist(lapply(keys, function(k) fit$estimates[[k]]))
      if (length(v) != kcs * nth) {
        stop("brms wants ", kcs * nth, " category-specific coefficients ",
             "but frmtmb holds ", length(v))
      }
      out[[nm]] <- matrix(v, nrow = kcs, ncol = nth, byrow = TRUE)
    } else if (grepl("^simo_", nm)) {
      j <- as.integer(sub(".*_(\\d+)$", "\\1", nm))
      zeta <- fit$estimates[[paste0("zeta", j)]]
      s <- exp(c(0, zeta))
      out[[nm]] <- s / sum(s)
    } else if (grepl("^(sd|z|L)_\\d+$", nm)) {
      i <- as.integer(sub("^[a-zA-Z]+_", "", nm))
      if (is.null(rtab)) {
        stop("Stan parameter ", nm, " needs the brms ranef table")
      }
      gp <- brms_group_pars(fit, sdat, rtab, i)
      if (startsWith(nm, "sd")) {
        out[[nm]] <- array(gp$sd, length(gp$sd))
        jac <- jac + gp$logJ
      } else if (startsWith(nm, "z")) {
        out[[nm]] <- gp$z
      } else {
        out[[nm]] <- gp$L
      }
    } else if (nm %in% names(links)) {
      # a dpar with no linear predictor: brms declares it on the natural
      # scale, frmtmb estimates it on the link scale
      fe <- brms_fe_of(fit, nm)
      out[[nm]] <- as.numeric(links[[nm]]$linkinv(fe[["(Intercept)"]]))
    } else {
      stop("no translation rule for Stan parameter ", nm)
    }
  }
  attr(out, "logJ") <- jac
  out
}

# Which entries of the unconstrained vector are the z blocks. The map
# from z to unconstrained space is the identity, so perturbing z and
# unconstraining again marks exactly its slots, whatever the declaration
# order and whatever else the program declares.
brms_z_index <- function(sf, pars) {
  zn <- grep("^z_\\d+$", names(pars), value = TRUE)
  if (!length(zn)) {
    return(integer())
  }
  bumped <- pars
  for (nm in zn) bumped[[nm]] <- pars[[nm]] + 1
  which(rstan::unconstrain_pars(sf, bumped) !=
          rstan::unconstrain_pars(sf, pars))
}

# The two checks of the plan, run at one point.
#
#   A  log_prob at frmtmb's estimates equals frmtmb's log density plus a
#      known constant.
#   B  the gradient there vanishes, which is what catches a wrong map:
#      a mistranslated parameter lands off brms's optimum.
#   C  when joint = TRUE the comparison is against minus RTMB's inner
#      objective plus the map's log-Jacobian, and B is asserted on the z
#      block only, because the outer gradient is not zero at that point.
brms_lp_check <- function(bform, family, data, fit, joint = FALSE,
                          const = 0, tol_grad = 1e-3, ...) {
  prior <- brms_flat_prior(bform, data = data, family = family, ...)
  code <- brms::make_stancode(bform, data = data, family = family,
                              prior = prior, ...)
  sdat <- brms_standata(bform, data = data, family = family,
                        prior = prior, ...)
  rtab <- brms_ranef_table(bform, data, family, prior, ...)
  mod <- brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))

  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  upars <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, upars, adjust_transform = FALSE,
                        gradient = FALSE)

  if (joint) {
    ours <- -fit$obj$env$f(fit$obj$env$last.par.best) + attr(pars, "logJ")
  } else {
    ours <- as.numeric(logLik(fit))
  }
  tol <- 1e-6 * max(1, abs(ours))
  testthat::expect_lt(abs(lp - ours - const), tol)

  grad <- rstan::grad_log_prob(sf, upars, adjust_transform = FALSE)
  if (joint) {
    idx <- brms_z_index(sf, pars)
    testthat::expect_gt(length(idx), 0)
    grad <- grad[idx]
  }
  testthat::expect_lt(max(abs(grad)), tol_grad)

  invisible(list(lp = lp, ours = ours, measured_const = lp - ours,
                 max_grad = max(abs(grad)), pars = pars, sdat = sdat,
                 sf = sf))
}

# Round trip through Stan's own constraint machinery. This is the unit
# test on the map alone: a simplex that does not sum to one, or an L
# that is not a valid Cholesky correlation factor, comes back changed
# even though log_prob would still return a number.
expect_par_roundtrip <- function(sf, pars, tol = 1e-10) {
  back <- rstan::constrain_pars(sf, rstan::unconstrain_pars(sf, pars))
  for (nm in names(pars)) {
    testthat::expect_true(nm %in% names(back),
                          label = paste0("constrain_pars returned ", nm))
    testthat::expect_lt(max(abs(as.numeric(back[[nm]]) -
                                  as.numeric(pars[[nm]]))), tol)
  }
  invisible(TRUE)
}
