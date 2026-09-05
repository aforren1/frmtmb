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

# Numeric tier: Stan compilation, minutes per model. This calls
# skip_unless_brms(), which calls skip_on_cran(), so outside R CMD check
# BOTH variables are needed to opt in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
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
# brms's column names so that the order is brms's, not frmtmb's. The
# suffix arrives in brms's spelling and goes through brms_lp_parts(),
# which is the identity for a univariate fit and turns "sigma_y1" into
# frmtmb's "y1_sigma" for a multivariate one.
brms_fe_of <- function(fit, sfx) {
  fe <- fixef(fit)
  key <- brms_lp_parts(fit, sfx)$frm
  if (is.null(fe[[key]])) {
    stop("frmtmb fit has no linear predictor named ", key)
  }
  fe[[key]]
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
  resps <- as.character(rows$resp)
  resps[is.na(resps)] <- ""
  list(id = id, group = group, coefs = rows$coef, dpars = rows$dpar,
       nlpars = rows$nlpar, resps = resps,
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
brms_frm_coef <- function(cols, coef, dpar, resp = "") {
  want <- brms_coef_to_frm(coef)
  suffix <- paste0(if (nzchar(dpar)) dpar else "mu", ":", want)
  if (nzchar(resp)) {
    # a multivariate block repeats the coefficient name once per
    # RESPONSE as well as once per dpar, so the bare suffix match below
    # is ambiguous and the response has to qualify it
    hit <- cols[cols == paste0(resp, ".", suffix)]
    if (length(hit) == 1L) {
      return(hit)
    }
  }
  if (want %in% cols) {
    return(want)
  }
  hit <- cols[cols == suffix | endsWith(cols, paste0(".", suffix))]
  if (length(hit) != 1L) {
    stop("frmtmb has no unique group-level column for ", coef,
         " under linear predictor ",
         if (nzchar(resp)) paste0(resp, "_", if (nzchar(dpar)) dpar else "mu")
         else if (nzchar(dpar)) dpar else "mu")
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
                 resp = info$resps,
                 MoreArgs = list(cols = colnames(sigma)),
                 USE.NAMES = FALSE)
  sigma <- sigma[want, want, drop = FALSE]
  sd <- sqrt(diag(sigma))
  corr <- sigma / tcrossprod(sd)
  lmat <- t(chol(corr))
  lmat[upper.tri(lmat)] <- 0
  r <- brms_ranef_block(fit, info$group)
  rwant <- mapply(brms_frm_coef, coef = info$coefs, dpar = info$dpars,
                  resp = info$resps,
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
    } else if (grepl("^(Ymi|Yl)_", nm)) {
      # mi(): brms declares Ymi_<resp> for the missing entries of a plain
      # mi() response and Yl_<resp> for EVERY entry of a measurement-error
      # one, and merges each with the observed data inside the program.
      # frmtmb keeps both kinds in one `miss` vector, per response in
      # formula order, and fit$frame$mi_map is the row-to-slot map. The
      # map to brms is the identity, so this block adds nothing to the
      # log-Jacobian.
      resp <- sub("^(Ymi|Yl)_", "", nm)
      mm <- fit$frame[["mi_map"]][[resp]]
      if (is.null(mm)) {
        stop("brms declares ", nm, " but the frmtmb fit has no mi() ",
             "block for response ", resp)
      }
      want <- if (startsWith(nm, "Ymi")) {
        as.integer(sdat[[paste0("Jmi_", resp)]])
      } else {
        seq_len(sdat[[paste0("N_", resp)]])
      }
      if (!identical(as.integer(mm$rows), as.integer(want))) {
        stop("brms holds ", length(want), " latent values for ", resp,
             " and frmtmb ", length(mm$rows), ", or in another order")
      }
      v <- unname(fit$estimates[["miss"]][mm$idx])
      out[[nm]] <- array(v, length(v))
    } else if (grepl("^bs(_(.+))?$", nm)) {
      # the unpenalized part of a smooth. brms puts it in Xs, which it
      # does NOT center, and frmtmb names the same columns
      # "<term>.fx<j>" among the population-level effects. The two
      # spellings of the null space span the same columns on a
      # different scale, so the coefficients go through the change of
      # basis; `bs` carries no prior, so any invertible map will do.
      sfx <- brms_dpar_of(sub("^bs_?", "", nm))
      fe <- brms_fe_of(fit, sfx)
      v <- fe[grepl("[.]fx[0-9]+$", names(fe))]
      ks <- sdat[[if (identical(sfx, "mu")) "Ks" else paste0("Ks_", sfx)]]
      if (length(v) != ks) {
        stop("brms wants ", ks, " unpenalized spline coefficients but ",
             "frmtmb holds ", length(v))
      }
      xb <- as.matrix(sdat[[if (identical(sfx, "mu")) "Xs" else {
        paste0("Xs_", sfx)
      }]])
      xf <- brms_lp_of(fit, sfx)$X
      xf <- as.matrix(xf)[, names(v), drop = FALSE]
      v <- brms_basis_map(xb, xf, orthogonal = FALSE) %*% unname(v)
      out[[nm]] <- array(as.numeric(v), length(v))
    } else if (grepl("^sds_", nm)) {
      # one smoothing SD per basis of the term. frmtmb makes each basis
      # its own homogeneous-diagonal block, so the SD is the square root
      # of any diagonal entry of that block's covariance.
      ip <- brms_idx_parts(nm, "sds")
      bks <- brms_smooth_term(fit, ip$dpar, ip$idx[[1]])
      sdv <- vapply(bks, function(b) {
        sqrt(brms_block_cov(fit, b)[1, 1])
      }, numeric(1))
      out[[nm]] <- array(unname(sdv), length(sdv))
      jac <- jac + sum(vapply(bks, `[[`, 0, "dim") * log(sdv))
    } else if (grepl("^zs_", nm)) {
      # brms builds s = sds[j] * zs, so zs is the wiggly part divided by
      # its own SD, as in the z rule for a group with one coefficient,
      # but first the coefficients are put in brms's basis: brms calls
      # mgcv::smoothCon() with diagonal.penalty = TRUE and frmtmb does
      # not, so for s() the two bases are the same columns in the
      # opposite order. The map has to be orthogonal or the two i.i.d.
      # priors are not the same prior, and brms_basis_map() checks that.
      ip <- brms_idx_parts(nm, "zs")
      bk <- brms_smooth_term(fit, ip$dpar, ip$idx[[1]])[[ip$idx[[2]]]]
      sdv <- sqrt(brms_block_cov(fit, bk)[1, 1])
      amat <- brms_basis_map(
        as.matrix(sdat[[paste0("Zs_", ip$idx[[1]], "_", ip$idx[[2]])]]),
        brms_smooth_z(fit, ip$dpar, bk))
      v <- as.numeric(amat %*% brms_block_b(fit, bk)) / sdv
      out[[nm]] <- array(v, length(v))
    } else if (grepl("^(sdgp|lscale)_", nm)) {
      pre <- sub("_.*$", "", nm)
      ip <- brms_idx_parts(nm, pre)
      bk <- brms_gp_block(fit, ip$dpar, ip$idx[[1]])
      th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
      kgp <- sdat[[paste0("Kgp_", ip$idx[[1]])]]
      if (!identical(as.integer(kgp), 1L)) {
        stop("brms splits GP term ", ip$idx[[1]], " into ", kgp,
             " sub-GPs; the by= spelling has no rule yet")
      }
      if (identical(pre, "sdgp")) {
        out[[nm]] <- array(exp(th[[1]]), 1)
      } else {
        dgp <- as.integer(sdat[[paste0("Dgp_", ip$idx[[1]])]])
        rho <- if (isTRUE(bk[["gp_iso"]])) rep(exp(th[[2]]), dgp) else {
          exp(th[-1])
        }
        if (identical(bk[["covstruct"]], "gp")) {
          # An exact gp() measures distances on the DATA scale in frmtmb
          # and on brms's unit-maximum-distance scale in Stan, so the two
          # length-scales differ by dmax. The HSGP form rescales its
          # inputs the same way brms does and needs no correction.
          rho <- rho / sdat[[paste0("dmax_", ip$idx[[1]])]]
        }
        out[[nm]] <- matrix(rho, nrow = 1, ncol = dgp)
      }
    } else if (grepl("^zgp_", nm)) {
      # brms builds the GP as L z with L the Cholesky factor of the
      # kernel, so z inverts it. For the HSGP the "kernel" is the
      # diagonal of spectral densities and this reduces to f / sd.
      #
      # For the EXACT form the two kernels are not the same matrix: see
      # the nugget divergence in dev/brms-likelihood-tests.md. Nothing
      # here can repair that, and the exact spelling is asserted
      # structurally instead of being run through the identity.
      ip <- brms_idx_parts(nm, "zgp")
      bk <- brms_gp_block(fit, ip$dpar, ip$idx[[1]])
      lmat <- t(chol(brms_block_cov(fit, bk)))
      v <- solve(lmat, brms_block_b(fit, bk))
      out[[nm]] <- array(as.numeric(v), length(v))
      jac <- jac + sum(log(diag(lmat)))
    } else if (nm %in% c("ar", "ma", "cosy", "Lcortime")) {
      # Residual autocorrelation. frmtmb keeps these on unconstrained
      # scales in `thetaac` and has one internal function that puts them
      # back on brms's scale under brms's own names, so this reads that
      # rather than repeating the Levinson recursion.
      ac <- brms_autocor_of(fit)
      th <- fit$estimates[["thetaac"]][ac[["theta_idx"]]]
      if (identical(nm, "Lcortime")) {
        # the unstructured factor IS frmtmb's parameterization: a lower
        # triangular matrix with unit row norms is a Cholesky
        # correlation factor already
        out[[nm]] <- as.matrix(us_chol_L(th, ac[["d"]]))
      } else {
        nat <- autocor_natural(th, ac)
        v <- nat[startsWith(names(nat), nm)]
        if (!length(v)) {
          stop("frmtmb's autocorrelation block has no ", nm, " parameter")
        }
        out[[nm]] <- if (identical(nm, "cosy")) {
          as.numeric(v)
        } else {
          array(unname(v), length(v))
        }
      }
    } else if (nm %in% c("sdcar", "car", "rcar", "zcar")) {
      # CAR fields. brms declares the proper form (escar) on the field
      # itself and the intrinsic one (icar) on a standardized field, so
      # only the second carries a Jacobian.
      bk <- brms_car_block(fit)
      th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
      ty <- bk[["aux_car"]][["type"]]
      if (!ty %in% c("escar", "icar")) {
        stop("brms parameterizes a ", ty, " field on a different set of ",
             "latent variables than frmtmb; see ",
             "dev/brms-likelihood-tests.md")
      }
      if (identical(nm, "sdcar")) {
        out[[nm]] <- exp(th[[1]])
      } else if (identical(nm, "car")) {
        out[[nm]] <- car_rho(th[[2]])
      } else if (identical(nm, "rcar")) {
        v <- brms_block_b(fit, bk)
        out[[nm]] <- array(v, length(v))
      } else {
        v <- brms_block_b(fit, bk) / exp(th[[1]])
        out[[nm]] <- array(v, length(v))
        jac <- jac + length(v) * th[[1]]
      }
    } else if (identical(nm, "Lrescor")) {
      out[[nm]] <- brms_lrescor(fit, sdat)
    } else if (!is.null(brms_link_of(fit, brms_lp_parts(fit, nm)))) {
      # a dpar with no linear predictor: brms declares it on the natural
      # scale, frmtmb estimates it on the link scale. The link comes
      # through brms_link_of() rather than off `links` directly because
      # family(fit)$links is empty for a multivariate fit, where each
      # response carries its own family.
      fe <- brms_fe_of(fit, nm)
      link <- brms_link_of(fit, brms_lp_parts(fit, nm))
      out[[nm]] <- as.numeric(link$linkinv(fe[["(Intercept)"]])) *
        brms_sigma_scale(fit, code, nm)
    } else {
      stop("no translation rule for Stan parameter ", nm)
    }
  }
  attr(out, "logJ") <- jac
  out
}

# The parameters that stand where frmtmb keeps an INNER parameter: the
# ones check C asserts a zero gradient on. Every one of them is declared
# unbounded by brms, so the map to unconstrained space is the identity
# and the bump below marks exactly its own slots.
#
#   z_<i>       group-level effects
#   zs_<i>_<j>  the penalized part of a smooth
#   zgp_<i>     the latent variables of a Gaussian process
#   Ymi_<r>     the missing entries of an mi() response
#   Yl_<r>      the latent values of a measurement-error response
#   rcar, zcar  a CAR field, on its own scale for the proper form and
#               standardized for the intrinsic one
brms_inner_pat <-
  "^(z_\\d+|zs_\\d+_\\d+|zgp_\\d+|Ymi_.+|Yl_.+|rcar|zcar)$"

# Which entries of the unconstrained vector those blocks occupy, found
# by perturbing them and diffing, so it needs no knowledge of the
# declaration order or of any other parameter's unconstrained size.
brms_inner_index <- function(sf, pars, pattern = brms_inner_pat) {
  zn <- grep(pattern, names(pars), value = TRUE)
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
    idx <- brms_inner_index(sf, pars)
    testthat::expect_gt(length(idx), 0)
    grad <- grad[idx]
  }
  testthat::expect_lt(max(abs(grad)), tol_grad)

  # A row that drifts should report a number, not only a fail, and the
  # plan's results table is built from these two. Printing is opt-in so
  # that an ordinary run stays quiet:
  #   options(frmtmb.brms_lp_report = TRUE)
  if (isTRUE(getOption("frmtmb.brms_lp_report", FALSE))) {
    cat(sprintf("LPCHECK const %.10g grad %.3g ours %.10g\n",
                lp - ours, max(abs(grad)), ours))
  }
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

# ---------------------------------------------------------------------
# Multivariate responses (rows 8 and 9).
#
# brms suffixes every parameter of a response with the response name and
# puts the dpar FIRST, so y1's sigma is `sigma_y1` and its distributional
# coefficients are `b_sigma_y1`. frmtmb keys its linear predictors the
# other way round, `y1_sigma`, so the two spellings are the same pair of
# names in the opposite order and one of them has to be rewritten.
# ---------------------------------------------------------------------

# The responses of a multivariate fit, in frmtmb's order. A univariate
# fit has none, which is what makes every rule below collapse to the
# spelling it had before this section existed.
brms_mv_resps <- function(fit) {
  if (!inherits(fit$bform, "frmtmb_mvformula")) {
    return(character())
  }
  names(fit$spec$responses)
}

# One brms parameter suffix, split into the response and the dpar it
# names, plus the key frmtmb files that linear predictor under. The
# longest matching response wins, because a response named `y` and one
# named `y_hi` would otherwise both match the suffix `sigma_y_hi`.
brms_lp_parts <- function(fit, sfx) {
  dp <- brms_dpar_of(sfx)
  resps <- brms_mv_resps(fit)
  hit <- resps[sfx == resps | endsWith(sfx, paste0("_", resps))]
  if (!length(hit)) {
    return(list(resp = "", dpar = dp, frm = dp))
  }
  r <- hit[[which.max(nchar(hit))]]
  dpar <- if (identical(sfx, r)) "mu" else sub(paste0("_", r, "$"), "", sfx)
  list(resp = r, dpar = dpar, frm = paste0(r, "_", dpar))
}

# The link frmtmb put on one dpar, or NULL when that dpar is not one of
# the family's. A univariate fit carries its links on the family object;
# for a multivariate fit family(fit)$links is EMPTY and there is one
# family per response instead, so the response has to choose it.
brms_link_of <- function(fit, parts) {
  fam <- family(fit)
  links <- if (nzchar(parts$resp)) fam[[parts$resp]]$links else fam$links
  links[[parts$dpar]]
}

# brms's response order, which is the row order of its `Y` array and so
# the order Lrescor's rows and columns are in. It is the order the
# responses appear in standata, not an alphabetical one.
brms_resp_order <- function(sdat) {
  sub("^Y_", "", grep("^Y_", names(sdat), value = TRUE))
}

# The residual-correlation Cholesky factor brms declares as Lrescor.
# frmtmb keeps the same matrix behind rescor_matrix(), which reads its
# `thetar` segment through the unit-row-norm parameterization, so the
# only work here is the factorization and putting the responses in
# brms's order.
brms_lrescor <- function(fit, sdat) {
  cmat <- rescor_matrix(fit)
  if (is.null(cmat)) {
    stop("brms declares Lrescor but the frmtmb fit has no rescor")
  }
  ord <- brms_resp_order(sdat)
  if (!all(ord %in% rownames(cmat))) {
    stop("brms orders the responses ", paste(ord, collapse = ", "),
         " but frmtmb has ", paste(rownames(cmat), collapse = ", "))
  }
  lmat <- t(chol(cmat[ord, ord, drop = FALSE]))
  lmat[upper.tri(lmat)] <- 0
  unname(lmat)
}

# ---------------------------------------------------------------------
# Latent blocks that are not group-level effects (rows 4, 10 and 11).
#
# frmtmb keeps a smooth's wiggly part, a GP's function values and an
# mi() response's imputations in the same inner parameter vector it uses
# for random effects, so all three are check C shapes: brms declares a
# standardized version of each and multiplies it up inside the program.
# ---------------------------------------------------------------------

# The dpar or nlpar a frmtmb block belongs to, with mu spelled the same
# way on both sides.
brms_block_dpar <- function(bk) {
  dp <- bk[["dpar"]]
  if (is.null(dp) || !nzchar(dp)) "mu" else dp
}

# brms indexes these parameters by position: sds_1 is the first smooth
# term of mu, zs_1_2 its second basis, sdgp_2 the second GP term. The
# name carries an optional dpar between the prefix and the numbers.
brms_idx_parts <- function(nm, prefix) {
  rest <- sub(paste0("^", prefix, "_"), "", nm)
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  nums <- suppressWarnings(as.integer(parts))
  k <- which(!is.na(nums))
  if (!length(k)) {
    stop("Stan parameter ", nm, " carries no index")
  }
  list(dpar = if (min(k) > 1L) {
         paste(parts[seq_len(min(k) - 1L)], collapse = "_")
       } else "mu",
       idx = nums[k])
}

# frmtmb's blocks of one covariance structure under one linear
# predictor, in frmtmb's order, which is the formula order brms also
# numbers by.
brms_blocks_of <- function(fit, cs, dpar = "mu") {
  bks <- fit$frame[["re_blocks"]]
  keep <- vapply(bks, function(b) {
    b[["covstruct"]] %in% cs && identical(brms_block_dpar(b), dpar)
  }, logical(1))
  bks[keep]
}

# One smooth TERM's blocks. brms gives a term one sds_<i> of length
# nb_<i> and one zs_<i>_<j> per basis, while frmtmb makes each basis its
# own block and labels them all with the term, so the term's blocks are
# the run that shares a label.
brms_smooth_term <- function(fit, dpar, i) {
  bks <- brms_blocks_of(fit, "smooth", dpar)
  if (!length(bks)) {
    stop("frmtmb fit has no smooth blocks under linear predictor ", dpar)
  }
  labs <- vapply(bks, `[[`, "", "term_label")
  terms <- split(bks, factor(labs, levels = unique(labs)))
  if (i > length(terms)) {
    stop("brms declares smooth term ", i, " but frmtmb has ",
         length(terms), " under linear predictor ", dpar)
  }
  terms[[i]]
}

# frmtmb's fitted values for one block, straight out of the inner
# parameter vector. There is no level reindexing here: the basis order
# of a smooth, and the unique-position order of a GP, are brms's own and
# are checked element for element in test-brms-agreement.R.
brms_block_b <- function(fit, bk) {
  as.numeric(fit$estimates[["b"]][bk[["b_idx"]]])
}

# The block's covariance at frmtmb's estimates, from the same registry
# entry the objective uses. VarCorr() reports only the marginal SD for a
# smooth or a GP block, so it cannot serve here.
brms_block_cov <- function(fit, bk) {
  th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
  as.matrix(covstruct_registry[[bk[["covstruct"]]]]$vcov(th, bk))
}

# One GP term's frmtmb block. brms numbers gp() terms 1..n per linear
# predictor and frmtmb keeps them in the same formula order; the exact
# and Hilbert-space forms are one sequence, since brms declares the same
# three parameters for both.
brms_gp_block <- function(fit, dpar, i) {
  bks <- brms_blocks_of(fit, c("gp", "hsgp"), dpar)
  if (i > length(bks)) {
    stop("brms declares GP term ", i, " but frmtmb has ", length(bks),
         " under linear predictor ", dpar)
  }
  bks[[i]]
}

# ---------------------------------------------------------------------
# Residual autocorrelation and CAR fields (rows 18 and 19).
# ---------------------------------------------------------------------

# The fit's single residual-autocorrelation block. brms declares `ar`,
# `cosy` and `Lcortime` once per model, with no response suffix in the
# univariate case, so a second block would make the name ambiguous.
brms_autocor_of <- function(fit) {
  acs <- fit$frame[["autocor"]]
  if (length(acs) != 1L) {
    stop("brms names its autocorrelation parameters once per model, so ",
         "this rule needs exactly one frmtmb block and the fit has ",
         length(acs))
  }
  acs[[1]]
}

# The fit's single CAR block, likewise: brms declares `rcar` and
# `sdcar` without an index.
brms_car_block <- function(fit) {
  bks <- fit$frame[["re_blocks"]]
  keep <- vapply(bks, function(b) identical(b[["covstruct"]], "car"),
                 logical(1))
  if (sum(keep) != 1L) {
    stop("brms declares one unindexed rcar, so this rule needs exactly ",
         "one frmtmb car() block and the fit has ", sum(keep))
  }
  bks[keep][[1]]
}

# The graph Laplacian brms's ICAR density is written on, built from the
# edge list in ITS standata rather than from frmtmb's adjacency matrix,
# so the constants below are a reference-free check on frmtmb's K.
brms_car_laplacian <- function(sdat) {
  n <- as.integer(sdat[["Nloc"]])
  w <- matrix(0, n, n)
  e1 <- as.integer(sdat[["edges1"]])
  e2 <- as.integer(sdat[["edges2"]])
  w[cbind(e1, e2)] <- 1
  w[cbind(e2, e1)] <- 1
  diag(rowSums(w)) - w
}

# What brms's CAR densities leave out.
#
# brms writes both of them unnormalized in the field, and frmtmb keeps a
# proper density, so the difference is a closed form in the DATA alone
# and is admitted the same way the mo() Dirichlet is. Every quantity
# here comes from brms's own standata.
#
#   escar   brms drops the -Nloc/2 log(2 pi) and the 0.5 log det D of a
#           proper CAR, leaving 0.5(Nloc log tau + sum log(1 - car e_i)
#           - tau q).
#   icar    brms drops the ICAR normalizer entirely and adds a
#           normalized soft sum-to-zero term, normal(sum(zcar) | 0,
#           0.001 Nloc), whose precision is exactly the rank-one term
#           frmtmb folds into K.
brms_car_const <- function(sdat, type) {
  n <- as.integer(sdat[["Nloc"]])
  if (identical(type, "escar")) {
    return(0.5 * (n * log(2 * pi) - sum(log(as.numeric(sdat[["Nneigh"]])))))
  }
  s <- 0.001 * n
  kmat <- brms_car_laplacian(sdat) + matrix(1 / s^2, n, n)
  0.5 * (n - 1) * log(2 * pi) - log(s) -
    0.5 * as.numeric(determinant(kmat, logarithm = TRUE)$modulus)
}

# frmtmb's linear-predictor frame for one brms suffix. The frame keys
# these "<resp>.<dpar>" and brms names the dpar alone in the univariate
# case, so the response is matched by suffix.
brms_lp_of <- function(fit, sfx) {
  parts <- brms_lp_parts(fit, sfx)
  lps <- fit$frame[["linpreds"]]
  hit <- names(lps)[endsWith(names(lps), paste0(".", parts$dpar)) &
                      (!nzchar(parts$resp) |
                         startsWith(names(lps), paste0(parts$resp, ".")))]
  if (length(hit) != 1L) {
    stop("frmtmb has no unique linear-predictor frame for ", sfx)
  }
  lps[[hit]]
}

# frmtmb's random-effect design columns for one block.
brms_smooth_z <- function(fit, dpar, bk) {
  as.matrix(brms_lp_of(fit, dpar)$Z)[, bk[["c_idx"]], drop = FALSE]
}

# The change of basis from frmtmb's design columns to brms's, read off
# the two design matrices: `zf = zb %*% A`, so a coefficient vector in
# frmtmb's basis becomes `A %*% v` in brms's.
#
# For a PENALIZED block the map must be orthogonal, because both
# packages put the same i.i.d. prior on the coefficients and only an
# orthogonal change of basis leaves that prior alone. For the
# unpenalized columns there is no prior and any invertible map is fine.
# A non-orthogonal penalized map would be a real divergence and is an
# error naming it rather than a silently different model.
brms_basis_map <- function(zb, zf, orthogonal = TRUE, tol = 1e-8) {
  zb <- as.matrix(zb)
  zf <- as.matrix(zf)
  if (!identical(dim(zb), dim(zf))) {
    stop("brms has a ", ncol(zb), "-column basis where frmtmb has ",
         ncol(zf))
  }
  amat <- solve(crossprod(zb), crossprod(zb, zf))
  if (max(abs(zf - zb %*% amat)) > tol * max(1, max(abs(zf)))) {
    stop("the two bases do not span the same space")
  }
  if (orthogonal &&
        max(abs(crossprod(amat) - diag(ncol(amat)))) > 1e-8) {
    stop("the change of basis between frmtmb's design and brms's is ",
         "not orthogonal, so the two i.i.d. priors are not the same ",
         "prior and the models differ")
  }
  amat
}

# brms's sigma is not always frmtmb's sigma.
#
# For ar(cov = TRUE) the generated cholesky_cor_ar1() returns the
# Cholesky factor of T / (1 - ar^2), not of the correlation matrix T,
# and normal_time_hom_lpdf multiplies it by sigma. So brms's sigma is
# the INNOVATION standard deviation of the AR process and frmtmb's is
# the marginal residual SD: the same model, two scales, related by
# sqrt(1 - ar^2). cosy() and unstr() carry proper correlation matrices
# and need no correction.
#
# The rule keys on the function brms actually emitted rather than on the
# family or the formula, so a program that stops dividing stops being
# corrected.
brms_sigma_scale <- function(fit, code, nm) {
  if (!identical(brms_lp_parts(fit, nm)$dpar, "sigma")) {
    return(1)
  }
  if (!grepl("cholesky_cor_ar1", code, fixed = TRUE)) {
    return(1)
  }
  ac <- brms_autocor_of(fit)
  nat <- autocor_natural(fit$estimates[["thetaac"]][ac[["theta_idx"]]], ac)
  sqrt(1 - nat[["ar[1]"]]^2)
}
