# Can L1 see every latent block of a pair at the points it reads?
#
#   Rscript dev/translator-audit.R <id> <outdir>
#
# From the worktree root; reads FRMTMB_EQUIV_LIB as dev/equiv-measure.R.
# Writes <outdir>/<id>.tsv, one row per latent block (the Stan blocks
# check C reads its gradient on: z, zs, zgp, Ymi, Yl, zme, zerr, rcar,
# zcar). Records nothing in the survey.
#
# The test: at each of the pair's three L1 points, the block's translated
# values are reversed (levels, basis coefficients or rows), and then each
# part of it that has its own SD (a matrix block's rows, a gr(by = )
# block's by-groups) is reversed alone (eq_audit_mutations()), which keeps
# its standard-normal prior exactly and moves only what the likelihood
# reads. If brms's log density moves by less than the L1 tolerance
# (1e-6 * max(1, |lp|)) at all three points, a wrong map of that block
# would pass L1: the block is BLIND there. A block's SD near 0 is the
# usual cause (review 2026-11-09, M1). The SDs the translator gives brms
# (sd_, sds_, sdgp_, sderr, sdme_) are printed beside each pair.
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- args[2]
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[smp$id == id, ]
key <- paste(row$feature_a, row$feature_b, sep = " x ")
# Not a sampled pair: a gr(by = ) block with one by-group's SD collapsed,
# built so that the by-group reversal has a case to fire on (no sampled
# pair has gr(by = )). Group A's levels carry a random intercept of SD
# 0.8, group B's none. Measured 2026-11-10: it stops in
# rstan::unconstrain_pars() ("sd_1 ... dims declared=(1,2); dims
# found=(1)"), because the suite's translator has no rule for gr(by = )'s
# sd matrix, so no L1 pair can hold such a block yet; the by-group
# reversal waits for that rule.
if (identical(id, "bycheck")) {
  key <- "audit gr(by = ) check"
  E$eq_specs[[key]] <- function() {
    d <- E$eq_dat(seed = 5001)
    d$f2 <- factor(ifelse(as.integer(d$g) <= 5L, "A", "B"))
    set.seed(5002)
    u <- stats::rnorm(nlevels(d$g), 0, 0.8)
    d$y <- 0.3 + 0.5 * d$x + rnorm(nrow(d)) +
      ifelse(d$f2 == "A", u[as.integer(d$g)], 0)
    list(bf = "bf(y ~ x + (1 | gr(g, by = f2)))", data = d)
  }
}
out_p <- file.path(outdir, paste0(id, ".tsv"))
put <- function(df) {
  utils::write.table(df, out_p, sep = "\t", row.names = FALSE, quote = TRUE)
}
none <- function(why) {
  put(data.frame(id = id, pair = key, block = "", n = NA, points = NA,
                 max_delta = NA, tol = NA, blind = NA, sds = "",
                 note = why))
  cat("audit", id, ":", why, "\n")
  quit(save = "no")
}

E$key <- key
res <- local(envir = E, {
  spec <- eq_specs[[key]]()
  if (length(spec$data2)) spec$brms_args$data2 <- spec$data2
  if (is.null(spec$prior) && !is.null(spec$prior_brms)) {
    spec$prior <- spec$prior_brms
  }
  if (!is.null(spec$skip) || !is.null(spec$nobrms)) {
    return(list(why = "no brms spelling or skipped"))
  }
  spec_l1 <- spec
  if (!is.null(spec$l1_without)) spec_l1$frm_args[spec$l1_without] <- NULL
  fits <- tryCatch({
    if (!is.null(spec$l1_points)) {
      lapply(0:2, function(k) eq_frm_fit(spec_l1, ctrl = list(
        optimizer = eq_jitter_optimizer(spec$l1_points$scale, k),
        restarts = 0, check_se = "ignore")))
    } else {
      f0 <- suppressWarnings(eq_frm_fit(spec_l1))
      c(list(f0), lapply(c(1, -1), function(s) {
        suppressWarnings(eq_frm_fit(spec_l1, ctrl = list(
          optimizer = eq_step_optimizer(s), restarts = 0,
          check_se = "ignore")))
      }))
    }
  }, error = function(e) e)
  if (inherits(fits, "error")) {
    return(list(why = paste("frm() failed:", conditionMessage(fits))))
  }
  bform <- eq_eval(spec$bf_brms %||% spec$bf, "brms")
  ft <- eq_fam_text(spec, "brms")
  family <- if (is.na(ft)) NULL else eq_eval(ft, "brms")
  prior <- eq_brms_prior(spec, bform, family, spec$data)
  data <- brms_mi_observed(fits[[1]], spec$data)
  args <- eq_brms_args(spec, data, prior)
  code <- do.call(brms::make_stancode, c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  rtab <- tryCatch(do.call(brms_ranef_table,
                           c(list(bform, data, family, prior),
                             spec$brms_args)), error = function(e) NULL)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  lp_of <- function(p) {
    rstan::log_prob(sf, rstan::unconstrain_pars(sf, p),
                    adjust_transform = FALSE, gradient = FALSE)
  }
  rev_block <- function(v) {
    if (is.matrix(v)) {
      if (ncol(v) > 1L) v[, ncol(v):1, drop = FALSE] else
        v[nrow(v):1, , drop = FALSE]
    } else {
      array(rev(as.numeric(v)), length(v))
    }
  }
  # A block's parts that carry their own SD, each reversed on its own:
  # the whole block; each row of a matrix block (one coefficient's z over
  # the levels, or one me() term's latent values), since a block with one
  # live and one collapsed row reads as seen when reversed whole (review
  # 2026-11-09, short check, F1: p087); and, for a gr(by = ) block, each
  # row within each by-group's levels (Jby_<i>), whose SDs differ by
  # group. Every reversal permutes standard-normal entries among
  # themselves, so the prior stays exact.
  eq_audit_mutations <- function(nm, v, sdat) {
    out <- list(list(label = nm, value = rev_block(v)))
    if (!is.matrix(v) || ncol(v) < 2L) return(out)
    for (r in seq_len(nrow(v))) {
      if (nrow(v) > 1L) {
        w <- v
        w[r, ] <- rev(v[r, ])
        out[[length(out) + 1L]] <- list(label = sprintf("%s[row %d]", nm, r),
                                         value = w)
      }
      jby <- if (grepl("^z_[0-9]+$", nm)) {
        sdat[[sub("^z_", "Jby_", nm)]]
      }
      if (length(jby) == ncol(v)) {
        for (b in sort(unique(jby))) {
          idx <- which(jby == b)
          w <- v
          w[r, idx] <- rev(v[r, idx])
          out[[length(out) + 1L]] <- list(
            label = sprintf("%s[row %d, by %d]", nm, r, b), value = w)
        }
      }
    }
    out
  }
  per <- list()
  sds <- ""
  for (k in seq_along(fits)) {
    pars <- tryCatch(stan_pars_from_fit(fits[[k]], sdat, code, rtab),
                     error = function(e) e)
    if (inherits(pars, "error")) {
      return(list(why = paste("translator:", conditionMessage(pars))))
    }
    if (k == 1L) {
      sn <- grep("^(sd|sds|sdgp|sdme)_|^sderr", names(pars), value = TRUE)
      sds <- paste(vapply(sn, function(n) paste0(n, "=", paste(
        signif(as.numeric(pars[[n]]), 3), collapse = "/")), ""),
        collapse = " ")
      # a correlation at +/-1 leaves a row of z with no weight even when
      # every SD is large: L's diagonal below the first entry goes to 0
      ln <- grep("^L_[0-9]+$", names(pars), value = TRUE)
      if (length(ln)) {
        sds <- paste(sds, paste(vapply(ln, function(n) paste0(
          "diag(", n, ")=", paste(signif(diag(as.matrix(pars[[n]])), 3),
                                  collapse = "/")), ""), collapse = " "))
      }
    }
    lp0 <- lp_of(pars)
    tol <- 1e-6 * max(1, abs(lp0))
    for (nm in grep(brms_inner_pat, names(pars), value = TRUE)) {
      for (mu in eq_audit_mutations(nm, pars[[nm]], sdat)) {
        v <- pars[[nm]]
        w <- mu$value
        d <- if (length(v) < 2L || isTRUE(all.equal(as.numeric(w),
                                                     as.numeric(v),
                                                     tolerance = 0))) {
          NA_real_
        } else {
          p2 <- pars
          p2[[nm]] <- w
          lp_of(p2) - lp0
        }
        per[[mu$label]] <- rbind(per[[mu$label]], c(delta = d, tol = tol))
      }
    }
  }
  list(per = per, sds = sds, points = length(fits),
       at = if (is.null(spec$l1_points)) "estimate and +/-0.05" else
         paste("fixed points, scale", spec$l1_points$scale))
})
if (!is.null(res$why)) none(res$why)
if (!length(res$per)) none(paste("no latent block; points:", res$at))
df <- do.call(rbind, lapply(names(res$per), function(nm) {
  m <- res$per[[nm]]
  d <- m[, "delta"]
  data.frame(id = id, pair = key, block = nm, n = nrow(m),
             points = res$at,
             max_delta = if (all(is.na(d))) NA else max(abs(d), na.rm = TRUE),
             tol = max(m[, "tol"]),
             blind = if (all(is.na(d))) NA else all(abs(d) < m[, "tol"]),
             sds = res$sds,
             note = if (all(is.na(d))) "one value or symmetric: untestable"
               else "", stringsAsFactors = FALSE)
}))
put(df)
cat("audit", id, ": blocks", nrow(df), "blind", sum(df$blind %in% TRUE), "\n")
