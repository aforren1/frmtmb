# The review's -83.81 and -83.74 (finding 3 of the equiv review): p094's
# L1 constant with frmtmb's nugget replaced by a smaller one, read through
# the translator before the order fix and after it, on brms's own program.
# The last row puts brms's exact nugget in frmtmb's kernel, so the two
# models are one and what is left is double precision on a kernel
# conditioned about 1e13.
#
#   Rscript dev/gpgap-review83.R
source("dev/gpgap-lib.R")
local(envir = E, {
  ns <- asNamespace("frmtmb")
  spec <- eq_specs[["lognormal x gp"]]()
  set_nug <- function(v) {
    unlockBinding("gp_nugget", ns)
    assign("gp_nugget", v, envir = ns)
    lockBinding("gp_nugget", ns)
  }
  set_corr <- function(fn) {
    unlockBinding("gp_corr", ns)
    assign("gp_corr", fn, envir = ns)
    lockBinding("gp_corr", ns)
  }
  bform <- eq_eval(spec$bf, "brms")
  family <- eq_eval(eq_fam_text(spec, "brms"), "brms")
  prior <- eq_brms_prior(spec, bform, family, spec$data)
  args <- eq_brms_args(spec, spec$data, prior)
  code <- do.call(brms::make_stancode, c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  sf <- gpgap_sf(code, sdat)
  new_order <- brms_gp_order
  read <- function(fit) {
    e <- fit$obj$env
    ours0 <- -e$f(e$last.par.best)
    one <- function(ord) {
      brms_gp_order <<- ord
      environment(stan_pars_from_fit) <<- environment(new_order)
      p <- stan_pars_from_fit(fit, sdat, code)
      gpgap_lp(sf, p) - (ours0 + attr(p, "logJ"))
    }
    old <- function(sdat, zname) {
      seq_len(NROW(sdat[[sub("^zgp_", "Xgp_", zname)]]))
    }
    r <- c(old = one(old), new = one(new_order))
    brms_gp_order <<- new_order
    r
  }
  orig <- get("gp_corr", envir = ns)
  for (v in c(1e-6, 1e-9, 1e-12)) {
    set_nug(v)
    fit <- suppressWarnings(eq_frm_fit(spec))
    r <- read(fit)
    cat(sprintf("GP83 relative nugget %.0e: conv %d, old order %.4f, brms order %.6f\n",
                v, fit$opt$convergence, r[["old"]], r[["new"]]))
  }
  set_nug(1e-6)
  set_corr(function(theta, blk) {
    Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
    exp(-Q) + diag(1, nrow(Q)) * (1e-12 * exp(-2 * theta[1]))
  })
  fit <- suppressWarnings(eq_frm_fit(spec))
  r <- read(fit)
  cat(sprintf("GP83 brms's nugget in frmtmb's kernel: conv %d, old order %.4f, brms order %.3e\n",
              fit$opt$convergence, r[["old"]], r[["new"]]))
  set_corr(orig)
})
