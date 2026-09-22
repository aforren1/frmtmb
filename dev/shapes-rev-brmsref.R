# Reviewer: brms is the external judge. Fit the SAME data in brms and
# record (a) the structure, dimnames, row order and class of every
# method item 2.6f moved, and (b) the posterior sd of the same slopes
# marginaleffects reports, which is the judge for priority 1.
#
# One R process; each model compiles once. Writes an RDS the comparison
# step reads, so nothing is compared inside a process that has brms and
# frmtmb on the same path.

Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"

# rev_data() only needs base R, so it is sourced without frmtmb
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
dd <- rev_data()
dd$y2 <- with(dd, rnorm(nrow(dd), 0.3 - 0.2 * x, 1))

ctl <- list(chains = 2, iter = 2000, refresh = 0, seed = 20260918,
            silent = 2, backend = "rstan")
run <- function(...) do.call(brm, c(list(...), ctl, list(data = dd)))

models <- list(
  gaussian = quote(run(bf(y ~ x + f) + gaussian())),
  binomial = quote(run(bf(bin ~ x + z) + bernoulli())),
  mixed    = quote(run(bf(ymix ~ x + (1 | g)) + gaussian())),
  ordinal  = quote(run(bf(ord ~ x) + cumulative())),
  distreg  = quote(run(bf(y ~ x + f, sigma ~ z) + gaussian())),
  multivar = quote(run(bf(mvbind(y, y2) ~ x) + gaussian()))
)

grab <- function(fit) {
  s <- summary(fit)
  safe <- function(e) tryCatch(suppressWarnings(e),
                               error = function(c)
                                 structure(list(msg = conditionMessage(c)),
                                           class = "revErr"))
  list(
    fixef        = safe(fixef(fit)),
    fixef_dimnm  = safe(dimnames(fixef(fit))),
    fixef_class  = safe(class(fixef(fit))),
    vcov         = safe(vcov(fit)),
    vcov_dimnm   = safe(dimnames(vcov(fit))),
    vcov_corr    = safe(dimnames(vcov(fit, correlation = TRUE))),
    summary_names= safe(names(s)),
    summary_class= safe(class(s)),
    fixed        = safe(s$fixed),
    fixed_cols   = safe(colnames(s$fixed)),
    fixed_rows   = safe(rownames(s$fixed)),
    spec_pars    = safe(s$spec_pars),
    random       = safe(s$random),
    cor_pars     = safe(s$cor_pars),
    ngrps        = safe(ngrps(fit)),
    ngrps_class  = safe(class(ngrps(fit))),
    fitted_dim   = safe(dim(fitted(fit))),
    fitted_dimnm = safe(dimnames(fitted(fit))),
    fitted_class = safe(class(fitted(fit))),
    resid_dim    = safe(dim(residuals(fit))),
    resid_dimnm  = safe(dimnames(residuals(fit))),
    predict_dim  = safe(dim(predict(fit))),
    predict_dimnm= safe(dimnames(predict(fit))),
    predict_nosum= safe(dim(predict(fit, summary = FALSE))),
    predict_nd3  = safe(dim(predict(fit, ndraws = 3))),
    variables    = safe(variables(fit)),
    nsamples_warn= safe(tryCatch(
      withCallingHandlers({ nsamples(fit); "no warning" },
        warning = function(w) { assign("W", conditionMessage(w),
                                       envir = globalenv())
          invokeRestart("muffleWarning") }),
      error = function(e) conditionMessage(e))),
    nsamples_msg = safe(get0("W", envir = globalenv())),
    parnames_msg = safe(tryCatch(
      withCallingHandlers({ parnames(fit); "ok" },
        warning = function(w) conditionMessage(w)),
      error = function(e) paste0("ERR: ", conditionMessage(e)))),
    post_samp_msg= safe(tryCatch(
      withCallingHandlers({ posterior_samples(fit); "ok" },
        warning = function(w) conditionMessage(w)),
      error = function(e) paste0("ERR: ", conditionMessage(e)))),
    nsamples_val = safe(suppressWarnings(nsamples(fit))),
    slopes       = safe({
      sl <- marginaleffects::avg_slopes(fit)
      data.frame(term = as.character(sl$term),
                 estimate = sl$estimate,
                 sd = if ("conf.low" %in% names(sl))
                   (sl$conf.high - sl$conf.low) / (2 * qnorm(0.975)) else NA)
    }),
    slopes_raw   = safe(as.data.frame(marginaleffects::avg_slopes(fit)))
  )
}

out <- list()
for (nm in names(models)) {
  cat("=== ", nm, " ===\n"); flush.console()
  t0 <- proc.time()
  fit <- tryCatch(eval(models[[nm]]), error = function(e) {
    cat("  FIT ERROR: ", conditionMessage(e), "\n"); NULL })
  cat("  elapsed ", round((proc.time() - t0)[3], 1), " s\n")
  if (is.null(fit)) next
  if (exists("W", envir = globalenv())) rm("W", envir = globalenv())
  out[[nm]] <- grab(fit)
  saveRDS(out, file.path(TREE, "dev/shapes-rev-brmsref.rds"))
}
cat("wrote dev/shapes-rev-brmsref.rds with ", length(out), " models\n")
