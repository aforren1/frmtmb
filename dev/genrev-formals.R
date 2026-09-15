LIB <- commandArgs(trailingOnly = TRUE)[1]
BASE <- commandArgs(trailingOnly = TRUE)[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
owners <- list(
  as_draws = "posterior", as_draws_array = "posterior",
  as_draws_df = "posterior", as_draws_list = "posterior",
  as_draws_matrix = "posterior", as_draws_rvars = "posterior",
  nchains = "posterior", ndraws = "posterior", niterations = "posterior",
  nvariables = "posterior", variables = "posterior", loo = "loo",
  loo_compare = "loo", waic = "loo", bayes_R2 = "rstantools",
  prior_summary = "rstantools", pp_check = "bayesplot",
  conditional_effects = "brms", expose_functions = "brms",
  hypothesis = "brms", posterior_summary = "brms", LOO = "brms",
  WAIC = "brms", ngrps = "brms", fixef = "nlme", ranef = "nlme",
  VarCorr = "nlme", refit = "generics")
extra <- list(ngrps = "lme4", refit = "lme4")
fmt <- function(f) {
  if (is.function(f)) paste(names(formals(f)), collapse = ",") else "<absent>"
}
gv <- function(p, n) tryCatch(getExportedValue(p, n), error = function(e) NULL)
say <- function(...) cat(sprintf(...))
fixns <- loadNamespace("frmtmb")
basef <- new.env()
lazyLoad(file.path(BASE, "frmtmb", "R", "frmtmb"), envir = basef)
say("%-22s %-11s %-26s %-26s %-26s\n", "name", "owner", "OWNER formals",
    "FIX core formals", "BASE core formals")
for (n in names(owners)) {
  o <- gv(owners[[n]], n)
  b <- if (exists(n, envir = basef, inherits = FALSE)) get(n, envir = basef) else NULL
  fx <- tryCatch(get(n, envir = fixns), error = function(e) NULL)
  flag <- if (!is.null(o) && !is.null(fx) && !identical(fmt(o), fmt(fx)))
    "  <== FIX != OWNER" else ""
  chg <- if (!is.null(b) && !identical(fmt(b), fmt(fx))) "  <== CHANGED" else ""
  say("%-22s %-11s %-26s %-26s %-26s%s%s\n", n, owners[[n]], fmt(o), fmt(fx),
      fmt(b), flag, chg)
  if (!is.null(extra[[n]])) {
    say("%-22s %-11s %-26s\n", "", extra[[n]], fmt(gv(extra[[n]], n)))
  }
}
say("\n== frmtmb.sample methods against the generic they now sit on ==\n")
sns <- loadNamespace("frmtmb.sample")
for (n in names(owners)) {
  mn <- paste0(n, ".frmtmb_draws")
  m <- tryCatch(get(mn, envir = sns), error = function(e) NULL)
  if (is.null(m)) next
  o <- gv(owners[[n]], n)
  if (is.null(o)) next
  gfor <- names(formals(o))
  mfor <- names(formals(m))
  miss <- setdiff(setdiff(gfor, "..."), mfor)
  dots <- ("..." %in% gfor) && !("..." %in% mfor)
  say("%-30s gen(%-10s)=%-18s meth=%-34s %s\n", mn, owners[[n]],
      paste(gfor, collapse = ","), paste(mfor, collapse = ","),
      if (length(miss) || dots)
        paste("MISMATCH:", paste(c(miss, if (dots) "..."), collapse = ","))
      else "ok")
}
cat("GENREVDONE\n")
