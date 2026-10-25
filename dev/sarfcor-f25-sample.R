# Fix round 1, F2 and F5, frmtmb.sample: the residual terms' names and
# brms's flat default on a bounded parameter, read from core's
# frmtmb_autocor_pars(). Run on the worktree and the stage-A copy:
#   Rscript dev/sarfcor-f25-sample.R <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
edit("extensions/frmtmb.sample/R/draws-brms.R", list(
  c("#' `cosy`, `cortime__`, `lagsar` and `errorsar` (frmtmb.spatial's
#' `sar()`), and `me()`'s `sdme_` and `corme__`. The sampler",
    "#' `cosy`, `cortime__`, the parameters of a residual term of core's
#' registry under the names it gives them (`lagsar`, `errorsar`;
#' `frmtmb_autocor_pars()`), and `me()`'s `sdme_` and `corme__`. The
#' sampler"),
  c("  # brms's names, as its regex_autocor_pars() lists them
  pat <- paste0(\"^(sd_|cor_|sds_|rescor__|cortime__|sdme_|corme__|\",
                \"ar[[]|ma[[]|cosy$|lagsar$|errorsar$)\")
  nm <- setdiff(grep(pat, names(ev), value = TRUE), have)",
    "  pat <- paste0(\"^(sd_|cor_|sds_|rescor__|cortime__|sdme_|corme__|\",
                \"ar[[]|ma[[]|cosy$)\")
  # a registry term names its own parameters; core lists them, so this
  # package keeps no list of those names
  nm <- setdiff(union(grep(pat, names(ev), value = TRUE),
                      intersect(frmtmb_autocor_pars(fit)$name, names(ev))),
                have)")))
edit("extensions/frmtmb.sample/R/sample.R", list(
  c("  # and on an unstr() residual correlation, class cortime, per response",
    "  # brms's (flat) on a residual term's bounded parameter (sar()'s
  # lagsar), between the bounds brms gives it: a uniform density on the
  # natural parameter, which core carries to the internal one with the
  # log-Jacobian of the term's map, so rho is flat on its interval and
  # not on the logit scale the term parameterizes it on
  ap <- frmtmb_autocor_pars(fit)
  for (i in which(is.finite(ap$lower) & is.finite(ap$upper))) {
    p <- set_prior(sprintf(\"uniform(%.17g, %.17g)\", ap$lower[i],
                           ap$upper[i]), class = ap$name[i],
                   resp = if (multi) ap$resp[i] else \"\")
    p[[1L]]$flat_interval <- TRUE
    add(p)
  }
  # and on an unstr() residual correlation, class cortime, per response"),
  c("                          if (identical(sp$class, \"xi\")) \"\"
                          else if (isTRUE(s$natural)) \"  [natural scale]\"",
    "                          if (identical(sp$class, \"xi\")) \"\"
                          else if (isTRUE(s$flat_interval))
                            \"  [flat on its interval, brms's (flat)]\"
                          else if (isTRUE(s$natural)) \"  [natural scale]\"")))
cat("edited", root, "\n")
