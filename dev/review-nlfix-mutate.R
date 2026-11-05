# Mutants of the lane's source, all behind REV_MUT, in the source copy only.
root <- commandArgs(trailingOnly = TRUE)[1]
sb <- function(path, old, new) {
  p <- file.path(root, path)
  s <- readChar(p, file.info(p)$size, useBytes = TRUE)
  crlf <- grepl("\r\n", s, fixed = TRUE)
  if (crlf) { old <- gsub("\n", "\r\n", old, fixed = TRUE); new <- gsub("\n", "\r\n", new, fixed = TRUE) }
  if (!grepl(old, s, fixed = TRUE, useBytes = TRUE)) stop("not found in ", path, ": ", old)
  s <- sub(old, new, s, fixed = TRUE, useBytes = TRUE)
  writeBin(charToRaw(s), p)
}
M <- 'Sys.getenv("REV_MUT")'
sb("R/frame.R", 'dp[["name"]] %in% (resp$family[["primary_dpars"]] %||% "mu")) {\n        # thresholds replace',
   paste0('(if (', M, ' == "is_primary") is_primary else dp[["name"]] %in% (resp$family[["primary_dpars"]] %||% "mu"))) {\n        # thresholds replace'))
sb("R/par-template.R", 'nl_check_body_read <- function(frame, nll, tpl) {\n',
   paste0('nl_check_body_read <- function(frame, nll, tpl) {\n  if (', M, ' == "guard_off") return(invisible(NULL))\n'))
sb("R/par-template.R", 'if (is.finite(v0) && is.finite(v1) && is.finite(v2) && v0 == v1 &&',
   paste0('if (', M, ' != "guard_obj_off" && is.finite(v0) && is.finite(v1) && is.finite(v2) && v0 == v1 &&'))
sb("R/par-template.R", 'if (has_terms && is.null(q[["Z"]])', paste0('if (', M, ' != "guard_design_off" && has_terms && is.null(q[["Z"]])'))
sb("R/par-template.R", 'nl_domain_start <- function(frame, tpl, taken, announce) {\n',
   paste0('nl_domain_start <- function(frame, tpl, taken, announce) {\n  if (', M, ' == "domain_off") return(tpl)\n'))
wrong <- function(tag) paste0('(if (', M, ' == "', tag, '") 1e-14 * sin(x * 1e14) else 0 * x)')
sb("R/ad-env.R", 'Phi = function(x) RTMB::pnorm(x),', paste0('Phi = function(x) RTMB::pnorm(x) + ', wrong("d_Phi"), ','))
sb("R/ad-env.R", 'inv_cloglog = function(x) -expm1(-exp(x)),', paste0('inv_cloglog = function(x) -expm1(-exp(x)) + ', wrong("d_cloglog"), ','))
sb("R/ad-env.R", 'RTMB::logspace_add(x + 0 * y, y + 0 * x)\n', paste0('RTMB::logspace_add(x + 0 * y, y + 0 * x) + ', wrong("d_lse"), '\n'))
sb("R/ad-env.R", 'RTMB::logspace_sub(x + 0 * y, y + 0 * x)\n', paste0('RTMB::logspace_sub(x + 0 * y, y + 0 * x) + ', wrong("d_lde"), '\n'))
sb("R/ad-env.R", 'hypot = function(x, y) sqrt(x * x + y * y),', paste0('hypot = function(x, y) sqrt(x * x + y * y) + ', wrong("d_hypot"), ','))
sb("R/ad-env.R", 'inv_logit = function(x) RTMB::plogis(x),', paste0('inv_logit = function(x) RTMB::plogis(x) + ', wrong("d_invlogit"), ','))
sb("R/ad-env.R", 'log1p_exp = function(x) RTMB::logspace_add(0 * x, x),', paste0('log1p_exp = function(x) if (', M, ' == "plain_l1pe") log1p(exp(x)) else RTMB::logspace_add(0 * x, x),'))
sb("R/ad-env.R", '1 + s / 2 - s * s / 2\n', paste0('if (', M, ' == "step_zero") 0.5 + s / 2 else 1 + s / 2 - s * s / 2\n'))
sb("R/predict.R", 'kept <- function(v) {\n    if (is.null(newdata)) {', paste0('kept <- function(v) {\n    if (is.null(newdata) && ', M, ' != "kept_off") {'))
sb("R/priors.R", 'if (!is.null(lp[["nl_body"]]) &&\n          !is.null(ob <- ord_lp_block(frame, spec, lp)) &&',
   paste0('if (', M, ' != "thres_rows" && !is.null(lp[["nl_body"]]) &&\n          !is.null(ob <- ord_lp_block(frame, spec, lp)) &&'))
sb("R/parse.R", 'if (!is.null(bform[["autocor"]])) {\n    dpars[[primaries[1L]]]', paste0('if (', M, ' != "ac_drop" && !is.null(bform[["autocor"]])) {\n    dpars[[primaries[1L]]]'))
sb("R/confint.R", 'for (nm in c("center", "cmc", "autocor")) {', paste0('for (nm in c("center", "cmc", if (', M, ' != "ac_refit") "autocor")) {'))
sb("R/parse.R", 'for (p in setdiff(primaries, names(nl_bodies))) {', paste0('for (p in if (', M, ' == "mm_first") primaries[1L] else setdiff(primaries, names(nl_bodies))) {'))
sb("R/ad-env.R", 'f <- get0(nm, envir = parent %||% baseenv(), mode = "function")\n    is.null(f) ||', paste0('f <- get0(nm, envir = parent %||% baseenv(), mode = "function")\n    ', M, ' == "shadow_off" || is.null(f) ||'))
# the fix this review asks for in sim_newdata_dpars(), as a mutant: none. Left to the lane.
cat("mutated\n")
