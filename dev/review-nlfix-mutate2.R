# Mutants of fix round 1's source, all behind REV_MUT, in the source copy only.
root <- commandArgs(trailingOnly = TRUE)[1]
sb <- function(path, old, new) {
  p <- file.path(root, path)
  s <- readChar(p, file.info(p)$size, useBytes = TRUE)
  if (grepl("\r\n", s, fixed = TRUE)) { old <- gsub("\n", "\r\n", old, fixed = TRUE); new <- gsub("\n", "\r\n", new, fixed = TRUE) }
  if (!grepl(old, s, fixed = TRUE, useBytes = TRUE)) stop("not found in ", path, ": ", old)
  writeBin(charToRaw(sub(old, new, s, fixed = TRUE, useBytes = TRUE)), p)
}
M <- 'Sys.getenv("REV_MUT")'
wrong <- function(tag) paste0('(if (', M, ' == "', tag, '") 1e-14 * sin(x * 1e14) else 0 * x)')
sb("R/ad-env.R", 'Phi = function(x) RTMB::pnorm(x),', paste0('Phi = function(x) RTMB::pnorm(x) + ', wrong("d_Phi"), ','))
sb("R/ad-env.R", 'exp(RTMB::logspace_sub(0 * u, -u))\n', paste0('exp(RTMB::logspace_sub(0 * u, -u)) + ', wrong("d_cloglog"), '\n'))
sb("R/ad-env.R", 'RTMB::logspace_add(x + 0 * y, y + 0 * x)\n', paste0('RTMB::logspace_add(x + 0 * y, y + 0 * x) + ', wrong("d_lse"), '\n'))
sb("R/ad-env.R", 'RTMB::logspace_sub(x + 0 * y, y + 0 * x)\n', paste0('RTMB::logspace_sub(x + 0 * y, y + 0 * x) + ', wrong("d_lde"), '\n'))
sb("R/ad-env.R", 'm * sqrt(1 + r * r)\n', paste0('m * sqrt(1 + r * r) + ', wrong("d_hypot"), '\n'))
sb("R/confint.R", 'for (nm in c("center", "cmc", "autocor")) {', paste0('for (nm in c("center", "cmc", if (', M, ' != "ac_refit") "autocor")) {'))
sb("R/ad-env.R", 'nl_pow <- function(x, y) {\n', paste0('nl_pow <- function(x, y) {\n  if (', M, ' == "pow_off") return(x^y)\n'))
sb("R/par-template.R", 'nl_fit_off_saddle <- function(run, frame, start, prior) {\n', paste0('nl_fit_off_saddle <- function(run, frame, start, prior) {\n  if (', M, ' == "saddle_off") return(run(start, TRUE))\n'))
sb("R/par-template.R", 'if (use_curvature && nl_coef_curved(frame, est, bc$coefs)) next', paste0('if (', M, ' != "curv_off" && use_curvature && nl_coef_curved(frame, est, bc$coefs)) next'))
sb("R/par-template.R", 'flagged <- !identical(fit$opt$convergence, 0L) ||', paste0('flagged <- ', M, ' != "flag_off" && !identical(fit$opt$convergence, 0L) ||'))
sb("R/par-template.R", 'nl_restart_credible <- function(r, ll) {\n', paste0('nl_restart_credible <- function(r, ll) {\n  if (', M, ' == "cred_off") return(TRUE)\n'))
sb("R/simulate-newdata.R", 'eta <- if (!is.null(pd$lp)) {', paste0('eta <- if (', M, ' != "sim_body" && !is.null(pd$lp)) {'))
cat("mutated\n")
