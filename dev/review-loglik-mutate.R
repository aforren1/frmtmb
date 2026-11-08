# Review of lane loglik: apply env-switched mutants to the copy in
# C:/Users/adf44/source/r/rev-loglik-src (never to the worktree).
# REV_MUT = nojac | twice | nooffset | halfnorm | remlnoref picks one at
# run time, so one install carries all of them; unset is the lane code.
src <- "C:/Users/adf44/source/r/rev-loglik-src"
sub1 <- function(path, old, new) {
  p <- file.path(src, path)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop(path, ": ", n, " matches of ", old)
  s <- sub(old, new, s, fixed = TRUE)
  writeLines(s, p, useBytes = TRUE)
}
sub1("R/methods-fit.R",
     '  -as.numeric(neg_log_prior_fn(ent)(fit[["estimates"]]))\n}',
     paste0(
'  mut <- Sys.getenv("REV_MUT")
  if (mut == "nojac") {
    oo <- options(frmtmb.rev_nojac = TRUE)
    on.exit(options(oo), add = TRUE)
  }
  if (mut == "nooffset") ent <- lapply(ent, function(e) { e$offset <- NULL; e })
  val <- -as.numeric(neg_log_prior_fn(ent)(fit[["estimates"]]))
  if (mut == "halfnorm") {
    val <- val + log(2) * sum(vapply(ent, function(e)
      identical(e$scale, "sd") * length(e$idx), 0))
  }
  val
}'))
sub1("R/methods-fit.R", '    if (isTRUE(fit$REML)) "beta")',
     '    if (isTRUE(fit$REML) && Sys.getenv("REV_MUT") != "remlnoref") "beta")')
sub1("R/methods-fit.R", '  structure(-object$opt$objective - lp,',
     '  structure(-object$opt$objective - lp * (if (Sys.getenv("REV_MUT") == "twice") 2 else 1),')
sub1("R/priors.R", '  prior_base_logdens(x, dist) + jac\n}',
     '  prior_base_logdens(x, dist) + (if (isTRUE(getOption("frmtmb.rev_nojac"))) 0 else jac)\n}')
sub1("R/priors.R", '             ac_trans_logjac(x, dist$trans))',
     '             (if (isTRUE(getOption("frmtmb.rev_nojac"))) 0 else ac_trans_logjac(x, dist$trans)))')
cat("mutants applied\n")
