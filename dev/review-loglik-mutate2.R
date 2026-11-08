# Review of lane loglik, short check of fix round 1: refresh the mutant
# copy in C:/Users/adf44/source/r/rev-loglik-src from the worktree's
# current package files, then apply the env-switched mutants again.
# REV_MUT = nojac | twice | nooffset | halfnorm | remlnoref | discpen |
# nochart; unset is the lane code. Never touches the worktree.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-loglik"
src <- "C:/Users/adf44/source/r/rev-loglik-src"
unlink(src, recursive = TRUE)
dir.create(src)
for (p in c("DESCRIPTION", "NAMESPACE", "NEWS.md", "R", "man", "inst",
            "data", "tests")) {
  file.copy(file.path(wt, p), src, recursive = TRUE)
}
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
  if (mut %in% c("nojac", "nochart")) {
    oo <- options(frmtmb.rev_nojac = (mut == "nojac"),
                  frmtmb.rev_nochart = (mut == "nochart"))
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
sub1("R/fit-end.R", '  ll <- -fit$opt$objective - lp',
     '  ll <- -fit$opt$objective - (if (Sys.getenv("REV_MUT") == "discpen") 0 else lp)')
sub1("R/priors.R", '  prior_base_logdens(x, dist) + jac\n}',
     '  prior_base_logdens(x, dist) + (if (isTRUE(getOption("frmtmb.rev_nojac"))) 0 else jac)\n}')
sub1("R/priors.R", '             ac_trans_logjac(x, dist$trans))',
     '             (if (isTRUE(getOption("frmtmb.rev_nojac"))) 0 else ac_trans_logjac(x, dist$trans)))')
# the simplex charts: density of w with no Jacobian of the chart
sub1("R/priors.R", '    return(dist$lognorm + sum(a * (v - lse)))',
     '    return(dist$lognorm + sum((a - isTRUE(getOption("frmtmb.rev_nochart"))) * (v - lse)))')
sub1("R/priors.R",
     '  out <- dist$lognorm - (D - 1) * log(1 + r2) + (D - 2) * log(2)',
     '  nch <- isTRUE(getOption("frmtmb.rev_nochart"))
  out <- if (nch) dist$lognorm else dist$lognorm - (D - 1) * log(1 + r2) + (D - 2) * log(2)')
sub1("R/priors.R", '    out <- out + (2 * a[k] - 1) * log(abs(u[k]))',
     '    out <- out + (if (nch) 2 * (a[k] - 1) else 2 * a[k] - 1) * log(abs(u[k]))')
cat("mutants applied\n")
