# Review of lane arstep: compare the snapshots of review-arstep-state.R.
# Within the lane process: default vs nosecond (a lost second descent must
# leave the fit as if it had not run), default vs nosettle (does the
# contamination mutant show), default vs atopt (a won second descent: do the
# readers read the kept optimum), default vs again (determinism).
# Across processes: lane default vs base default on rows the second lost.
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
L <- readRDS(file.path(RV_OUT, "state-lane.rds"))
B <- readRDS(file.path(RV_OUT, "state-base.rds"))
strip <- function(o) { o$gp_second <- NULL; o }
fields <- c("opt", "est", "last.par", "last.par.best", "value.best", "vcov",
            "summ", "pred", "fitted", "ll", "start", "warns", "sdr_par")
cmp <- function(a, b) {
  a$opt <- strip(a$opt); b$opt <- strip(b$opt)
  vapply(fields, function(f) {
    x <- a[[f]]; y <- b[[f]]
    if (identical(x, y)) return("same")
    if (is.list(x) || is.character(x) || is.null(x) || is.null(y)) {
      if (f == "opt") {
        d <- setdiff(union(names(x), names(y)),
                     names(x)[vapply(names(x), function(k) identical(x[[k]], y[[k]]), NA)])
        return(paste0("DIFF[", paste(d, collapse = "/"), "]"))
      }
      if (f == "est") {
        num <- unlist(x); num2 <- unlist(y)
        if (length(num) == length(num2)) {
          return(sprintf("DIFF %.1e", max(abs(num - num2) / pmax(1, abs(num2)))))
        }
      }
      return("DIFF")
    }
    if (length(x) != length(y)) return("DIFF(len)")
    sprintf("DIFF %.1e", max(abs(as.numeric(x) - as.numeric(y)) /
                              pmax(1, abs(as.numeric(y))), na.rm = TRUE))
  }, "")
}
for (k in names(L)) {
  l <- L[[k]]
  won <- isTRUE(l$default$opt$gp_second[["gain"]] > 0)
  cat(sprintf("== %s second %s (gain %.4g); alt left in env: %s\n", k,
              if (won) "WON" else "lost", l$default$opt$gp_second[["gain"]],
              l$default$alt_left))
  show <- function(lbl, r) {
    cat(sprintf("SC %-18s %-22s %s\n", k, lbl,
                paste(names(r), r, sep = "=", collapse = " ")))
  }
  show("lane default~again", cmp(l$default, l$again))
  show("base default~again", cmp(B[[k]]$default, B[[k]]$again))
  show("default~nosecond", cmp(l$default, l$nosecond))
  show("default~nosettle", cmp(l$default, l$nosettle))
  show("default~atopt", cmp(l$default, l$atopt))
  show("lane~base (x-proc)", cmp(l$default, B[[k]]$default))
  # the readers at the kept optimum: last.par.best outer == opt$par
  env_ok <- function(s) {
    lp <- s$last.par.best
    r <- grepl("^b$", names(lp))
    outer <- lp[!r]
    c(best_at_opt = isTRUE(all(unname(outer) == unname(s$opt$par))),
      value_at_opt = isTRUE(s$value.best == s$opt$objective),
      last_eq_best = identical(s$last.par, s$last.par.best),
      sdr_at_opt = isTRUE(all(unname(s$sdr_par) == unname(s$opt$par))))
  }
  for (arm in c("default", "nosettle", "nosecond")) {
    e <- env_ok(l[[arm]])
    cat(sprintf("SE %-18s %-9s %s\n", k, arm, paste(names(e), e, sep = "=", collapse = " ")))
  }
  cat(sprintf("SO %-18s opt evals default %s nosecond %s | nonfinite_trials %s %s | iterations %s %s\n", k,
              l$default$opt$evals, l$nosecond$opt$evals,
              l$default$opt$nonfinite_trials %||% NA, l$nosecond$opt$nonfinite_trials %||% NA,
              l$default$opt$iterations, l$nosecond$opt$iterations))
}
cat("DONE\n")
