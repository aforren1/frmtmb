# Lane ci691: does nlminb (PORT) return a FINITE trial it rejected as
# `par`, with `objective` the best value, so that fn(par) differs from
# objective? nlminb_best_par() repairs this only when that last trial
# was not finite. Runs nlminb on the dry-run objective of the cs()
# ordinal mixture of test-optima.R, seeds 1 to 20, each from the default
# start and from 4 starts moved by a few ulps, and classifies each run:
#   same     fn(par) is the reported objective
#   nonfin   par is the last trial and it was not finite (repaired today)
#   finrej   par is the last trial, finite, and not the reported value
#   other    fn(par) differs and par is not the last trial
#   Rscript dev/ci691-port-probe.R <lib>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
optima_csmix_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
tab <- list()
for (s in 1:20) {
  u <- suppressWarnings(frm(bf(y ~ x + cs(z)),
                            family = mixture(cumulative(), sratio()),
                            data = optima_csmix_data(s),
                            dry_run = "objective"))
  obj <- u$obj
  for (k in 0:4) {
    st <- obj$par
    if (k > 0) {
      set.seed(k)
      st <- st * (1 + sample(c(-4:-1, 1:4), length(st), TRUE) *
                    .Machine$double.eps)
    }
    fnw <- frmtmb:::nlminb_trial_fn(obj$fn)
    r <- suppressWarnings(stats::nlminb(st, fnw$fn, obj$gr,
      control = list(eval.max = 1000, iter.max = 1000)))
    last <- fnw$last()
    best <- fnw$best()
    v <- obj$fn(r$par)
    is_last <- identical(unname(as.numeric(r$par)),
                         unname(as.numeric(last$par)))
    cls <- if (identical(v, r$objective)) "same" else if (is_last &&
      !is.finite(last$val)) "nonfin" else if (is_last) "finrej" else "other"
    tab[[length(tab) + 1L]] <- data.frame(
      seed = s, k = k, code = r$convergence, class = cls,
      diff = v - r$objective, best_is_obj = identical(best$val, r$objective),
      fn_best = identical(obj$fn(best$par), best$val))
  }
}
tab <- do.call(rbind, tab)
print(tab[tab$class != "same", ], digits = 4)
cat("\nruns", nrow(tab), "\n")
print(table(tab$class, tab$code))
cat("best value is the reported objective in", sum(tab$best_is_obj), "of",
    nrow(tab), "\n")
cat("fn(best par) reproduces best value in", sum(tab$fn_best), "of",
    nrow(tab), "\n")
