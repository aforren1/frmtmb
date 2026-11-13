# Review of lane ordstart: which refit paths rerun the jittered starts,
# and what each costs, on cum_sratio seed 4 of the 40-fit data (a step
# function at 0.81.0 that the lane's starts move to a proper maximum).
# Counts mixture_starts_pick() calls (one per mixture_starts() run that
# reached its pick) and times each path.
#   Rscript ordstart-review-refits.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
ns <- asNamespace("frmtmb")
rec <- new.env()
rec$picks <- 0L
suppressMessages(trace("mixture_starts_pick", where = ns, print = FALSE,
  exit = quote(assign("picks", get("picks", envir = .GlobalEnv$rec) + 1L,
                      envir = .GlobalEnv$rec))))
assign("rec", rec, envir = .GlobalEnv)
q <- function(e) suppressWarnings(suppressMessages(e))
probe <- function(label, expr) {
  rec$picks <- 0L
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(q(eval.parent(substitute(expr))), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("PATH %-34s picks %3d elapsed %7.2f %s\n", label, rec$picks,
              el, if (inherits(r, "error")) paste("ERROR", conditionMessage(r))
              else ""))
  invisible(r)
}
d <- mk(4, mix = TRUE)
fam <- mixture(cumulative(), sratio())
f0 <- probe("frm starts = 0", frm(bf(y ~ x + z), family = fam, data = d,
                                  control = frmtmb_control(starts = 0)))
f <- probe("frm default", frm(bf(y ~ x + z), family = fam, data = d))
cat("logLik default", as.numeric(logLik(f)), "starts = 0",
    as.numeric(logLik(f0)), "\n")
set.seed(1)
ynew <- simulate(f, nsim = 1)[[1]]
r1 <- probe("refit(newresp)", refit(f, ynew))
cat("refit logLik", as.numeric(logLik(r1)), "\n")
probe("frm_bootstrap(nsim = 5)", frm_bootstrap(f, nsim = 5, seed = 1))
probe("frm_allfit()", frm_allfit(f))
probe("frm_allfit(start_from_mle = FALSE)",
      frm_allfit(f, start_from_mle = FALSE))
d2 <- d[-(1:10), ]
probe("update(data = )", update(f, data = d2))
probe("update(formula. = ~ . - z)", update(f, formula. = ~ . - z))
probe("confint(method = 'profile', parm = 1)",
      confint(f, parm = names(f$opt$par)[1], method = "profile"))
probe("diagnose()", diagnose(f))
probe("summary()", summary(f))
cat("DONE\n")
