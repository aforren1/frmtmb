# Reviewer of lane optrap, consolidation check: an R_PROFILE_USER that
# logs every newton_polish() call a test file makes: random effects or
# not, the `ill` flag, the number of parameters, the steps, the gain,
# the objective before and after, and the seconds. Lines go to
# $OPTRAP_REV4_LOG/<pid>.txt; logging failures are swallowed so the fit
# is unchanged.
local({
  logdir <- Sys.getenv("OPTRAP_REV4_LOG")
  if (!nzchar(logdir)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    if (!exists("newton_polish", envir = ns, inherits = FALSE)) return()
    out_file <- file.path(logdir, paste0(Sys.getpid(), ".txt"))
    a <- commandArgs(trailingOnly = TRUE)
    tag <- if (length(a) >= 3) paste0(a[2], "/", basename(a[3])) else "?"
    orig <- get("newton_polish", envir = ns)
    w <- function(obj, opt, bounds, par_units = NULL, control = NULL,
                  ill = FALSE) {
      t0 <- proc.time()[[3]]
      res <- orig(obj, opt, bounds, par_units, control, ill = ill)
      el <- proc.time()[[3]] - t0
      tryCatch(cat(sprintf("POLISH\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%.3f\n",
                           tag, length(obj$env$random) > 0L, isTRUE(ill),
                           length(opt$par), as.integer(opt$convergence %||% NA),
                           format(res$polish[["steps"]] %||% 0),
                           format(opt$objective, digits = 15),
                           format(res$objective, digits = 15), el),
                   file = out_file, append = TRUE),
               error = function(e) NULL)
      res
    }
    environment(w) <- environment()
    unlockBinding("newton_polish", ns)
    assign("newton_polish", w, envir = ns)
    lockBinding("newton_polish", ns)
  })
})
