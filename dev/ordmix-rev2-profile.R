# Reviewer of lane ordmix, re-check: an R profile that wraps the INSTALLED
# lane build's check_convergence() in its namespace, without rebuilding
# anything, and appends one line per fit to $ORDMIX_REV2_LOG/<pid>.txt:
# whether the gradient at the optimum is finite, whether the new
# non-finite-gradient warning was raised, and the family. The test file
# comes from $ORDMIX_REV2_FILE. Logging failures are swallowed so the
# fit is unchanged.
local({
  logdir <- Sys.getenv("ORDMIX_REV2_LOG")
  if (!nzchar(logdir)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    orig <- get("check_convergence", envir = ns)
    wrap <- function(fit, control) {
      g <- tryCatch(drop(fit$obj$gr(fit$opt$par)),
                    error = function(e) NULL)
      out <- orig(fit, control)
      tryCatch({
        nf <- !is.null(g) && length(g) > 0 && any(!is.finite(g))
        w <- any(grepl("gradient at the reported optimum is not finite",
                       out$warnings, fixed = TRUE))
        fam <- paste(vapply(fit$spec$responses, function(r) {
          paste(as.character(r$family[["family"]]), collapse = "/")
        }, ""), collapse = "+")
        bad <- if (nf) paste(unique(names(fit$opt$par)[!is.finite(g)]),
                             collapse = ",") else ""
        cat(sprintf("CC file=%s nonfinite=%s warned=%s conv=%s fam=%s bad=%s\n",
                    Sys.getenv("ORDMIX_REV2_FILE"), nf, w,
                    fit$opt$convergence, gsub(" ", "", fam), bad),
            file = file.path(logdir, paste0(Sys.getpid(), ".txt")),
            append = TRUE)
      }, error = function(e) NULL)
      out
    }
    environment(wrap) <- environment()
    unlockBinding("check_convergence", ns)
    assign("check_convergence", wrap, envir = ns)
    lockBinding("check_convergence", ns)
  })
})
