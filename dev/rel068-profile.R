# The ordmix review's instrumented profile (dev/ordmix-rev2-profile.R),
# for the 0.68.0 merged tree, with one more hook: lane fixes'
# nl_flat_message() is wrapped too, so that a fit where both the flat
# warning and the non-finite-gradient warning fire can be seen. Used as
# R_PROFILE_USER by dev/rel068-prof-par.sh. Each line goes to
# $ORDMIX_REV2_LOG/<pid>.txt: "FLAT" when nl_flat_message() returns a
# message, then "CC" from check_convergence() for the same fit (frm()
# runs the flat check first). Logging failures are swallowed so the fit
# is unchanged.
local({
  logdir <- Sys.getenv("ORDMIX_REV2_LOG")
  if (!nzchar(logdir)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    out_file <- file.path(logdir, paste0(Sys.getpid(), ".txt"))
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
            file = out_file, append = TRUE)
      }, error = function(e) NULL)
      out
    }
    environment(wrap) <- environment()
    orig_flat <- get("nl_flat_message", envir = ns)
    wrap_flat <- function(obj, opt, frame) {
      m <- orig_flat(obj, opt, frame)
      if (!is.null(m)) {
        tryCatch(cat(sprintf("FLAT file=%s\n",
                             Sys.getenv("ORDMIX_REV2_FILE")),
                     file = out_file, append = TRUE),
                 error = function(e) NULL)
      }
      m
    }
    environment(wrap_flat) <- environment()
    for (nm in c("check_convergence", "nl_flat_message")) {
      unlockBinding(nm, ns)
      assign(nm, if (nm == "check_convergence") wrap else wrap_flat,
             envir = ns)
      lockBinding(nm, ns)
    }
  })
})
