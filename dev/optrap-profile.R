# Lane optrap: an R_PROFILE_USER that logs every fit a test file makes,
# so that two builds' runs of the suites can be compared fit by fit
# (dev/optrap-fitdiff.R): which fits change their optimum, and by how
# much. A fit is one call of fit_assembled(); its key is a digest of
# the response, the outer parameter names, REML, quadrature and the
# call, numbered by occurrence within the file, so the same fit made
# by the same test line matches across builds. An error is logged and
# re-raised. Lines go to $OPTRAP_PROF_LOG/<pid>.txt; logging failures
# are swallowed so the fit is unchanged.
local({
  logdir <- Sys.getenv("OPTRAP_PROF_LOG")
  if (!nzchar(logdir)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    out_file <- file.path(logdir, paste0(Sys.getpid(), ".txt"))
    a <- commandArgs(trailingOnly = TRUE)
    tag <- if (length(a) >= 3) paste0(a[2], "/", basename(a[3])) else "?"
    o_fa <- get("fit_assembled", envir = ns)
    w_fa <- function(spec, frame, bform, cl, ...) {
      key <- tryCatch(digest::digest(list(frame[["y"]], deparse(cl),
                                          list(...)[["REML"]],
                                          isTRUE(list(...)[["quadrature"]]))),
                      error = function(e) "nokey")
      fit <- tryCatch(o_fa(spec, frame, bform, cl, ...), error = function(e) {
        tryCatch(cat(sprintf("FIT\t%s\t%s\tERROR\t%s\n", tag, key,
                             substr(gsub("[\r\n\t]+", " ",
                                         conditionMessage(e)), 1, 150)),
                     file = out_file, append = TRUE),
                 error = function(e2) NULL)
        stop(e)
      })
      tryCatch({
        op <- fit$opt
        unfitted <- is.null(op)
        cat(sprintf("FIT\t%s\t%s\tOK\t%s\t%s\t%d\t%s\t%s\t%s\n", tag, key,
                    if (unfitted) "NA" else
                      format(op$objective, digits = 15),
                    if (unfitted) "NA" else op$convergence,
                    if (unfitted) 0L else length(op$par),
                    paste(names(fit$opt$par), collapse = ","),
                    if (unfitted) "unfitted" else
                      substr(op$message %||% "", 1, 40),
                    tryCatch(substr(gsub("[\t\r\n ]+", " ", paste(
                      c(deparse(bform[["formula"]] %||% cl[["formula"]]),
                        vapply(spec$responses, function(r) {
                          as.character(r$family[["family"]])[1]
                        }, "")), collapse = " ")), 1, 200),
                      error = function(e) "")),
            file = out_file, append = TRUE)
      }, error = function(e) NULL)
      fit
    }
    environment(w_fa) <- environment()
    unlockBinding("fit_assembled", ns)
    assign("fit_assembled", w_fa, envir = ns)
    lockBinding("fit_assembled", ns)
  })
})
