# Lane optstart: an R_PROFILE_USER that logs every fit a test file makes
# (dev/optrap-profile.R, with the warnings each fit raised), so that two
# builds' runs of the suites can be compared fit by fit
# (dev/optstart-fitdiff.R). A fit is one call of fit_assembled(); its
# key is a digest of the response, the call, REML and quadrature,
# numbered by occurrence within the file. The warnings are those raised
# inside that call, recorded and passed on unchanged. Lines go to
# $OPTSTART_PROF_LOG/<pid>.txt; logging failures are swallowed so the
# fit is unchanged.
local({
  logdir <- Sys.getenv("OPTSTART_PROF_LOG")
  if (!nzchar(logdir)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    out_file <- file.path(logdir, paste0(Sys.getpid(), ".txt"))
    a <- commandArgs(trailingOnly = TRUE)
    tag <- if (length(a) >= 3) paste0(a[2], "/", basename(a[3])) else "?"
    o_fa <- get("fit_assembled", envir = ns)
    clean <- function(x, k) substr(gsub("[\r\n\t]+", " ", x), 1, k)
    w_fa <- function(spec, frame, bform, cl, ...) {
      key <- tryCatch(digest::digest(list(frame[["y"]], deparse(cl),
                                          list(...)[["REML"]],
                                          isTRUE(list(...)[["quadrature"]]))),
                      error = function(e) "nokey")
      ws <- character(0)
      fit <- tryCatch(withCallingHandlers(
        o_fa(spec, frame, bform, cl, ...),
        warning = function(w) {
          ws <<- c(ws, clean(conditionMessage(w), 90))
        }), error = function(e) {
        tryCatch(cat(sprintf("FIT\t%s\t%s\tERROR\t%s\n", tag, key,
                             clean(conditionMessage(e), 150)),
                     file = out_file, append = TRUE),
                 error = function(e2) NULL)
        stop(e)
      })
      tryCatch({
        op <- fit$opt
        unfitted <- is.null(op)
        cat(sprintf("FIT\t%s\t%s\tOK\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n", tag, key,
                    if (unfitted) "NA" else
                      format(op$objective, digits = 15),
                    if (unfitted) "NA" else op$convergence,
                    if (unfitted) 0L else length(op$par),
                    if (unfitted) "unfitted" else
                      clean(op$message %||% "", 40),
                    if (length(ws)) paste(ws, collapse = " || ") else "-",
                    paste(names(fit$opt$par), collapse = ","),
                    tryCatch(clean(paste(
                      c(deparse(bform[["formula"]] %||% cl[["formula"]]),
                        vapply(spec$responses, function(r) {
                          as.character(r$family[["family"]])[1]
                        }, "")), collapse = " "), 200),
                      error = function(e) "")),
            file = out_file, append = TRUE)
      }, error = function(e) NULL)
      # Item R1: where the standard errors come from the finite-difference
      # Hessian in natural units (par_units other than 1), that Hessian
      # (fit_outer_hessian(), the build's own route) against a Richardson
      # Jacobian of the AD gradient. HESS lines: npar, the largest rise
      # tau = 0.5 * H_ii * (1e-3 u_i)^2 of the route's first pass, the
      # largest relative SE error, and whether each matrix is positive
      # definite. Only fits with at most 40 outer parameters whose
      # gradient takes under 0.2 s; the objective's state is restored.
      tryCatch({
        u <- fit$par_units
        p <- fit$opt$par
        if (!is.null(u) && any(u != 1) && length(p) <= 40 &&
              !inherits(fit, "frmtmb_unfitted")) {
          obj <- fit$obj
          saved <- get("obj_state_save", ns)(obj)
          on.exit(get("obj_state_restore", ns)(obj, saved), add = TRUE)
          t1 <- system.time(g0 <- obj$gr(p))[["elapsed"]]
          if (t1 < 0.2 && all(is.finite(g0))) {
            Hf <- get("fit_outer_hessian", ns)(fit)$H
            Hr <- numDeriv::jacobian(function(x) as.numeric(obj$gr(x)), p)
            Hr <- (Hr + t(Hr)) / 2
            pd <- function(H) {
              all(is.finite(H)) && !inherits(try(chol(H), silent = TRUE),
                                              "try-error")
            }
            se <- function(H) {
              tryCatch(sqrt(diag(solve(H))), error = function(e) NA,
                       warning = function(e) NA)
            }
            sf <- se(Hf)
            sr <- se(Hr)
            ok <- is.finite(sf) & is.finite(sr) & sr > 0
            err <- if (any(ok)) max(abs(sf[ok] / sr[ok] - 1)) else NA
            tau <- max(0.5 * diag(Hf) * (1e-3 * u)^2, na.rm = TRUE)
            cat(sprintf("HESS\t%s\t%s\t%d\t%.3g\t%.3g\t%s\t%s\t%d\n", tag, key,
                        length(p), tau, err, pd(Hf), pd(Hr), sum(!ok)),
                file = out_file, append = TRUE)
          }
        }
      }, error = function(e) NULL)
      fit
    }
    environment(w_fa) <- environment()
    unlockBinding("fit_assembled", ns)
    assign("fit_assembled", w_fa, envir = ns)
    lockBinding("fit_assembled", ns)
  })
})
