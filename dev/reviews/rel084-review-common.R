# rel084 release review: shared setup. REVLIB = r23 (default) or r22.
which_lib <- Sys.getenv("REVLIB", "r23")
.libPaths(c(paste0("C:/Users/adf44/source/r/rellib-", which_lib),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")

# every condition a call raises, by class, with the value or the error
conds_of <- function(expr) {
  said <- list()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      said[[length(said) + 1L]] <<- list(k = "W", cls = class(w),
                                         msg = conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      said[[length(said) + 1L]] <<- list(k = "M", cls = class(m),
                                         msg = conditionMessage(m))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      said[[length(said) + 1L]] <<- list(k = "E", cls = class(e),
                                         msg = conditionMessage(e))
      NULL
    })
  list(v = v, said = said)
}
cls_line <- function(r) {
  if (!length(r$said)) return("-")
  paste(vapply(r$said, function(s) {
    paste0(s$k, ":", s$cls[[1]])
  }, ""), collapse = " ")
}
has_cls <- function(r, cl, k = NULL) {
  any(vapply(r$said, function(s) {
    cl %in% s$cls && (is.null(k) || s$k == k)
  }, NA))
}
msg_of <- function(r, cl) {
  for (s in r$said) if (cl %in% s$cls) return(s$msg)
  NA_character_
}
