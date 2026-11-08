# Reviewer's shared setup. REV_ARM = "after" (lane lib first), "before"
# (rellib-r18 only), or a path to a mutant library put first.
rev_arm <- Sys.getenv("REV_ARM", "after")
rev_libs <- c(
  if (rev_arm == "after") "C:/Users/adf44/source/r/wt-bigx-lib"
  else if (rev_arm != "before") rev_arm,
  "C:/Users/adf44/source/r/rellib-r18",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(rev_libs)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", rev_arm, "frmtmb", format(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")
has_basis_code <- exists("opt_basis", envir = asNamespace("frmtmb"))
cat("opt_basis present:", has_basis_code, "\n")

# fit with warnings and messages caught
rfit <- function(...) {
  w <- character(0); wc <- character(0); msg <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(...), error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 160))
      wc <<- c(wc, paste(setdiff(class(x), c("warning", "condition")),
                         collapse = "/"))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msg <<- c(msg, substr(conditionMessage(m), 1, 160))
      invokeRestart("muffleMessage")
    })
  list(f = f, w = w, wc = wc, msg = msg, err = inherits(f, "error"))
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

beta_of <- function(f) f$estimates[["beta"]]
ll_of <- function(f) as.numeric(logLik(f))
se_beta <- function(f) {
  v <- tryCatch(quiet(vcov(f)), error = function(e) NULL)
  if (is.null(v)) return(NA)
  sqrt(pmax(diag(as.matrix(v)), 0))
}
basis_on <- function(f) {
  !is.null(tryCatch(f$obj$env$.frmtmb_opt_basis, error = function(e) NULL))
}
short_w <- function(r) {
  if (!length(r$w)) "none" else paste(unique(paste0("[", r$wc, "] ",
                                                    substr(r$w, 1, 70))),
                                      collapse = " | ")
}
