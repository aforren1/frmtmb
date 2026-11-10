# Lane code8 shared setup. CODE8_ARM = "before" (rellib-r19 only) or
# "after" (the lane library first).
code8_arm <- Sys.getenv("CODE8_ARM", "after")
.libPaths(c(
  if (code8_arm == "after") "C:/Users/adf44/source/r/wt-code8-lib",
  "C:/Users/adf44/source/r/rellib-r19",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", code8_arm, "frmtmb", format(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")
fns <- asNamespace("frmtmb")

quiet <- function(expr) suppressMessages(suppressWarnings(expr))
ll_of <- function(f) as.numeric(logLik(f))

# fit with every warning and message caught, with its classes
rfit <- function(...) {
  w <- character(0); wc <- character(0)
  m <- character(0); mc <- character(0)
  cls <- function(x) paste(setdiff(class(x), c("warning", "message",
                                               "condition", "simpleWarning",
                                               "simpleMessage")),
                           collapse = "/")
  f <- withCallingHandlers(
    tryCatch(frm(...), error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); wc <<- c(wc, cls(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x)); mc <<- c(mc, cls(x))
      invokeRestart("muffleMessage")
    })
  list(f = f, w = w, wc = wc, m = m, mc = mc, err = inherits(f, "error"))
}

# What the stop is, measured after the fit with the package's own
# instruments: the Newton gain of conv_headroom() and whether every
# standard error of the fixed parameters is finite.
stop_info <- function(r) {
  if (r$err) {
    return(data.frame(code = NA, msg = substr(conditionMessage(r$f), 1, 60),
                      gain = NA, se_fin = NA, n_lost = NA))
  }
  f <- r$f
  g <- tryCatch(drop(f$obj$gr(f$opt$par)), error = function(e) NULL)
  room <- tryCatch(fns$conv_headroom(f, g), error = function(e) NA_real_)
  sdr <- tryCatch(quiet(fns$sdr_of(f)), error = function(e) NULL)
  se <- if (is.null(sdr)) NA else sqrt(diag(as.matrix(sdr$cov.fixed)))
  data.frame(code = f$opt$convergence,
             msg = substr(f$opt$message %||% "", 1, 30),
             gain = room,
             se_fin = !anyNA(se) && all(is.finite(se)),
             n_lost = length(sdr$se_lost))
}

cond_txt <- function(r) {
  a <- if (length(r$w)) paste0("W[", r$wc, "] ", substr(r$w, 1, 60))
  b <- if (length(r$m)) paste0("M[", r$mc, "] ", substr(r$m, 1, 60))
  x <- gsub("[\t\n]", " ", c(a, b))
  if (!length(x)) "none" else paste(unique(x), collapse = " || ")
}

`%||%` <- function(a, b) if (is.null(a)) b else a
