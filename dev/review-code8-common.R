# Review of lane code8: shared setup. RV_ARM = "before" (rellib-r19) or
# "after" (the lane library first).
rv_arm <- Sys.getenv("RV_ARM", "after")
.libPaths(c(
  if (rv_arm == "after") "C:/Users/adf44/source/r/wt-code8-lib",
  "C:/Users/adf44/source/r/rellib-r19",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", rv_arm, format(packageVersion("frmtmb")), find.package("frmtmb"),
    "ritz:", exists("se_ridge_ritz", ns, inherits = FALSE),
    "tol:", if (exists("conv_at_max_tol", ns)) ns$conv_at_max_tol, "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
rfit <- function(...) {
  w <- character(0); wc <- character(0)
  m <- character(0); mc <- character(0)
  cls <- function(x) paste(setdiff(class(x), c("warning", "message",
    "condition", "simpleWarning", "simpleMessage", "frmtmb_warning",
    "frmtmb_message")), collapse = "/")
  t0 <- proc.time()[["elapsed"]]
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
  list(f = f, w = w, wc = wc, m = m, mc = mc, err = inherits(f, "error"),
       time = proc.time()[["elapsed"]] - t0)
}
cond_txt <- function(r) {
  a <- if (length(r$w)) paste0("W[", r$wc, "] ", substr(r$w, 1, 50))
  b <- if (length(r$m)) paste0("M[", r$mc, "] ", substr(r$m, 1, 50))
  x <- gsub("[\t\n]", " ", c(a, b))
  if (!length(x)) "none" else paste(unique(x), collapse = " || ")
}
se_of <- function(f) {
  sdr <- quiet(ns$sdr_of(f))
  se <- sqrt(diag(as.matrix(sdr$cov.fixed)))
  names(se) <- ns$outer_par_names(f)
  list(se = se, lost = sdr$se_lost)
}
dir.create("dev/review-code8-out", showWarnings = FALSE)
