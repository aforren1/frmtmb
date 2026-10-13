# How many deprecation warnings parnames() gives on a fit and on draws,
# by load order and call route. Each case runs in a fresh Rscript.
# Usage: Rscript dev/moparity-parnames.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args) && nzchar(args[[1]])) args[[1]] else NULL
libs <- c(lib, "C:/Users/adf44/source/r/rellib-r7",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
run_child <- function(code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste0(deparse(libs), collapse = "")),
               code), f)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  c("--vanilla", shQuote(f)),
                                  stdout = TRUE, stderr = TRUE))
  paste(out, collapse = "\n")
}
setup_fit <- c(
  "set.seed(1); d <- data.frame(x = rnorm(40)); d$y <- rnorm(40, d$x)",
  "fit <- frmtmb::frm(frmtmb::bf(y ~ x), family = gaussian(), data = d)")
setup_draws <- c(setup_fit,
  "tpl <- fit$frame$par_template",
  "est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))",
  "M <- matrix(rep(est, each = 3), 3,",
  "  dimnames = list(NULL, frmtmb::brms_par_labels(fit)))",
  "obj <- structure(list(stanfit = NULL, draws = cbind(",
  "  frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0), fit = fit),",
  "  class = 'frmtmb_draws')")
count <- function(call) c(
  "n <- 0L",
  sprintf("r <- tryCatch(withCallingHandlers(%s, warning = function(w) {",
          call),
  "  n <<- n + 1L; invokeRestart('muffleWarning') }),",
  "  error = function(e) { cat('ERR:', conditionMessage(e), '\\n'); NULL })",
  "cat('NWARN:', n, '\\n')")
orders <- list(
  frmtmb_only = "suppressMessages(library(frmtmb))",
  brms_then_frmtmb = c("suppressMessages(library(brms))",
                       "suppressMessages(library(frmtmb))"),
  frmtmb_then_brms = c("suppressMessages(library(frmtmb))",
                       "suppressMessages(library(brms))"),
  frmtmb_brms_ns = c("suppressMessages(library(frmtmb))",
                     "invisible(loadNamespace('brms'))"))
routes <- c(bare = "parnames(obj)", frmtmb = "frmtmb::parnames(obj)",
            brms = "brms::parnames(obj)",
            saved = "pn(obj)", unloaded = "pb(obj)")
for (kind in c("fit", "draws")) {
  for (on in names(orders)) {
    ord <- orders[[on]]
    if (kind == "draws") ord <- c(ord,
                                  "suppressMessages(library(frmtmb.sample))")
    for (rn in names(routes)) {
      if (rn == "brms" && on == "frmtmb_only") next
      pre <- if (rn == "saved") {
        # frmtmb's binding captured before brms loads, called after
        c("suppressMessages(library(frmtmb))", "pn <- frmtmb::parnames",
          "invisible(loadNamespace('brms'))", ord)
      } else ord
      if (rn == "unloaded") {
        # brms's generic captured, brms unloaded, then called
        pre <- c(ord, "pb <- brms::parnames",
                 "invisible(unloadNamespace('brms'))")
      }
      if (rn %in% c("saved", "unloaded") && on != "frmtmb_only") next
      setup <- if (kind == "fit") c(setup_fit, "obj <- fit") else setup_draws
      out <- run_child(c(pre, setup, count(routes[[rn]])))
      nw <- sub(".*NWARN: *([0-9]+).*", "\\1", out)
      if (!grepl("NWARN", out)) nw <- paste("FAILED:", out)
      err <- if (grepl("ERR:", out)) {
        paste0(" (the call errored: ",
               trimws(sub("(?s).*ERR:([^\n]*).*", "\\1", out, perl = TRUE)),
               ")")
      } else ""
      cat(sprintf("%-6s %-18s %-7s NWARN %s%s\n", kind, on, rn, nw, err))
    }
  }
}
