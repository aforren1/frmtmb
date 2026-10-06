# Punch round 1, B2: instrument a copy of the lane source so a full-suite
# run records every fit that reaches check_convergence(), every firing
# of the non-finite-gradient warning and of the degenerate-component
# warning, and the reach of every ordinal-mixture fit. One log line
# each, appended to $FRMTMB_P1_LOG/<pid>.txt with the test file the
# process runs (dev/ordmix-runtest.R's second argument). The copy is
# installed into dev/ordmix-p1-check/instr-lib only, and put first on
# the library path by dev/ordmix-runtest.R's "mut:" arm; the behavior
# is the lane build's, plus the appends.
src <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-p1-check/instr-src"
lib <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-p1-check/instr-lib"
patch <- function(file, anchor, add, after = TRUE) {
  p <- file.path(src, "R", file)
  s <- readLines(p)
  i <- which(s == anchor)
  if (length(i) != 1L) stop("anchor not unique in ", file, ": ", anchor)
  s <- if (after) append(s, add, i) else append(s, add, i - 1L)
  writeLines(s, p)
}
writeLines(c(
  "p1_log <- function(kind, ...) {",
  "  dir <- Sys.getenv(\"FRMTMB_P1_LOG\")",
  "  if (!nzchar(dir)) return(invisible())",
  "  a <- commandArgs(TRUE)",
  "  tf <- if (length(a) >= 2L) basename(a[2]) else \"?\"",
  "  # a logging failure must not change what the fit does",
  "  tryCatch(cat(kind, tf, paste(..., collapse = \" \"), \"\\n\",",
  "      file = file.path(dir, paste0(Sys.getpid(), \".txt\")),",
  "      append = TRUE), error = function(e) NULL)",
  "  invisible()",
  "}"), file.path(src, "R", "zz-p1-log.R"))
fams <- "paste(vapply(fit$spec$responses, function(r) paste(as.character(r$family[[\"family\"]] %||% \"?\"), collapse = \"/\"), \"\"), collapse = \"+\")"
patch("fit.R", "  if (inherits(gvec, \"try-error\")) gvec <- NULL",
      paste0("  p1_log(\"FIT\", ", fams, ", \"gradfinite=\", ",
             "is.null(gvec) || all(is.finite(gvec)))"))
patch("fit.R", "    bad <- names(fit$opt$par)[!is.finite(gvec)]",
      paste0("    p1_log(\"GRADNF\", ", fams, ", \"pars=\", ",
             "paste(unique(bad), collapse = \",\"), \"conv=\", ",
             "fit$opt$convergence, \"formula=\", ",
             "gsub(\"[[:space:]]+\", \" \", paste(deparse(",
             "fit$spec$responses[[1]]$formula), collapse = \"\")))"))
patch("families.R", "  bad <- dg$reach > 50 | dg$collapsed",
      paste0("  p1_log(\"MIXFIT\", fit$spec$responses[[resp]]$family",
             "[[\"family\"]], \"reach=\", signif(max(dg$reach), 4), ",
             "\"collapsed=\", any(dg$collapsed), \"fires=\", ",
             "any(dg$reach > 50 | dg$collapsed))"))
patch("families.R", "  frm_warning(fam[[\"family\"]], \": \", paste(what, collapse = \"; \"),",
      paste0("  p1_log(\"DEGEN\", fam[[\"family\"]], paste(what, ",
             "collapse = \"; \"))"), after = FALSE)
rcmd <- file.path(R.home("bin"), "R.exe")
st <- system2(rcmd, c("CMD", "INSTALL", paste0("--library=", lib),
                      "--no-multiarch", "--no-test-load", shQuote(src)))
if (st != 0) stop("install failed")
cat("instrumented build installed in", lib, "\n")
