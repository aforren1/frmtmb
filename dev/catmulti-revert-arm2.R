# Fix round 2 of lane catmulti: test-refcat.R against the ROUND-1
# loglik_df(), which rebuilt the reference model from the fit's model
# frame. The round-1 bodies are put back into the loaded namespace for
# this session only (nothing installed).
#   Rscript dev/catmulti-revert-arm2.R
# Output: dev/catmulti-out/revert-arm2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
put <- function(nm, f) {
  environment(f) <- ns
  unlockBinding(nm, ns)
  assign(nm, f, envir = ns)
  lockBinding(nm, ns)
}
# round 1, verbatim but for the cache, which a fresh fit has empty
put("loglik_df", function(object) {
  n <- n_outer_est(object)
  fams <- lapply(object$spec$responses, `[[`, "family")
  na_ref <- vapply(fams, function(f) is.function(f[["refcat_first"]]), NA)
  if (!any(na_ref)) return(n)
  tryCatch({
    fam <- fams[na_ref][[1L]]
    if (length(object$spec$responses) > 1L) frm_stop("mv", call. = FALSE)
    own <- intersect(names(object$bform[["pforms"]] %||% list()),
                     fam[["primary_dpars"]])
    if (length(own)) frm_stop("only partly free", call. = FALSE)
    bf2 <- object$bform
    bf2$family <- fam[["refcat_first"]]()
    fr2 <- assemble_frame(parse_spec(bf2), object$frame[["data_frame"]],
                          data2 = object$data2 %||% list())
    as.integer(n - (frame_n_outer(object$frame, object) -
                      frame_n_outer(fr2, object)))
  }, frmtmb_error = function(e) {
    structure(NA_integer_, df_why = conditionMessage(e))
  })
})
# round 1's SE sentence came with every lost SE of an NA fit
f <- get("se_lost_message", ns)
src <- paste(deparse(f), collapse = "\n")
old <- "is.function(r$family[[\"refcat_first\"]]) &&"
stopifnot(grepl(old, src, fixed = TRUE))
src <- sub(old, "is.function(r$family[[\"refcat_first\"]]) || ", src,
           fixed = TRUE)
put("se_lost_message", eval(parse(text = src)))
r <- as.data.frame(test_file("tests/testthat/test-refcat.R",
                             package = "frmtmb", env = test_env("frmtmb"),
                             reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("REVERTED test-refcat.R pass=%d fail=%d err=%d\n",
            sum(r$passed), sum(r$failed), sum(r$error)))
for (b in r[r$failed > 0 | r$error, "test"]) cat("  not passing:", b, "\n")
