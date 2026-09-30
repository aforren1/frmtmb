# Reviewer of lane defects, recheck: mutants of the punch-round fixes,
# applied in the loaded namespace (no install), one test file each.
#   Rscript dev/defects-rev-mutant2.R <mutant> <pkg> <test file>
a <- commandArgs(trailingOnly = TRUE)
mut <- a[1]; p <- a[2]; tf <- a[3]
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p,
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(testthat); library(frmtmb)})
mutants <- list(
  M9 = list("frmtmb", "frame_observed_y", "y[miss] <- NA", "NULL"),
  M10 = list("frmtmb.sample", "draws_observed_y", "y[miss] <- NA", "NULL"),
  M11 = list("frmtmb", "fitted.frmtmb_fit",
             "paste0(\"eta\", seq_len(dim(out)[3L]))", "NULL"),
  M12 = list("frmtmb", "assemble_frame",
             "(check_response || identical(resp$family[[\"type\"]], \"ordinal\"))",
             "check_response")
)
m <- mutants[[mut]]
suppressMessages(loadNamespace(m[[1]]))
ns <- asNamespace(m[[1]])
f <- get(m[[2]], envir = ns)
src <- deparse(f, width.cutoff = 500L)
stopifnot("mutation site not found" = any(grepl(m[[3]], src, fixed = TRUE)))
g <- eval(parse(text = gsub(m[[3]], m[[4]], src, fixed = TRUE)), envir = ns)
environment(g) <- ns
utils::assignInNamespace(m[[2]], g, ns = m[[1]])
if (grepl("[.]frmtmb_fit$", m[[2]])) {
  registerS3method(sub("[.]frmtmb_fit$", "", m[[2]]), "frmtmb_fit", g,
                   envir = ns)
}
suppressMessages(library(p, character.only = TRUE))
res <- as.data.frame(test_file(tf, package = p, env = test_env(p),
                               reporter = "silent", stop_on_failure = FALSE))
cat("MUTANT", mut, m[[2]], basename(tf), "pass=", sum(res$passed),
    "fail=", sum(res$failed), "err=", sum(res$error), "\n")
bad <- res[res$failed > 0 | res$error, "test"]
if (length(bad)) cat("  caught by:", paste(bad, collapse = " | "), "\n")
