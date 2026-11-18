# Q6 of the codesize lane: measure the candidates to move out of core
# under dev/core-boundary.md guardrail 3: lines in R/, the core code that
# calls in, exports, the test files that name the candidate's functions
# and their time under covr, and how many release commits changed it.
#   Rscript dev/codesize-candidates.R
# Reads dev/codesize-out/parsed.rds, functions.tsv and, when present,
# dev/codesize-cov/results.tsv (per-test-file seconds under covr: an order
# of magnitude only).
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
funs <- read.delim(file.path(out, "functions.tsv"), stringsAsFactors = FALSE)
funs <- funs[!is.na(funs$name), ]
rt <- file.path(dirname(out), "codesize-cov/results.tsv")
tt <- if (file.exists(rt)) read.delim(rt, stringsAsFactors = FALSE) else NULL

cands <- list(
  `spectral (frequency-domain models)` = list(files = "spectral"),
  `importance sampling of the Laplace approximation` =
    list(files = "importance"),
  `multiple imputation (frm_multiple, Rubin's rules)` =
    list(files = "multiple"),
  `downstream glue (emmeans, marginaleffects, insight)` =
    list(files = c("interop", "insight")),
  `case-deletion influence, allFit` = list(files = c("influence", "allfit")),
  `extreme-value families (gev, frechet, discrete_weibull)` =
    list(files = "families-extreme"),
  `curve inference` = list(files = c("curve", "curve-deriv", "curve-feature",
                                     "curve-cov")),
  `nested quadrature` = list(files = "quadrature-nested"),
  `cluster sandwich` = list(files = "sandwich"),
  `brms post-fit surface (standata, expose_functions, ...)` =
    list(files = c("brms-postfit", "brms-utilities")),
  `compat table` = list(files = "compat"))

ref_sym <- function(pd) {
  pd <- pd[pd$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL", "SPECIAL"), ]
  pd
}
core_tok <- lapply(P$core, function(x) ref_sym(x$pd))
ext_tok <- lapply(P$ext, function(x) if (!is.null(x$pd)) ref_sym(x$pd))
since <- c("0.79.0" = "b9d9d8f0", "0.78.0" = "d5d0056d")
cat("BEGIN codesize-candidates\n")
for (nm in names(cands)) {
  fl <- cands[[nm]]$files
  fx <- funs[sub("[.]R$", "", funs$file) %in% fl, ]
  own <- fx$name
  lines <- sum(vapply(fl, function(f) length(P$core[[f]]$lines), 1L))
  # callers in other core files: functions whose body names one of ours
  callers <- character()
  for (f in setdiff(names(P$core), fl)) {
    x <- P$core[[f]]; s <- core_tok[[f]]
    hit <- s$line1[s$text %in% own]
    if (!length(hit)) next
    t <- x$top
    who <- vapply(hit, function(l) {
      i <- which(t$first <= l & t$last >= l)[1]
      if (is.na(i) || is.na(t$name[i])) paste0(f, ".R:(top level)") else
        paste0(t$name[i], " (", f, ".R)")
    }, "")
    callers <- c(callers, who)
  }
  callers <- sort(unique(callers))
  tf <- names(ext_tok)[vapply(ext_tok, function(s)
    !is.null(s) && any(s$text %in% own), TRUE)]
  tcore <- basename(tf[grepl("^tests/testthat/test-", tf)])
  text <- unique(sub("^extensions/([^/]+)/.*$", "\\1",
                     tf[grepl("^extensions/", tf)]))
  ttime <- if (!is.null(tt)) sum(tt$sec[tt$file %in% tcore], na.rm = TRUE)
  ch <- vapply(since, function(c0) length(system2("git", c(
    "-C", P$wt, "log", "--oneline", paste0(c0, "..c6383527"), "--",
    paste0("R/", fl, ".R")), stdout = TRUE)), 1L)
  cat(sprintf("\n## %s\nfiles: %s; lines %d, code %d, functions %d, exported %d\n",
              nm, paste0(fl, ".R", collapse = " "), lines, sum(fx$code),
              nrow(fx), sum(fx$exported)))
  cat(sprintf("commits changing it since 0.78.0: %d (since 0.79.0: %d)\n",
              ch[["0.78.0"]], ch[["0.79.0"]]))
  cat(sprintf("core test files naming its functions: %d (covr seconds %s); extensions naming them: %s\n",
              length(tcore), if (is.null(ttime)) "NA" else format(ttime),
              if (length(text)) paste(text, collapse = " ") else "none"))
  cat(sprintf("core callers outside it (%d): %s\n", length(callers),
              paste(head(callers, 25), collapse = "; ")))
}
cat("END codesize-candidates\n")
