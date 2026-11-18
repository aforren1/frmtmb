# Parse core R/, its tests, the extensions and inst/ once, for the
# codesize scripts (dev/codesize-*.R). Writes dev/codesize-out/parsed.rds.
#   Rscript dev/codesize-parse.R [tree]
# `tree` defaults to the lane worktree. Read-only on the tree.
a <- commandArgs(trailingOnly = TRUE)
wt <- if (length(a)) a[1] else "C:/Users/adf44/source/r/frmtmb-wt-codesize"
out <- file.path("C:/Users/adf44/source/r/frmtmb-wt-codesize",
                 "dev/codesize-out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Subsystems of the lane brief, plus the rest of core in named groups.
sub_map <- c(
  families = "families families-cdf families-extreme families-simplex links
    links-brms thres structure",
  `fit-end checks and SE check` = "se-check se-reparam fit-end
    fit-open-ends fit-restart weights-unbounded ald-laplace mixture-starts",
  priors = "priors",
  `prediction and draws` = "predict predict-brms simulate-new
    simulate-newdata draws-generics integrated-design sampling-api",
  `compat and interop` = "compat interop insight brms-vocabulary
    brms-shapes brms-names brms-utilities brms-postfit sugar generic-owners
    extension-home loo-generics",
  `parse and frame` = "parse frame bf re-formula subset refit-data",
  `conditional effects` = "conditional-effects ce-levels
    conditional-smooths",
  `autocor/covstruct` = "autocor autocor-latent autocor-seam autocor-loglik
    covstruct covstruct-seam gr-by",
  `fit engine and starts` = "fit objective par-template autoscale start-ls
    gp-start ad-env quadrature-nested allfit fit-end",
  `inference extras` = "confint importance sandwich bootstrap influence
    multiple joint-sparse extra-deriv",
  `curve inference` = "curve curve-deriv curve-feature curve-cov",
  `model terms (me, gp, ps, spectral)` = "me ps gp-by spectral",
  `post-fit methods` = "methods-fit diagnostics loo scales",
  `infrastructure` = "utils conditions zzz registries hazard-containers",
  `docs-only` = "data frmtmb-package srr-stats-standards")
sub_of <- unlist(lapply(names(sub_map), function(s) {
  f <- strsplit(trimws(gsub("\\s+", " ", sub_map[[s]])), " ")[[1]]
  setNames(rep(s, length(f)), f)
}))
# fit-end.R is listed in both; the checks win
sub_of <- sub_of[!duplicated(names(sub_of))]

parse_one <- function(f) {
  p <- tryCatch(parse(f, keep.source = TRUE), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  pd <- utils::getParseData(p, includeText = TRUE)
  pd <- pd[pd$terminal, c("line1", "col1", "line2", "token", "text")]
  pd <- pd[order(pd$line1, pd$col1), ]
  sr <- attr(p, "srcref")
  top <- data.frame(
    first = vapply(sr, function(s) s[1], 1L),
    last = vapply(sr, function(s) s[3], 1L))
  top$name <- NA_character_
  top$is_fun <- FALSE
  for (i in seq_along(p)) {
    e <- p[[i]]
    if (is.call(e) && length(e) == 3 && is.name(e[[1]]) &&
        as.character(e[[1]]) %in% c("<-", "=", "<<-")) {
      lhs <- e[[2]]
      nm <- if (is.name(lhs) || is.character(lhs)) as.character(lhs) else NA
      top$name[i] <- nm
      rhs <- e[[3]]
      top$is_fun[i] <- is.call(rhs) && identical(rhs[[1]], as.name("function"))
    }
  }
  list(lines = readLines(f, warn = FALSE), pd = pd, top = top)
}

rf <- sort(list.files(file.path(wt, "R"), "[.]R$", full.names = TRUE))
core <- setNames(lapply(rf, parse_one), sub("[.]R$", "", basename(rf)))
stopifnot(!vapply(core, is.null, TRUE))
ext_files <- c(
  list.files(file.path(wt, "tests"), "[.]R$", recursive = TRUE,
             full.names = TRUE),
  list.files(file.path(wt, "extensions"), "[.]R$", recursive = TRUE,
             full.names = TRUE),
  list.files(file.path(wt, "inst"), "[.]R$", recursive = TRUE,
             full.names = TRUE))
ext_files <- ext_files[!grepl("/(revdep|docs)/", ext_files)]
ext <- lapply(ext_files, function(f) {
  x <- parse_one(f)
  if (is.null(x)) list(lines = readLines(f, warn = FALSE), pd = NULL) else x
})
names(ext) <- sub(paste0("^", wt, "/"), "", ext_files)
ns <- parse(file.path(wt, "NAMESPACE"))
exports <- unique(unlist(lapply(ns, function(e)
  if (identical(e[[1]], as.name("export"))) vapply(as.list(e)[-1],
                                               as.character, "") )))
s3 <- do.call(rbind, lapply(ns, function(e) {
  if (!identical(e[[1]], as.name("S3method"))) return(NULL)
  x <- vapply(as.list(e)[-1], function(z)
    if (is.character(z)) z else paste(deparse(z), collapse = ""), "")
  x <- sub("^.*::", "", x)
  data.frame(generic = x[1], class = x[2],
             method = if (length(x) > 2) x[3] else paste0(x[1], ".", x[2]))
}))
miss <- setdiff(names(core), names(sub_of))
if (length(miss)) stop("files with no subsystem: ", paste(miss, collapse = " "))
saveRDS(list(core = core, ext = ext, exports = exports, s3 = s3,
             sub_of = sub_of, wt = wt), file.path(out, "parsed.rds"))
cat("parsed", length(core), "core files,", length(ext), "other files,",
    length(exports), "exports,", nrow(s3), "S3 registrations\n")
cat("other files that did not parse:",
    sum(vapply(ext, function(x) is.null(x$pd), TRUE)), "\n")
