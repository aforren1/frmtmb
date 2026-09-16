## Is there a FOURTH hop? Systematic instead of one generic at a time.
## For every S3 generic frmtmb registers a method on, scan the R base
## packages and the installed model-interop packages for calls that pass
## a NAMED argument the frmtmb method has no formal for and does not
## .allow. Those are the calls the refusal now breaks.
## Usage: Rscript dev/asrev-hopscan.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1],
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

# what each method really accepts: its formals, plus whatever the
# .allow in its own body names literally
allowed_of <- function(fn) {
  fo <- setdiff(names(formals(fn)), "...")
  b <- body(fn)
  extra <- character()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(as.character(e[[1L]])[1L], "frm_check_dots")) {
        a <- as.list(e)[-1L]
        al <- a[[".allow"]]
        if (!is.null(al)) {
          extra <<- c(extra, all.names(al), as.character(unlist(
            lapply(as.list(al), function(z)
              if (is.character(z)) z else character()))))
        }
      }
      for (i in seq_along(e)) if (!is.null(e[[i]])) try(walk(e[[i]]),
                                                        silent = TRUE)
    }
  }
  try(walk(b), silent = TRUE)
  unique(c(fo, extra))
}

meths <- list()
for (pkg in c("frmtmb", "frmtmb.sample")) {
  reg <- parseNamespaceFile(pkg,
                            dirname(system.file(package = pkg)))$S3methods
  for (i in seq_len(nrow(reg))) {
    gen <- reg[i, 1]; fun <- if (is.na(reg[i, 3]))
      paste(reg[i, 1], reg[i, 2], sep = ".") else reg[i, 3]
    o <- tryCatch(get(fun, envir = asNamespace(pkg)),
                  error = function(e) NULL)
    if (!is.function(o) || !"..." %in% names(formals(o))) next
    meths[[gen]] <- unique(c(meths[[gen]], allowed_of(o)))
  }
}
cat("generics with a dots-refusing frmtmb method:", length(meths), "\n\n")

scan_pkgs <- c("stats", "base", "utils", "graphics", "methods",
               "emmeans", "insight", "marginaleffects", "bayesplot",
               "loo", "posterior", "brms", "lme4", "nlme", "mgcv",
               "testthat", "knitr")
scan_pkgs <- Filter(function(p) requireNamespace(p, quietly = TRUE),
                    scan_pkgs)
found <- 0L
for (pk in scan_pkgs) {
  ns <- asNamespace(pk)
  for (nm in ls(ns, all.names = TRUE)) {
    o <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
    if (!is.function(o)) next
    src <- tryCatch(deparse(body(o)), error = function(e) character())
    if (!length(src)) next
    for (gen in names(meths)) {
      # a generic whose name is not an identifier ([[, $, [) cannot be
      # matched this way and is scanned by its literal text instead
      if (!grepl("^[A-Za-z._][A-Za-z._0-9]*$", gen)) next
      pat <- paste0("(^|[^._[:alnum:]])", gsub("\\.", "[.]", gen), "\\(")
      ln <- grep(pat, src, value = TRUE)
      if (!length(ln)) next
      for (l in ln) {
        a <- regmatches(l, gregexpr("[A-Za-z._][A-Za-z._0-9]* *=", l))[[1]]
        a <- trimws(sub("=$", "", trimws(a)))
        extra <- setdiff(a, c(meths[[gen]], "..."))
        # names that are almost certainly arguments of some OTHER call
        # on the same deparsed line are unavoidable noise; report and
        # let the reader judge
        if (length(extra)) {
          found <- found + 1L
          cat(sprintf("%-10s %-28s %-14s would-refuse: %s\n    %s\n",
                      pk, nm, gen, paste(extra, collapse = ", "),
                      trimws(substr(l, 1, 96))))
        }
      }
    }
  }
}
cat("\ncandidate sites:", found, "\n")
