## Insert `frm_check_dots(...)` as the first statement of every
## registered S3 method that takes `...` and never reads it.
##
## Run: Rscript dev/argspell-insert.R <pkgdir>
## The target list comes from the INSTALLED base package, so the script
## reports what it did not find in the sources rather than assuming.
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- if (length(args)) args[1] else "."
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
suppressMessages(library(pkg, character.only = TRUE))

# Generics another package owns and calls itself with arguments of its
# own choosing. Refusing an unknown name there breaks the caller.
# Measured, not assumed: see dev/argspell-exempt2.R.
exempt <- c(
  "emm_basis.frmtmb_fit", "get_predict.frmtmb_fit",
  "get_vcov.frmtmb_fit", "get_varcov.frmtmb_fit",
  "get_parameters.frmtmb_fit", "find_formula.frmtmb_fit",
  "find_random.frmtmb_fit",
  "get_coef.frmtmb_fit", "set_coef.frmtmb_fit")

ns <- asNamespace(pkg)
reg <- parseNamespaceFile(pkg, dirname(system.file(package = pkg)))$S3methods
nm <- ifelse(is.na(reg[, 3]), paste(reg[, 1], reg[, 2], sep = "."), reg[, 3])
nm <- setdiff(unique(nm), exempt)
targets <- Filter(function(f) {
  o <- tryCatch(get(f, envir = ns), error = function(e) NULL)
  if (is.null(o) || !is.function(o)) return(FALSE)
  if (!"..." %in% names(formals(o))) return(FALSE)
  !any(grepl("...", deparse(body(o)), fixed = TRUE))
}, nm)
cat("targets:", length(targets), "\n")

files <- list.files(file.path(pkgdir, "R"), pattern = "[.]R$",
                    full.names = TRUE)
found <- character()
noblock <- character()
for (fp in files) {
  txt <- readLines(fp, warn = FALSE)
  ex <- parse(text = paste(txt, collapse = "\n"), keep.source = TRUE)
  ins <- list()
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !identical(as.character(e[[1L]]), "<-")) next
    if (!is.name(e[[2L]])) next
    fn <- as.character(e[[2L]])
    if (!fn %in% targets) next
    rhs <- e[[3L]]
    if (!is.call(rhs) || !identical(as.character(rhs[[1L]]), "function")) next
    b <- rhs[[3L]]
    # A `{` call carries one srcref per element and element ONE is the
    # brace itself, so the first statement is srl[[2]]. Taking srl[[1]]
    # put the insert one line above the brace, which lands inside a
    # multi-line argument list; every such file then failed to parse.
    srl <- attr(b, "srcref")
    ok <- is.call(b) && identical(as.character(b[[1L]]), "{") &&
      !is.null(srl) && length(srl) >= 2L
    if (ok) {
      s1 <- as.integer(srl[[2L]])
      # a one-line body puts the first statement on the brace's own
      # line, where there is no line to insert before
      ok <- s1[1L] > as.integer(srl[[1L]])[1L]
    }
    if (!ok) {
      noblock <- c(noblock, fn)
      next
    }
    found <- c(found, fn)
    ins[[length(ins) + 1L]] <- list(
      line = s1[1L] - 1L,
      text = paste0(strrep(" ", s1[5L] - 1L), "frm_check_dots(...)"))
  }
  if (!length(ins)) next
  for (k in order(vapply(ins, function(z) z[["line"]], 0L),
                  decreasing = TRUE)) {
    at <- ins[[k]]$line
    txt <- append(txt, ins[[k]]$text, after = at)
  }
  con <- file(fp, open = "wb")
  writeLines(txt, con, sep = "\n")
  close(con)
  cat("patched", basename(fp), length(ins), "\n")
}
cat("\ninserted into", length(found), "of", length(targets), "targets\n")
miss <- setdiff(targets, found)
if (length(miss)) cat("NOT FOUND in sources:\n  ",
                      paste(miss, collapse = "\n  "), "\n")
if (length(noblock)) cat("NO BRACED BODY (do by hand):\n  ",
                         paste(noblock, collapse = "\n  "), "\n")
