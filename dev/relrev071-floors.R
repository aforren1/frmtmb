# Release reviewer of 0.71.0: which core symbols each extension's R code
# names (frmtmb:: and frmtmb:::, importFrom(frmtmb, ...)) are missing
# from a given core build, and which of them changed signature.
#   Rscript dev/relrev071-floors.R
ext <- c("frmtmb.sample", "frmtmb.eam", "frmtmb.latent", "frmtmb.spline")
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/"
cores <- c(r9 = "C:/Users/adf44/source/r/rellib-r9",
           r8 = "C:/Users/adf44/source/r/rellib-r8",
           r6 = "C:/Users/adf44/source/r/rellib-r6")
ns_syms <- lapply(cores, function(l) {
  r <- callr_like <- NULL
  e <- new.env()
  f <- file.path(l, "frmtmb", "R", "frmtmb")
  lazyLoad(f, envir = e)
  v <- read.dcf(file.path(l, "frmtmb", "DESCRIPTION"), "Version")[1]
  list(v = v, env = e)
})
for (p in ext) {
  files <- list.files(file.path(root, p, "R"), "\\.R$", full.names = TRUE)
  src <- unlist(lapply(files, readLines, warn = FALSE))
  src <- sub("#.*$", "", src)
  m <- unlist(regmatches(src, gregexpr("frmtmb:::?[A-Za-z_.][A-Za-z0-9_.]*",
                                       src)))
  syms <- unique(sub("^frmtmb:::?", "", m))
  nsf <- readLines(file.path(root, p, "NAMESPACE"))
  imp <- regmatches(nsf, regexpr("^importFrom\\(frmtmb,.*\\)", nsf))
  imp <- sub("^importFrom\\(frmtmb, *([^)]*)\\)$", "\\1", imp)
  syms <- unique(c(syms, gsub("\"", "", imp)))
  if (any(grepl("^import\\(frmtmb\\)", nsf))) cat(p, "imports all of frmtmb\n")
  dv <- read.dcf(file.path(root, p, "DESCRIPTION"), c("Version", "Depends",
                                                         "Imports"))
  cat("\n==", p, dv[1, "Version"], "| needs:",
      regmatches(paste(dv[1, 2:3], collapse = " "),
                 regexpr("frmtmb \\([^)]*\\)",
                         paste(dv[1, 2:3], collapse = " "))),
      "| core symbols named:", length(syms), "\n")
  for (k in names(cores)) {
    e <- ns_syms[[k]]$env
    miss <- syms[!vapply(syms, exists, NA, envir = e, inherits = FALSE)]
    cat(sprintf("  missing from %s (frmtmb %s): %s\n", k, ns_syms[[k]]$v,
                if (length(miss)) paste(sort(miss), collapse = ", ") else
                  "none"))
  }
  e9 <- ns_syms$r9$env; e8 <- ns_syms$r8$env
  ch <- syms[vapply(syms, function(s) {
    exists(s, envir = e9, inherits = FALSE) &&
      exists(s, envir = e8, inherits = FALSE) &&
      is.function(get(s, e9)) &&
      !identical(names(formals(get(s, e9))), names(formals(get(s, e8))))
  }, NA)]
  cat("  formals changed 0.70.0 -> 0.71.0:",
      if (length(ch)) paste(ch, collapse = ", ") else "none", "\n")
}
