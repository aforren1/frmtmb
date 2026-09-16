## Reviewer: recount the frmtmb_draws methods and re-run the positional
## criterion, independently of dev/brmsmatch-beyond.R.
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(brms))

ns <- asNamespace("frmtmb.sample")
cat("frmtmb.sample ", format(packageVersion("frmtmb.sample")),
    " from ", dirname(system.file(package = "frmtmb.sample")), "\n")
cat("brms          ", format(packageVersion("brms")), "\n\n")

## --- count 1: the NAMESPACE file, which is the ground truth ----------
nsf <- file.path(system.file(package = "frmtmb.sample"), "NAMESPACE")
ln <- readLines(nsf)
s3 <- grep("^S3method\\(", ln, value = TRUE)
dr <- grep("frmtmb_draws", s3, value = TRUE)
## S3method(gen,cls) or S3method(pkg::gen,cls)
gen_of <- function(x) {
  x <- sub("^S3method\\(", "", sub("\\)\\s*$", "", x))
  p <- strsplit(x, ",")[[1L]]
  sub("^.*::", "", trimws(p[[1L]]))
}
gens <- vapply(dr, gen_of, "")
cat("S3method() lines mentioning frmtmb_draws: ", length(dr), "\n")
cat("distinct generic names:                   ",
    length(unique(gens)), "\n")
dup <- names(which(table(gens) > 1L))
cat("registered under >1 owner:                ",
    paste(sort(dup), collapse = " "), "\n")

## --- count 2: the methods table, as brmsmatch-beyond does ------------
tb <- get(".__S3MethodsTable__.", envir = ns, inherits = FALSE)
own <- grep("\\.frmtmb_draws$", ls(tb), value = TRUE)
cat("in this package's own S3 methods table:   ", length(own), "\n\n")

## --- the positional criterion over every generic ---------------------
## truncate each side at its own `...`, compare name for name.
trunc_at_dots <- function(f) {
  n <- names(formals(f))
  i <- match("...", n)
  if (is.na(i)) n else if (i == 1L) character(0) else n[seq_len(i - 1L)]
}
resolve_ours <- function(g) {
  nm <- paste0(g, ".frmtmb_draws")
  f <- tryCatch(get(nm, envir = tb, inherits = FALSE),
                error = function(e) NULL)
  if (!is.null(f)) return(f)
  ## registered into another package's table
  for (p in loadedNamespaces()) {
    t2 <- tryCatch(get(".__S3MethodsTable__.", envir = asNamespace(p),
                       inherits = FALSE), error = function(e) NULL)
    if (is.null(t2)) next
    f <- tryCatch(get(nm, envir = t2, inherits = FALSE),
                  error = function(e) NULL)
    if (!is.null(f)) return(f)
  }
  tryCatch(getS3method(g, "frmtmb_draws"), error = function(e) NULL)
}
resolve_brms <- function(g) {
  bt <- get(".__S3MethodsTable__.", envir = asNamespace("brms"),
            inherits = FALSE)
  for (cl in c("brmsfit", "brmsfit_multiple")) {
    f <- tryCatch(get(paste0(g, ".", cl), envir = bt, inherits = FALSE),
                  error = function(e) NULL)
    if (!is.null(f)) return(f)
    f <- tryCatch(get(paste0(g, ".", cl), envir = asNamespace("brms")),
                  error = function(e) NULL)
    if (!is.null(f)) return(f)
  }
  NULL
}

rows <- list(); nb <- 0L; ndiv <- 0L; nag <- 0L
for (g in sort(unique(gens))) {
  fo <- resolve_ours(g); fb <- resolve_brms(g)
  if (is.null(fb)) { nb <- nb + 1L; next }
  if (is.null(fo)) { rows[[length(rows) + 1L]] <-
      c(g, "NA", "OURS NOT RESOLVED", ""); next }
  a <- trunc_at_dots(fb); b <- trunc_at_dots(fo)
  k <- min(length(a), length(b))
  d <- which(a[seq_len(k)] != b[seq_len(k)])
  if (length(d)) {
    ndiv <- ndiv + 1L
    rows[[length(rows) + 1L]] <-
      c(g, d[[1L]], a[[d[[1L]]]], b[[d[[1L]]]])
  } else nag <- nag + 1L
}
cat("== positional criterion over every registered generic ==\n")
cat(sprintf("%-22s %-4s %-18s %s\n", "method", "pos", "brms's",
            "this package's"))
for (r in rows) cat(sprintf("%-22s %-4s %-18s %s\n", r[1], r[2], r[3],
                            r[4]))
cat("\ndiverge: ", ndiv, "   agree as far as both go: ", nag,
    "   no brmsfit method: ", nb, "\n\n")

## --- the full signature, both sides, for every comparable generic ----
cat("== full signatures, brms then ours, for every generic ==\n")
for (g in sort(unique(gens))) {
  fo <- resolve_ours(g); fb <- resolve_brms(g)
  if (is.null(fb) || is.null(fo)) {
    cat(sprintf("%-22s  %s\n", g,
                if (is.null(fb)) "(no brmsfit method)" else
                  "(ours not resolved)"))
    next
  }
  cat(sprintf("%-22s B: %s\n", g,
              paste(names(formals(fb)), collapse = " ")))
  cat(sprintf("%-22s O: %s\n", "",
              paste(names(formals(fo)), collapse = " ")))
}
cat("DONE\n")
