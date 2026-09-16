## The same positional audit, over EVERY `frmtmb_draws` method this
## package registers, not only the 28 OWNED generic names the review
## covered. It says whether the ten were the whole of the defect.
##
## The method table of frmtmb.sample's own namespace is NOT the place to
## look: it holds only the 28 methods on names this package also defines
## a generic for. `parseNamespaceFile()` lists every registration.
##
##   Rscript dev/brmsmatch-beyond.R
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
for (p in c("brms", "bayesplot", "bridgesampling", "coda", "loo",
            "posterior", "rstantools", "graphics", "stats", "nlme",
            "lme4", "gratia")) {
  q(requireNamespace(p, quietly = TRUE))
}

ns <- parseNamespaceFile("frmtmb.sample",
                         dirname(find.package("frmtmb.sample")))
m <- ns$S3methods
gen <- unique(m[m[, 2] == "frmtmb_draws", 1])

cut_dots <- function(f) {
  a <- names(formals(f))
  a[seq_len(match("...", a, nomatch = length(a) + 1L) - 1L)]
}
pkgs <- c("brms", "bayesplot", "bridgesampling", "coda", "loo",
          "posterior", "rstantools", "graphics", "stats", "nlme",
          "lme4", "gratia", "frmtmb", "frmtmb.sample")
find_method <- function(nm, cls) {
  for (p in pkgs) {
    t <- tryCatch(get(".__S3MethodsTable__.", envir = asNamespace(p),
                      inherits = FALSE), error = function(e) NULL)
    if (is.null(t)) next
    f <- tryCatch(get(paste0(nm, ".", cls), envir = t,
                      inherits = FALSE), error = function(e) NULL)
    if (is.function(f)) return(f)
  }
  for (p in pkgs) {
    f <- tryCatch(get(paste0(nm, ".", cls), envir = asNamespace(p),
                      inherits = FALSE), error = function(e) NULL)
    if (is.function(f)) return(f)
  }
  NULL
}

cat("frmtmb_draws generics this package registers a method on: ",
    length(gen), "\n\n")
cat(sprintf("%-22s %-4s %-18s %-18s\n", "method", "pos", "brms's",
            "this package's"))
bad <- 0L; good <- 0L; nocmp <- 0L; badnames <- character()
nocmpnames <- character()
for (nm in sort(gen)) {
  b <- find_method(nm, "brmsfit")
  f <- find_method(nm, "frmtmb_draws")
  if (!is.function(b) || !is.function(f)) {
    nocmp <- nocmp + 1L
    nocmpnames <- c(nocmpnames, nm)
    next
  }
  bn <- cut_dots(b)
  fn <- cut_dots(f)
  kk <- min(length(bn), length(fn))
  d <- if (kk) which(bn[seq_len(kk)] != fn[seq_len(kk)]) else integer()
  if (!length(d)) { good <- good + 1L; next }
  bad <- bad + 1L
  badnames <- c(badnames, nm)
  cat(sprintf("%-22s %-4d %-18s %-18s\n", nm, d[[1L]], bn[[d[[1L]]]],
              fn[[d[[1L]]]]))
}
cat("\ndiverge: ", bad, "   agree as far as both go: ", good,
    "   no brmsfit method to compare: ", nocmp, "\n")
cat("diverging: ", paste(badnames, collapse = " "), "\n")
cat("not comparable: ", paste(nocmpnames, collapse = " "), "\n")

## What "agree as far as both go" CANNOT see.
##
## The criterion truncates each side at ITS OWN `...`, so a method whose
## arguments run out before brms's scores as agreeing however many
## positional slots brms still has. Those slots are not empty: a
## positional call lands in our `...` and is ignored, which is an answer
## to a different question with nothing said. The count below is the
## size of that blind spot, and it is why six is the count the criterion
## can see and not the remaining count.
cat("\n== the blind spot: our arguments run out before brms's ==\n")
cat(sprintf("%-22s %-5s %-6s %s\n", "method", "ours", "brms's",
            "brms slots past our end"))
short <- 0L
shortnames <- character()
for (nm in sort(gen)) {
  b <- find_method(nm, "brmsfit")
  f <- find_method(nm, "frmtmb_draws")
  if (!is.function(b) || !is.function(f)) next
  bn <- cut_dots(b)
  fn <- cut_dots(f)
  if (length(bn) <= length(fn)) next
  short <- short + 1L
  shortnames <- c(shortnames, nm)
  extra <- bn[(length(fn) + 1L):length(bn)]
  cat(sprintf("%-22s %-5d %-6d %s\n", nm, length(fn), length(bn),
              paste(utils::head(extra, 5L), collapse = " ")))
}
cat("\nscored 'agree' but shorter than brms's: ", short, "\n")
cat("names: ", paste(shortnames, collapse = " "), "\n")
cat("DONE\n")
