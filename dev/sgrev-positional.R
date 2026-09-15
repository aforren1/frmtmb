# Item 4, second half, audited in full rather than by example. For
# every one of the 28 names, line up frmtmb.sample's frmtmb_draws
# METHOD against brms's brmsfit METHOD and report the FIRST argument
# position at which a positional caller would be answering a
# different question. brms is the tiebreaker, so brms's method is the
# standard, not the owner's generic.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
for (p in c("brms", "bayesplot", "bridgesampling", "coda", "gratia",
            "loo", "posterior", "rstantools", "frmtmb.sample"))
  q(requireNamespace(p, quietly = TRUE))
GEN <- readRDS("dev/sgrev-out/names.rds")
tbl <- function(p) tryCatch(get(".__S3MethodsTable__.",
  envir = asNamespace(p), inherits = FALSE), error = function(e) NULL)
ft <- tbl("frmtmb.sample")
find_brms <- function(nm) {
  for (p in c("brms", "bayesplot", "bridgesampling", "coda", "loo",
              "posterior", "rstantools")) {
    t <- tbl(p); if (is.null(t)) next
    m <- tryCatch(get(paste0(nm, ".brmsfit"), envir = t,
                      inherits = FALSE), error = function(e) NULL)
    if (is.function(m)) return(m)
  }
  NULL
}
cat("== positional-argument agreement with brms, all 28 ==\n")
cat("position 1 is the object. A row is listed when the first",
    "disagreeing\nposition is one a caller can actually pass",
    "positionally, that is, before\nbrms's own `...`.\n\n")
cat(sprintf("%-20s %-4s %-18s %-18s\n", "generic", "pos",
            "brms's method", "frmtmb.sample's"))
bad <- 0; good <- 0; nocmp <- 0
for (nm in GEN$names) {
  b <- find_brms(nm)
  f <- tryCatch(get(paste0(nm, ".frmtmb_draws"), envir = ft,
                    inherits = FALSE), error = function(e) NULL)
  if (!is.function(b) || !is.function(f)) { nocmp <- nocmp + 1; next }
  bn <- names(formals(b)); fn <- names(formals(f))
  bn <- bn[seq_len(match("...", bn, nomatch = length(bn) + 1L) - 1L)]
  fn <- fn[seq_len(match("...", fn, nomatch = length(fn) + 1L) - 1L)]
  k <- min(length(bn), length(fn))
  pos <- NA
  if (k > 0) {
    d <- which(bn[seq_len(k)] != fn[seq_len(k)])
    if (length(d)) pos <- d[1]
  }
  if (is.na(pos)) { good <- good + 1; next }
  bad <- bad + 1
  cat(sprintf("%-20s %-4d %-18s %-18s\n", nm, pos, bn[pos], fn[pos]))
}
cat("\ndisagree at some positional slot: ", bad,
    "   agree as far as both go: ", good,
    "   no brms method to compare: ", nocmp, "\n")

cat("\n== the one that is a silent wrong answer, not an error ==\n")
cat("brms's idiom posterior_epred(fit, newdata, NA) means",
    "re_formula = NA,\nwhich drops the random effects. Here position",
    "3 is `resp`.\n")
q(library(frmtmb)); q(library(frmtmb.sample))
ds <- readRDS("dev/stan-cache/sgrev-draws.rds")
nd <- ds$fit$frame$data[1:5, , drop = FALSE]
set.seed(1)
a <- tryCatch(posterior_epred(ds, nd, NA, ndraws = 20),
              error = function(e) conditionMessage(e))
set.seed(1)
b <- tryCatch(posterior_epred(ds, nd, re_formula = NA, ndraws = 20),
              error = function(e) conditionMessage(e))
set.seed(1)
cc <- tryCatch(posterior_epred(ds, nd, ndraws = 20),
               error = function(e) conditionMessage(e))
show <- function(lbl, x) cat(sprintf("  %-42s %s\n", lbl,
  if (is.character(x)) paste("ERROR:", substr(x, 1, 50)) else
    paste0("[", paste(format(x[1, ], digits = 6), collapse = " "), "]")))
show("posterior_epred(ds, nd, NA)  (positional)", a)
show("posterior_epred(ds, nd, re_formula = NA)", b)
show("posterior_epred(ds, nd)  (with the REs)", cc)
if (!is.character(a) && !is.character(cc))
  cat("  positional NA gives the SAME answer as no third argument: ",
      identical(a, cc), "\n")
if (!is.character(b) && !is.character(cc))
  cat("  re_formula = NA differs from the default, as it should: ",
      !identical(b, cc), "\n")
cat("DONE\n")
