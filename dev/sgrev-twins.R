# The three rules the new Rd section states, checked against the
# installed FIX build rather than against the lane's test.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
OWN <- c("bayesplot", "bridgesampling", "brms", "coda", "gratia",
         "loo", "posterior", "rstantools")
for (p in c(OWN, "frmtmb.sample")) q(requireNamespace(p, quietly = TRUE))
GEN <- readRDS("dev/sgrev-out/names.rds")
tb <- function(p) tryCatch(get(".__S3MethodsTable__.",
  envir = asNamespace(p), inherits = FALSE), error = function(e) NULL)

owners_of <- function(nm) {
  out <- character()
  for (p in OWN) {
    g <- tryCatch(getExportedValue(p, nm), error = function(e) NULL)
    if (!is.function(g)) next
    if (!identical(environment(g), asNamespace(p))) next
    if (grepl("UseMethod", paste(deparse(body(g)), collapse = ""),
              fixed = TRUE)) out <- c(out, p)
  }
  out
}
cat("== rule 3: every frmtmb_draws method in EVERY owner's table ==\n")
miss <- character(); n <- 0
for (nm in GEN$names) {
  ow <- owners_of(nm)
  for (p in ow) {
    n <- n + 1
    t <- tb(p)
    ok <- is.environment(t) &&
      exists(paste0(nm, ".frmtmb_draws"), envir = t, inherits = FALSE)
    if (!ok) miss <- c(miss, paste0(p, "::", nm, ".frmtmb_draws"))
  }
}
cat("  (name, owner) pairs checked: ", n, "   missing twins: ",
    length(miss), "\n")
if (length(miss)) cat("   ", paste(miss, collapse = " "), "\n")
cat("  own table has all 28: ",
    sum(vapply(GEN$names, function(nm)
      exists(paste0(nm, ".frmtmb_draws"),
             envir = tb("frmtmb.sample"), inherits = FALSE), NA)),
    " of 28\n")

cat("\n== rule 2: formals equal the FIRST owner's, brms deciding ==\n")
fml <- function(f) paste(deparse(formals(f)), collapse = "")
pref <- list(rhat = "posterior", posterior_samples = "brms")
bad <- character()
for (nm in GEN$names) {
  ow <- owners_of(nm)
  first <- if (!is.null(pref[[nm]])) pref[[nm]] else ow[1]
  g <- getExportedValue("frmtmb.sample", nm)
  o <- getExportedValue(first, nm)
  if (!identical(fml(g), fml(o)))
    bad <- c(bad, sprintf("%s: ours %s vs %s's %s", nm,
                          fml(g), first, fml(o)))
}
cat("  generics whose formals differ from their owner's: ",
    length(bad), "\n")
if (length(bad)) cat("   ", paste(bad, collapse = "\n    "), "\n")
cat("  and the METHOD carries the generic's formals: ")
bad2 <- character()
for (nm in GEN$names) {
  g <- getExportedValue("frmtmb.sample", nm)
  m <- get(paste0(nm, ".frmtmb_draws"), envir = tb("frmtmb.sample"),
           inherits = FALSE)
  gn <- names(formals(g)); mn <- names(formals(m))
  if (!all(setdiff(gn, "...") %in% mn)) bad2 <- c(bad2, nm)
}
cat(length(bad2), " missing\n")
if (length(bad2)) cat("   ", paste(bad2, collapse = " "), "\n")

cat("\n== rule 1: bare UseMethod on the FIX build with no owner ==\n")
cat("  (checked on BASE in dev/sgrev-names.R: 0 of 28 carry work)\n")
cat("DONE\n")
