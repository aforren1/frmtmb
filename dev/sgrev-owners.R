# Ownership and formals, measured against BRMS rather than against the
# owner: brms is the project's tiebreaker.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(p) suppressMessages(suppressWarnings(
  requireNamespace(p, quietly = TRUE)))
for (p in c("brms", "bayesplot", "bridgesampling", "coda", "gratia",
            "loo", "posterior", "rstantools", "bbmle", "rstan"))
  q(p)

GEN <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/sgrev-out/names.rds")
nms <- GEN$names
fml <- function(f) if (is.function(f))
  sub("^function ", "", sub("[[:space:]]*NULL$", "",
      paste(deparse(args(f)), collapse = " "))) else NA_character_

tbl_of <- function(p) tryCatch(get(".__S3MethodsTable__.",
  envir = asNamespace(p), inherits = FALSE), error = function(e) NULL)

cat("== gratia posterior_samples methods ==\n")
t <- tbl_of("gratia")
cat(paste(grep("^posterior_samples", ls(t, all.names = TRUE),
                value = TRUE), collapse = " "), "\n")
cat("gratia version: ", as.character(packageVersion("gratia")), "\n")
cat("gratia exports posterior_samples: ",
    "posterior_samples" %in% getNamespaceExports("gratia"), "\n\n")

cat("== bbmle parnames ==\n")
f <- tryCatch(getExportedValue("bbmle", "parnames"), error = function(e) NULL)
cat("is a function: ", is.function(f), "  body has UseMethod: ",
    is.function(f) && grepl("UseMethod",
      paste(deparse(body(f)), collapse = ""), fixed = TRUE), "\n")
cat("bbmle parnames formals: ", fml(f), "\n")
cat("bbmle also has parnames<- : ",
    "parnames<-" %in% getNamespaceExports("bbmle"), "\n\n")

cat("== per name: brms's GENERIC (as brms sees it), brms's ",
    "brmsfit METHOD, each definer's generic, frmtmb.sample's now ==\n",
    sep = "")
cat(sprintf("%-20s %-52s %s\n", "name", "what", "formals"))
for (nm in nms) {
  cat("--", nm, "\n")
  # what brms itself calls: the value of the name INSIDE brms's ns
  bg <- tryCatch(get(nm, envir = asNamespace("brms")),
                 error = function(e) NULL)
  bgenv <- if (is.function(bg)) environmentName(environment(bg)) else NA
  cat(sprintf("   %-50s %s\n",
      paste0("brms's own ", nm, "() comes from ", bgenv), fml(bg)))
  # brms's method for brmsfit, wherever it is registered
  for (p in c("brms", "bayesplot", "bridgesampling", "coda", "loo",
              "posterior", "rstantools", "gratia")) {
    t <- tbl_of(p); if (is.null(t)) next
    m <- tryCatch(get(paste0(nm, ".brmsfit"), envir = t,
                      inherits = FALSE), error = function(e) NULL)
    if (is.function(m))
      cat(sprintf("   %-50s %s\n",
          paste0("brmsfit method in ", p, "'s table, from ",
                 environmentName(environment(m))), fml(m)))
  }
  # every package that DEFINES a generic of this name
  for (p in c("brms", "bayesplot", "bridgesampling", "coda", "gratia",
              "loo", "posterior", "rstantools", "bbmle", "rstan")) {
    g <- tryCatch(getExportedValue(p, nm), error = function(e) NULL)
    if (!is.function(g)) next
    if (!identical(environment(g), asNamespace(p))) next
    isgen <- grepl("UseMethod", paste(deparse(body(g)), collapse = ""),
                   fixed = TRUE)
    cat(sprintf("   %-50s %s\n",
        paste0("DEFINED by ", p, if (isgen) " (generic)" else
               " (NOT a generic)"), fml(g)))
  }
  fx <- tryCatch(getExportedValue("frmtmb.sample", nm),
                 error = function(e) NULL)
  cat(sprintf("   %-50s %s\n", "frmtmb.sample FIX generic", fml(fx)))
  mx <- tryCatch(get(paste0(nm, ".frmtmb_draws"),
                     envir = tbl_of("frmtmb.sample"), inherits = FALSE),
                 error = function(e) NULL)
  cat(sprintf("   %-50s %s\n", "frmtmb.sample frmtmb_draws METHOD",
              fml(mx)))
}
cat("DONE\n")
