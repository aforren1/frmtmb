# Follow-up to dev/sgrev-rhat.R. Two questions the first run raised:
#   a) rhat(ds) and every other draws accessor do not name the same
#      variables;
#   b) is rhat the ONLY diagnostic that delegates to bayesplot on the
#      stanfit where brms uses posterior on the draws?
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
q(library(posterior)); q(library(bayesplot)); q(library(rstan))
q(requireNamespace("brms", quietly = TRUE))

ds <- readRDS("dev/stan-cache/sgrev-draws.rds")
cat("== a) the variable names ==\n")
cat("variables(ds)        : ", paste(variables(ds), collapse = " "), "\n")
cat("names(rhat(ds))      : ", paste(names(rhat(ds)), collapse = " "), "\n")
cat("names(neff_ratio(ds)): ", paste(names(neff_ratio(ds)),
                                     collapse = " "), "\n")
cat("colnames(as.matrix(ds)): ",
    paste(colnames(as.matrix(ds)), collapse = " "), "\n")
cat("in rhat() but not variables(): ",
    paste(setdiff(names(rhat(ds)), variables(ds)), collapse = " "), "\n")
cat("in variables() but not rhat(): ",
    paste(setdiff(variables(ds), names(rhat(ds))), collapse = " "), "\n\n")

cat("== b) every draws method whose body reaches bayesplot ==\n")
t <- get(".__S3MethodsTable__.", envir = asNamespace("frmtmb.sample"),
         inherits = FALSE)
for (k in sort(ls(t, all.names = TRUE))) {
  f <- get(k, envir = t, inherits = FALSE)
  if (!is.function(f)) next
  b <- paste(deparse(body(f)), collapse = " ")
  if (grepl("bayesplot", b, fixed = TRUE))
    cat("  ", k, ": ", substr(b, 1, 120), "\n", sep = "")
}

cat("\n== brms's own diagnostics, for the same names ==\n")
pt <- get(".__S3MethodsTable__.", envir = asNamespace("posterior"),
          inherits = FALSE)
bt <- get(".__S3MethodsTable__.", envir = asNamespace("bayesplot"),
          inherits = FALSE)
for (nm in c("neff_ratio", "nuts_params", "log_posterior")) {
  m <- NULL
  for (tb in list(pt, bt)) {
    mm <- tryCatch(get(paste0(nm, ".brmsfit"), envir = tb,
                       inherits = FALSE), error = function(e) NULL)
    if (is.function(mm)) m <- mm
  }
  cat("-- brms's ", nm, ".brmsfit:\n", sep = "")
  if (is.function(m)) cat(paste(deparse(body(m)), collapse = "\n"), "\n")
  mb <- tryCatch(get(paste0(nm, ".stanfit"), envir = bt,
                     inherits = FALSE), error = function(e) NULL)
  cat("-- bayesplot's ", nm, ".stanfit:\n", sep = "")
  if (is.function(mb)) cat(paste(deparse(body(mb)), collapse = "\n"), "\n")
}

cat("\n== neff_ratio numbers ==\n")
a <- neff_ratio(ds)
arr <- as_draws_array(ds)
sd <- posterior::summarise_draws(arr, ess_bulk = posterior::ess_bulk)
pb <- setNames(sd$ess_bulk / posterior::ndraws(arr), sd$variable)
cm <- intersect(names(a), names(pb))
cat(sprintf("%-14s %14s %14s %10s\n", "variable", "neff_ratio(ds)",
            "posterior bulk", "rel.diff"))
for (v in cm) cat(sprintf("%-14s %14.8f %14.8f %10.3e\n", v, a[[v]],
                          pb[[v]], abs(a[[v]] - pb[[v]]) / pb[[v]]))
cat("max relative difference: ",
    format(max(abs(a[cm] - pb[cm]) / pb[cm]), digits = 6), "\n")
cat("DONE\n")
