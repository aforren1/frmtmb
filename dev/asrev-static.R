## REVIEW 2.5e: static claims. brms formals (risk 2), the frame attack on
## frm_check_dots, the exemption list, and the swallow detector's holes.
LIB <- "C:/Users/adf44/source/r/asrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

cat("==== A. brms 2.23.0 formals, risk 2 ====\n")
cat("brms", as.character(packageVersion("brms")), "\n")
bns <- asNamespace("brms")
for (g in c("posterior_epred", "posterior_linpred", "posterior_predict",
            "predictive_error", "predictive_interval", "pp_check",
            "fitted", "predict", "residuals")) {
  m <- tryCatch(getFromNamespace(paste0(g, ".brmsfit"), "brms"),
                error = function(e) NULL)
  if (is.null(m)) { cat(sprintf("%-22s <no brmsfit method>\n", g)); next }
  fo <- names(formals(m))
  cat(sprintf("%-22s re_formula=%-5s re.form=%-5s | %s\n", g,
              "re_formula" %in% fo, "re.form" %in% fo,
              paste(fo, collapse = ", ")))
}
cat("\nprepare_predictions formals:\n")
pp <- tryCatch(getFromNamespace("prepare_predictions", "brms"),
               error = function(e) NULL)
if (!is.null(pp)) cat(paste(names(formals(pp)), collapse = ", "), "\n")
ppd <- tryCatch(getFromNamespace("prepare_predictions.brmsfit", "brms"),
                error = function(e) NULL)
if (!is.null(ppd)) {
  cat("prepare_predictions.brmsfit:\n")
  cat(paste(names(formals(ppd)), collapse = ", "), "\n")
}
cat("\npp_check.brmsfit body mentions prepare_predictions? ",
    any(grepl("prepare_predictions",
              deparse(body(getFromNamespace("pp_check.brmsfit", "brms"))),
              fixed = TRUE)), "\n")
cat("pp_check.brmsfit body: re.form mentioned? ",
    any(grepl("re.form",
              deparse(body(getFromNamespace("pp_check.brmsfit", "brms"))),
              fixed = TRUE)), "\n")
cat("predictive_interval.brmsfit body:\n")
print(body(getFromNamespace("predictive_interval.brmsfit", "brms")))

cat("\n==== B. frame attack on frm_check_dots ====\n")
fcd <- getFromNamespace("frm_check_dots", "frmtmb")
cfn <- getFromNamespace("calling_fun_name", "frmtmb")

# B1 plain call
m1 <- function(object, ...) { fcd(...); "ran" }
cat("B1 direct   :", tryCatch(m1(1, bogus = 2), error = conditionMessage), "\n")

# B2 do.call with a character name
cat("B2 do.call(chr):",
    tryCatch(do.call("m1", list(1, bogus = 2)), error = conditionMessage), "\n")
# B3 do.call with the function OBJECT
cat("B3 do.call(fn) :",
    tryCatch(do.call(m1, list(1, bogus = 2)), error = conditionMessage), "\n")

# B4 S3 dispatch through UseMethod
gen <- function(object, ...) UseMethod("gen")
gen.asrev <- function(object, ...) { fcd(...); "ran" }
x <- structure(list(), class = c("asrev", "asrev_base"))
cat("B4 UseMethod:", tryCatch(gen(x, bogus = 2), error = conditionMessage), "\n")

# B5 NextMethod
gen.asrev2 <- function(object, ...) NextMethod()
gen.asrev_base <- function(object, ...) { fcd(...); "ran" }
y <- structure(list(), class = c("asrev2", "asrev_base"))
cat("B5 NextMethod:", tryCatch(gen(y, bogus = 2), error = conditionMessage), "\n")

# B6 nested helper: a method that delegates the check to a helper
helper <- function(...) fcd(...)
m6 <- function(object, realarg = 1, ...) { helper(...); "ran" }
cat("B6 nested   :", tryCatch(m6(1, bogus = 2), error = conditionMessage), "\n")

# B7 anonymous / lapply frame
cat("B7 lapply   :",
    tryCatch((function(object, ...) fcd(...))(1, bogus = 2),
             error = conditionMessage), "\n")

# B8 real package methods through do.call and through the generic
cat("\n-- real methods --\n")

cat("\n==== C. the exemption list ====\n")
ns <- asNamespace("frmtmb")
exempt <- c(
  "emm_basis.frmtmb_fit", "recover_data.frmtmb_fit",
  "get_coef.frmtmb_fit", "get_predict.frmtmb_fit", "get_vcov.frmtmb_fit",
  "set_coef.frmtmb_fit", "get_varcov.frmtmb_fit",
  "get_parameters.frmtmb_fit", "find_formula.frmtmb_fit",
  "find_random.frmtmb_fit", "find_statistic.frmtmb_fit",
  "link_function.frmtmb_fit", "link_inverse.frmtmb_fit")
ar_refusers <- c("stop", "fit_no_draws", "multiple_no_draws")
ar_refuses_always <- function(fn) {
  b <- body(fn)
  head_is_refusal <- function(e) {
    is.call(e) && as.character(e[[1L]])[1L] %in% ar_refusers
  }
  if (!is.call(b) || !identical(as.character(b[[1L]]), "{")) {
    return(head_is_refusal(b))
  }
  length(b) == 2L && head_is_refusal(b[[2L]])
}
ar_swallows <- function(fn) {
  fo <- names(formals(fn))
  if (!"..." %in% fo) return(FALSE)
  if (any(grepl("...", deparse(body(fn)), fixed = TRUE))) return(FALSE)
  !ar_refuses_always(fn)
}
for (e in exempt) {
  o <- tryCatch(get(e, envir = ns), error = function(x) NULL)
  if (is.null(o)) { cat(sprintf("%-28s MISSING\n", e)); next }
  cat(sprintf("%-28s swallows=%-5s formals: %s\n", e, ar_swallows(o),
              paste(names(formals(o)), collapse = ", ")))
}

cat("\n==== D. swallow-detector holes ====\n")
## D1: a literal '...' inside a STRING in the body defeats the fixed grep
d1 <- function(object, ...) { message("blah ... blah"); object }
cat("D1 string-ellipsis swallower detected as swallower? ",
    ar_swallows(d1), "  (want TRUE)\n")
## D2: a comment is dropped by deparse, so it cannot fool it
d2 <- function(object, ...) { object }
cat("D2 plain swallower: ", ar_swallows(d2), " (want TRUE)\n")
## D3: unconditional refusal for an UNRELATED reason
d3 <- function(object, ...) stop("not implemented yet")
cat("D3 unrelated unconditional stop: swallows=", ar_swallows(d3),
    " (exempted by shape)\n")
## D4: a refusal whose head is a helper NOT in ar_refusers
d4 <- function(object, ...) frmtmb_abort("no")
cat("D4 stop via another helper: swallows=", ar_swallows(d4), "\n")

cat("\n-- do any real method bodies contain a literal '...' in a string? --\n")
scan_pkg <- function(pkg) {
  n <- asNamespace(pkg)
  s3 <- parseNamespaceFile(pkg, dirname(system.file(package = pkg)))
  reg <- s3$S3methods
  nm <- ifelse(is.na(reg[, 3]), paste(reg[, 1], reg[, 2], sep = "."), reg[, 3])
  nm <- unique(nm)
  for (f in nm) {
    o <- tryCatch(get(f, envir = n), error = function(e) NULL)
    if (is.null(o) || !is.function(o)) next
    fo <- names(formals(o))
    if (!"..." %in% fo) next
    d <- deparse(body(o))
    # lines that mention ... only inside a quoted string
    lines <- grep("...", d, fixed = TRUE, value = TRUE)
    strong <- grep("(\\.\\.\\.[,)= ]|\\(\\.\\.\\.\\)|list\\(\\.\\.\\.|\\.\\.\\.length|\\.\\.\\.names|\\.\\.\\.elt)",
                   lines)
    if (length(lines) && !length(strong)) {
      cat(sprintf("  %s::%s  ONLY-IN-STRING: %s\n", pkg, f,
                  trimws(lines[1L])))
    }
  }
}
scan_pkg("frmtmb"); scan_pkg("frmtmb.sample")

cat("\n==== E. fit_no_draws / multiple_no_draws are unconditional ====\n")
for (f in c("fit_no_draws", "multiple_no_draws")) {
  o <- tryCatch(get(f, envir = ns), error = function(e) NULL)
  if (is.null(o)) { cat(f, "MISSING in frmtmb\n"); next }
  cat("--", f, "--\n"); print(o)
}

cat("\n==== F. the rawNamespace export seam ====\n")
cat("frm_check_dots exported from frmtmb? ",
    "frm_check_dots" %in% getNamespaceExports("frmtmb"), "\n")
cat("visible via frmtmb:: ? ",
    tryCatch(is.function(get("frm_check_dots", asNamespace("frmtmb"))),
             error = function(e) FALSE), "\n")

cat("\nDONE\n")
