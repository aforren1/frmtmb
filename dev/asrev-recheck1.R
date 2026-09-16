## BLOCKER 1 recheck: the parse-tree detector, both ways, and fix 3.
LIB <- "C:/Users/adf44/source/r/asrev-lib2"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

ar_dots_names <- c("...", "..1", "..2", "..3",
                   "...length", "...names", "...elt")
ar_refusers <- c("stop", "fit_no_draws", "multiple_no_draws")
ar_refuses_always <- function(fn) {
  b <- body(fn)
  hd <- function(e) is.call(e) && as.character(e[[1L]])[1L] %in% ar_refusers
  if (!is.call(b) || !identical(as.character(b[[1L]]), "{")) return(hd(b))
  length(b) == 2L && hd(b[[2L]])
}
ar_swallows <- function(fn) {
  fo <- names(formals(fn))
  if (!"..." %in% fo) return(FALSE)
  if (any(all.names(body(fn)) %in% ar_dots_names)) return(FALSE)
  !ar_refuses_always(fn)
}

cat("==== A. the detector, the OTHER way: legitimate dots readers ====\n")
cases <- list(
  "list(...) at top level"        = function(object, ...) list(...),
  "list(...) nested two deep"     = function(object, ...)
    do.call(base::c, list(list(object), list(...))),
  "...elt(1) nested"              = function(object, ...)
    if (...length()) identity(...elt(1L)) else object,
  "...elt only, no ...length"     = function(object, ...) ...elt(1L),
  "...names() only"               = function(object, ...) ...names(),
  "..1 deep inside a call"        = function(object, ...)
    base::identity(base::identity(..1)),
  "..4 (past the enumerated set)" = function(object, ...) ..4,
  "..11 (two digits)"             = function(object, ...) ..11,
  "match.call()$... "             = function(object, ...)
    match.call(expand.dots = FALSE)$...,
  "x[[i, ...]]"                   = function(x, i, ...) unclass(x)[[i, ...]],
  "NextMethod() only"             = function(object, ...) NextMethod(),
  "sys.call() only"               = function(object, ...) sys.call(),
  "string '...' only"             = function(object, ...) {
    message("a ... b"); object },
  "cat('... ') only"              = function(object, ...) {
    cat("... ", 1, " more\n", sep = ""); object },
  "plain swallower"               = function(object, ...) object
)
want <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
          FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
bad <- 0L
for (i in seq_along(cases)) {
  got <- ar_swallows(cases[[i]])
  flag <- if (identical(got, want[i])) "  " else "<<"
  if (!identical(got, want[i])) bad <- bad + 1L
  cat(sprintf("%s %-32s swallows=%-5s want=%-5s\n", flag, names(cases)[i],
              got, want[i]))
}
cat("mismatches:", bad, "\n")
cat("NOTE: NextMethod() and sys.call() forward the dots at RUNTIME but\n")
cat("name no dots symbol, so the detector calls them swallowers. That is\n")
cat("conservative, not unsafe: it can only over-report.\n")

cat("\n==== B. no method still swallows ====\n")
methods_of <- function(pkg) {
  reg <- parseNamespaceFile(pkg,
                            dirname(system.file(package = pkg)))$S3methods
  unique(ifelse(is.na(reg[, 3]), paste(reg[, 1], reg[, 2], sep = "."),
                reg[, 3]))
}
exempt <- c("emm_basis.frmtmb_fit", "get_predict.frmtmb_fit",
            "get_vcov.frmtmb_fit", "get_varcov.frmtmb_fit",
            "get_parameters.frmtmb_fit", "find_formula.frmtmb_fit",
            "find_random.frmtmb_fit", "get_coef.frmtmb_fit",
            "set_coef.frmtmb_fit")
for (pkg in c("frmtmb", "frmtmb.sample")) {
  ns <- asNamespace(pkg)
  bad <- Filter(function(f) {
    o <- tryCatch(get(f, envir = ns), error = function(e) NULL)
    !is.null(o) && is.function(o) && ar_swallows(o)
  }, setdiff(methods_of(pkg), exempt))
  cat(pkg, ": still swallowing ", length(bad), "  ",
      paste(bad, collapse = ", "), "\n", sep = "")
}

cat("\n==== C. print.frmtmb_par_template is guarded ====\n")
set.seed(3)
dd <- data.frame(x = rnorm(60)); dd$y <- 1 + 0.6 * dd$x + rnorm(60)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
pt <- par_template(fit)
cat("print(pt, nosucharg = 1) ->",
    tryCatch({ utils::capture.output(print(pt, nosucharg = 1)); "NO ERROR" },
             error = conditionMessage), "\n")
cat("print(pt, n = 2) still works ->",
    tryCatch({ utils::capture.output(print(pt, n = 2)); "ok" },
             error = conditionMessage), "\n")

cat("\n==== D. fix 3: predictive_interval and pp_check, ours ====\n")
for (nm in c("predictive_interval.frmtmb_draws", "pp_check.frmtmb_draws",
             "pp_check.frmtmb_fit")) {
  pk <- if (grepl("draws", nm)) "frmtmb.sample" else "frmtmb"
  cat(sprintf("  %-34s %s\n", nm,
              paste(names(formals(getFromNamespace(nm, pk))),
                    collapse = ", ")))
}
cat("\n  dual set (ours re.form == brms re.form):\n")
for (g in c("posterior_epred", "posterior_linpred", "posterior_predict",
            "predictive_error", "predictive_interval", "pp_check")) {
  ours <- tryCatch("re.form" %in%
                     names(formals(getFromNamespace(paste0(g,
                       ".frmtmb_draws"), "frmtmb.sample"))),
                   error = function(e) NA)
  bm <- tryCatch("re.form" %in%
                   names(formals(getFromNamespace(paste0(g, ".brmsfit"),
                                                  "brms"))),
                 error = function(e) NA)
  cat(sprintf("    %-20s ours %-5s brms-declared %-5s\n", g, ours, bm))
}
cat("\nDONE\n")
