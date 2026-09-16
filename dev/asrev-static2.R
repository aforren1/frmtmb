## REVIEW 2.5e, second static pass: the swallow detector against a
## parse-tree detector, the 38 unguarded methods, the exemption's
## necessity, and the conditional_effects() message the guard broke.
.libPaths(c("C:/Users/adf44/source/r/asrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

ar_refusers <- c("stop", "fit_no_draws", "multiple_no_draws")
ar_refuses_always <- function(fn) {
  b <- body(fn)
  hd <- function(e) is.call(e) && as.character(e[[1L]])[1L] %in% ar_refusers
  if (!is.call(b) || !identical(as.character(b[[1L]]), "{")) return(hd(b))
  length(b) == 2L && hd(b[[2L]])
}
# the test's detector
ar_swallows <- function(fn) {
  fo <- names(formals(fn))
  if (!"..." %in% fo) return(FALSE)
  if (any(grepl("...", deparse(body(fn)), fixed = TRUE))) return(FALSE)
  !ar_refuses_always(fn)
}
# a detector that cannot be fooled by a string: the PARSE TREE
dots_names <- c("...", "..1", "..2", "..3", "...length", "...names",
                "...elt")
parse_swallows <- function(fn) {
  fo <- names(formals(fn))
  if (!"..." %in% fo) return(FALSE)
  an <- tryCatch(all.names(body(fn)), error = function(e) character())
  if (any(an %in% dots_names)) return(FALSE)
  !ar_refuses_always(fn)
}

methods_of <- function(pkg) {
  reg <- parseNamespaceFile(pkg, dirname(system.file(package = pkg)))$S3methods
  nm <- ifelse(is.na(reg[, 3]), paste(reg[, 1], reg[, 2], sep = "."), reg[, 3])
  unique(nm)
}
cat("==== detector agreement ====\n")
for (pkg in c("frmtmb", "frmtmb.sample")) {
  ns <- asNamespace(pkg)
  dis <- character()
  for (f in methods_of(pkg)) {
    o <- tryCatch(get(f, envir = ns), error = function(e) NULL)
    if (!is.function(o)) next
    if (!identical(ar_swallows(o), parse_swallows(o))) {
      dis <- c(dis, sprintf("%s (grep=%s parse=%s)", f,
                            ar_swallows(o), parse_swallows(o)))
    }
  }
  cat(pkg, ": disagreements ", length(dis), "\n", sep = "")
  if (length(dis)) cat(paste0("  ", dis), sep = "\n")
}

cat("\n==== methods exempted BY SHAPE (unconditional refusal) ====\n")
for (pkg in c("frmtmb", "frmtmb.sample")) {
  ns <- asNamespace(pkg)
  hits <- character()
  for (f in methods_of(pkg)) {
    o <- tryCatch(get(f, envir = ns), error = function(e) NULL)
    if (!is.function(o) || !"..." %in% names(formals(o))) next
    if (any(grepl("...", deparse(body(o)), fixed = TRUE))) next
    if (ar_refuses_always(o)) hits <- c(hits, f)
  }
  cat(pkg, ": ", length(hits), "\n", sep = "")
  cat(paste(strwrap(paste(hits, collapse = ", "), 76), collapse = "\n"), "\n")
}

cat("\n==== the conditional_effects() message the guard made wrong ====\n")
set.seed(3)
dd <- data.frame(x = rnorm(60))
dd$y <- 1 + 0.6 * dd$x + rnorm(60)
dl <- list(dd, dd)
fm <- tryCatch(frm_multiple(bf(y ~ x) + gaussian(), data = dl),
               error = function(e) e)
if (inherits(fm, "error")) {
  cat("frm_multiple failed:", conditionMessage(fm), "\n")
} else {
  cat("conditional_effects(fm, \"x\"):",
      tryCatch({ conditional_effects(fm, "x"); "NO ERROR" },
               error = conditionMessage), "\n")
  cat("conditional_effects(fm, effects = \"x\"):",
      tryCatch({ conditional_effects(fm, effects = "x"); "NO ERROR" },
               error = conditionMessage), "\n")
}

cat("\n==== the exemption's necessity, constructed ====\n")
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
if (requireNamespace("insight", quietly = TRUE)) {
  for (lab in c("component", "effects")) {
    cl <- list(quote(insight::get_parameters), quote(fit))
    cl[[lab]] <- "all"
    cat(sprintf("  get_parameters(fit, %s = 'all') -> %s\n", lab,
                tryCatch({ eval(as.call(cl)); "ok (swallowed)" },
                         error = conditionMessage)))
  }
  cat("  n_parameters(fit, component = 'conditional') ->",
      tryCatch({ insight::n_parameters(fit, component = "conditional");
                 "ok" }, error = conditionMessage), "\n")
}
if (requireNamespace("emmeans", quietly = TRUE)) {
  cat("  emm_basis dots: ",
      paste(names(formals(getFromNamespace("emm_basis.frmtmb_fit",
                                           "frmtmb"))), collapse = ", "),
      "\n")
}
cat("\nDONE\n")
