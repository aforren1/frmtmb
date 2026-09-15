# Deep ownership audit for Phase 2.5.
#
# The seeded dev/generics-audit.R answers "who exports this name".  That
# is not enough to pick an adoption target: a package that exports a
# name may DEFINE the generic or may import it from somewhere else, and
# an S3 method registers into the namespace where its GENERIC is
# defined.  So the question that decides the design is "where does
# brms's <gen>.brmsfit method actually live".
#
# Uses parseNamespaceFile(), never grep: nlme declares its exports in a
# multi-name block that a one-name-per-line grep cannot see.
LIB <- "C:/Users/adf44/source/r/rellib-r3"
UL  <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
SYS <- file.path(R.home(), "library")
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib", UL, SYS))

frm_gens <- c(
  "as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
  "as_draws_matrix", "as_draws_rvars", "bayes_R2", "conditional_effects",
  "expose_functions", "fixef", "hypothesis", "latent_probs", "loo",
  "LOO", "loo_compare", "nchains", "ndraws", "ngrps", "niterations",
  "nvariables", "par_template", "posterior_summary", "pp_check",
  "prior_summary", "ranef", "refit", "VarCorr", "variables", "waic",
  "WAIC"
)

cands <- c("brms", "posterior", "loo", "rstantools", "bayesplot",
           "nlme", "generics", "lme4", "glmmTMB", "cmdstanr")

# Does pkg DEFINE the generic, or import it?  parseNamespaceFile gives
# both the export list and the importFrom map.
info <- function(pkg) {
  ns <- tryCatch(parseNamespaceFile(pkg, dirname(find.package(pkg))),
                 error = function(e) NULL)
  if (is.null(ns)) return(NULL)
  imp <- list()
  for (i in ns$imports) {
    if (is.list(i) && length(i) == 2L) {
      from <- i[[1]]
      nms <- i[[2]]
      for (n in nms) imp[[n]] <- from
    }
  }
  list(exports = ns$exports, imports = imp,
       s3 = if (is.null(ns$S3methods)) matrix(character(), 0, 3)
            else ns$S3methods)
}

ii <- lapply(cands, function(p)
  tryCatch(info(p), error = function(e) NULL))
names(ii) <- cands

cat("### which package DEFINES each generic (export and not imported)\n\n")
cat(sprintf("%-20s %s\n", "generic", "defines / (imports from)"))
cat(strrep("-", 78), "\n")
for (g in frm_gens) {
  bits <- character()
  for (p in cands) {
    x <- ii[[p]]
    if (is.null(x)) next
    if (!(g %in% x$exports)) next
    src <- x$imports[[g]]
    bits <- c(bits, if (is.null(src)) sprintf("%s:DEF", p)
                    else sprintf("%s<-%s", p, src))
  }
  cat(sprintf("%-20s %s\n", g,
              if (length(bits)) paste(bits, collapse = "  ") else "-"))
}

cat("\n### where brms's <gen>.brmsfit method is registered\n\n")
suppressMessages(loadNamespace("brms"))
tab <- function(pkg) {
  ns <- tryCatch(asNamespace(pkg), error = function(e) NULL)
  if (is.null(ns)) return(character())
  if (!exists(".__S3MethodsTable__.", envir = ns, inherits = FALSE))
    return(character())
  ls(get(".__S3MethodsTable__.", envir = ns, inherits = FALSE))
}
loaded <- loadedNamespaces()
cat(sprintf("%-20s %-28s %s\n", "generic", "brmsfit method lives in",
            "generic env"))
cat(strrep("-", 78), "\n")
for (g in frm_gens) {
  m <- paste0(g, ".brmsfit")
  where <- loaded[vapply(loaded, function(p) m %in% tab(p), logical(1))]
  ge <- tryCatch(environmentName(environment(
    get(g, envir = asNamespace("brms")))), error = function(e) "-")
  cat(sprintf("%-20s %-28s %s\n", g,
              if (length(where)) paste(where, collapse = ",") else "(none)",
              ge))
}

cat("\n### formals of each owner's generic vs frmtmb's\n\n")
fr <- tryCatch(asNamespace("frmtmb"), error = function(e) NULL)
sig <- function(f) paste(names(formals(f)), collapse = ", ")
for (g in frm_gens) {
  ours <- tryCatch(sig(get(g, envir = fr)), error = function(e) NA)
  row <- sprintf("%-20s frmtmb(%s)", g, ours)
  for (p in cands) {
    x <- ii[[p]]
    if (is.null(x) || !(g %in% x$exports)) next
    o <- tryCatch(sig(getExportedValue(p, g)), error = function(e) NA)
    if (!is.na(o) && !identical(o, ours))
      row <- paste0(row, sprintf("  != %s(%s)", p, o))
  }
  cat(row, "\n")
}

cat("\nDONE\n")
