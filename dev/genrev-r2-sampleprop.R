# genrev round 2: frmtmb.sample attached, brms only LOADED.  For every
# name frmtmb.sample re-exports from core that is a borrowed generic,
# report what UseMethod() would run on a brmsfit: the class method, a
# .default (and whose), nothing, or the name is not reachable at all.
# The name list is read from frmtmb.sample's own NAMESPACE, not typed.
LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
own <- function(f) {
  e <- environment(f); if (is.null(e)) return("base")
  n <- environmentName(topenv(e)); if (nzchar(n)) n else "<anon>"
}
nsi <- parseNamespaceFile("frmtmb.sample", dirname(system.file(package = "frmtmb.sample")))
imp <- Filter(function(x) is.list(x) && x[[1]] == "frmtmb", nsi$imports)
fromcore <- unique(unlist(lapply(imp, `[[`, 2)))
borrowed <- c("as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
              "as_draws_matrix", "as_draws_rvars", "nchains", "ndraws",
              "niterations", "nvariables", "variables", "loo", "loo_compare",
              "waic", "bayes_R2", "prior_summary", "pp_check",
              "conditional_effects", "expose_functions", "hypothesis",
              "posterior_summary", "LOO", "WAIC", "ngrps", "fixef", "ranef",
              "VarCorr", "refit")
gens <- intersect(borrowed, intersect(fromcore, nsi$exports))
q(library(frmtmb.sample))
q(loadNamespace("brms"))
res <- character()
for (g in gens) {
  f <- tryCatch(get(g, envir = globalenv()), error = function(e) NULL)
  if (!is.function(f)) { res[g] <- "absent"; next }
  genv <- environment(f)
  tb <- get(".__S3MethodsTable__.", envir = genv)
  if (exists(paste0(g, ".brmsfit"), envir = tb, inherits = FALSE)) {
    res[g] <- "class"
  } else if (exists(paste0(g, ".default"), envir = tb, inherits = FALSE)) {
    res[g] <- paste0("default(", own(get(paste0(g, ".default"), envir = tb)), ")")
  } else res[g] <- "none"
}
cat(sprintf("LIB %s\n", LIB))
cat(sprintf("re-exported borrowed generics: %d\n", length(gens)))
cat(sprintf("class method reached: %d  default: %d  none: %d  absent: %d\n",
            sum(res == "class"), sum(startsWith(res, "default")),
            sum(res == "none"), sum(res == "absent")))
cat("not class:", paste(sprintf("%s=%s", names(res)[res != "class"],
                                  res[res != "class"]), collapse = ", "), "\n")
cat("GENREVDONE\n")
