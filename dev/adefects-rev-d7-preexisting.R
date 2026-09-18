source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
# Is the "fallback table loses the sibling's method when the owner was
# loaded first" property NEW to log_lik, or does it already hold for the
# generics frmtmb.sample already re-exported?
loadNamespace("rstantools"); loadNamespace("loo"); loadNamespace("posterior")
suppressMessages(library(frmtmb)); suppressMessages(library(frmtmb.sample))
tb <- function(p) if (isNamespaceLoaded(p))
  ls(get(".__S3MethodsTable__.", envir = asNamespace(p))) else character(0)
fr <- tb("frmtmb"); sa <- tb("frmtmb.sample")
gens <- c("log_lik", "loo", "waic", "bayes_R2", "posterior_epred",
          "posterior_predict", "as_draws_df", "variables", "ndraws",
          "hypothesis", "conditional_effects", "pp_check", "LOO", "WAIC")
for (g in gens) {
  cat(sprintf("%-20s frmtmb:%-28s sample:%s\n", g,
              paste(grep(paste0("^", g, "\\."), fr, value = TRUE),
                    collapse = ","),
              paste(grep(paste0("^", g, "\\."), sa, value = TRUE),
                    collapse = ",")))
}
cat("--- after unloading the owners ---\n")
for (p in c("rstantools", "loo", "posterior")) {
  try(unloadNamespace(p), silent = TRUE)
}
ds <- structure(list(), class = "frmtmb_draws")
for (g in gens) {
  r <- tryCatch(do.call(g, list(ds)), error = function(e) e)
  m <- if (inherits(r, "condition")) conditionMessage(r) else "VALUE"
  cat(sprintf("%-20s %s\n", g,
              if (grepl("no applicable method", m, fixed = TRUE))
                "NO_METHOD" else
              if (grepl("could not find function", m, fixed = TRUE))
                "NO_GENERIC" else "REACHED"))
}
