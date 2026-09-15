.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
suppressMessages(library(frmtmb))

fake <- structure(list(), class = "brmsfit")

owner <- function(nm) {
  f <- tryCatch(get(nm), error = function(e) NULL)
  if (is.null(f)) return("absent")
  e <- environment(f)
  if (is.null(e)) return("base/primitive")
  n <- environmentName(topenv(e))
  if (!nzchar(n)) "?" else n
}

gens <- c("conditional_effects", "hypothesis", "pp_check", "bayes_R2",
          "posterior_summary", "prior_summary", "expose_functions",
          "fixef", "ranef", "VarCorr",
          "loo", "loo_compare", "waic",
          "as_draws_df", "variables", "ndraws", "nchains")

cat(sprintf("%-20s %-12s %-10s %-10s\n",
            "generic", "resolves-to", "brms-meth", "dispatch"))
for (g in gens) {
  own <- owner(g)
  has <- !is.null(tryCatch(getS3method(g, "brmsfit", optional = TRUE),
                           error = function(e) NULL))
  disp <- tryCatch({
    do.call(g, list(fake))
    "ok"
  }, error = function(e) {
    m <- conditionMessage(e)
    if (grepl("no applicable method", m)) "NO METHOD" else "dispatched"
  })
  cat(sprintf("%-20s %-12s %-10s %-10s\n", g, own, has, disp))
}

cat("\n--- owner of the generic BEFORE frmtmb is attached ---\n")
cat("run separately; see second block\n")
