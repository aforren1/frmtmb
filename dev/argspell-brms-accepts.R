## Does brms ACCEPT `re.form` on each method, as opposed to DECLARING
## it? The review falsified the formals-deep reading: the whole body of
## `predictive_interval.brmsfit` is `posterior_predict(object, ...)`,
## and `posterior_predict.brmsfit` has a `re.form` formal, so brms
## honors `predictive_interval(x, re.form = NA)` one frame down.
##
## Run: Rscript dev/argspell-brms-accepts.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", as.character(utils::packageVersion("brms")), "\n\n")

forwards_to <- c("posterior_predict", "posterior_epred",
                 "posterior_linpred", "prepare_predictions")

declares <- function(gen) {
  m <- tryCatch(getFromNamespace(paste0(gen, ".brmsfit"), "brms"),
                error = function(e) NULL)
  !is.null(m) && "re.form" %in% names(formals(m))
}

gens <- c("posterior_epred", "posterior_linpred", "posterior_predict",
          "predictive_error", "predictive_interval", "pp_check")
for (g in gens) {
  m <- getFromNamespace(paste0(g, ".brmsfit"), "brms")
  nms <- all.names(body(m))
  callees <- intersect(nms, forwards_to)
  reach <- Filter(declares, callees)
  cat(sprintf("%-20s declares %-5s body calls {%s} declaring {%s}\n",
              g, declares(g), paste(callees, collapse = ", "),
              paste(reach, collapse = ", ")))
}

cat("\nbody of predictive_interval.brmsfit:\n")
print(body(getFromNamespace("predictive_interval.brmsfit", "brms")))
cat("\nbody of pp_check.brmsfit, dots forwarding lines:\n")
b <- deparse(body(getFromNamespace("pp_check.brmsfit", "brms")))
writeLines(grep("[.][.][.]", b, value = TRUE))
