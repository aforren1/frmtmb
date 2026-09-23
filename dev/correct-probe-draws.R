# Lane wt-correct, item 2 on draws: pp_check() of frmtmb.sample's draws
# for every type bayesplot exposes, on the gaussian design of
# dev/correct-probe.R. Usage: Rscript dev/correct-probe-draws.R <arm>
arm <- commandArgs(trailingOnly = TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
if (identical(arm, "base")) .libPaths(.libPaths()[-1L])
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-data.R")
d <- correct_data_gauss()
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 1)))
ty <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
ns <- asNamespace("bayesplot")
for (t in ty) {
  fa <- names(formals(get(paste0("ppc_", t), ns)))
  args <- list(ds, type = t, ndraws = 20)
  if ("group" %in% fa) args$group <- "g"
  if ("x" %in% fa) args$x <- "x"
  r <- tryCatch({
    p <- suppressWarnings(suppressMessages(do.call(pp_check, args)))
    invisible(suppressWarnings(suppressMessages(ggplot2::ggplot_build(p))))
    "OK"
  }, error = function(e) {
    paste0("ERROR [", class(e)[1], "]: ",
           substr(gsub("\n", " ", conditionMessage(e)), 1, 120))
  })
  cat(sprintf("%-28s %s\n", t, r))
}
cat("violin:", tryCatch({pp_check(ds, type = "violin"); "OK"},
                        error = function(e) conditionMessage(e)), "\n")
