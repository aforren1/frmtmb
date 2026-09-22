# Reviewer, priority 1: did any standard error MOVE?
#
# Runs on base and on lane and writes an RDS of everything numeric that
# a standard error can be read from. The comparison is a separate step
# so that nothing is compared inside a process that has only one build
# on its path.
#
#   Rscript dev/shapes-rev-se.R base
#   Rscript dev/shapes-rev-se.R lane

arg <- commandArgs(trailingOnly = TRUE)
which_lib <- if (identical(arg[1], "base")) "base" else "lane"
lib <- if (which_lib == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
cat("lib:", lib, "frmtmb", as.character(packageVersion("frmtmb")), "\n")

dd <- rev_data()
fits <- rev_fits(dd)
out <- list()

safe <- function(e) tryCatch(suppressWarnings(e), error = function(c)
  structure(list(msg = conditionMessage(c)), class = "revErr"))

for (nm in names(fits)) {
  f <- fits[[nm]]
  r <- list()
  # the covariance itself, both spellings
  r$vcov <- safe(as.matrix(vcov(f)))
  r$vcov_full <- safe(as.matrix(vcov(f, full = TRUE)))
  r$vcov_est <- safe(if (which_lib == "lane")
    as.matrix(frmtmb::vcov_estimated(f)) else as.matrix(vcov(f)))
  # everything that reads a covariance
  r$confint <- safe(as.matrix(confint(f)))
  r$summary_coef <- safe(as.matrix(summary(f)$coefficients))
  r$hypothesis <- safe({
    h <- hypothesis(f, "x = 0")
    as.matrix(h$hypothesis[, sapply(h$hypothesis, is.numeric)])
  })
  r$sandwich <- safe(as.matrix(vcov(f, cluster = ~ g)))
  r$se_fixed <- safe({
    v <- if (which_lib == "lane") frmtmb::vcov_estimated(f) else vcov(f)
    setNames(sqrt(diag(as.matrix(v))), rownames(as.matrix(v)))
  })
  # interop seams
  r$get_varcov <- safe(as.matrix(insight::get_varcov(f)))
  r$get_params <- safe(as.data.frame(insight::get_parameters(f)))
  r$emm <- safe({
    e <- emmeans::emmeans(f, "x", at = list(x = c(-1, 0, 1)))
    as.matrix(as.data.frame(summary(e))[, c("emmean", "SE")])
  })
  r$slopes <- safe({
    s <- marginaleffects::avg_slopes(f)
    as.matrix(s[, c("estimate", "std.error")])
  })
  r$slopes_term <- safe({
    s <- marginaleffects::avg_slopes(f)
    as.character(s$term)
  })
  r$mfx_pred <- safe({
    p <- marginaleffects::predictions(
      f, newdata = marginaleffects::datagrid(model = f, x = c(-1, 0, 1)))
    as.matrix(as.data.frame(p)[, c("estimate", "std.error")])
  })
  r$influence <- safe({
    i <- influence(f)
    list(dfbetas = as.matrix(dfbetas(i)),
         cooks = as.numeric(cooks.distance(i)))
  })
  out[[nm]] <- r
}

# importance, which reads standard errors through its own path
out$importance <- safe(as.data.frame(frm_importance(fits$gaussian)))

saveRDS(out, file.path(TREE, paste0("dev/shapes-rev-se-", which_lib,
                                    ".rds")))
cat("wrote dev/shapes-rev-se-", which_lib, ".rds\n", sep = "")
