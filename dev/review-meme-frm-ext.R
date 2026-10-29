# Review of lane meme: frmtmb on the extended shapes, both routes.
#   Rscript dev/review-meme-frm-ext.R  -> dev/review-meme-out/frm-ext.rds
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-meme"
source(file.path(wt, "dev/review-meme-shapes.R"))
d <- rv_data()
sh <- rv_shapes(bf, mvbf, set_mecor, set_rescor)
tabdf <- function(t) {
  t <- as.data.frame(t)
  t[, intersect(c("prior", "class", "coef", "group", "resp", "dpar",
                  "nlpar", "lb", "ub", "source"), names(t))]
}
mkprior <- function(pr) {
  if (!length(pr)) return(NULL)
  Reduce(`+`, lapply(pr, function(r) {
    set_prior(r[4], class = r[1], coef = r[2], group = r[3])
  }))
}
err <- function(e) paste("ERROR:", conditionMessage(e))
out <- list()
for (nm in names(sh)) {
  cat("====", nm, "\n")
  s <- sh[[nm]]
  r <- list()
  r$default <- tryCatch(tabdf(default_prior(s$f, data = d)), error = err)
  r$default_sample <- tryCatch(tabdf(default_prior(s$f, data = d,
                                                   route = "sample")),
                               error = err)
  print(r$default_sample)
  pr <- mkprior(s$pr)
  if (!is.null(pr)) {
    r$validate <- tryCatch(tabdf(validate_prior(pr, s$f, data = d)),
                           error = err)
    print(r$validate)
    fit <- tryCatch(suppressWarnings(frm(s$f, data = d, prior = pr)),
                    error = err)
    if (inherits(fit, "frmtmb_fit")) {
      r$fit_variables <- variables(fit)
      r$fit_prior_summary <- tabdf(prior_summary(fit))
      sm <- summary(fit)
      r$fit_summary_me <- if (!is.null(sm$me)) rownames(sm$me)
      cat("fit variables:", r$fit_variables, "\n", fill = 78)
    } else {
      r$fit_error <- fit
      cat("FIT", fit, "\n")
    }
    ds <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
      s$f, data = d, prior = pr, sample_prior = "yes", chains = 1,
      iter = 300, seed = 1, refresh = 0))), error = err)
    if (!is.character(ds)) {
      r$variables <- variables(ds)
      r$prior_draws <- names(prior_draws(ds))
      r$prior_summary <- tabdf(prior_summary(ds))
      cat("draws variables:", r$variables, "\n", fill = 78)
      cat("prior_draws:", r$prior_draws, "\n", fill = 78)
      print(r$prior_summary)
    } else {
      r$sample_error <- ds
      cat("SAMPLE", ds, "\n")
    }
  } else {
    r$fit_error <- tryCatch({frm(s$f, data = d); "fitted"}, error = err)
    cat("no-prior fit:", r$fit_error, "\n")
  }
  out[[nm]] <- r
}
saveRDS(out, file.path(wt, "dev/review-meme-out/frm-ext.rds"))
cat("DONE\n")
