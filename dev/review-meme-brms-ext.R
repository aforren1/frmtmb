# Review of lane meme: brms 2.23.0 on the extended shapes.
#   Rscript dev/review-meme-brms-ext.R  -> dev/review-meme-out/brms-ext.rds
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
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
  do.call(c, lapply(pr, function(r) {
    set_prior(r[4], class = r[1], coef = r[2], group = r[3])
  }))
}
out <- list()
for (nm in names(sh)) {
  cat("====", nm, "\n")
  s <- sh[[nm]]
  r <- list()
  r$default <- tryCatch(tabdf(default_prior(s$f, data = d)),
                        error = function(e) conditionMessage(e))
  print(r$default)
  pr <- mkprior(s$pr)
  if (!is.null(pr)) {
    r$validate <- tryCatch(tabdf(validate_prior(pr, s$f, data = d)),
                           error = function(e) conditionMessage(e))
    print(r$validate)
  }
  r$stancode <- tryCatch(as.character(stancode(s$f, data = d, prior = pr,
                                                sample_prior = "yes")),
                         error = function(e) conditionMessage(e))
  r$standata <- tryCatch(names(standata(s$f, data = d, prior = pr)),
                         error = function(e) conditionMessage(e))
  if (!is.null(pr)) {
    fit <- tryCatch(suppressWarnings(suppressMessages(brm(
      s$f, data = d, prior = pr, sample_prior = "yes", chains = 1,
      iter = 300, seed = 1, refresh = 0))),
      error = function(e) conditionMessage(e))
    if (inherits(fit, "brmsfit")) {
      r$variables <- variables(fit)
      r$prior_draws <- names(prior_draws(fit))
      r$prior_summary <- tabdf(prior_summary(fit))
      sm <- summary(fit)
      r$summary <- lapply(sm[vapply(sm, is.data.frame, NA)], rownames)
      r$summary$random <- lapply(sm$random, rownames)
      cat("variables:", grep("^(Xme|lp__|lprior|Ymi)", r$variables,
                             value = TRUE, invert = TRUE), "\n",
          fill = 78)
      cat("prior_draws:", r$prior_draws, "\n", fill = 78)
      print(r$prior_summary)
    } else {
      r$fit_error <- fit
      cat("FIT ERROR:", fit, "\n")
    }
  }
  out[[nm]] <- r
}
saveRDS(out, file.path(wt, "dev/review-meme-out/brms-ext.rds"))
cat("DONE\n")
