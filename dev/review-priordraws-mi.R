# Reviewer: where brms 2.23.0 puts the mi() measurement term of a
# response with known noise under sample_prior = "only", and what
# frm_sample() does with that model.
#   Rscript dev/review-priordraws-mi.R
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + d$x + rnorm(40)
d$sdy <- 0.5
sc <- brms::stancode(brms::bf(y | mi(sdy) ~ x), data = d,
                     prior = brms::set_prior("normal(0, 1)", class = "b"),
                     sample_prior = "only")
ln <- strsplit(sc, "\n")[[1]]
i <- grep("model \\{", ln)
cat(ln[i:length(ln)], sep = "\n")
r <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
  bf(y | mi(sdy) ~ x), data = d, family = gaussian(),
  prior = set_prior("normal(0, 1)", class = "b"), sample_prior = "only",
  chains = 1, iter = 200, seed = 1, refresh = 0))), error = function(e) e)
cat("frm_sample only:", if (inherits(r, "error")) conditionMessage(r) else "sampled", "\n")
