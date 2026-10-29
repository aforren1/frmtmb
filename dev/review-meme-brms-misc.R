# Review of lane meme: brms 2.23.0's Stan code for prior_Lme at d = 3
# and for a negative lower bound on sdme. No compile.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
sc <- stancode(bf(y ~ me(x, sx) + me(z, sz) + me(w, sw)), data = d,
               prior = c(set_prior("normal(0,1)", class = "meanme"),
                         set_prior("exponential(1)", class = "sdme"),
                         set_prior("lkj(2)", class = "corme")),
               sample_prior = "yes")
cat(grep("prior_(Lme|sdme|meanme)", strsplit(sc, "\n")[[1]], value = TRUE),
    sep = "\n")
sc2 <- stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                prior = set_prior("exponential(1)", class = "sdme", lb = -1))
cat(grep("sdme_1", strsplit(sc2, "\n")[[1]], value = TRUE), sep = "\n")
sc3 <- stancode(bf(y ~ x + (1 | g)), data = d,
                prior = set_prior("exponential(1)", class = "sd", lb = -1))
cat(grep("sd_1", strsplit(sc3, "\n")[[1]], value = TRUE), sep = "\n")
