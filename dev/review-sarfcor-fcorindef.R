# Reviewer: does brms's own Stan program for fcor() take an indefinite M?
# The cached fcor model of the lane's lpcheck, data from standata().
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/stan-cache")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib", "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-sarfcor/tests/testthat/helper-brms.R", envir = env)
set.seed(1); n <- 30
V <- 0.6^abs(outer(seq_len(n), seq_len(n), "-"))
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 5), each = 5)))
d$y <- 1 + 0.5 * d$x + as.vector(t(chol(V)) %*% rnorm(n)) * 0.8
form <- brms::bf(y ~ x + fcor(V))
prior <- env$brms_flat_prior(form, data = d, family = gaussian(), data2 = list(V = V))
code <- brms::make_stancode(form, data = d, family = gaussian(), prior = prior, data2 = list(V = V))
mod <- env$brms_stan_model(code)
Vi <- V; Vi[1, 2] <- Vi[2, 1] <- 1.3
cat("min eigen of Vi:", min(eigen(Vi)$values), "\n")
sd <- env$brms_standata(form, data = d, family = gaussian(), prior = prior, data2 = list(V = Vi))
r <- tryCatch({s <- rstan::sampling(mod, data = sd, chains = 1, iter = 20, refresh = 0)
  paste("draws:", nrow(as.matrix(s)))}, error = function(e) paste("ERROR", conditionMessage(e)))
cat("brms Stan program with indefinite M:", r, "\n")
