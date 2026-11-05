# Fix round 2: which families read "discrete"; the rlddm change on a copy; a mixture's type.
source("dev/review-nlfix-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("frmtmb")
fams <- list(poisson = poisson(), bernoulli = bernoulli(), binomial = binomial(), negbinomial = negbinomial(), beta_binomial = beta_binomial(), hurdle_poisson = hurdle_poisson(), hurdle_negbinomial = hurdle_negbinomial(), zip = zero_inflated_poisson(), zib = zero_inflated_binomial(), zibb = zero_inflated_beta_binomial(), geometric = geometric(), compois = compois(), discrete_weibull = discrete_weibull(), gaussian = gaussian(), hurdle_gamma = hurdle_gamma(), hurdle_lognormal = hurdle_lognormal(), zero_inflated_beta = zero_inflated_beta(), categorical = categorical(), multinomial = multinomial(), dirichlet_multinomial = dirichlet_multinomial(), cumulative = cumulative(),
  mix_pp = mixture(poisson(), poisson()), mix_gg = mixture(gaussian(), gaussian()))
for (nm in names(fams)) cat(sprintf("%-22s %s\n", nm, as.character(fams[[nm]][["type"]] %||% "<null>")))
suppressPackageStartupMessages({ library(frmtmb.learn); library(frmtmb.latent); library(frmtmb.eam); library(frmtmb.survival) })
cat("extension families:\n")
ex <- list(rlddm = rlddm(subject = id, trial = trial), bandit = bandit2arm_delta(subject = id, trial = trial), hmm_p = frmtmb.latent::hmm(2, poisson()), hmm_g = frmtmb.latent::hmm(2, gaussian()), wiener = wiener(), gng = wiener_gng(deadline = 1), rdm = rdm(2), cox = cox())
for (nm in names(ex)) cat(sprintf("%-22s %s\n", nm, as.character(ex[[nm]][["type"]] %||% "<null>")))
set.seed(13010); n <- 200; d <- data.frame(x = rnorm(n)); d$y2 <- c(rep(2L, 100), 5L + rbinom(100, 1, 0.01))
g <- tt(frm(y2 ~ 1, family = frmtmb.latent::hmm(2, poisson()), data = transform(d, id = factor(rep(1:20, each = 10)), tt = rep(1:10, 20))))
cat("hmm poisson:", if (iserr(g)) substr(g$msg, 1, 120) else as.numeric(logLik(g)), "\n")
