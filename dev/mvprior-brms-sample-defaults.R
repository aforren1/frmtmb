# Lane wt-mvprior: brms 2.23.0's default_prior() on the data of the
# several-location blocks in frmtmb.sample's test-default-priors-brms.R
# (dpb_several_data(), copied here verbatim).
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages(library(brms))
dpb_several_data <- function() {
  set.seed(2309)
  n <- 80
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(1:8, length.out = n)))
  d$cat <- factor(c("a", "b", "c")[1 + (d$x + rlogis(n) > 0) +
                                      (d$z + rlogis(n) > 0.5)])
  d$ym <- ifelse(rbinom(n, 1, 0.5) == 1, 3 + d$x, -1 + d$x) +
    rnorm(n, 0, 0.5)
  d
}
d <- dpb_several_data()
show <- function(tab) {
  tab <- as.data.frame(tab)
  print(tab[tab$prior != "", c("prior", "class", "group", "dpar")],
        row.names = FALSE)
}
show(default_prior(cat ~ x + (1 | g), d, family = categorical()))
show(default_prior(ym ~ x + (1 | g), d,
                   family = mixture(gaussian(), gaussian())))
