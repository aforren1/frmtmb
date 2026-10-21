# Reviewer: what brms 2.23.0 itself accepts and writes for the two
# families (stancode()/default_prior(), no compile).
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(6, 10),
                 e = exp(runif(n)))
dd$y <- rnbinom(n, mu = 3, size = 2)
dd$N <- sample(3:9, n, TRUE)
Y <- t(sapply(dd$N, function(k) rmultinom(1, k, c(0.3, 0.3, 0.4))))
colnames(Y) <- c("a", "b", "c")
dd$Y <- Y
dd$w <- rnorm(n)
nb2 <- brms::brmsfamily("negbinomial2")
lik <- function(code) {
  l <- strsplit(code, "\n")[[1]]
  l[grepl("target \\+=|_lpmf|_lpdf|lprior", l)]
}
probe("nb2 sigma ~ z stancode",
      lik(stancode(bf(y ~ x, sigma ~ z), dd, family = nb2)))
probe("nb2 y | rate(e) stancode", lik(stancode(y | rate(e) ~ x, dd,
                                              family = nb2)))
probe("nb2 cens stancode", {
  dd$cc <- rep(0:1, length.out = n)
  lik(stancode(y | cens(cc) ~ x, dd, family = nb2))
})
probe("nb2 trunc stancode", lik(stancode(y | trunc(lb = 0, ub = 30) ~ x, dd,
                                         family = nb2)))
probe("mixture(nb2, poisson) stancode",
      lik(stancode(y ~ x, dd, family = mixture(nb2, poisson()))))
probe("mvbf nb2 + gaussian stancode",
      lik(stancode(bf(y ~ x, family = nb2) + bf(w ~ x) + set_rescor(FALSE),
                   dd)))
probe("default_prior nb2 y ~ x", default_prior(y ~ x, dd, family = nb2))
probe("default_prior nb2 identity sigma",
      default_prior(y ~ x, dd, family = brms::brmsfamily(
        "negbinomial2", link_sigma = "identity")))
probe("default_prior nb2 softplus sigma",
      default_prior(y ~ x, dd, family = brms::brmsfamily(
        "negbinomial2", link_sigma = "softplus")))
probe("dm phi ~ z stancode",
      lik(stancode(bf(Y | trials(N) ~ x, phi ~ z), dd,
                   family = dirichlet_multinomial())))
probe("dm lpmf body", {
  l <- strsplit(stancode(Y | trials(N) ~ x, dd,
                         family = dirichlet_multinomial()), "\n")[[1]]
  i <- grep("dirichlet_multinomial_logit2_lpmf", l)[1]
  l[i + 0:12]
})
probe("dm without trials()", lik(stancode(Y ~ x, dd,
                                          family = dirichlet_multinomial())))
probe("dm mvbf + gaussian",
      lik(stancode(bf(Y | trials(N) ~ x, family = dirichlet_multinomial()) +
                     bf(w ~ x) + set_rescor(FALSE), dd)))
probe("dm refcat NA", lik(stancode(Y | trials(N) ~ x, dd,
                                   family = dirichlet_multinomial(refcat = NA))))
probe("dm mixture", lik(stancode(Y | trials(N) ~ x, dd,
                                 family = mixture(dirichlet_multinomial(),
                                                  multinomial()))))
probe("dm default_prior", default_prior(Y | trials(N) ~ x, dd,
                                        family = dirichlet_multinomial()))
probe("dm cens", {
  dd$cc <- 0
  lik(stancode(Y | trials(N) + cens(cc) ~ x, dd,
               family = dirichlet_multinomial()))
})
probe("dm se()", lik(stancode(Y | trials(N) + weights(N) ~ x, dd,
                              family = dirichlet_multinomial())))
probe("dm trials not equal row sums", {
  d2 <- dd
  d2$N[1] <- d2$N[1] + 1
  sd <- standata(Y | trials(N) ~ x, d2, family = dirichlet_multinomial())
  c(sd$trials[1], sum(sd$Y[1, ]))
})
probe("brmsfamily negbinomial2 links", nb2$link_sigma)
probe("brms .family_negbinomial2", brms:::.family_negbinomial2())
probe("brms .family_dirichlet_multinomial",
      brms:::.family_dirichlet_multinomial())
