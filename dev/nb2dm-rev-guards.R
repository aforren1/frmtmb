# Reviewer: each new refusal with its condition present and absent.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
ns <- asNamespace("frmtmb")
nb <- get("fam_negbinomial", ns)()
fin <- nb
fin[["family_finalize"]] <- function(fam, y, aterms) fam
id <- function(d) d
probe("fam_reparam: inner with family_finalize (present)",
      get("fam_reparam", ns)(fin, "t", c("mu", "shape"),
                             list(mu = "log", shape = "log"), id, list()))
probe("fam_reparam: inner without (absent)",
      get("fam_reparam", ns)(nb, "t", c("mu", "shape"),
                             list(mu = "log", shape = "log"), id,
                             list())$family)
set.seed(1)
d <- data.frame(x = rnorm(60), N = 10L)
d$yb <- rbinom(60, 10, 0.4)
d$Y <- t(sapply(1:60, function(i) rmultinom(1, 10, c(.2, .3, .5))))
probe("binomial without trials: bernoulli hint kept",
      frm(yb ~ x, family = binomial(), data = d))
probe("dirichlet_multinomial without trials: no hint",
      frm(Y ~ x, family = dirichlet_multinomial(), data = d))
probe("mvbf(binomial, dm), neither with trials: hint kept?",
      frm(mvbf(bf(yb ~ x, family = binomial()),
               bf(Y ~ x, family = dirichlet_multinomial())), data = d))
probe("dm link 'logit' accepted, 'probit' refused",
      list(dirichlet_multinomial("logit")$link,
           tryCatch(dirichlet_multinomial("probit"),
                    error = conditionMessage)))
probe("dm(3) number as link: K hint",
      tryCatch(dirichlet_multinomial(3), error = conditionMessage))
probe("dm with link_phi identity (allowed by brms?)",
      dirichlet_multinomial(link_phi = "identity")$link_phi)
probe("brms dm link_phi identity",
      brms::dirichlet_multinomial(link_phi = "identity")$link_phi)
probe("nb2 dharma (not polytomous) works",
      class(dharma_residuals(frm(yb ~ x, family = brmsfamily("negbinomial2"),
                                 data = d), nsim = 10)))
probe("dm dharma refused", dharma_residuals(frm(Y | trials(N) ~ x,
  family = dirichlet_multinomial(), data = d), nsim = 10))
