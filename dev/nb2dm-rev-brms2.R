# Reviewer: brms 2.23.0 standata() without trials(), has_trials(), and the
# mixture / mvbf default priors frmtmb.sample must match.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
suppressMessages(library(brms))
suppressMessages(library(frmtmb.sample))
cat("frmtmb.sample", dirname(system.file(package = "frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(6, 10))
dd$y <- rnbinom(n, mu = 3, size = 2)
dd$N <- sample(3:9, n, TRUE)
Y <- t(sapply(dd$N, function(k) rmultinom(1, k, c(0.3, 0.3, 0.4))))
colnames(Y) <- c("a", "b", "c")
dd$Y <- Y
dd$w <- rnorm(n)
probe("brms has_trials(dm)", brms:::has_trials(dirichlet_multinomial()))
probe("brms standata dm without trials",
      names(standata(Y ~ x, dd, family = brms::dirichlet_multinomial())))
probe("brms standata multinomial without trials",
      names(standata(Y ~ x, dd, family = brms::multinomial())))
rows <- function(tab) {
  tab <- as.data.frame(tab)
  stats::setNames(tab$prior, paste(tab$class, tab$coef, tab$group, tab$resp,
                                   tab$dpar, sep = "|"))
}
cmp <- function(fr, br, fam_f, fam_b) {
  ours <- rows(frmtmb::default_prior(fr, dd, family = fam_f,
                                     route = "sample"))
  theirs <- rows(brms::default_prior(br, dd, family = fam_b))
  dens <- theirs[nzchar(theirs)]
  data.frame(key = names(dens), brms = unname(dens),
             frmtmb = unname(ours[names(dens)]),
             same = unname(ours[names(dens)]) == unname(dens))
}
probe("default priors: mixture(nb2, poisson)",
      cmp(y ~ x, y ~ x,
          frmtmb::mixture(frmtmb::brmsfamily("negbinomial2"),
                          stats::poisson()),
          brms::mixture(brms::brmsfamily("negbinomial2"), poisson())))
probe("default priors: mixture(negbinomial, poisson)",
      cmp(y ~ x, y ~ x,
          frmtmb::mixture(frmtmb::negbinomial(), stats::poisson()),
          brms::mixture(brms::negbinomial(), poisson())))
probe("default priors: mvbf(nb2, gaussian)", {
  ours <- rows(frmtmb::default_prior(
    frmtmb::mvbf(frmtmb::bf(y ~ x, family = frmtmb::brmsfamily("negbinomial2")),
                 frmtmb::bf(w ~ x, family = stats::gaussian())),
    dd, route = "sample"))
  theirs <- rows(brms::default_prior(
    brms::bf(y ~ x, family = brms::brmsfamily("negbinomial2")) +
      brms::bf(w ~ x) + brms::set_rescor(FALSE), dd))
  dens <- theirs[nzchar(theirs)]
  data.frame(key = names(dens), brms = unname(dens),
             frmtmb = unname(ours[names(dens)]))
})
probe("default priors: mvbf(dm, gaussian)", {
  ours <- rows(frmtmb::default_prior(
    frmtmb::mvbf(frmtmb::bf(Y | trials(N) ~ x,
                            family = frmtmb::dirichlet_multinomial()),
                 frmtmb::bf(w ~ x, family = stats::gaussian())),
    dd, route = "sample"))
  theirs <- rows(brms::default_prior(
    brms::bf(Y | trials(N) ~ x, family = brms::dirichlet_multinomial()) +
      brms::bf(w ~ x) + brms::set_rescor(FALSE), dd))
  dens <- theirs[nzchar(theirs)]
  data.frame(key = names(dens), brms = unname(dens),
             frmtmb = unname(ours[names(dens)]))
})
probe("default priors: nb2 identity link_sigma",
      cmp(y ~ x, y ~ x,
          frmtmb::brmsfamily("negbinomial2", link_sigma = "identity"),
          brms::brmsfamily("negbinomial2", link_sigma = "identity")))
probe("default priors: nb2 softplus link_sigma",
      cmp(y ~ x, y ~ x,
          frmtmb::brmsfamily("negbinomial2", link_sigma = "softplus"),
          brms::brmsfamily("negbinomial2", link_sigma = "softplus")))
probe("default priors: dm softplus link_phi",
      cmp(Y | trials(N) ~ x, Y | trials(N) ~ x,
          frmtmb::dirichlet_multinomial(link_phi = "softplus"),
          brms::dirichlet_multinomial(link_phi = "softplus")))
probe("default priors: dm refcat NA",
      cmp(Y | trials(N) ~ x, Y | trials(N) ~ x,
          frmtmb::dirichlet_multinomial(refcat = NA),
          brms::dirichlet_multinomial(refcat = NA)))
