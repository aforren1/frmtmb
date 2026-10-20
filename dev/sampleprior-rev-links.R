# Reviewer of lane sampleprior: non-default dpar links. A softplus sigma
# without and with a formula (the Jacobian and the scale), and an
# identity-link nu, whose brms bound is 1 (brms: real<lower=1> nu).
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
lt <- function(v, df, mu, s) stats::dt((v - mu) / s, df, log = TRUE) - log(s)
mx <- mean(x)
sp <- function(e) log1p(exp(e))
cmp <- function(lab, f, bfun, pts) {
  od <- cap_obj(f[[1]], f[[2]], f[[3]], reparameterize = FALSE)
  of <- cap_obj(f[[1]], f[[2]], f[[3]], prior = "flat",
                reparameterize = FALSE)
  fl <- vapply(pts, function(q) logprior_at(od, of, q), 1)
  bl <- vapply(pts, bfun, 1)
  cat(sprintf("%-26s max |d frm - d brms| %.3g\n", lab,
              max(abs((fl - fl[1]) - (bl - bl[1])))))
  invisible(list(od = od, of = of))
}
# brms: sigma without a formula is on the identity, scale mad(y) = 9.9;
# with sigma ~ x on softplus, def_scale_prior on a link outside its list
# keeps 2.5 (stancode in dev/sampleprior-rev-log/links-stancode.txt)
for (fb in list(list(y ~ x, "nof"), list(brms::bf(y ~ x, sigma ~ x), "x"))) {
  sc <- brms::stancode(fb[[1]], data = D$wide,
                       family = brms::brmsfamily("gaussian",
                                                 link_sigma = "softplus"))
  s <- strsplit(as.character(sc), "\n")[[1]]
  cat("brms softplus", fb[[2]], ":", grep("lprior \\+=", s, value = TRUE),
      "\n")
}
fam <- frmtmb::brmsfamily("gaussian", link_sigma = "softplus")
set.seed(3)
p <- c(27, 2.4, 2.3)
pts <- lapply(1:5, function(k) p + rnorm(3, 0, 0.4))
cmp("softplus sigma, no formula", list(y ~ x, D$wide, fam),
    function(q) lt(q[1] + q[2] * mx, 3, 33.4, 9.9) +
      lt(sp(q[3]), 3, 0, 9.9) + log(stats::plogis(q[3])), pts)
p <- c(27, 2.4, 2.3, 0)
pts <- lapply(1:5, function(k) p + rnorm(4, 0, 0.4))
cmp("softplus sigma ~ x", list(bf(y ~ x, sigma ~ x), D$wide, fam),
    function(q) lt(q[1] + q[2] * mx, 3, 33.4, 9.9) +
      lt(q[3] + q[4] * mx, 3, 0, 2.5), pts)
# identity nu: inside brms's support the density is gamma(2, 0.1) with
# no Jacobian; below 1 brms has none
famn <- frmtmb::brmsfamily("student", link_nu = "identity")
p <- c(2.4, 0.6, 0.2, 4)
pts <- lapply(1:5, function(k) { q <- p + rnorm(4, 0, 0.3); q[4] <- 2 + k; q })
r <- cmp("identity nu, nu > 1", list(y ~ x, D$real, famn),
         function(q) lt(q[1] + q[2] * mx, 3, 4.1, 2.5) +
           lt(exp(q[3]), 3, 0, 2.5) + q[3] +
           stats::dgamma(q[4], 2, 0.1, log = TRUE), pts)
for (nu in c(3, 1.01, 0.9, 0.5, -0.5)) {
  q <- p; q[4] <- nu
  cat(sprintf("identity nu = %5.2f: default-prior nll %s, flat nll %s\n", nu,
              format(r$od$env$f(q)), format(r$of$env$f(q))))
}
sc <- brms::stancode(y ~ x, data = D$real,
                     family = brms::brmsfamily("student", link_nu = "identity"))
s <- strsplit(as.character(sc), "\n")[[1]]
cat("brms identity nu:", grep("nu;|nu \\|", s, value = TRUE), "\n")
