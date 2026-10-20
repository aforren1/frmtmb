# Reviewer of lane sampleprior: rtmb-pitfalls item 22 in the DEFAULT
# densities. beta(1, 1) and gamma(1, rate) carry (shape - 1) * log(x),
# a constant-zero operand when the shape is 1, which the tape folds to
# 0: on doubles the density is NaN outside the support, on the tape it
# is finite. Harmless under a link that maps into the support; with an
# identity link the default no longer rejects a value brms's bound does.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
core <- asNamespace("frmtmb")
for (s in c("beta(1, 1)", "gamma(1, 1)", "gamma(2, 1)")) {
  dd <- core$parse_prior_dist(s)
  v <- if (grepl("beta", s)) c(-0.2, 1.2) else -0.5
  tp <- RTMB::MakeTape(function(z) core$prior_base_logdens(z, dd), 0)
  cat(sprintf("%-11s at %s: doubles %s, tape %s\n", s,
              paste(v, collapse = ", "),
              paste(format(sapply(v, core$prior_base_logdens, dist = dd)),
                    collapse = ", "),
              paste(format(sapply(v, function(z) tp(z))), collapse = ", ")))
}
# the sampler's own objective, with link_zi = "identity"
fam <- tryCatch(frmtmb::brmsfamily("zero_inflated_poisson",
                                   link_zi = "identity"),
                error = function(e) NULL)
if (!is.null(fam)) {
  od <- cap_obj(y ~ x, D$zcnt, fam)
  of <- cap_obj(y ~ x, D$zcnt, fam, prior = "flat")
  print(od$env$par)
  k <- which(names(od$env$par) == "betad")[1]
  for (zi in c(0.3, 0.01, -0.01, -0.05, 1.05)) {
    q <- od$env$par; q[k] <- zi
    cat(sprintf("zip identity zi = %5.2f: default nll %s, flat nll %s\n", zi,
                format(od$env$f(q)), format(of$env$f(q))))
  }
  sc <- strsplit(as.character(brms::stancode(y ~ x, data = D$zcnt,
    family = brms::brmsfamily("zero_inflated_poisson", link_zi = "identity"))),
    "\n")[[1]]
  cat("brms:", grep("zi;|zi \\|", sc, value = TRUE), "\n")
}
# wiener bs under an identity link, gamma(1, 1)
famw <- tryCatch(frmtmb.eam::wiener(link_bs = "identity"),
                 error = function(e) conditionMessage(e))
if (inherits(famw, "family") || is.list(famw)) {
  od <- cap_obj(rt | dec(upper) ~ x, D$wien, famw)
  of <- cap_obj(rt | dec(upper) ~ x, D$wien, famw, prior = "flat")
  print(od$env$par)
  k <- which(names(od$env$par) == "betad")[1]
  for (bs in c(1.4, 0.5, -0.5)) {
    q <- od$env$par; q[k] <- bs
    cat(sprintf("wiener identity bs = %4.1f: default nll %s, flat nll %s\n",
                bs, format(od$env$f(q)), format(of$env$f(q))))
  }
} else cat("wiener(link_bs = 'identity'):", famw, "\n")
# a user's beta(1, 1) on a coefficient: no link, no likelihood guard
od <- cap_obj(y ~ x, D$unit, frmtmb::brmsfamily("beta"),
              prior = set_prior("beta(1, 1)", class = "b", coef = "x"))
for (b in c(0.5, 2, -1)) {
  q <- od$env$par; q[2] <- b
  cat(sprintf("beta(1, 1) on b_x, b_x = %4.1f: sampler nll %s\n", b,
              format(od$env$f(q))))
}
