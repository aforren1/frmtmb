# Reviewer of lane sampleprior: a three-part mixture with a formula on
# one weight only (theta2 ~ x). brms then has no joint simplex; what does
# it put on theta3, and what does frm_sample() add?
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$mix
bfm <- brms::bf(y ~ 1, theta2 ~ x)
bfam <- brms::mixture(gaussian, gaussian, gaussian)
print(brms::default_prior(bfm, data = d, family = bfam))
sc <- strsplit(as.character(brms::stancode(bfm, data = d, family = bfam)),
               "\n")[[1]]
cat(grep("theta|lprior", sc, value = TRUE), sep = "\n")
ffam <- frmtmb::mixture(gaussian(), gaussian(), gaussian())
print(default_prior(bf(y ~ 1, theta2 ~ x), data = d, family = ffam,
                    route = "sample"))
od <- cap_obj(bf(y ~ 1, theta2 ~ x), d, ffam, reparameterize = FALSE)
of <- cap_obj(bf(y ~ 1, theta2 ~ x), d, ffam, prior = "flat",
              reparameterize = FALSE)
print(od$env$par)
p <- od$env$par
for (k in which(names(p) == "betad")) {
  q1 <- p; q2 <- p; q2[k] <- q2[k] + 0.7
  cat("betad[", k, "] +0.7: d logprior",
      logprior_at(od, of, q2) - logprior_at(od, of, q1), "\n")
}
