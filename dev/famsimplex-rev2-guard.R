# Reviewer round 2: the location-constant refusal, guard-absent side.
#   Rscript famsimplex-rev2-guard.R lane|base
# Output: dev/famsimplex-rev2-guard-<arm>.log (diffed by the reviewer)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.latent)
})
cat("arm:", arm, " lib:", find.package("frmtmb"), find.package("frmtmb.eam"),
    find.package("frmtmb.latent"), "\n")
src <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R"
if (file.exists(src)) source(src)
run <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(suppressMessages(suppressWarnings(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (inherits(r, "frmtmb_fit")) {
    cat(sprintf("logLik %.10f  npar %d\n", as.numeric(logLik(r)),
                length(r$opt$par)))
    cat("par:", sprintf("%.8f", r$opt$par), "\n")
  } else cat(substr(r, 1, 400), "\n")
}
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
cl <- rbinom(n, 1, 0.4)
d$y <- ifelse(cl == 1, 2 + 1 * d$x, -1 + 0.5 * d$z) + rnorm(n, 0, 0.7)
E <- cbind(0, 0.3 + 0.5 * d$z, -0.2 + 0.6 * d$x)
P <- exp(E) / rowSums(exp(E))
d$f <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample(3, 1, prob = p))])
d$y1 <- 1 + d$x + rnorm(n)
d$y2 <- 2 - d$x + rnorm(n, 0, 2)
d$yo <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), ordered_result = TRUE)
d$ynl <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.3)
mix <- mixture(gaussian(), gaussian())
run("mixture theta1 = 0.5", frm(bf(y ~ x, theta1 = 0.5), family = mix, data = d))
run("mixture theta1 ~ x", frm(bf(y ~ x, theta1 ~ x), family = mix, data = d))
run("mixture sigma2 = 1", frm(bf(y ~ x, sigma2 = 1), family = mix, data = d))
run("mixture sigma1 = 'sigma2'",
    frm(bf(y ~ x, sigma1 = "sigma2"), family = mix, data = d))
run("mixture mu2 = 'mu1'", frm(bf(y ~ x, mu2 = "mu1"), family = mix, data = d))
run("mixture mu1 = 0, mu2 = 1, theta1 = 0.5 (names)",
    frm(bf(y ~ x, mu1 = 0, mu2 = 1, theta1 = 0.5), family = mix, data = d))
run("mixture mu2 ~ 1", frm(bf(y ~ x, mu2 ~ 1), family = mix, data = d))
run("categorical muc ~ 1", frm(bf(f ~ x, muc ~ 1), family = categorical(),
                               data = d))
run("mvbf gaussian sigma = 1 + gaussian",
    frm(mvbf(bf(y1 ~ x, sigma = 1), bf(y2 ~ x)), data = d))
run("mvbf categorical + gaussian sigma = 2",
    frm(mvbf(bf(f ~ x, family = categorical()), bf(y2 ~ x, sigma = 2)),
        data = d))
run("mvbf categorical muc = 0 (location)",
    frm(mvbf(bf(f ~ x, muc = 0, family = categorical()), bf(y2 ~ x)),
        data = d))
run("nl with sigma constant",
    frm(bf(ynl ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE, sigma = 0.3),
        data = d, family = gaussian(),
        control = frmtmb_control(start = list())))
run("nl with constant nlpar",
    frm(bf(ynl ~ a * exp(b * x), a ~ 1, b = 0.3, nl = TRUE), data = d,
        family = gaussian()))
run("cumulative disc = 1", frm(bf(yo ~ x, disc = 1), family = cumulative(),
                               data = d))
run("gaussian sigma = 1", frm(bf(y1 ~ x, sigma = 1), family = gaussian(),
                              data = d))
if (exists("fs_data_dirichlet")) {
  dd <- fs_data_dirichlet(seed = 11)
  run("dirichlet phi = 5", frm(bf(Y ~ x, phi = 5), family = dirichlet(),
                               data = dd))
  dl <- fs_data_logistic_normal(seed = 13)
  run("logistic_normal sigmac = 1",
      frm(bf(Y ~ x, sigmac = 1), family = logistic_normal(), data = dl))
}
# frmtmb.eam
dw <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.25)
run("wiener bias = 0.5", frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                             family = wiener(), data = dw))
run("wiener bias = 0.5, ndt = 0.2",
    frm(bf(rt | dec(upper) ~ 1, bias = 0.5, ndt = 0.2), family = wiener(),
        data = dw))
dlb <- lba_simulate(300, v = c(2.2, 1.4, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
run("lba(3)", frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = dlb))
run("lba(3) A = 0.5", frm(bf(rt | vint(choice) ~ 1, A = 0.5), family = lba(3),
                          data = dlb))
run("lba(3) ndt = 0.2", frm(bf(rt | vint(choice) ~ 1, ndt = 0.2),
                            family = lba(3), data = dlb))
run("lba(3) v3 = 0.8 (location)",
    frm(bf(rt | vint(choice) ~ 1, v3 = 0.8), family = lba(3), data = dlb))
run("lba(3) v1 = 2.2 (location)",
    frm(bf(rt | vint(choice) ~ 1, v1 = 2.2), family = lba(3), data = dlb))
drd <- rdm_simulate(300, v = c(2.2, 1.4), A = 0.5, k = 0.5, ndt = 0.2)
run("rdm(2)", frm(bf(rt | vint(choice) ~ 1), family = rdm(2), data = drd))
run("rdm(2) A = 0.5", frm(bf(rt | vint(choice) ~ 1, A = 0.5),
                          family = rdm(2), data = drd))
run("rdm(2) v2 = 1.4 (location)", frm(bf(rt | vint(choice) ~ 1, v2 = 1.4),
                                      family = rdm(2), data = drd))
# frmtmb.latent
set.seed(4501)
dh <- do.call(rbind, lapply(1:10, function(id) {
  s <- integer(20)
  s[1] <- 1L
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)
  for (t in 2:20) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
  data.frame(id = id, t = 1:20, y = rnorm(20, c(0, 3)[s], 0.6))
}))
hf <- hmm(K = 2, gaussian(), time = t, group = id)
run("hmm", frm(bf(y ~ 1), family = hf, data = dh))
run("hmm sigma1 = 0.6", frm(bf(y ~ 1, sigma1 = 0.6), family = hf, data = dh))
run("hmm mu1 = 0 (location)", frm(bf(y ~ 1, mu1 = 0), family = hf, data = dh))
run("hmm mu2 = 3 (location)", frm(bf(y ~ 1, mu2 = 3), family = hf, data = dh))
cat("\nhmm primaries:", hf$primary_dpars, "\n")
cat("lba primaries:", lba(3)$primary_dpars, "\n")
