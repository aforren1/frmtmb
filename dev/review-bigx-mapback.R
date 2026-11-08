# Everything read from a fit on a shifted design must be in natural
# coordinates: compare with the same model on the centered column,
# transformed exactly. Poisson (1|g), m = 1e6, data seed 2.
# REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
m <- 1e6
set.seed(2)
n <- 150L; G <- 15L
g <- factor(rep(seq_len(G), length.out = n))
z <- rnorm(n)
eta <- 0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]
d <- data.frame(y = rpois(n, exp(eta)), z = z, w = z + m, g = g)
fz <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
rw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d)
fw <- rw$f
cat("warn:", short_w(rw), " basis:", basis_on(fw), "\n")
A <- matrix(c(1, 0, -m, 1), 2)  # beta_w = A beta_z
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
say <- function(lbl, v) cat(sprintf("%-38s %s\n", lbl, format(v, digits = 4)))

bz <- quiet(fixef(fz))[, 1]; bw <- quiet(fixef(fw))[, 1]
say("fixef: rel diff vs A b_z", rel(bw, drop(A %*% bz)))
Vz <- quiet(vcov(fz)); Vw <- quiet(vcov(fw))
say("vcov: rel diff vs A V_z A'", rel(Vw, A %*% Vz %*% t(A)))
say("vcov slope rel diff", rel(Vw[2, 2], Vz[2, 2]))
cw <- quiet(confint(fw, method = "wald")); cz <- quiet(confint(fz, method = "wald"))
print(cw); print(cz)
pw <- tryCatch(withCallingHandlers(confint(fw, parm = "w", method = "profile"), warning = function(x) { cat("  profile w warning:", conditionMessage(x), "
"); invokeRestart("muffleWarning") }),
               error = function(e) e)
pz <- tryCatch(quiet(confint(fz, parm = "z", method = "profile")),
               error = function(e) e)
cat("profile w:\n"); print(pw); cat("profile z:\n"); print(pz)
nd <- data.frame(z = c(-1, 0, 1.5), g = factor(c(1, 2, 3)))
nd$w <- nd$z + m
say("predict(newdata, summary) mean rel", {
  set.seed(1); a <- quiet(predict(fw, newdata = nd, ndraws = 4000))[, 1]
  set.seed(1); b <- quiet(predict(fz, newdata = nd, ndraws = 4000))[, 1]
  rel(a, b)
})
say("fitted() rel", rel(quiet(fitted(fw))[, 1], quiet(fitted(fz))[, 1]))
say("fitted(newdata, linear) rel",
    rel(quiet(fitted(fw, newdata = nd, scale = "linear"))[, 1],
        quiet(fitted(fz, newdata = nd, scale = "linear"))[, 1]))
say("ranef rel", rel(as.numeric(quiet(ranef(fw))$g)[1:15], as.numeric(quiet(ranef(fz))$g)[1:15]))
sw <- capture.output(quiet(print(summary(fw))))
cat(sw[grepl("^\\s*(Intercept|w)\\b", sw)], sep = "\n")
say("logLik diff", ll_of(fw) - ll_of(fz))

# update(): a refit on new data, and the same data
d2 <- d; set.seed(9); d2$y <- rpois(n, exp(eta))
uw <- rfit_u <- tryCatch(quiet(update(fw, newdata = d2)), error = function(e) e)
uz <- quiet(update(fz, newdata = d2))
if (inherits(uw, "error")) cat("update error:", conditionMessage(uw), "\n") else {
  say("update(newdata) slope shift / SE",
      (quiet(fixef(uw))["w", 1] - quiet(fixef(uz))["z", 1]) /
        quiet(fixef(uz))["z", 2])
  say("update(newdata) dll", ll_of(uw) - ll_of(uz))
  say("update basis on", basis_on(uw))
}
# start echo: the fit's own start and a start from its estimates
st <- tryCatch(fw$start, error = function(e) NULL)
cat("fit$start:\n"); str(st)
fs <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
           start = list(beta = unname(bw)))
say("refit from its own estimates: dll", ll_of(fs$f) - ll_of(fw))
say("  its iterations", fs$f$opt$iterations)
cat("  warn:", short_w(fs), "\n")
# frm_allfit
af <- tryCatch(quiet(frm_allfit(fw)), error = function(e) e)
if (inherits(af, "error")) cat("allfit error:", conditionMessage(af), "\n") else {
  cat("allfit:\n")
  for (nm in names(af)) {
    o <- af[[nm]]
    if (inherits(o, "frmtmb_fit")) {
      cat(sprintf("  %-28s slope shift %+.2e SE dll %+.2e conv %s\n", nm,
                  (quiet(fixef(o))["w", 1] - bz[["z"]] ) /
                    quiet(fixef(fz))["z", 2],
                  ll_of(o) - ll_of(fz), o$opt$convergence))
    } else cat("  ", nm, class(o)[1], "\n")
  }
}
# saved fit
tf <- tempfile(fileext = ".rds")
saveRDS(fw, tf); fr <- readRDS(tf)
say("readRDS vcov rel", rel(quiet(vcov(fr)), Vw))
say("readRDS fitted(newdata) rel",
    rel(quiet(fitted(fr, newdata = nd, scale = "linear"))[, 1],
        quiet(fitted(fw, newdata = nd, scale = "linear"))[, 1]))
say("saved size MB", file.size(tf) / 2^20)
# priors on class b: slope prior identical under the shift (intercept
# flat in both), so the MAP fits must agree
for (sdp in c(1, 0.1)) {
  pr <- set_prior(sprintf("normal(0, %g)", sdp), class = "b")
  pz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d,
                  prior = pr))
  pw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d, prior = pr)
  cat(sprintf("prior normal(0,%g) on b: slope w %.6f z %.6f shift %+.2e SE; warn %s\n",
              sdp, quiet(fixef(pw$f))["w", 1], quiet(fixef(pz))["z", 1],
              (quiet(fixef(pw$f))["w", 1] - quiet(fixef(pz))["z", 1]) /
                quiet(fixef(pz))["z", 2], short_w(pw)))
}
# a prior on the intercept: brms puts it on the centered intercept
pr <- set_prior("normal(0, 1)", class = "Intercept")
pz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d, prior = pr))
pw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d, prior = pr)
cat(sprintf("prior normal(0,1) on Intercept: slope w %.6f z %.6f; warn %s; err %s\n",
            if (pw$err) NA else quiet(fixef(pw$f))["w", 1],
            quiet(fixef(pz))["z", 1], short_w(pw),
            if (pw$err) conditionMessage(pw$f) else "none"))
