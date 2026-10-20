# Survey of every post-fit path that refits or re-evaluates the formula,
# on models with transformed terms. Arg: base | lane.
# Log: dev/drop1raw-out/survey-<arm>.txt (seed 5, n = 200).
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
res <- function(lab, expr) {
  r <- tryCatch({
    v <- q(expr)
    paste("ok", if (is.character(v)) paste(v, collapse = " | ") else "")
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 120)))
  cat(sprintf("%-46s %s\n", lab, r))
}
mk <- function() {
  set.seed(5)
  n <- 200
  d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                  g = factor(sample(letters[1:8], n, TRUE)),
                  w = runif(n, 0.5, 2))
  d$z <- 1 + d$x + log(d$x2) + rnorm(8)[d$g] + rnorm(n)
  d$cc <- ifelse(seq_len(n) %% 5 == 0, "interval", "none")
  d$z2 <- d$z + 1
  d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
  d
}
d <- mk()
rows <- function(t) rownames(t)

cat("\n== drop1() ==\n")
fs <- list(
  log = z ~ x + log(x2), scale = z ~ scale(x) + x2, I = z ~ I(x^2) + x2,
  ns = z ~ x + splines::ns(x2, 3), polyx = z ~ poly(x, 2) + x2,
  inter = z ~ x * log(x2), re = z ~ x + log(x2) + (1 | g),
  resp = log(z + 10) ~ x + x2, plain = z ~ x + x2,
  cens = z | cens(cc, z2) ~ x + x2, weights = z | weights(w) ~ x + log(x2),
  pois = cnt ~ x + log(x2))
for (nm in names(fs)) {
  fam <- if (nm == "pois") poisson() else gaussian()
  fit <- q(frm(bf(fs[[nm]]), family = fam, data = d))
  res(paste("drop1", nm, deparse1(fs[[nm]])), rows(drop1(fit)))
}
fit <- q(frm(bf(z ~ x + log(x2)), data = d))
res("drop1 scope = ~ x", rows(drop1(fit, scope = ~ x)))
res("drop1 scope = 'log(x2)'", rows(drop1(fit, scope = "log(x2)")))

cat("\n== other paths, data present ==\n")
fit <- q(frm(bf(z ~ x + log(x2) + x3), data = d))
res("update(. ~ . - x)", format(nobs(update(fit, . ~ . - x))))
res("update(. ~ . + x3:x)", format(nobs(update(fit, . ~ . + x3:x))))
res("add1(~ . + I(x3^2))", rows(add1(fit, ~ . + I(x3^2))))
res("step()", deparse1(formula(step(fit, trace = 0))))
fr <- q(frm(bf(z ~ x + log(x2) + (1 | g)), data = d, REML = TRUE))
fr2 <- q(frm(bf(z ~ log(x2) + (1 | g)), data = d, REML = TRUE))
res("anova(refit = TRUE)", format(nrow(anova(fr, fr2, refit = TRUE))))
res("influence() default data", format(nrow(influence(
  q(frm(bf(z ~ x + log(x2) + (1 | g)), data = d)), groups = "g")$fixed)))
res("confint(method = 'boot')", format(nrow(confint(fit, method = "boot",
                                                    nsim = 4, seed = 1))))
res("confint(method = 'profile')", format(nrow(confint(fit, parm = "b_x",
                                                       method = "profile"))))
res("frm_allfit()", format(length(frm_allfit(fit))))
res("emmeans(~ x)", format(nrow(as.data.frame(emmeans::emmeans(fit, ~ x)))))
res("marginaleffects avg_slopes()", format(nrow(
  marginaleffects::avg_slopes(fit))))
res("conditional_effects()", paste(names(conditional_effects(fit)),
                                   collapse = ","))
res("insight::get_data()", paste(names(insight::get_data(fit)),
                                 collapse = ","))
res("vcov_cluster(~ g) (g not in model)",
    format(nrow(vcov_cluster(fit, cluster = ~ g))))
dl <- list(d, transform(d, z = z + rnorm(nrow(d), 0, 0.1)))
res("frm_multiple() + anova", {
  m1 <- frm_multiple(bf(z ~ x + log(x2)), data = dl)
  m0 <- frm_multiple(bf(z ~ log(x2)), data = dl)
  format(nrow(anova(m1, m0)))
})

cat("\n== data removed from the workspace ==\n")
fit <- q(frm(bf(z ~ x + log(x2) + x3), data = d))
fitp <- q(frm(bf(z ~ x + x3), data = d))
fitr <- q(frm(bf(z ~ x + log(x2) + (1 | g)), data = d))
saveRDS(list(fit = fit, fitp = fitp, fitr = fitr),
        sprintf("dev/drop1raw-out/fits-%s.rds", arm))
d_keep <- d
rm(d)
res("drop1() log", rows(drop1(fit)))
res("drop1() plain", rows(drop1(fitp)))
res("update(. ~ . - x)", format(nobs(update(fit, . ~ . - x))))
res("update() plain (. ~ . - x)", format(nobs(update(fitp, . ~ . - x))))
res("add1(~ . + I(x3^2))", rows(add1(fit, ~ . + I(x3^2))))
res("step()", deparse1(formula(step(fit, trace = 0))))
res("influence() default data", format(nrow(influence(fitr,
                                                      groups = "g")$fixed)))
res("emmeans(~ x)", format(nrow(as.data.frame(emmeans::emmeans(fit, ~ x)))))
res("marginaleffects avg_slopes()", format(nrow(
  marginaleffects::avg_slopes(fit))))
res("conditional_effects()", paste(names(conditional_effects(fit)),
                                   collapse = ","))
res("vcov_cluster(~ g) (g not in model)",
    format(nrow(vcov_cluster(fit, cluster = ~ g))))
res("vcov_cluster(~ g) (g in model)",
    format(nrow(vcov_cluster(fitr, cluster = ~ g))))
d <- d_keep

cat("\n== NA rows: drop1() rows and the frozen basis ==\n")
dn <- d
dn$x3[c(3, 17, 40)] <- NA
fit <- q(frm(bf(z ~ splines::ns(x2, 3) + x3), data = dn))
res("drop1() with NA in x3", {
  t <- drop1(fit, test = "Chisq")
  sprintf("%s=%.10f", rownames(t), t$AIC)
})
l <- lm(z ~ splines::ns(x2, 3) + x3, data = dn)
res("lm drop1 AIC + const", {
  t <- drop1(l)
  sprintf("%s=%.10f", rownames(t),
          t$AIC + nobs(l) * (log(2 * pi) + 1) + 2)
})
