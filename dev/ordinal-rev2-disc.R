# Reviewer re-check (punch round 1), lane ordinal: B2. (1) A modeled disc,
# a disc fixed by the user at 2, and one fixed by the user at 1, on every
# reader. (2) The downstream readers of fixef(flatten = TRUE) and of the
# coefficient vector on a flexible fit, base against lane with identical().
# Usage: Rscript ... <base|lane>. Data seed 20261013, sampler seed 5.
# Output: dev/ordinal-rev2-log-disc-<arm>.txt and -<arm>.rds
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)
  library(emmeans)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261013)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(letters[1:6], n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x +
  rnorm(6, 0, 0.4)[as.integer(d$g)]
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
grab <- function(expr) tryCatch(suppressWarnings(expr), error = function(e)
  paste("ERROR:", conditionMessage(e)))
show <- function(lab, x) { cat("\n##", lab, "\n"); print(x) }

if (arm == "lane") {
  fits <- list(
    modeled = frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d),
    fixed2 = grab(frm(bf(y ~ x, disc = 2), family = sratio(), data = d)),
    fixed1 = grab(frm(bf(y ~ x, disc = 1), family = sratio(), data = d)),
    default = frm(y ~ x, family = sratio(), data = d))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    cat("\n==========", nm, "\n")
    if (is.character(f)) { cat(f, "\n"); next }
    s <- utils::capture.output(print(summary(f)))
    show("summary: Links / Fixed dpar / disc lines",
         grep("Links|Fixed dpar|disc", s, value = TRUE))
    show("summary()$fixed_dpars", summary(f)$fixed_dpars)
    show("fixef rows", rownames(fixef(f)))
    show("fixef(flatten = TRUE)", fixef(f, flatten = TRUE))
    show("coef names", names(coef(f)))
    show("confint rows", rownames(grab(confint(f))))
    show("variables", variables(f))
    show("hypothesis disc_z", grab(hypothesis(f, "disc_z = 0")$hypothesis[, 1:3]))
    show("emmeans dpar disc", grab(as.data.frame(emmeans(f, ~ x, dpar = "disc"))))
    show("frm_linpred(dpar = disc) head", grab(head(frm_linpred(f, dpar = "disc"))))
    ds <- grab(suppressMessages(frm_sample(f, chains = 1, iter = 100,
                                           refresh = 0, seed = 5)))
    show("draws variables", grab(variables(ds)))
    show("draws fixef rows", grab(rownames(fixef(ds))))
    show("draws summary disc lines",
         grab(grep("disc|Links", utils::capture.output(print(summary(ds))),
                   value = TRUE)))
  }
}

# downstream readers on a flexible fit, both arms, compared with identical()
f <- frm(y ~ x + z, family = cumulative(), data = d)
f0 <- frm(y ~ x, family = cumulative(), data = d)
fr <- frm(y ~ x + (1 | g), family = sratio(), data = d)
out <- list(
  fixef_flat = fixef(f, flatten = TRUE),
  boot = grab(frm_bootstrap(f, nsim = 4, seed = 3)$t),
  allfit = grab({a <- frm_allfit(f); utils::capture.output(print(a))}),
  influence = grab({i <- influence(fr); list(dfbeta = dfbeta(i),
                                            cooks = cooks.distance(i))}),
  cooks_fit = grab(cooks.distance(f)),
  drop1 = grab(as.data.frame(drop1(f, test = "Chisq"))),
  anova = grab(as.data.frame(anova(f0, f))),
  vcov_names = rownames(vcov(f)),
  confint = grab(confint(f)))
saveRDS(out, file.path(wt, paste0("dev/ordinal-rev2-disc-", arm, ".rds")))
if (arm == "lane") {
  show("allFit print (lane)", out$allfit)
  show("dfbeta colnames (lane)", colnames(out$influence$dfbeta))
}
