# Review of lane names: odds and ends on the lane build.
# Usage: Rscript dev/reviews/names-review-misc.R <base|lane>
a <- commandArgs(TRUE)
arm <- a[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
q <- function(e) tryCatch(suppressWarnings(suppressMessages(e)),
                          error = function(err) paste("ERROR:",
                                                      conditionMessage(err)))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(sprintf("g%02d", 1:10), each = 20)),
                f = factor(rep(c("fa", "fb"), length.out = n)))
d$y <- 1 + 0.5 * d$x + rnorm(10, 0, 0.5)[d$g] + rnorm(n)
f_nl <- q(frm(bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1 + (1 | g),
                 nl = TRUE), data = d))
cat("\n== nl coef names/cols\n"); print(lapply(coef(f_nl), colnames))
cat("== nl ranef cols\n"); print(lapply(unclass(ranef(f_nl)), colnames))
cat("== nl variables\n"); print(variables(f_nl))
f_d <- q(frm(bf(y ~ x + (1 | g), sigma ~ (1 | g)), data = d))
cat("== dist VarCorr names/rows\n")
print(lapply(VarCorr(f_d), function(v) rownames(v$sd)))
cat("== dist summary Multilevel block\n")
s <- capture.output(summary(f_d))
writeLines(s[grep("Multilevel", s):(grep("Regression|Population", s)[1] - 1)])
cat("== 0 + Intercept dist variables\n")
f0 <- q(frm(bf(y ~ 0 + Intercept + x, sigma ~ z), data = d))
print(if (is.character(f0)) f0 else variables(f0))
cat("== mo ML variables and hypothesis\n")
d$k <- sample(1:4, n, TRUE)
fm <- q(frm(y ~ mo(k) + (1 | g), data = d))
print(variables(fm))
print(q(hypothesis(fm, "simo_mok1[2] > simo_mok1[1]", class = NULL)$hypothesis))
cat("== two-cat mo\n")
d$k2 <- sample(1:2, n, TRUE)
fm2 <- q(frm(y ~ mo(k2), data = d))
print(variables(fm2))
cat("== mo in sigma\n")
fm3 <- q(frm(bf(y ~ 1, sigma ~ mo(k)), data = d))
print(if (is.character(fm3)) fm3 else variables(fm3))
cat("== ranef groups= and pars= on merged\n")
print(lapply(unclass(q(ranef(f_d, pars = "Intercept"))), colnames))
cat("== [[ by label, condVar attr\n")
r <- ranef(f_d, condVar = TRUE)
print(dim(attr(r$g, "condSD")))
print(dim(attr(r[["sigma: 1 | g"]], "condSD")))
cat("== frm_multiple ranef and pooled varcorr labels (||)\n")
dl <- lapply(1:2, function(i) { dd <- d; dd$y <- dd$y + rnorm(n, 0, .01); dd })
fmu <- q(frm_multiple(bf(y ~ x + (1 + f || g)), data = dl))
if (!is.character(fmu)) {
  print(names(ranef(fmu)))
  print(unique(fmu$pooled_varcorr$grp))
  print(unique(confint_varcorr(fmu$fits[[1]])$block))
}
cat("== frm_multiple summary Formula\n")
fmu2 <- q(frm_multiple(bf(y ~ x + (1 | g), sigma ~ z), data = dl))
if (!is.character(fmu2)) {
  s <- capture.output(summary(fmu2))
  writeLines(s[1:6])
}
