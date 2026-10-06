# Small repros behind the defects and the brms checks in
# dev/vigport-findings.md, one line of output each.
#
#   Rscript dev/vigport-repros.R [lib]
#
# brms is reached with brms:: and never attached; its side runs through
# brms::stancode() or brms::standata(), so nothing compiles.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")
show <- function(tag, e) cat(sprintf("%-52s %s\n", tag, tryCatch(
  {force(e); "ACCEPTED"}, error = function(err)
    paste("REFUSED:", substr(gsub("\n", " ", conditionMessage(err)), 1, 140)))))
set.seed(1)

## R1. stancode() / standata() / pp_mixture() on a fit have no refusal.
## With the brms namespace loaded (any brms:: call loads it) the first
## two reach brms's default method and say "Data must be specified".
d <- data.frame(x = rnorm(40)); d$y <- d$x + rnorm(40)
f <- frm(y ~ x, data = d)
show("R1 stancode(fit), brms not loaded", stancode(f))
show("R1 pp_mixture(fit)", pp_mixture(f))
invisible(loadNamespace("brms"))
show("R1 stancode(fit), brms loaded", stancode(f))
show("R1 standata(fit), brms loaded", standata(f))

## R2. plot() of a fit refuses brms's plot.brmsfit() arguments as
## unknown, and suggests `x` for `N`.
kidney <- brms::kidney
fk <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal())
grDevices::pdf(NULL)
show("R2 plot(fit, N = 2)", plot(fk, N = 2, ask = FALSE))
show("R2 plot(fit, variable = '^b', regex = TRUE)",
     plot(fk, variable = "^b", regex = TRUE))
invisible(grDevices::dev.off())

## R3. update() keeps a class-wide cor prior after the update removes
## the last correlation, and refuses; brms's update() refits. brms also
## refuses the prior when it is written on such a model directly.
pr <- c(set_prior("normal(0,5)", class = "b"),
        set_prior("cauchy(0,2)", class = "sd"),
        set_prior("lkj(2)", class = "cor"))
f1 <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal(), prior = pr)
show("R3 update(fit, drop the correlation)",
     update(f1, formula. = ~ . - (1 + age | patient) + (1 | patient)))
show("R3 brms: lkj prior on (1 | patient), direct",
     brms::stancode(time | cens(censored) ~ age * sex + disease +
                      (1 | patient), data = kidney,
                    family = brms::lognormal(),
                    prior = brms::set_prior("lkj(2)", class = "cor")))

## R4. fixef() has no method for a frm_multiple() result, and the
## pooled table names the intercept "(Intercept)".
data("nhanes", package = "mice")
imp <- mice::mice(nhanes, m = 3, print = FALSE, seed = 1)
fm <- frm_multiple(bmi ~ age * chl, data = imp)
show("R4 fixef(frm_multiple)", fixef(fm))
cat(sprintf("%-52s %s\n", "R4 rownames(x$pooled)",
            paste(rownames(fm$pooled), collapse = ", ")))

## R5. hurdle_cumulative() fits equidistant and sum-to-zero thresholds,
## which vignettes/brms-migration.Rmd (lines 42 and 277) says it does
## not take.
dd <- data.frame(x = rnorm(400))
z <- cut(dd$x + rlogis(400), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
dd$y <- ifelse(runif(400) < 0.3, 0L, z)
for (th in c("equidistant", "sum_to_zero")) {
  show(paste0("R5 hurdle_cumulative(threshold = '", th, "')"),
       frm(y ~ x, data = dd, family = hurdle_cumulative(threshold = th)))
}

## R6. Two spellings frmtmb refuses "as brms does", checked on brms.
loss <- brms::loss
nl <- brms::bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
               ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE)
pn <- brms::prior(normal(5000, 1000), nlpar = "ult") +
  brms::prior(normal(1, 2), nlpar = "omega") +
  brms::prior(normal(45, 10), nlpar = "theta")
show("R6 brms: nonlinear, class sd without nlpar",
     brms::stancode(nl, data = loss,
                    prior = pn + brms::set_prior("cauchy(0,2)", class = "sd")))
show("R6 brms: y ~ x, class Intercept, dpar sigma",
     brms::stancode(y ~ x, data = d,
                    prior = brms::set_prior("normal(0, 1)",
                                            class = "Intercept",
                                            dpar = "sigma")))
show("R6 brms: a hypothesis with no relation",
     brms::hypothesis(data.frame(x = rnorm(50)), "x"))

## R7. brms 2.23.0 drops brm(threshold =) into its dots: the program has
## no delta. The family argument is the spelling that works in both.
inhaler <- brms::inhaler
sc <- brms::stancode(rating ~ period + carry + cs(treat) + (1 | subject),
                     data = inhaler, family = brms::sratio(),
                     threshold = "equidistant")
cat(sprintf("%-52s %s\n", "R7 brms: delta in program, threshold = argument",
            grepl("delta", sc)))
show("R7 frmtmb: sratio(threshold = 'equidistant')",
     frm(rating ~ period + carry + cs(treat) + (1 | subject),
         data = inhaler, family = sratio(threshold = "equidistant")))
