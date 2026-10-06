# Punch round 1: what conditional_effects(fit5, "income:age") holds on
# a seed where printing it fails (seed 7 of dev/vigport-pr1-mo.R).
#   Rscript dev/vigport-pr1-mo7.R [lib] [seed]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
s <- if (length(args) > 1) as.integer(args[2]) else 7L
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "seed", s, "\n")
set.seed(s)
income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(income_options, 100, TRUE),
                 levels = income_options, ordered = TRUE)
mean_ls <- c(30, 60, 70, 75)
ls <- mean_ls[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)
dat$age <- rnorm(100, mean = 40, sd = 10)
fit5 <- frm(ls ~ mo(income) * age, data = dat)
cat("convergence", fit5$opt$convergence, "\n")
print(fixef(fit5))
ce <- conditional_effects(fit5, "income:age")
d <- ce[[1]]
print(d[, intersect(c("income", "age", "effect2__", "estimate__", "se__",
                      "lower__", "upper__"), names(d))])
cat("non-finite estimate/se/lower/upper:",
    sum(!is.finite(d$estimate__)), sum(!is.finite(d$se__)),
    sum(!is.finite(d$lower__)), sum(!is.finite(d$upper__)), "\n")
print(tryCatch(vcov(fit5), error = function(e) conditionMessage(e)))
