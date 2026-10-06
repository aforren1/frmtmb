# Punch round 1: the reviewer's unseeded spell pass failed one row the
# seeded passes run, brms_monotonic.14.2,
# conditional_effects(fit5, "income:age"), with "need finite 'ylim'
# values" when the result is printed. The data are random, so this
# sweeps seeds over brms_monotonic's own data code and counts failures.
#
#   Rscript dev/vigport-pr1-mo.R [lib] [nseeds]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
ns <- if (length(args) > 1) as.integer(args[2]) else 200L
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")
grDevices::pdf(NULL)
rows <- list()
for (s in seq_len(ns)) {
  set.seed(s)
  income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(income_options, 100, TRUE),
                   levels = income_options, ordered = TRUE)
  mean_ls <- c(30, 60, 70, 75)
  ls <- mean_ls[income] + rnorm(100, sd = 7)
  dat <- data.frame(income, ls)
  dat$age <- rnorm(100, mean = 40, sd = 10)
  w5 <- character()
  fit5 <- withCallingHandlers(frm(ls ~ mo(income) * age, data = dat),
                              warning = function(w) {
                                w5 <<- c(w5, conditionMessage(w))
                                invokeRestart("muffleWarning")
                              })
  fit1 <- suppressWarnings(frm(ls ~ mo(income), data = dat))
  nan5 <- any(!is.finite(fixef(fit5)[, "Est.Error"]))
  nan1 <- any(!is.finite(fixef(fit1)[, "Est.Error"]))
  r <- tryCatch({
    ce <- conditional_effects(fit5, "income:age")
    print(ce)
    "OK"
  }, error = function(e) conditionMessage(e))
  rows[[s]] <- data.frame(seed = s, ce_ok = r == "OK", msg = r,
                          nan_se_fit5 = nan5, nan_se_fit1 = nan1,
                          conv5 = fit5$opt$convergence,
                          warned5 = length(w5) > 0,
                          warn5 = if (length(w5)) substr(w5[1], 1, 90) else "")
}
X <- do.call(rbind, rows)
cat(sprintf("conditional_effects(fit5) print failed on %d of %d seeds\n",
            sum(!X$ce_ok), ns))
cat(sprintf("fit5 has a non-finite standard error on %d seeds; fit1 on %d\n",
            sum(X$nan_se_fit5), sum(X$nan_se_fit1)))
cat("failures with a non-finite fit5 SE:", sum(!X$ce_ok & X$nan_se_fit5),
    " without:", sum(!X$ce_ok & !X$nan_se_fit5), "\n")
cat("fit5 optimizer code nonzero on", sum(X$conv5 != 0), "seeds\n")
cat("fit5 warned on", sum(X$warned5), "seeds; of the non-finite-SE seeds,",
    sum(X$warned5 & X$nan_se_fit5), "warned\n")
cat("first warnings on non-finite-SE seeds:\n")
print(utils::head(table(X$warn5[X$nan_se_fit5]), 5))
cat("failing seeds:", paste(utils::head(X$seed[!X$ce_ok], 20),
                            collapse = ", "), "\n")
