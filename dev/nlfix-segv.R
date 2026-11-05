# The segfault of spec ord_hurdle_disc: which step?
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(frmtmb))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
d <- E$eq_dat("ordinal", seed = 4952)
d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
cat("frm linear\n")
f0 <- suppressWarnings(frm(bf(y ~ x, disc ~ 0 + z, hu ~ z),
                           family = hurdle_cumulative(), data = d))
print(logLik(f0))
cat("frm nl\n")
f1 <- suppressWarnings(frm(bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, hu ~ z, nl = TRUE),
                           family = hurdle_cumulative(), data = d))
print(logLik(f1))
cat("brms stancode, linear\n")
sc <- suppressMessages(brms::stancode(brms::bf(y ~ x, disc ~ 0 + z, hu ~ z), data = d,
                     family = brms::hurdle_cumulative()))
cat(nchar(sc), "\n")
cat("brms stancode, nl\n")
sc <- suppressMessages(brms::stancode(brms::bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, hu ~ z, nl = TRUE),
                     data = d, family = brms::hurdle_cumulative()))
cat(nchar(sc), "\n")
cat(grep("hurdle_cumulative|target \\+=", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
