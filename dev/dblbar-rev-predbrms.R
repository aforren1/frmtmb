# Reviewer, lane dblbar: frmtmb's predictions of (1 + f || g) on the rows
# dev/dblbar-rev-brms-refo.R gives brms, for a side-by-side reading;
# plus frmtmb.sample posterior_epred() at old levels, || against three
# blocks. Data seed 20261013; sampler seed 11.
# Output: dev/dblbar-rev-log/predbrms.txt
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
sink(file.path(wt, "dev/dblbar-rev-log/predbrms.txt"), split = TRUE)
dd <- dblbar_data()
dd$fb <- as.numeric(dd$f == "b")
dd$fc <- as.numeric(dd$f == "c")
A <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
B <- frm(bf(y ~ f + (1 | g) + (0 + fb | g) + (0 + fc | g)),
         family = gaussian(), data = dd)
nd_old <- dd[c(1, 5, 9, 40), ]
s <- function(m) paste(sprintf("%.3f", m[, 1]), collapse = " ")
cat("full       ", s(fitted(A, newdata = nd_old)), "\n")
cat("(1 | g)    ", s(fitted(A, newdata = nd_old, re_formula = ~ (1 | g))),
    "\n")
cat("NA         ", s(fitted(A, newdata = nd_old, re_formula = NA)), "\n")
nd_new <- nd_old
nd_new$g <- c("NEW1", "NEW2", "NEW1", "NEW3")
fg <- fitted(A, newdata = nd_new, allow_new_levels = TRUE,
             sample_new_levels = "gaussian")
cat("new gauss  ", s(fg), " Err", sprintf("%.3f", fg[, 2]), "\n")
dA <- frm_sample(A, chains = 4, iter = 1500, seed = 11, refresh = 0)
dB <- frm_sample(B, chains = 4, iter = 1500, seed = 11, refresh = 0)
ea <- posterior_epred(dA, newdata = nd_old)
eb <- posterior_epred(dB, newdata = nd_old)
cat("draws epred old A", sprintf("%.3f", colMeans(ea)), "sd",
    sprintf("%.3f", apply(ea, 2, sd)), "\n")
cat("draws epred old B", sprintf("%.3f", colMeans(eb)), "sd",
    sprintf("%.3f", apply(eb, 2, sd)), "\n")
ra <- posterior_epred(dA, newdata = nd_old, re_formula = ~ (1 | g))
rb <- posterior_epred(dB, newdata = nd_old, re_formula = ~ (1 | g))
cat("draws epred (1|g) A", sprintf("%.3f", colMeans(ra)), "\n")
cat("draws epred (1|g) B", sprintf("%.3f", colMeans(rb)), "\n")
set.seed(9)
pa <- posterior_predict(dA, newdata = nd_old)
cat("posterior_predict A dim", dim(pa), "means", sprintf("%.3f",
                                                          colMeans(pa)), "\n")
sink()
