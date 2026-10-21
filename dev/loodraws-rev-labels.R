# Reviewer of lane loodraws: frmtmb's draw labels for the random-intercept
# and mvbf models, against brms's variable names.
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
B <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/loodraws-rev-out/brmsfits.rds")
dd <- B$dd
fit <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
print(frmtmb::brms_par_labels(fit))
print(str(fit$frame$par_template))
print(head(colnames(as.matrix(B$fits$ri)), 20))
fm <- suppressWarnings(frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE), data = dd))
print(frmtmb::brms_par_labels(fm))
print(grep("^b_|^sigma", colnames(as.matrix(B$fits$mv)), value = TRUE))
x <- frmtmb.sample:::draws_to_natural
print(args(x))
