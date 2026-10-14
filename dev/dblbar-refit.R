# Lane dblbar: refit() on a (1 + f || g) fit keeps the one diag block
# and gives the fit of the same model on the new response. Data seed
# 20261013, new response seed 5. Output: dev/dblbar-log/refit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
dd <- dblbar_data()
fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
set.seed(5)
y2 <- dd$y + stats::rnorm(nrow(dd), 0, 0.3)
r <- refit(fit, newresp = y2)
d2 <- dd
d2$y <- y2
direct <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = d2)
sink(file.path(wt, "dev/dblbar-log/refit.txt"), split = TRUE)
cat("blocks:", vapply(r$frame$re_blocks, function(b) {
  paste0(b$covstruct, "[", paste(b$cnms, collapse = ","), "]")
}, ""), "\n")
cat("names:", grep("^sd_|^cor_", variables(r), value = TRUE), "\n")
cat("max relative difference refit vs direct fit:",
    format(max(abs(unlist(r$estimates) - unlist(direct$estimates)) /
                 pmax(abs(unlist(direct$estimates)), 1e-8)), digits = 3),
    "\n")
cat("logLik refit", format(as.numeric(logLik(r)), digits = 12), "direct",
    format(as.numeric(logLik(direct)), digits = 12), "\n")
sink()
