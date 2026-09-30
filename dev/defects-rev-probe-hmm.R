# Reviewer of lane defects: frmtmb.latent's hmm() with se(), column
# present and column missing, on the lane core (latent from rellib-r3).
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("frmtmb from", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
set.seed(1)
dd <- data.frame(id = rep(1:10, each = 12), t = rep(1:12, 10),
                 y = rnorm(120), sdv = 1)
fam <- function() hmm(K = 2, gaussian(), time = t, group = id)
m <- function(e) tryCatch(e, error = function(err) conditionMessage(err))
cat("present:", substr(m(frm(bf(y | se(sdv) ~ 1), family = fam(), data = dd)), 1, 160), "\n")
cat("missing:", substr(m(frm(bf(y | se(nope) ~ 1), family = fam(), data = dd)), 1, 160), "\n")
