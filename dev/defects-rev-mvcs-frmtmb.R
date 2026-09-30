# Reviewer of lane defects, recheck of B2: frmtmb's multivariate
# fitted(scale = "linear") against brms's (dev/defects-rev-mvcs-brms.R).
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/defects-rev-mvcs-data.R")
d <- mvcs_data()
br <- readRDS("dev/defects-rev-log/mvcs-brms.rds")
for (nm in names(mvcs_models)) {
  cat("=====", nm, "\n")
  f <- frm(eval(mvcs_models[[nm]]), data = d)
  a <- tryCatch(fitted(f, scale = "linear"),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(a)) { cat(a, "\n"); next }
  b <- br[[nm]]$fitted
  cat("frmtmb dim:", dim(a), " layers:", dimnames(a)[[3]], "\n")
  cat("brms   dim:", dim(b), " layers:", dimnames(b)[[3]], "\n")
  cat("identical dimnames:", identical(dimnames(a), dimnames(b)), "\n")
  if (identical(dim(a), dim(b))) {
    e <- a[, "Estimate", ] - b[, "Estimate", ]
    sc <- apply(abs(b[, "Estimate", ]), 2, max)
    cat("max |frmtmb - brms| Estimate per layer / layer max|brms|:",
        round(apply(abs(e), 2, max) / sc, 3), "\n")
    cat("Est.Error ratio frmtmb/brms, median per layer:",
        round(apply(a[, "Est.Error", ] / b[, "Est.Error", ], 2, median), 3), "\n")
  }
  fe <- fixef(f)
  cat("frmtmb fixef:\n"); print(round(fe[, 1], 3))
  cat("brms posterior means:\n"); print(round(br[[nm]]$means, 3))
  if (nm == "one_cs") {
    u <- fitted(f, scale = "linear", resp = "yo")
    cat("resp = 'yo' equals the stacked eta layers:",
        isTRUE(all.equal(unclass(u), unclass(a[, , c("eta1", "eta2")]),
                         check.attributes = FALSE)), "\n")
    bx <- fe["yo_x", 1]; c1 <- fe["yo_z[1]", 1]; c2 <- fe["yo_z[2]", 1]
    cat("max |eta1 - (b x + cs1 z)|:", max(abs(a[, "Estimate", "eta1"] - (bx * d$x + c1 * d$z))),
        " eta2:", max(abs(a[, "Estimate", "eta2"] - (bx * d$x + c2 * d$z))), "\n")
  }
}
