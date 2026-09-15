## Item 3 of the review: is the bias real, and is it a WRONG ANSWER or
## a DIFFERENT QUESTION?
##
## Every rung of a replicate saw the SAME data at the SAME n, so the
## rungs can be differenced seed by seed. That is far more powerful
## than comparing five means against 0.5, and it also settles the
## attribution: a finite-sample-n explanation cannot bias two rungs and
## leave three alone when all five have the same n.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("dev/cohrev-read.R")

rungs <- c("cond", "smooth", "id", "idcond", "full")
m <- read_kv("dev/coh-recovery-main.tsv")
tb <- table(m$seed)
d <- m[as.character(m$seed) %in% names(tb)[tb == 5L], ]
w <- reshape(d[, c("seed", "rung", "est")], idvar = "seed",
             timevar = "rung", direction = "wide")
names(w) <- sub("^est\\.", "", names(w))

pair <- function(a, b) {
  x <- w[[a]] - w[[b]]
  sprintf("%-16s %+.5f  se %.5f  t %6.2f", paste0(a, " - ", b),
          mean(x), sd(x) / sqrt(length(x)),
          mean(x) / (sd(x) / sqrt(length(x))))
}
cat("==== paired differences between rungs, 148 replicates ====\n")
for (p in list(c("cond", "full"), c("smooth", "full"),
               c("id", "full"), c("idcond", "full"),
               c("cond", "smooth"), c("smooth", "id"),
               c("id", "idcond"))) {
  cat(pair(p[1L], p[2L]), "\n")
}
cat("\nagainst the truth, 0.5, unpaired:\n")
for (rg in rungs) {
  x <- w[[rg]] - 0.5
  cat(sprintf("%-8s %+.5f  se %.5f  t %6.2f\n", rg, mean(x),
              sd(x) / sqrt(length(x)), mean(x) / (sd(x) / sqrt(length(x)))))
}

## How much variation in the linear predictor does each rung leave
## unmodeled? Under a nonlinear link an omitted, balanced source of
## variation attenuates the contrast: the estimator is consistent for a
## MARGINAL contrast rather than the conditional 0.5 the truth carries.
## The logit/probit approximation puts that attenuation at
## 1 / sqrt(1 + 0.346 V) with V the omitted variance.
fr <- seq_len(60L) / 60
bump <- 0.8 * exp(-0.5 * ((fr - 0.35) / 0.12)^2)
vb <- mean((bump - mean(bump))^2)
cat(sprintf("\nvariance of the frequency bump over the grid: %.5f\n", vb))
vv <- c(cond = vb + 0.35^2 + 0.20^2, smooth = 0.35^2 + 0.20^2,
        id = 0.20^2, idcond = 0, full = 0)
cat(sprintf("%-8s omitted V %7.4f  predicted contrast %.4f  predicted bias %+.5f  observed %+.5f\n",
            names(vv), vv, 0.5 / sqrt(1 + 0.346 * vv),
            0.5 / sqrt(1 + 0.346 * vv) - 0.5,
            vapply(names(vv), function(r) mean(w[[r]]) - 0.5,
                   numeric(1))))
cat("\nobserved MINUS the correct model's own offset:\n")
cat(sprintf("%-8s %+.5f  against predicted %+.5f\n", names(vv),
            vapply(names(vv), function(r) mean(w[[r]] - w$full),
                   numeric(1)),
            0.5 / sqrt(1 + 0.346 * vv) - 0.5))
