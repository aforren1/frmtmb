# Reviewer, item 2: predict() and frmtmb.sample's posterior_predict()
# on rows whose cs() thresholds cross, cumulative() (new at 0.68.0) and
# hurdle_cumulative() (lane ordmix). Uses the fit of dev/relrev-cs.R.
#   Rscript dev/relrev-cs2.R > dev/relrev-log/cs2.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
S <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev-log/cs-fit.rds")
fit <- S$fit; nd <- S$nd; crossing <- S$crossing
cat("crossing rows:", which(crossing), "\n")
set.seed(11)
ps <- predict(fit, newdata = nd, summary = FALSE)
cat("predict(summary = FALSE) class", class(ps), "dim", dim(ps), "\n")
cat("NA replicates per row:", colSums(is.na(ps)), "of", nrow(ps), "\n")
set.seed(11)
pp <- predict(fit, newdata = nd)
cat("predict() proportions on the crossing rows (from the non-NA replicates):\n")
print(round(pp[crossing, , drop = FALSE], 4))
cat("\n== hurdle_cumulative (lane ordmix), the same construction\n")
d <- fit$data
d$yh <- ifelse(runif(nrow(d)) < 0.2, 0L, d$y)
fh <- suppressWarnings(frm(yh ~ cs(x), family = hurdle_cumulative(), data = d))
lp <- Filter(function(l) identical(l$dpar, "mu"), fh$frame$linpreds)[[1L]]
b <- fh$estimates[[lp$cs[[1L]]$par]]
tau <- frmtmb:::ord_threshold_values(family(fh), fh$estimates$tau_raw)
xc <- (tau[2] - tau[1]) / (b[2] - b[1])
ndh <- data.frame(x = c(0, 3 * xc))
cat("hurdle crossing x", xc, "\n")
set.seed(11)
psh <- predict(fh, newdata = ndh, summary = FALSE)
cat("hurdle NA replicates per row:", colSums(is.na(psh)), "of", nrow(psh), "\n")
print(round(predict(fh, newdata = ndh), 4))
cat("hurdle fitted NaN per row:", rowSums(is.nan(fitted(fh, newdata = ndh)[, "Estimate", ])), "\n")

cat("\n== frmtmb.sample on the cumulative() cs() fit\n")
dr <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400, warmup = 200,
                                  seed = 7, refresh = 0))
pe <- posterior_epred(dr, newdata = nd)
cat("posterior_epred class", class(pe), "dim", dim(pe), "\n")
cat("epred NaN draws per row (any category):",
    apply(pe, 2, function(m) sum(apply(is.nan(m), 1, any))), "\n")
pq <- posterior_predict(dr, newdata = nd)
cat("posterior_predict class", class(pq), "typeof", typeof(pq), "dim", dim(pq), "\n")
cat("posterior_predict NA draws per row:", colSums(is.na(pq)), "of", nrow(pq), "\n")
# per draw: does row r cross at that draw?
th <- as.matrix(dr)
cat("draw columns:", paste(head(colnames(th), 12), collapse = " "), "\n")
