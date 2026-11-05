# Generate the counts the review quotes, from the logs. And the two consumers my
# first script called wrongly: predict(newdata) and car() with dimnames.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
root <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out"
rd <- function(d) {
  f <- list.files(file.path(root, d), "^[^-]*-[0-9]+[.]tsv$", full.names = TRUE)
  do.call(rbind, lapply(f, function(x) {
    if (grepl("\tERROR\t", readLines(x, 1L))) {
      cat("ERROR in", x, "\n")
      return(NULL)
    }
    read.delim(x, quote = "", colClasses = c(warn = "character"))
  }))
}
a1 <- rd("accuracy")
a2 <- rd("accuracy2")
com <- c("design", "seed", "rounds", "n_warn", "max_z_est", "max_z_lap", "newton",
         "se_ratio_min", "se_ratio_max", "ess_min")
m <- rd("multid")
scale_d <- c("scale1e5", "pois1e5", "year")
all <- rbind(a1[com], a2[!a2$design %in% scale_d, com], m[com])
cat("fits", nrow(all), "designs", length(unique(all$design)), "\n")
print(table(all$design))
print(table(all$rounds))
bad <- all$design == "re_sigma_small" & all$seed == 1
cat("excl re_sigma_small 1: max z", max(all$max_z_est[!bad], na.rm = TRUE),
    "newton", max(all$newton[!bad], na.rm = TRUE),
    "lap max z", max(all$max_z_lap[!bad], na.rm = TRUE), "\n")
cat("SE ratio range", range(all$se_ratio_min[!bad], na.rm = TRUE),
    range(all$se_ratio_max[!bad], na.rm = TRUE), "\n")
print(all[which(all$se_ratio_max > 1.08 | all$se_ratio_min < 0.92), c("design", "seed", "se_ratio_min", "se_ratio_max", "n_warn")])
cat("warned", sum(all$n_warn > 0), "\n")
print(all[all$n_warn > 0 & all$ess_min >= 0.25, c("design", "seed", "n_warn", "ess_min")])
aw <- do.call(rbind, lapply(c(list.files(file.path(root, "accuracy"), "away", full.names = TRUE),
                              list.files(file.path(root, "accuracy2"), "away", full.names = TRUE)),
                            read.delim))
aw <- aw[!aw$design %in% scale_d, ]
cat("away points", nrow(aw), " finite ad_vs_fd", sum(is.finite(aw$ad_vs_fd)),
    " max", max(aw$ad_vs_fd, na.rm = TRUE), " median", median(aw$ad_vs_fd, na.rm = TRUE), "\n")
ub <- do.call(rbind, lapply(list.files(file.path(root, "unbiased"), "tsv$", full.names = TRUE), read.delim))
print(format(ub[ub$pt %in% c("anchor", "sd_up", "sd_dn", "b1_3se"), c("design", "pt", "err_logdet", "err_joint", "mean_ratio_min", "mean_ratio_max", "max_abs_z")], digits = 3), row.names = FALSE)
cat("max err_logdet", max(ub$err_logdet), "max err_joint", max(ub$err_joint), "\n")
mc <- do.call(rbind, lapply(list.files(file.path(root, "mcse"), "tsv$", full.names = TRUE), read.delim))
options(width = 180)
print(format(mc[order(mc$design, mc$draws), c("design", "draws", "nok", "sd_over_mcse", "cover95", "warned", "mean_err", "max_abs_z", "ess_min")], digits = 3), row.names = FALSE)

set.seed(21)
G <- 20
g <- factor(rep(seq_len(G), each = 6))
x <- rnorm(120)
d <- data.frame(y = rbinom(120, 1, plogis(0.2 + 0.6 * x + rnorm(G, 0, 0.9)[g])), x, g)
fi <- frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d, importance = 200L)
cat("\npredict newdata:\n")
print(head(predict(fi, newdata = d[1:3, ]), 3))
print(head(frm_linpred(fi, newdata = d[1:3, ]), 3))
tf <- tempfile(fileext = ".rds")
saveRDS(fi, tf)
b <- readRDS(tf)
cat("saved fit predict and linpred agree:",
    isTRUE(all.equal(frm_linpred(b, newdata = d[1:3, ]), frm_linpred(fi, newdata = d[1:3, ]))), "\n")
W <- matrix(0, 8, 8)
W[cbind(1:7, 2:8)] <- 1
W <- W + t(W)
dimnames(W) <- list(as.character(1:8), as.character(1:8))
ds <- data.frame(y = rpois(80, 2), x = rnorm(80), loc = factor(rep(as.character(1:8), 10)))
cat("car:", tryCatch({
  frm(bf(y ~ x + car(W, gr = loc)), family = poisson(), data = ds, data2 = list(W = W),
      importance = 64L)
  "fit"
}, error = function(e) substr(conditionMessage(e), 1, 90)), "\n")
