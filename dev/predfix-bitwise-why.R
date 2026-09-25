# Which elements of each moved fit differ, and what the base arm held
# there. Rscript dev/predfix-bitwise-why.R base1 lane
a <- commandArgs(trailingOnly = TRUE)
dir <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
A <- readRDS(paste0(dir, "bitwise-", a[1], ".rds"))
B <- readRDS(paste0(dir, "bitwise-", a[2], ".rds"))
els <- c("error", "ll", "par", "conv", "se_full", "se_fix", "sd_re")
desc <- function(v) {
  if (inherits(v, "err")) paste("ERROR:", substr(unclass(v), 1, 50)) else
    if (is.null(v)) "NULL" else paste(length(v), "numbers")
}
rows <- list()
for (k in names(A)) {
  d <- els[!vapply(els, function(e) identical(A[[k]][[e]], B[[k]][[e]]), NA)]
  if (!length(d)) next
  rows[[k]] <- data.frame(fit = k, differ = paste(d, collapse = ","),
    base_se = desc(A[[k]]$se_fix), lane_se = desc(B[[k]]$se_fix))
}
r <- do.call(rbind, rows)
rownames(r) <- NULL
options(width = 200)
print(r)
mv <- r$fit[!grepl("^TARGET", r$fit)]
cat("\nbattery fits moved:", length(mv), "; of them ll, par and conv all",
    "identical:", sum(vapply(mv, function(k) identical(A[[k]][c("ll",
    "par", "conv")], B[[k]][c("ll", "par", "conv")]), NA)),
    "; base vcov() an error on:", sum(vapply(mv, function(k)
    inherits(A[[k]]$se_fix, "err"), NA)), "\n")
