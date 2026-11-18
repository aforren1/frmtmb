# rel084-kt copy of dev/nlband-fix1.R: sources dev/rel084-kt/nlband-common.R (rellib-r23). No other change.
# Lane nlband, fix round 1 (review dev/reviews/2026-11-17-nlband.md):
#   SEP   separated bernoulli band data (y = 1 for x > 2), seeds 1 to 10
#         (item 1); a fit must warn or be refused, not kept silent
#   STUD  student band data (y = 1 / (x - 2) + 0.2 t3), seeds 1 to 3, the
#         default against the fit from (1, -2) (item 2)
#   NEAR  a pole 0.003 below the first row above the hole, seeds 1 to 3
#         (item 4, filed)
#   BAND  the band data, seeds 1 to 3 (5115 still fixed)
# NLBAND_ARM=before|after Rscript dev/nlband-fix1.R [sections]
args <- commandArgs(trailingOnly = TRUE)
secs <- if (length(args)) strsplit(args[1], ",")[[1]] else
  c("sep", "stud", "near", "band")
source("dev/rel084-kt/nlband-common.R")
show <- function(tag, r) {
  if (iserr(r)) return(cat(tag, "| ERROR", substr(r$msg, 1, 150), "\n"))
  se <- tryCatch(sqrt(diag(suppressWarnings(vcov(r)))),
                 error = function(e) NA)
  cat(sprintf("%s | ll %.8g code %d beta %s se %s | warn %s msg %s\n",
              tag, llv(r), r$opt$convergence,
              paste(signif(r$opt$par, 6), collapse = " "),
              paste(signif(se, 4), collapse = " "), wsum(r), msum(r)))
  for (m in attr(r, "m")) if (grepl("pole", m)) cat("   MSG", m, "\n")
}
if ("sep" %in% secs) for (s in 1:10) {
  d <- band_data(s)
  d$y <- as.integer(d$x > 2)
  show(paste("SEP seed", s), tt(frm(band_form, data = d,
                                    family = bernoulli())))
}
if ("stud" %in% secs) for (s in 1:3) {
  set.seed(s)
  x <- band_data(s)$x
  d <- data.frame(x = x, y = 1 / (x - 2) + 0.2 * rt(length(x), 3))
  show(paste("STUD seed", s, "default"),
       tt(frm(band_form, data = d, family = student())))
  show(paste("STUD seed", s, "from (1,-2)"),
       tt(frm(band_form, data = d, family = student(),
              start = list(beta = c(1, -2)),
              control = frmtmb_control(nl_starts = FALSE))))
}
if ("near" %in% secs) for (s in 1:3) {
  d <- band_data(s)
  p <- min(d$x[d$x > 2]) - 0.003
  set.seed(s + 700)
  d$y <- 1 / (d$x - p) + rnorm(nrow(d), 0, 0.2)
  show(paste("NEAR seed", s, "default"), tt(frm(band_form, data = d)))
  show(paste("NEAR seed", s, "from pole"),
       tt(frm(band_form, data = d, start = list(beta = c(1, -p)),
              control = frmtmb_control(nl_starts = FALSE))))
}
if ("band" %in% secs) for (s in 1:3) {
  show(paste("BAND seed", s), tt(frm(band_form, data = band_data(s))))
}
cat("DONE\n")
