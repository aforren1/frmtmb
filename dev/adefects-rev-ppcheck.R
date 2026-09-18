source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260917); n <- 40
d <- data.frame(x = rnorm(n), g = factor(rep_len(paste0("g", 1:4), n)))
d$y <- rnorm(n, 1 + 0.5 * d$x, 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
types <- c("violin_grouped", "stat_grouped", "dens_overlay_grouped",
           "intervals_grouped", "ribbon_grouped", "freqpoly_grouped",
           "ecdf_overlay_grouped", "stat_freqpoly_grouped")
for (ty in types) {
  r <- tryCatch({
    p <- pp_check(fit, type = ty, group = "g", ndraws = 5)
    invisible(print(p))
    "OK (plotted)"
  }, error = function(e) paste0("ERR: ",
       substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1, 80)),
     warning = function(w) paste0("WARN: ",
       substr(conditionMessage(w), 1, 60)))
  cat(sprintf("  %-24s %s\n", ty, r))
}
cat("\n  ungrouped control:\n")
for (ty in c("dens_overlay", "stat", "violin")) {
  r <- tryCatch({ invisible(print(pp_check(fit, type = ty, ndraws = 5)));
                  "OK" },
    error = function(e) paste0("ERR: ",
      substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1, 70)))
  cat(sprintf("  %-24s %s\n", ty, r))
}
