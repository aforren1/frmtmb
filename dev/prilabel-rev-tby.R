# Review of lane prilabel: one Student-t gr(by = ) term, which frmtmb
# splits into one block per by-level, on both arms.
#   arg 1: "lane" or "base". Data seed 21 (as dev/prilabel-rev-misc.R).
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, paste0("dev/prilabel-rev-log/tby-", arm, ".txt")),
     split = TRUE)
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(21)
n <- 360
g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
d <- data.frame(g, x = rnorm(n), z = rnorm(n),
                hg = factor(ifelse(as.integer(g) <= 12L, "u", "v")))
d$y <- 1 + 0.5 * d$x + rnorm(24, 0, 0.8)[d$g] + rnorm(n, 0, 0.6)
for (f in list(y ~ x + (1 | gr(g, by = hg, dist = "student")),
               y ~ x + (1 + x | gr(g, by = hg, dist = "student")))) {
  r <- tryCatch(frm(bf(f), data = d), error = function(e) {
    paste("ERROR:", conditionMessage(e))
  })
  cat("\n--", deparse1(f), "\n")
  cat(if (is.character(r)) r else
    paste("fits, logLik", format(as.numeric(logLik(r)), digits = 10),
          "blocks", length(r$frame$re_blocks)), "\n")
}
sink()
