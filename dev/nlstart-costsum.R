# Summarize dev/nlstart-cost.R output. Usage: Rscript dev/nlstart-costsum.R <dir>
dir <- commandArgs(trailingOnly = TRUE)[1]
l <- unlist(lapply(list.files(dir, "[.]txt$", full.names = TRUE),
                   function(f) grep("^COST ", readLines(f), value = TRUE)))
x <- read.table(text = sub("^COST ", "", l), stringsAsFactors = FALSE,
                col.names = c("family", "seed", "slope", "body", "ll_off",
                              "ll_on", "ll_off2", "t_off", "t_on", "t_off2"))
cat("rows", nrow(x), "\n")
cat("seconds: off", sum(x$t_off), " on", sum(x$t_on), " off2", sum(x$t_off2),
    "\n")
cat("ratio on / off:", round(sum(x$t_on) / sum(x$t_off), 3),
    "  control off2 / off:", round(sum(x$t_off2) / sum(x$t_off), 3), "\n")
cat("by body, on / mean(off, off2):\n")
print(round(tapply(x$t_on, x$body, sum) /
              tapply((x$t_off + x$t_off2) / 2, x$body, sum), 3))
cat("by family, on / mean(off, off2):\n")
print(round(tapply(x$t_on, x$family, sum) /
              tapply((x$t_off + x$t_off2) / 2, x$family, sum), 3))
cat("logLik on - off: higher > 1e-4", sum(x$ll_on - x$ll_off > 1e-4,
                                          na.rm = TRUE),
    " lower < -1e-4", sum(x$ll_on - x$ll_off < -1e-4, na.rm = TRUE),
    "  off2 vs off differ > 1e-8:", sum(abs(x$ll_off2 - x$ll_off) > 1e-8,
                                        na.rm = TRUE), "\n")
