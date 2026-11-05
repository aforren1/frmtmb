# Short check 2, item 3: the candidate gate, per fit and on the lane's design labels.
args <- commandArgs(TRUE)
dir <- paste0("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/", args[[1]])
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE, comment.char = "")))
x <- x[is.finite(x$ll) & is.finite(x$x_ls), ]
ord <- x$design %in% c("ord4", "cum3", "cum5", "cum7", "f_cumpro")
x$wrong <- abs(x$ls - x$x_ls) > log(1.25) | abs(x$ar - x$x_ar) > .3
x$sh <- pmax(abs(x$b_ls), abs(x$b_ar))
x$sh[!is.finite(x$sh)] <- Inf
x$rs <- exp(x$ls - x$x_ls); x$dar <- x$ar - x$x_ar
x$cov <- abs(x$ls - log(x$s_true)) <= 1.96 * x$se_ls
des <- do.call(rbind, lapply(split(x, paste(x$set, x$design)), function(d) data.frame(set = d$set[1], design = d$design[1], bad = median(d$rs, na.rm = TRUE) < .9 | median(d$rs, na.rm = TRUE) > 1.15 | abs(median(d$dar, na.rm = TRUE)) > .15 | mean(d$cov, na.rm = TRUE) < .8)))
x <- merge(x, des, by = c("set", "design"))
ord <- x$design %in% c("ord4", "cum3", "cum5", "cum7", "f_cumpro")
tab <- function(fire, lab, truth, what) {
  o <- sapply(split(seq_len(nrow(x)), x$set), function(i) sprintf("%s %.3f/%.3f", x$set[i][1], mean(fire[i][truth[i]]), mean(fire[i][!truth[i]])))
  cat(sprintf("%-6s %-30s %s  ALL %.3f/%.3f  J=%.3f\n", what, lab, paste(o, collapse = "  "), mean(fire[truth]), mean(fire[!truth]), mean(fire[truth]) - mean(fire[!truth])))
}
for (tt in list(list(x$wrong, "fit"), list(x$bad, "design"))) {
  tab(x$warned == 1, "lane's warning", tt[[1]], tt[[2]])
  for (k in c(.1, .15, .2)) tab(x$sh > k | ord, sprintf("|shift| > %.2f or ordinal", k), tt[[1]], tt[[2]])
  tab(x$sh > .15, "|shift| > 0.15 alone", tt[[1]], tt[[2]])
}
g <- x$sh > .15 | ord
d <- do.call(rbind, lapply(split(x, paste(x$set, x$design)), function(q) data.frame(set = q$set[1], design = q$design[1], wrong = round(mean(q$wrong), 2), lane = round(mean(q$warned == 1), 2), gate = round(mean(q$sh > .15 | q$design %in% c("ord4", "cum3", "cum5", "cum7", "f_cumpro")), 2))))
cat("\nthe fresh set, by design\n")
print(d[d$set == "fresh", ][order(-d$wrong[d$set == "fresh"]), ], row.names = FALSE)
cat("\ndesigns the gate gets most wrong:\n")
print(d[d$wrong >= .5 & d$gate < .5, ], row.names = FALSE)
print(d[d$wrong <= .2 & d$gate > .25, ], row.names = FALSE)
cat("\nmean correction seconds", round(mean(x$t_tk, na.rm = TRUE), 3), " at n", round(mean(x$n)), "\n")
