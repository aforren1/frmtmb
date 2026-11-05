# Short check 2, item 3: per fit, not per design. Does a fit's own number say that
# THIS fit misreads? Truth: log sderr off the exact MLE by more than log(1.25), or
# ar off by more than 0.3. Scores: the lane's warning (a yes or no), and the
# first-order shift. Also fits whose exact MLE is at a bound are kept.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/tk"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE, comment.char = "")))
x <- x[is.finite(x$ll) & is.finite(x$x_ls), ]
x$wrong <- abs(x$ls - x$x_ls) > log(1.25) | abs(x$ar - x$x_ar) > .3
x$sh <- pmax(abs(x$b_ls), abs(x$b_ar))
x$sh[!is.finite(x$sh)] <- Inf
cat("fits", nrow(x), " wrong", sum(x$wrong), " shift NA", sum(!is.finite(x$sh)), "\n")
tab <- function(fire, lab) {
  o <- sapply(split(seq_len(nrow(x)), x$set), function(i) sprintf("%s %.3f/%.3f", x$set[i][1], mean(fire[i][x$wrong[i]]), mean(fire[i][!x$wrong[i]])))
  cat(sprintf("%-26s %s  ALL %.3f/%.3f  J=%.3f\n", lab, paste(o, collapse = "  "), mean(fire[x$wrong]), mean(fire[!x$wrong]), mean(fire[x$wrong]) - mean(fire[!x$wrong])))
}
tab(x$warned == 1, "lane's warning")
tab((x$warned == 1) | (x$margin < 2 & !is.na(x$margin)), "lane's + margin")
for (k in c(.03, .05, .075, .1, .15, .2, .3)) tab(x$sh > k, sprintf("|shift| > %.3f", k))
# a by-family view of the gate at one threshold
k <- .1
x$fire <- x$sh > k
d <- do.call(rbind, lapply(split(x, paste(x$set, x$design)), function(q) data.frame(set = q$set[1], design = q$design[1], fits = nrow(q), wrong = round(mean(q$wrong), 2), lane = round(mean(q$warned == 1), 2), gate = round(mean(q$fire), 2), med_shift = round(median(q$sh[is.finite(q$sh)]), 3))))
print(d[order(d$set, -d$wrong), ], row.names = FALSE)
cat("\ndisagreements at the design level (lane fires on a design that is rarely wrong, or the reverse):\n")
cat("wrong >= 0.5: lane < 0.5:", sum(d$wrong >= .5 & d$lane < .5), " gate < 0.5:", sum(d$wrong >= .5 & d$gate < .5), " of", sum(d$wrong >= .5), "\n")
cat("wrong <= 0.2: lane > 0.25:", sum(d$wrong <= .2 & d$lane > .25), " gate > 0.25:", sum(d$wrong <= .2 & d$gate > .25), " of", sum(d$wrong <= .2), "\n")
print(d[d$wrong <= .2 & d$lane > .25, c("set", "design", "wrong", "lane", "gate")], row.names = FALSE)
print(d[d$wrong >= .5 & d$lane < .5, c("set", "design", "wrong", "lane", "gate")], row.names = FALSE)
print(d[d$wrong >= .5 & d$gate < .5, c("set", "design", "wrong", "lane", "gate")], row.names = FALSE)
print(d[d$wrong <= .2 & d$gate > .25, c("set", "design", "wrong", "lane", "gate")], row.names = FALSE)
cat("\nfresh set alone, lane against gate, by design: see the table\n")
