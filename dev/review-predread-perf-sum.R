# Item 3: ratios of the timing arms, per round and over rounds, against the control.
O <- "dev/review-predread-out"
rd <- function(arm, r) {
  f <- file.path(O, sprintf("perf-%s-%d.tsv", arm, r))
  if (file.exists(f)) cbind(read.delim(f), round = r)
}
b <- do.call(rbind, lapply(1:3, rd, arm = "before"))
a <- do.call(rbind, lapply(1:3, rd, arm = "after"))
tmb <- tapply(b$tmin, b$reader, min)
tma <- tapply(a$tmin, a$reader, min)
rounds <- sapply(unique(a$round), function(r) {
  x <- a[a$round == r, ]
  y <- b[b$round == r, ]
  setNames(x$tmin / y$tmin[match(x$reader, y$reader)], x$reader)
})
out <- data.frame(reader = names(tma), before_s = signif(as.numeric(tmb[names(tma)]), 3),
                  after_s = signif(as.numeric(tma), 3),
                  ratio_of_mins = round(as.numeric(tma) / as.numeric(tmb[names(tma)]), 3))
print(cbind(out, round(rounds[out$reader, , drop = FALSE], 3)), row.names = FALSE)
cat("rounds:", length(unique(a$round)), "\n")
