# Item 2 summary (dev/fitedge-nl2-sweep.R): each arm's default against the best
# of the reference and both arms. Short: more than 1e-3 nats below it.
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
rd <- function(f) read.delim(file.path(o, f), stringsAsFactors = FALSE)
k <- c("grp", "body", "seed", "slope")
for (g in commandArgs(TRUE)) {
  b <- rd(sprintf("nl2-%s-fit-base.tsv", g)); l <- rd(sprintf("nl2-%s-fit-lane.tsv", g))
  r <- rd(sprintf("nl2-%s-ref-lane.tsv", g))
  m <- merge(merge(b, l, by = k, suffixes = c(".b", ".l")),
             setNames(r[, c(k, "ll")], c(k, "ref")), by = k)
  best <- pmax(m$ref, m$ll.b, m$ll.l, na.rm = TRUE)
  for (arm in c("b", "l")) {
    ll <- m[[paste0("ll.", arm)]]; cv <- m[[paste0("conv.", arm)]]
    nw <- m[[paste0("nwarn.", arm)]]
    sh <- is.na(ll) | ll < best - 1e-3
    cat(sprintf("NL2 %s %s rows %d short %d silent %d worst %.3f time %.1f s\n", g,
                c(b = "base", l = "lane")[arm], nrow(m), sum(sh),
                sum(sh & cv %in% 0 & nw %in% 0), min(ll - best, na.rm = TRUE),
                sum(m[[paste0("sec.", arm)]])))
  }
  cat(sprintf("NL2 %s lane above base %d, below base %d\n", g,
              sum(m$ll.l > m$ll.b + 1e-3, na.rm = TRUE),
              sum(m$ll.l < m$ll.b - 1e-3, na.rm = TRUE)))
  d <- which(abs(m$ll.l - m$ll.b) > 1e-3 | is.na(m$ll.l) != is.na(m$ll.b))
  for (i in d) {
    cat(sprintf("NL2ROW %s %s seed %d slope %.1f base %.4f (conv %s w %s) lane %.4f (conv %s w %s) ref %.4f\n",
                g, m$body[i], m$seed[i], m$slope[i], m$ll.b[i], m$conv.b[i], m$nwarn.b[i],
                m$ll.l[i], m$conv.l[i], m$nwarn.l[i], m$ref[i]))
  }
}
