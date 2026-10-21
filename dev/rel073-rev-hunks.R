# Release reviewer, 0.73.0: for every file two or more lanes changed,
# the expected count of each line (the 4854e425 count plus every lane's
# net change) against the release tree's count. A mismatch is a line
# lost, duplicated or resurrected by the merge, or a resolve edit.
#   Rscript dev/rel073-rev-hunks.R   (log dev/rel073-rev-log/hunks.txt)
rel <- "C:/Users/adf44/source/r/frmtmb-wt-release"
lanes <- c("loodraws", "cdfs", "nb2dm", "fuzzfix")
st <- lapply(lanes, function(l) {
  x <- readLines(file.path(rel, "dev/rel073-log/merge-work", paste0(l, ".status")))
  sub("^[A-Z]\t", "", x)
})
names(st) <- lanes
all <- table(unlist(st))
multi <- names(all)[all >= 2]
gen <- "^(man/|NAMESPACE$|extensions/[^/]+/NAMESPACE$|R/links-brms.R$|.*helper-brms-suite.R$)"
multi <- multi[!grepl(gen, multi)]
cat("files changed by two or more lanes, generated left out:", length(multi),
    "\n")
rd <- function(f) {
  if (!file.exists(f)) return(character())
  sub("\r$", "", readLines(f, warn = FALSE))
}
blob <- function(f) {
  sub("\r$", "", system2("git", c("-C", rel, "show", paste0("4854e425:", f)),
                         stdout = TRUE))
}
cnt <- function(x) table(factor(x, levels = unique(c(x))))
for (f in multi) {
  who <- lanes[vapply(st, function(s) f %in% s, TRUE)]
  b <- blob(f)
  keys <- unique(c(b, unlist(lapply(who, function(l)
    rd(file.path(dirname(rel), paste0("frmtmb-wt-", l), f)))),
    rd(file.path(rel, f))))
  tab <- function(x) as.vector(table(factor(x, levels = keys)))
  exp <- tab(b)
  for (l in who) {
    exp <- exp + tab(rd(file.path(dirname(rel), paste0("frmtmb-wt-", l), f))) -
      tab(b)
  }
  got <- tab(rd(file.path(rel, f)))
  bad <- which(exp != got)
  cat(sprintf("%-55s lanes %-26s lines %5d mismatched keys %d\n", f,
              paste(who, collapse = ","), length(rd(file.path(rel, f))),
              length(bad)))
  for (i in head(bad, 12)) {
    cat(sprintf("    expected %d got %d: %s\n", exp[i], got[i],
                substr(keys[i], 1, 110)))
  }
  if (length(bad) > 12) cat("    ...", length(bad) - 12, "more\n")
}
