# Reviewer copy of dev/setier-verdict-diff.R reading dev/setier-rev-log.
# Lane setier: are the standard-error verdicts of two suite runs the
# same? Per test file, the multiset of standard-error conditions raised
# (boundary messages, standard-error warnings, separation warnings), each
# cut to its parameter list, from the .cond logs of dev/setier-run1.R.
#   Rscript dev/setier-verdict-diff.R <run a> <run b>
a <- commandArgs(TRUE)
root <- "dev/setier-rev3-log"
key <- function(f) {
  if (!file.exists(f)) return(character(0))
  x <- readLines(f, warn = FALSE)
  x <- x[grepl(paste0("^M\tBoundary [(]singular[)] fit|",
                      "^W\tStandard errors are not available|",
                      "^W\tThe data separate the outcomes"), x)]
  # the part that names the parameters; numbers (counts of moved
  # observations, the optimizer's message) are left out
  x <- sub("^(M\tBoundary [(]singular[)] fit: [^)]*[)]).*", "\\1", x)
  x <- sub("(parameters[.] [^:]*):.*", "\\1", x)
  x <- sub("^(W\tThe data separate the outcomes [^:]*:[^;]*every).*",
           "\\1", x)
  sort(x)
}
fa <- list.files(file.path(root, a[1]), pattern = "[.]cond$")
fb <- list.files(file.path(root, a[2]), pattern = "[.]cond$")
files <- union(fa, fb)
nd <- 0L
nl <- 0L
for (f in sort(files)) {
  ka <- key(file.path(root, a[1], f))
  kb <- key(file.path(root, a[2], f))
  if (identical(ka, kb)) next
  nd <- nd + 1L
  only_a <- setdiff(ka, kb)
  only_b <- setdiff(kb, ka)
  # a multiset difference: how many conditions one run raised that the
  # other did not, counted with multiplicity
  u <- union(ka, kb)
  nl <- nl + sum(abs(table(factor(ka, u)) - table(factor(kb, u))))
  cat("==", sub("[.]cond$", "", f), ": ", length(ka), " vs ", length(kb),
      "\n", sep = "")
  for (x in only_a) cat("  only ", a[1], ": ", substr(x, 1, 150), "\n",
                        sep = "")
  for (x in only_b) cat("  only ", a[2], ": ", substr(x, 1, 150), "\n",
                        sep = "")
}
tot_a <- sum(vapply(fa, function(f) length(key(file.path(root, a[1], f))),
                    0L))
tot_b <- sum(vapply(fb, function(f) length(key(file.path(root, a[2], f))),
                    0L))
cat(sprintf("VERDICT %s vs %s: conditions %d vs %d; files differing %d; ",
            a[1], a[2], tot_a, tot_b, nd),
    sprintf("conditions not matched %d\n", nl))
