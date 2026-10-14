# Lane postfix: the workaround rows (dev/brms-vignettes/_workarounds.R)
# whose outcome differs between rellib-r8 and the lane build.
#   Rscript dev/postfix-wa-diff.R
rd <- function(arm) {
  utils::read.csv(file.path("dev/postfix-bv-out", arm, "workarounds.csv"),
                  stringsAsFactors = FALSE)
}
b <- rd("wa-before")
a <- rd("wa-after")
stopifnot(identical(b$label, a$label))
cat("workaround rows:", nrow(a), "; the brms spelling runs on",
    sum(b$ok), "before and", sum(a$ok), "after\n")
d <- which(b$ok != a$ok | b$msg != a$msg)
for (i in d) {
  cat(sprintf("%s\n    before ok=%s %s\n    after  ok=%s %s\n", a$label[i],
              b$ok[i], substr(b$msg[i], 1, 140), a$ok[i],
              substr(a$msg[i], 1, 140)))
}
