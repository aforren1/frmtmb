# Guardrail 3 for the mixture_mvn() move: the lines of R/ that move,
# read from the base tree (a git show of 19b162b9, not the edited
# worktree).
#
#   Rscript dev/mvnlatent-measure.R <base families.R> <base compat.R>
a <- commandArgs(trailingOnly = TRUE)
src <- readLines(a[1])
objs <- c("mvn_cov_models", "mvn_cov_spec", "mixture_mvn", "mvn_sim_rows")
ex <- parse(a[1], keep.source = TRUE)
sr <- lapply(ex, function(e) attr(e, "srcref"))
srs <- attr(ex, "srcref")
nm <- vapply(ex, function(e) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) {
    as.character(e[[2]])
  } else ""
}, "")
classify <- function(l) {
  ifelse(grepl("^\\s*$", l), "blank",
         ifelse(grepl("^\\s*#", l), "comment", "code"))
}
tot <- c(code = 0, comment = 0, blank = 0)
for (o in objs) {
  i <- which(nm == o)
  s <- srs[[i]][1]
  e <- srs[[i]][3]
  # the comment block right above the object (roxygen or plain)
  b <- s
  while (b > 1 && grepl("^\\s*#", src[b - 1])) b <- b - 1
  cl <- table(factor(classify(src[b:e]), c("code", "comment", "blank")))
  cat(sprintf("%-16s lines %d-%d: code %d, comment %d, blank %d\n", o, b,
              e, cl[["code"]], cl[["comment"]], cl[["blank"]]))
  tot <- tot + as.numeric(cl)
}
cat(sprintf("TOTAL families.R: code %d, comment %d, blank %d, all %d\n",
            tot[1], tot[2], tot[3], sum(tot)))
cp <- readLines(a[2])
hit <- grep("mixture_mvn", cp)
cat("compat.R lines naming mixture_mvn:", length(hit), "\n")
rows <- grep('^\\s*r[(].*"mixture_mvn"', cp)
cat("compat.R rules naming mixture_mvn:", length(rows), "\n")
fam_src <- length(src)
cat("families.R total lines:", fam_src, "\n")
allR <- Sys.glob(file.path(dirname(a[1]), "*.R"))
nR <- sum(vapply(allR, function(f) length(readLines(f)), 0))
cat("R/ total lines:", nR, "in", length(allR), "files\n")
cat(sprintf("share of R/ lines: %.2f%%\n", 100 * sum(tot) / nR))
