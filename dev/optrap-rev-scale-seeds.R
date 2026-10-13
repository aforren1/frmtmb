# Reviewer of lane optrap: the seeds behind dev/optrap-rev-scale-sum.R's
# shortfalls, per case.
#   Rscript dev/optrap-rev-scale-seeds.R <case> [<case> ...]
for (cs in commandArgs(TRUE)) {
  a <- read.delim(sprintf("dev/optrap-rev-log/scale-%s-lane.tsv", cs),
                  quote = "")
  b <- read.delim(sprintf("dev/optrap-rev-log/scale-%s-r7.tsv", cs),
                  quote = "")
  a$llu_r7 <- b$llu
  a$code_r7 <- b$code
  a$warn_r7 <- b$warn
  ref <- tapply(pmax(a$llu, a$ext, b$llu, na.rm = TRUE), a$seed, max,
                na.rm = TRUE)
  k <- as.character(a$seed)
  a$gap <- a$llu - ref[k]
  a$gapr7 <- a$llu_r7 - ref[k]
  a$extgap <- a$ext - ref[k]
  x <- a[a$gap < -1e-4 | a$gapr7 < -1e-4, ]
  print(x[, c("case", "seed", "scale", "code", "code_r7", "gap", "gapr7",
              "extgap", "nwarn")])
  cat("lane warnings:\n")
  print(table(substr(a$warn[nzchar(a$warn)], 1, 90)))
}
