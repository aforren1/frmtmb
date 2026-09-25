# The seeds driver 4 refits: targets minus what rds-check.txt says fitted.
x <- readLines("dev/phase3b-log/rds-check.txt")
tg <- list(cleft = 1:60, contfix = c(1:20, 61:80), contdl = 1:40,
           collapse = 1:30, cont = 1:30, cens = 1:40)
out <- character(0)
for (arm in names(tg)) {
  l <- grep(paste0("^", arm, " "), x, value = TRUE)
  miss <- suppressWarnings(as.integer(strsplit(sub(".*missing: *", "", l),
                                               " ")[[1]]))
  miss <- intersect(tg[[arm]], miss[!is.na(miss)])
  out <- c(out, paste(arm, miss))
}
writeLines(out, "dev/phase3b-log/recov5-todo.txt")
cat(length(out), "jobs\n")
