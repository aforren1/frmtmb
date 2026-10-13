# Reviewer of the 0.70.0 consolidation: is every lane's development NEWS
# bullet in the release NEWS once, and is every release bullet a lane's
# or the consolidation's? Bullets are compared as whitespace-collapsed
# text, with seflags' "BREAKING: " prefix removed.
#   Rscript dev/relrev070-news.R
R <- "C:/Users/adf44/source/r"
bullets <- function(lines) {
  lines <- sub("\r$", "", lines)
  out <- character(0)
  cur <- NULL
  flush <- function() {
    if (!is.null(cur)) out <<- c(out, gsub("[[:space:]]+", " ",
                                           paste(cur, collapse = " ")))
    cur <<- NULL
  }
  for (l in lines) {
    if (grepl("^[*-] ", l)) {
      flush()
      cur <- sub("^[*-] ", "", l)
    } else if (grepl("^#", l) || !nzchar(trimws(l))) {
      flush()
    } else if (!is.null(cur)) {
      cur <- c(cur, trimws(l))
    }
  }
  flush()
  trimws(out)
}
dev_section <- function(path) {
  x <- readLines(path, warn = FALSE)
  s <- grep("^# .*\\(development version\\)", x)
  if (!length(s)) return(character(0))
  e <- grep("^# ", x)
  e <- e[e > s[1]][1]
  x[(s[1] + 1):(e - 1)]
}
norm <- function(b) sub("^BREAKING: ", "", b)
rel <- readLines(file.path(R, "frmtmb-wt-release/NEWS.md"), warn = FALSE)
e <- grep("^# frmtmb 0[.]69[.]1", rel)
relb <- norm(bullets(rel[1:(e - 1)]))
cat("release 0.70.0 bullets:", length(relb), "\n")
seen <- character(0)
for (l in c("moparity", "optrap", "seflags")) {
  lb <- norm(bullets(dev_section(file.path(R, paste0("frmtmb-wt-", l),
                                           "NEWS.md"))))
  hits <- vapply(lb, function(b) sum(relb == b), 0L)
  cat(l, ": bullets", length(lb), "| found once", sum(hits == 1),
      "| missing", sum(hits == 0), "| more than once", sum(hits > 1), "\n")
  for (b in lb[hits != 1]) cat("  NOT ONCE (", sum(relb == b), "): ",
                               substr(b, 1, 160), "\n", sep = "")
  seen <- c(seen, lb)
}
extra <- relb[!relb %in% seen]
cat("release bullets from no lane:", length(extra), "\n")
for (b in extra) cat("  EXTRA: ", substr(b, 1, 200), "\n", sep = "")
# The BREAKING-marked seflags bullets: under which heading are they?
sec <- NA_character_
for (i in seq_len(e - 1)) {
  if (grepl("^## ", rel[i])) sec <- rel[i]
  if (grepl("^[*-] ", rel[i]) && grepl("alias|curved ridge", rel[i],
                                       ignore.case = TRUE)) {
    cat("line", i, "under", sec, ":", substr(rel[i], 1, 100), "\n")
  }
}
cat("BREAKING: prefixes left in 0.70.0:", sum(grepl("BREAKING", rel[1:(e - 1)])),
    "\n")
# frmtmb.sample
for (pkg in c("frmtmb.sample")) {
  rs <- readLines(file.path(R, "frmtmb-wt-release/extensions", pkg,
                            "NEWS.md"), warn = FALSE)
  e2 <- grep("^# ", rs)[2]
  rb <- bullets(rs[1:(e2 - 1)])
  lb <- bullets(dev_section(file.path(R, "frmtmb-wt-moparity/extensions",
                                      pkg, "NEWS.md")))
  cat(pkg, ": head", rs[1], "| release bullets", length(rb), "| lane",
      length(lb), "| lane found", sum(lb %in% rb), "\n")
  hd <- system2("git", c("-C", file.path(R, "frmtmb-wt-release"), "show",
                         paste0("2884ccbd:extensions/", pkg, "/NEWS.md")),
                stdout = TRUE)
  cat(pkg, ": released sections preserved:",
      identical(sub("\r$", "", hd), sub("\r$", "", rs[e2:length(rs)])), "\n")
}
re <- readLines(file.path(R, "frmtmb-wt-release/extensions/frmtmb.eam/NEWS.md"),
                warn = FALSE)
he <- system2("git", c("-C", file.path(R, "frmtmb-wt-release"), "show",
                       "2884ccbd:extensions/frmtmb.eam/NEWS.md"), stdout = TRUE)
e3 <- grep("^# ", re)[2]
cat("frmtmb.eam: head", re[1], "| released preserved:",
    identical(sub("\r$", "", he), sub("\r$", "", re[e3:length(re)])), "\n")
cat(re[1:(e3 - 1)], sep = "\n")
