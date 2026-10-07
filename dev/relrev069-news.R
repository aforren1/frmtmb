# Reviewer of the 0.69.0 consolidation: is every lane NEWS bullet in
# the release NEWS exactly once (whitespace normalized), and what does
# the release hold that no lane wrote?
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
norm <- function(x) gsub("[[:space:]]+", " ", trimws(x))
dev_section <- function(path, head) {
  x <- readLines(path, warn = FALSE)
  s <- grep(paste0("^# ", head, " "), x)
  i <- s[1]
  j <- s[2]
  x[(i + 1):(j - 1)]
}
bullets <- function(lines) {
  starts <- grep("^[*] ", lines)
  if (!length(starts)) return(character(0))
  ends <- c(starts[-1] - 1, length(lines))
  b <- mapply(function(a, z) {
    l <- lines[a:z]
    l <- l[!grepl("^#", l)]
    # a bullet ends at its first blank line
    k <- which(!nzchar(trimws(l)))
    if (length(k)) l <- l[seq_len(k[1] - 1)]
    norm(paste(l, collapse = " "))
  }, starts, ends)
  unname(b)
}
for (pk in c("frmtmb", "frmtmb.sample")) {
  sub <- if (pk == "frmtmb") "NEWS.md" else
    "extensions/frmtmb.sample/NEWS.md"
  relb <- bullets(dev_section(file.path(rel, sub), pk))
  cat("==", pk, ": release bullets", length(relb), "\n")
  seen <- integer(length(relb))
  for (l in c("ciharden", "surface", "setier", "optima")) {
    f <- file.path(root, paste0("frmtmb-wt-", l), sub)
    x <- readLines(f, warn = FALSE)
    if (!any(grepl(paste0("^# ", pk, " [(]development version[)]"), x))) {
      cat(l, ": no development section\n"); next
    }
    lb <- bullets(dev_section(f, pk))
    for (b in lb) {
      hit <- which(relb == b)
      seen[hit] <- seen[hit] + 1
      if (length(hit) != 1) {
        # near match: first 60 characters
        near <- which(substr(relb, 1, 60) == substr(b, 1, 60))
        cat(sprintf("  %s bullet found %d times exactly, %d by prefix: %s\n",
                    l, length(hit), length(near), substr(b, 1, 100)))
        if (length(near) == 1) {
          a <- strsplit(b, " ")[[1]]; r <- strsplit(relb[near], " ")[[1]]
          cat("    lane-only words:", setdiff(a, r), "\n")
          cat("    release-only words:", setdiff(r, a), "\n")
          seen[near] <- seen[near] + 1
        }
      }
    }
    cat(l, ":", length(lb), "bullets\n")
  }
  cat("release bullets matched by no lane bullet:\n")
  for (k in which(seen == 0)) cat("  ", substr(relb[k], 1, 140), "\n")
  cat("release bullets matched more than once:", sum(seen > 1), "\n")
}
