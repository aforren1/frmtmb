# Reviewer round 2, lane dblbar: NEWS sentence lengths, split at full
# stops only (the lane's dev/dblbar-wordcount.R also splits at ; and :).
# Headings are dropped; code spans count as one word.
n <- readLines("C:/Users/adf44/source/r/frmtmb-wt-dblbar/NEWS.md")
i <- grep("^# frmtmb 0.70.0", n)[1]
n <- n[seq_len(i - 1L)]
n <- n[!grepl("^#", n)]
txt <- paste(n, collapse = " ")
txt <- gsub("`[^`]*`", "CODE", txt)
txt <- gsub("[*]", "", txt)
s <- unlist(strsplit(txt, "(?<=[.])\\s+(?=[A-Z(]|CODE)", perl = TRUE))
s <- s[nzchar(trimws(s))]
w <- vapply(strsplit(trimws(s), "\\s+"), length, 1L)
cat("sentences:", length(s), " max words:", max(w), " over 25:",
    sum(w > 25), " over 20:", sum(w > 20), "\n")
for (k in which(w > 20)) cat(w[k], ":", substr(s[k], 1, 110), "\n")
