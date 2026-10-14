# Lane dblbar: sentences over 25 words in NEWS.md's development section
# (ASD-STE100). Code spans count as one word each.
n <- readLines("C:/Users/adf44/source/r/frmtmb-wt-dblbar/NEWS.md")
i <- grep("^# frmtmb 0.70.0", n)[1]
txt <- paste(n[seq_len(i - 1L)], collapse = " ")
txt <- gsub("`[^`]*`", "CODE", txt)
txt <- gsub("[*]", "", txt)
s <- unlist(strsplit(txt, "(?<=[.;:])\\s+(?=[A-Z(]|CODE)", perl = TRUE))
w <- vapply(strsplit(trimws(s), "\\s+"), length, 1L)
cat("sentences:", length(s), " over 25 words:", sum(w > 25), "\n")
for (k in which(w > 25)) cat(w[k], ":", substr(s[k], 1, 100), "\n")
