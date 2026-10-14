# Reviewer of lane optstart, round 2: sentences of the development
# NEWS section with more than 25 words (list items count as sentences,
# code spans count as one word).
x <- readLines("NEWS.md")
end <- grep("^# frmtmb 0.70.0", x)[1] - 1
x <- x[2:end]
x <- x[!grepl("^#", x)]
txt <- paste(gsub("^\\s*[-*]\\s+", "\n", x), collapse = " ")
txt <- gsub("`[^`]*`", "CODE", txt)
s <- unlist(strsplit(txt, "(?<=[.;:])\\s+(?=[A-Z*(])|\n", perl = TRUE))
s <- trimws(gsub("\\s+", " ", s))
s <- s[nzchar(s)]
n <- lengths(strsplit(s, " "))
cat(length(s), "sentences; over 25 words:", sum(n > 25), "\n\n")
for (i in which(n > 25)) cat(n[i], ":", s[i], "\n\n")
