# Lane illse2 (from dev/illse-rev2-news.R): sentences over 25 words in
# the development-version section of NEWS.md. Code spans count as one
# word; decimals do not split sentences; a list item starts a new
# sentence.
f <- "C:/Users/adf44/source/r/frmtmb-wt-illse2/NEWS.md"
x <- readLines(f)
e <- grep("^# ", x)[2] - 1
x <- x[2:e]
x <- x[!grepl("^#", x)]
x <- sub("^[[:space:]]*[*-][[:space:]]+", "\u00a7 ", x)
txt <- paste(x, collapse = " ")
txt <- gsub("`[^`]*`", "CODE", txt)
txt <- gsub("[0-9]+[.][0-9]+", "N", txt)
txt <- gsub("[*][*]", "", txt)
s <- trimws(unlist(strsplit(txt, "(?<=[.!?])[[:space:]]+(?=[A-Z0-9(\"])|\u00a7",
                            perl = TRUE)))
s <- s[nzchar(s)]
n <- vapply(strsplit(s, "[[:space:]]+"), function(w) sum(nzchar(w)), 1L)
cat("sentences", length(s), "| over 25 words", sum(n > 25),
    "| longest", max(n), "\n")
for (i in which(n > 25)) cat(sprintf("  [%d] %s\n", n[i], substr(s[i], 1, 150)))
