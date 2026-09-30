# Punch round 1: replace the development-version section of NEWS.md
# with dev/formrobust-news-p1.md. Record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/"
p <- paste0(wt, "NEWS.md")
x <- readLines(p)
i <- which(x == "# frmtmb 0.66.0")
stopifnot(length(i) == 1L)
new <- readLines(paste0(wt, "dev/formrobust-news-p1.md"))
con <- file(p, "wb")
writeLines(c(new, "", x[i:length(x)]), con, sep = "\r\n")
close(con)
