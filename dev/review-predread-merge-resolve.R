# Item 8: resolve the textual conflicts of the trial merge the only way that keeps
# all three lanes' intent, so the merged tree can be installed and run.
M <- "C:/Users/adf44/source/r/rev-predread-src/merge"
res <- function(f, pick) {
  p <- file.path(M, f)
  x <- readLines(p, warn = FALSE)
  out <- character()
  i <- 1L
  while (i <= length(x)) {
    if (startsWith(x[i], "<<<<<<< ")) {
      j1 <- i + which(startsWith(x[(i + 1L):length(x)], "||||||| "))[1L]
      j2 <- i + which(startsWith(x[(i + 1L):length(x)], "======="))[1L]
      j3 <- i + which(startsWith(x[(i + 1L):length(x)], ">>>>>>> "))[1L]
      ours <- if (j1 > i + 1L) x[(i + 1L):(j1 - 1L)] else character()
      theirs <- if (j3 > j2 + 1L) x[(j2 + 1L):(j3 - 1L)] else character()
      out <- c(out, pick(ours, theirs))
      i <- j3 + 1L
    } else {
      out <- c(out, x[i])
      i <- i + 1L
    }
  }
  writeLines(out, p)
}
# predict.R: latres's early returns, then predread's test
res("R/predict.R", function(ours, theirs) {
  lst <- length(theirs)
  c(theirs[-lst], ours)
})
# compat.R: nlfix's row text (the row is about nl)
res("R/compat.R", function(ours, theirs) theirs)
for (f in c("NEWS.md", "extensions/frmtmb.sample/NEWS.md")) {
  res(f, function(ours, theirs) c(ours, theirs))
}
left <- system2("grep", c("-l", shQuote("^<<<<<<< "), "-r", shQuote(file.path(M, "R")),
                          shQuote(file.path(M, "extensions/frmtmb.sample/R"))),
                stdout = TRUE)
cat("files with markers left:", length(left), "\n")
print(left)
