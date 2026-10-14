setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
o <- "dev/famuni-out/"
blk <- function(nm, x) c(paste0("<!-- BEGIN generated: ", nm, " -->"),
                         "```text", x, "```",
                         paste0("<!-- END generated: ", nm, " -->"))
res <- function(f) grep("^RESULT", readLines(file.path(o, f)), value = TRUE)
ck <- function() {
  fl <- "dev/famuni-check/frmtmb/check.log"
  if (!file.exists(fl)) return("not run")
  con <- file(fl, encoding = "UTF-16LE")
  x <- gsub("[^ -~]", "", readLines(con, warn = FALSE))
  close(con)
  s <- grep("^Status|NOTE$|WARNING$|ERROR$", x, value = TRUE)
  if (!length(s)) "still running" else s
}
blocks <- list(
  R2XIEDGE = function() readLines(file.path(o, "r2-xiedge-sum.txt")),
  R2TESTS = function() c(res("r2-t-extreme.txt"),
                         res("r2-t-test-gev-tmp-xi.R.txt"),
                         res("r2-t-test-default-priors-brms.R.txt"),
                         res("r2-t-wlinks.txt")),
  R2CHECK = function() c("[dev/famuni-check/frmtmb/check.log]", ck()))
f <- "dev/famuni-findings.md"
md <- readLines(f)
for (nm in names(blocks)) {
  b0 <- which(md == paste0("<!-- BEGIN generated: ", nm, " -->"))
  b1 <- which(md == paste0("<!-- END generated: ", nm, " -->"))
  i <- which(md == paste0("PASTE-", nm))
  new <- blk(nm, blocks[[nm]]())
  if (length(b0) == 1L) md <- c(md[seq_len(b0 - 1L)], new, md[-seq_len(b1)])
  else if (length(i) == 1L) md <- c(md[seq_len(i - 1L)], new, md[-seq_len(i)])
}
writeLines(md, f)
cat("filled\n")
