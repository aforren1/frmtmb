# Lane ceplot punch 1: paste the generated blocks of section 10 into
# dev/ceplot-findings.md, so no count in it is typed by hand.
#   Rscript dev/ceplot-findings-fill-p1.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
lg <- file.path(wt, "dev/ceplot-log")
f <- file.path(wt, "dev/ceplot-findings.md")
rd <- function(p) sub("\r$", "", readLines(file.path(lg, p), warn = FALSE))
gen <- function(...) {
  sub("\r$", "", system2(file.path(R.home("bin"), "Rscript"),
                         c(file.path(wt, "dev/ceplot-findings-gen.R"), ...),
                         stdout = TRUE))
}
mut <- unlist(lapply(list.files(file.path(lg, "p1mut"), full.names = TRUE),
                     function(p) {
  l <- sub("\r$", "", readLines(p, warn = FALSE))
  c(grep("^RESULT", l, value = TRUE), grep("^  fails:", l, value = TRUE))
}))
seen <- unlist(lapply(list.files(lg, "^p1-seen-failing-", full.names = TRUE),
                      function(p) {
  l <- sub("\r$", "", readLines(p, warn = FALSE))
  c(paste0(basename(p), ":"),
    paste0("  ", substr(grep("^── [0-9]+[.] (Failure|Error)", l,
                             value = TRUE), 1, 100)))
}))
check <- unlist(lapply(c("frmtmb", "frmtmb.sample"), function(p) {
  cl <- rd(file.path("..", "ceplot-check", p, paste0(p, ".Rcheck"),
                     "00check.log"))
  hit <- grep("^\\* checking .*\\.\\.\\. .*(NOTE|WARNING|ERROR)", cl)
  c(paste0(p, ": ", grep("^Status:", cl, value = TRUE)),
    if (length(hit)) paste0("  ", unlist(lapply(hit, function(i) {
      cl[i + 0:1]
    }))))
}))
noreg <- rd("p1-noreg-compare.txt")
blocks <- list(
  P1OLD = c(rd("p1-rev-oldlevels.txt"), rd("p1-rev-oldlevels2.txt")),
  P1NOREG = c(grep("DIFFERS", noreg, value = TRUE), tail(noreg, 1)),
  P1MUT = mut,
  P1SEEN = seen,
  P1FULL = gen("p1full"),
  P1GATED = gen("p1gated"),
  P1CHECK = check
)
s <- readLines(f, warn = FALSE)
for (k in names(blocks)) {
  i <- which(s == k)
  if (length(i) != 1L) stop("placeholder ", k, " found ", length(i), " times")
  s <- c(s[seq_len(i - 1L)], blocks[[k]], s[-seq_len(i)])
}
writeLines(s, f)
cat("filled\n")
