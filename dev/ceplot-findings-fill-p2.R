# Lane ceplot punch 2: paste the generated blocks of section 11 into
# dev/ceplot-findings.md.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
lg <- file.path(wt, "dev/ceplot-log")
f <- file.path(wt, "dev/ceplot-findings.md")
rd <- function(p) sub("\r$", "", readLines(file.path(lg, p), warn = FALSE))
gen <- function(...) {
  sub("\r$", "", system2(file.path(R.home("bin"), "Rscript"),
                         c(file.path(wt, "dev/ceplot-findings-gen.R"), ...),
                         stdout = TRUE))
}
check <- unlist(lapply(c("frmtmb", "frmtmb.sample"), function(p) {
  cl <- rd(file.path("..", "ceplot-check", p, paste0(p, ".Rcheck"),
                     "00check.log"))
  hit <- grep("^[*] checking .*[.][.][.] .*(NOTE|WARNING|ERROR)", cl)
  c(paste0(p, ": ", grep("^Status:", cl, value = TRUE)),
    if (length(hit)) paste0("  ", unlist(lapply(hit, function(i) {
      cl[i + 0:1]
    }))))
}))
blocks <- list(
  P2OLD = c(rd("p2-rev-oldlevels4.txt"), rd("p2-rev-oldlevels5.txt")),
  P2SUITE = gen("p2suite", "p2gated"),
  P2CHECK = check
)
s <- readLines(f, warn = FALSE)
for (k in names(blocks)) {
  i <- which(s == k)
  if (length(i) != 1L) stop("placeholder ", k, " found ", length(i), " times")
  s <- c(s[seq_len(i - 1L)], blocks[[k]], s[-seq_len(i)])
}
writeLines(s, f)
cat("filled\n")
