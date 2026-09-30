# Lane ceplot: paste the generated blocks into dev/ceplot-findings.md,
# so no count in it is typed by hand. Each placeholder line is replaced
# by the named log's content (or by the summarizer's output).
#   Rscript dev/ceplot-findings-fill.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
lg <- file.path(wt, "dev/ceplot-log")
f <- file.path(wt, "dev/ceplot-findings.md")
rd <- function(p) sub("\r$", "", readLines(file.path(lg, p), warn = FALSE))
gen <- function(...) {
  sub("\r$", "", system2(file.path(R.home("bin"), "Rscript"),
                         c(file.path(wt, "dev/ceplot-findings-gen.R"), ...),
                         stdout = TRUE))
}
flips <- rd("flips-lane.txt")
flips <- flips[grepl("^lane ", flips)]
flips <- sub("^lane ([^ ]+) : brms (brmsfit-methods:[0-9]+) ", "- \\1 \\2: ",
             flips)
blocks <- list(
  "CROSSED-BRMS" = rd("crossed-brms.txt"),
  "MMBY" = rd("mmby-lane.txt"),
  "OLDLEVELS" = rd("oldlevels.txt"),
  "SEENFAIL" = rd("seen-failing-summary.txt"),
  "FULL" = gen("full2"),
  "GATED" = gen("gated2"),
  "CHECK" = unlist(lapply(c("frmtmb", "frmtmb.sample"), function(p) {
    cl <- rd(file.path("..", "ceplot-check", p, paste0(p, ".Rcheck"),
                       "00check.log"))
    hit <- grep("^\\* checking .*\\.\\.\\. .*(NOTE|WARNING|ERROR)", cl)
    c(paste0(p, ": ", grep("^Status:", cl, value = TRUE)),
      paste0("  ", unlist(lapply(hit, function(i) cl[i + 0:1]))))
  })),
  "FLIPS" = c(flips, "",
              "On the base build both methods files hold every row as",
              "recorded (dev/ceplot-log/gated2.log, the two `base` lines).")
)
s <- readLines(f, warn = FALSE)
for (k in names(blocks)) {
  i <- which(s == k)
  if (length(i) != 1L) stop("placeholder ", k, " found ", length(i), " times")
  s <- c(s[seq_len(i - 1L)], blocks[[k]], s[-seq_len(i)])
}
writeLines(s, f)
cat("filled\n")
