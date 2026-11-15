# Lane famend: the counts the findings quote, generated from the logs.
#   Rscript dev/famend-summary.R
o <- "dev/famend-out"
rd <- function(f) readLines(file.path(o, f), warn = FALSE)
num <- function(x) as.numeric(x)

cat("== frechet probes (dev/famend-probe-frechet.R), lane build\n")
for (f in c("frechet-n300-lane.txt", "frechet-near-lane.txt")) {
  L <- grep("^seed ", rd(f), value = TRUE)
  nu1 <- num(sub(".*min\\(nu - 1\\) ([^ ]+) .*", "\\1", L))
  dll <- num(sub(".* dll ([^ ]+) gain.*", "\\1", L))
  gain <- num(sub(".* gain ([^ ]+) .*", "\\1", L))
  warned <- grepl("frechet: nu ran to", L)
  edge <- gain >= 0
  cat(sprintf(paste0("%s: fits %d; gain >= 0 (edge) %d, warned %d; edge ",
                     "nu - 1 %s, gain %s, slope in nu alone %s; interior ",
                     "%d: min nu - 1 %.3g, largest gain %.3g\n"),
              f, length(L), sum(edge), sum(warned),
              paste(signif(nu1[edge], 3), collapse = " "),
              paste(signif(gain[edge], 3), collapse = " "),
              paste(signif(dll[edge], 3), collapse = " "),
              sum(!edge), min(nu1[!edge]), max(gain[!edge])))
}

cat("\n== csmix (dev/famend-repro.R csmix)\n")
for (f in c("repro-base-3.txt", "repro-lane-3.txt")) {
  L <- grep("^csmix seed", rd(f), value = TRUE)
  code <- num(sub(".* code ([0-9]+) .*", "\\1", L))
  gap <- num(sub(".* min gap ([^ ]+) .*", "\\1", L))
  named <- grepl("named TRUE", L)
  cat(sprintf(paste0("%s: fits %d; code 1 %d (gaps %.3g to %.3g), named %d ",
                     "of them; code 0 %d (gaps %s), named %d\n"),
              f, length(L), sum(code == 1), min(gap[code == 1]),
              max(gap[code == 1]), sum(named & code == 1), sum(code == 0),
              paste(signif(gap[code == 0], 3), collapse = " "),
              sum(named & code == 0)))
}
