# Cut cox() out of core's compat registry (R/compat.R), saving the cut
# rules for frmtmb.survival's registration. Run once.
f <- "R/compat.R"
x <- readLines(f)
one <- function(pat) {
  i <- grep(pat, x, fixed = TRUE)
  if (length(i) != 1L) stop("expected one line for: ", pat)
  i
}
# the vocabulary and the two CDF groups
i <- one("\"sratio\", \"cratio\", \"acat\", \"categorical\", \"von_mises\", \"cox\",")
x[i] <- sub(" \"cox\",", "", x[i], fixed = TRUE)
i <- one("          \"inverse.gaussian\", \"cox\", \"frechet\", \"gen_extreme_value\",")
x[i] <- sub(" \"cox\",", "", x[i], fixed = TRUE)
i <- one("                     \"inverse.gaussian\", \"cox\", \"frechet\",")
x[i] <- sub(" \"cox\",", "", x[i], fixed = TRUE)
# the rules naming cox, each one r() line and its note
cut <- character()
for (pat in c("  r(\"cens()\", \"cox\", \"works\",",
              "  r(\"trunc()\", \"cox\", \"conditional\",",
              "  r(\"mvbf\", \"cox\", \"refused\",",
              "  r(\"simulate\", \"cox\", \"refused\",",
              "  r(\"predict\", \"cox\", \"conditional\",",
              "  r(\"fitted\", \"cox\", \"refused\",")) {
  i <- one(pat)
  stopifnot(grepl("\")$", x[i + 1L]))
  cut <- c(cut, x[i], x[i + 1L])
  x[c(i, i + 1L)] <- NA_character_
}
x <- x[!is.na(x)]
i <- one("  # that gap off the family registry and finds cox alone, which has its")
stopifnot(grepl("own rule below", x[i + 1L]))
x[i] <- "  # that gap off the family registry and finds none: cox(), which"
x[i + 1L] <- paste0("  # has no simulator, is in frmtmb.survival and brings its own",
                    " rule. The")
x[i + 2L] <- sub("^  # own rule below. The group this replaces also listed tweedie, compois",
                 "  # group this replaces also listed tweedie, compois", x[i + 2L])
writeLines(x, f)
writeLines(cut, "dev/splitsurv-out/cox-compat-cut.R")
cat(length(cut), "lines cut\n")
