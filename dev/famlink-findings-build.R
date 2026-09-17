# Assemble dev/famlink-findings.md from dev/famlink-findings.src.md, with
# every <<BLOCK:name>> replaced by generated output, so no count in the
# findings is typed. Run after the logs it reads exist:
#   Rscript dev/famlink-findings-build.R
rd <- function(f) readLines(f, warn = FALSE, encoding = "UTF-8")
from <- function(f, start) {
  x <- rd(f)
  i <- grep(start, x)
  if (!length(i)) stop(f, " has no line matching ", start)
  x[i[1]:length(x)]
}
fence <- function(x) c("```", x, "```")

consumers <- local({
  files <- c(list.files("R", "[.]R$", full.names = TRUE),
             Sys.glob("extensions/*/R/*.R"))
  pat <- '\\[\\["links"\\]\\]|\\$links\\b|"links"'
  hits <- unlist(lapply(files, function(f) {
    x <- rd(f)
    keep <- grepl(pat, x, perl = TRUE) & !grepl("^\\s*#", x)
    rep(f, sum(keep))
  }))
  tab <- sort(table(hits), decreasing = TRUE)
  c(sprintf("lines reading or writing `links`, comments excluded: %d in %d files",
            length(hits), length(tab)),
    sprintf("  %-40s %d", names(tab), as.integer(tab)))
})

suites <- unlist(lapply(
  c("dev/famlink-suite-core-log.txt",
    Sys.glob("dev/famlink-suite-frmtmb.*-log.txt")),
  function(f) {
    x <- rd(f)
    i <- grep("^---- GENERATED COUNTS", x)
    j <- grep("END GENERATED COUNTS", x)
    if (length(i) != 1L || length(j) != 1L) {
      return(paste(f, ": INCOMPLETE, no generated block"))
    }
    bad <- grep("CRASHED|fail +[1-9]|error +[1-9]", x, value = TRUE)
    c(x[(i + 1L):(j - 1L)], if (length(bad)) paste("  !", bad), "")
  }))

check <- if (file.exists("dev/famlink-check-log.txt")) {
  x <- rd("dev/famlink-check-log.txt")
  grep("^Status|WARNING|ERROR|NOTE|^\\* checking .* \\.\\.\\. (NOTE|WARNING|ERROR)",
       x, value = TRUE)
} else "not run"

blocks <- list(
  "defects-base" = fence(from("dev/famlink-defects-base-log.txt", "==== summary")),
  "defects-lane" = fence(from("dev/famlink-defects-lane-log.txt", "==== summary")),
  "consumers" = fence(consumers),
  "refusals-base" = fence(rd("dev/famlink-refusals-base-log.txt")),
  "refusals-lane" = fence(rd("dev/famlink-refusals-lane-log.txt")),
  "falsealarm" = fence(rd("dev/famlink-falsealarm-sum-log.txt")),
  "ledger" = rd("dev/famlink-ledger-sum-log.txt"),
  "newtest" = fence(c(
    paste("base:", grep("^FAMLINK", rd("dev/famlink-newtest-base-log.txt"), value = TRUE)),
    paste("lane:", grep("^FAMLINK", rd("dev/famlink-newtest-lane-log.txt"), value = TRUE)))),
  "suites" = fence(suites),
  "insightscan" = fence(rd("dev/famlink-punch-insightscan-log.txt")),
  "view" = fence(rd("dev/famlink-punch-view-log.txt")),
  "nlminbinf" = fence(rd("dev/famlink-1b-nlminb-inf-log.txt")),
  "nancheck" = fence(rd("dev/famlink-1b-nancheck-log.txt")),
  "nansurvey" = fence(from("dev/famlink-1b-nansurvey-after-log.txt",
                           "^---- GENERATED")),
  "invgauss-base" = fence(rd("dev/famlink-punch-invgauss-base-log.txt")),
  "invgauss-lane" = fence(rd("dev/famlink-punch-invgauss-lane-log.txt")),
  "check" = fence(check),
  "p2-invgauss" = fence(c(
    rd("dev/famlink-rev2-invgauss-sum-log.txt")[1:25],
    grep("^clean replicates", rd("dev/famlink-rev2-invgauss-sum-log.txt"),
         value = TRUE))),
  "p2-undercount-before" = fence(grep(
    "nlminb calls", rd("dev/famlink-p2-undercount-before-log.txt"),
    value = TRUE)),
  "p2-undercount-after" = fence(grep(
    "nlminb calls", rd("dev/famlink-p2-undercount-after-log.txt"),
    value = TRUE)),
  "p2-newtest" = fence(c(
    paste("unfixed:", grep("^FAMLINK|`actual`|`expected`",
                           rd("dev/famlink-p2-newtest-unfixed-log.txt"),
                           value = TRUE)),
    paste("fixed:  ", grep("^FAMLINK",
                           rd("dev/famlink-p2-newtest-fixed-log.txt"),
                           value = TRUE)))),
  "p2-scale" = fence(grep("^FAMLINK|frmtmb.* from",
                          rd("dev/famlink-p2-scale-small-log.txt"),
                          value = TRUE)),
  "p2-scalebern" = fence(rd("dev/famlink-p2-scale-bern-log.txt")),
  "p2-scan" = fence(rd("dev/famlink-p2-trials-scan-before-log.txt")),
  "p2-fuzzbinom" = fence(grep("[^[:space:]]",
                              rd("dev/famlink-p2-fuzzbinom-log.txt"),
                              value = TRUE)),
  "p2-fuzz" = fence(grep("^FAMLINK|frmtmb from",
                         rd("dev/famlink-p2-fuzz-log.txt"), value = TRUE)),
  "p2-trials" = fence(c(
    paste("before:", rd("dev/famlink-rev2-trials-cmp-log.txt")[1:4]),
    paste("after: ", rd("dev/famlink-p2-trials-cmp-log.txt")[1:4]))),
  "p2-gptest" = fence(unique(grep(
    "^FAMLINK|^Error [(]", rd("dev/famlink-p2-gptest-unfixed-log.txt"),
    value = TRUE))),
  "p2-dollar" = fence(rd("dev/famlink-p2-dollarcost-log.txt"))
)

src <- rd("dev/famlink-findings.src.md")
out <- character(0)
for (ln in src) {
  m <- regmatches(ln, regexec("^<<BLOCK:([a-z0-9-]+)>>$", ln))[[1]]
  if (length(m)) {
    b <- blocks[[m[2]]]
    if (is.null(b)) stop("no block ", m[2])
    out <- c(out, b)
  } else {
    out <- c(out, ln)
  }
}
writeLines(out, "dev/famlink-findings.md", useBytes = TRUE)
cat("wrote dev/famlink-findings.md,", length(out), "lines\n")
