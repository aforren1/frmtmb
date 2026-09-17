# Reviewer, lane wt-conditions: run ONE test file in ONE process on one
# library arm, and record message, call and class of every condition an
# expect_error(), expect_warning(), expect_message() or
# expect_condition() caught, keyed so the two arms can be joined.
#   Rscript dev/conditions-rev-sweep-onefile.R lane|base <file> <outdir>
av <- commandArgs(trailingOnly = TRUE)
arm <- av[1L]; file <- av[2L]; outdir <- av[3L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             mustWork = FALSE))
suppressMessages(library(testthat))
pkg <- regmatches(file, regexpr("frmtmb[.][a-z]+", file))
if (!length(pkg)) pkg <- "frmtmb"
suppressMessages(library(pkg, character.only = TRUE))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")

rev_env <- new.env()
rev_env$rows <- list()
rev_env$seen <- list()
one <- function(x) paste(gsub("\r", "", x), collapse = "\\n")
# the text try() stores, built by try()'s own code from the condition
try_text <- function(e) {
  call <- conditionCall(e)
  if (!is.null(call)) {
    dcall <- deparse(call, nlines = 1L)
    prefix <- paste("Error in", dcall, ": ")
    LONG <- 75L
    sm <- strsplit(conditionMessage(e), "\n")[[1L]]
    w <- 14L + nchar(dcall, type = "w") + nchar(sm[1L], type = "w")
    if (is.na(w)) w <- 14L + nchar(dcall, type = "b") +
      nchar(sm[1L], type = "b")
    if (w > LONG) prefix <- paste0(prefix, "\n  ")
  } else prefix <- "Error : "
  paste0(prefix, conditionMessage(e), "\n")
}
rev_record <- function(fn, obj, v) {
  if (!inherits(v, "condition")) return(invisible())
  key <- paste(fn, one(deparse(obj, width.cutoff = 200L)))
  n <- (rev_env$seen[[key]] %||% 0L) + 1L
  rev_env$seen[[key]] <- n
  rev_env$rows[[length(rev_env$rows) + 1L]] <- data.frame(
    file = file, key = paste0(key, " #", n),
    class = paste(class(v), collapse = "/"),
    message = one(conditionMessage(v)),
    call = one(deparse(conditionCall(v), width.cutoff = 500L)),
    trytext = one(try_text(v)),
    stringsAsFactors = FALSE)
}
for (fn in c("expect_error", "expect_warning", "expect_message",
             "expect_condition")) {
  for (w in list(asNamespace("testthat"),
                 as.environment("package:testthat"))) {
    suppressMessages(trace(fn, where = w, print = FALSE,
      exit = bquote(.(rev_record)(.(fn), substitute(object),
                                  returnValue(default = NULL)))))
  }
}

res <- test_file(file, reporter = "silent", package = pkg)
df <- as.data.frame(res)
cat(sprintf("BLOCKS %d PASS %d FAIL %d ERROR %d SKIP %d  %s\n",
            nrow(df), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), file))
rows <- do.call(rbind, rev_env$rows)
if (is.null(rows)) rows <- data.frame(file = character(0),
  key = character(0), class = character(0), message = character(0),
  call = character(0), trytext = character(0))
cat("CAUGHT", nrow(rows), "\n")
saveRDS(rows, file.path(outdir, paste0(gsub("/", "_", file), ".rds")))
cat("DONE\n")
