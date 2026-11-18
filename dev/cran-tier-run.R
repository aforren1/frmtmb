# Lane cran: run the proposed CRAN test tier in ONE session, the way
# R CMD check runs tests/testthat.R (NOT_CRAN unset), with every block
# on the moved list skipped, and report its time and results. With mode
# "trace", count the calls that reach each export of core as well.
#
#   Rscript dev/cran-tier-run.R <worktree> <moved.tsv> <time|trace> <out prefix>
#
# The skip is a test_that() placed in the test environment ahead of
# testthat's; a moved block is matched by file and description. This is
# a measurement of the proposal, not the mechanism the package would
# ship (see dev/cran-findings.md).
local({
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(!nzchar(Sys.getenv("NOT_CRAN")))
a <- commandArgs(trailingOnly = TRUE)
root <- a[1]
mv <- utils::read.delim(a[2], quote = "", comment.char = "")
mode <- a[3]
outp <- a[4]
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
key <- paste(mv$file, mv$desc, sep = "\r")

hits <- new.env()
if (identical(mode, "trace")) {
  ex <- sort(unique(parseNamespaceFile("frmtmb", dirname(root))$exports))
  for (e in ex) {
    # a re-export (fixef from nlme) has no binding of its own here
    if (!exists(e, envir = ns, inherits = FALSE) ||
        bindingIsActive(e, ns) || !is.function(get(e, ns))) next
    assign(e, 0L, envir = hits)
    tr <- bquote(assign(.(e), get(.(e), envir = .(hits)) + 1L,
                        envir = .(hits)))
    suppressMessages(trace(e, tracer = tr, where = ns, print = FALSE))
  }
}

env <- testthat::test_env("frmtmb")
cur <- new.env()
env$test_that <- function(desc, code) {
  if (paste(cur$file, desc, sep = "\r") %in% key) {
    return(invisible(testthat::test_that(desc,
                                         testthat::skip("cran tier: moved"))))
  }
  eval(substitute(testthat::test_that(desc, code)), parent.frame())
}
tdir <- file.path(root, "tests", "testthat")
files <- sort(list.files(tdir, "^test-.*[.]R$"))
t0 <- proc.time()
# test_dir() would source the helpers once but cannot tell the override
# which file runs, so the files go one by one with the helpers and the
# setup file sourced once into the shared environment, as test_check()
# does
th0 <- proc.time()
testthat::source_test_helpers(tdir, env = env)
testthat::source_test_setup(tdir, env = env)
th <- proc.time() - th0
per <- data.frame(file = files, real = NA_real_, cpu = NA_real_,
                  pass = NA, fail = NA, err = NA, skip = NA, warn = NA)
for (i in seq_along(files)) {
  cur$file <- files[i]
  ti <- proc.time()
  r <- tryCatch(
    as.data.frame(testthat::test_file(file.path(tdir, files[i]),
                                      reporter = "silent",
                                      package = "frmtmb", env = env,
                                      load_helpers = FALSE)),
    error = function(e) {
      cat("FILEERROR\t", files[i], "\t", conditionMessage(e), "\n", sep = "")
      NULL
    })
  d <- proc.time() - ti
  per$real[i] <- d[[3]]
  per$cpu[i] <- sum(d[c(1, 2)])
  if (!is.null(r)) {
    per$pass[i] <- sum(r$passed)
    per$fail[i] <- sum(r$failed)
    per$err[i] <- sum(r$error)
    per$skip[i] <- sum(r$skipped)
    per$warn[i] <- sum(r$warning)
  }
}
tt <- proc.time() - t0
cat("COUNTS\tpass=", sum(per$pass, na.rm = TRUE),
    "\tfail=", sum(per$fail, na.rm = TRUE),
    "\terr=", sum(per$err, na.rm = TRUE),
    "\tskipblocks=", sum(per$skip, na.rm = TRUE),
    "\twarn=", sum(per$warn, na.rm = TRUE),
    "\tfileerrors=", sum(is.na(per$pass)), "\n", sep = "")
utils::write.table(per[order(-per$real), ], paste0(outp, "-perfile.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("SESSION\tmode=", mode, "\tfiles=", length(files),
    "\thelpers_real=", sprintf("%.2f", th[[3]]),
    "\ttotal_real=", sprintf("%.2f", tt[[3]]),
    "\ttotal_cpu=", sprintf("%.2f", sum(tt[c(1, 2)])), "\n", sep = "")
if (identical(mode, "trace")) {
  n <- unlist(mget(ls(hits), envir = hits))
  utils::write.table(data.frame(export = names(n), calls = n),
                     paste0(outp, "-calls.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE)
  cat("TRACED\t", length(n), "\treached=", sum(n > 0), "\n", sep = "")
}
})
