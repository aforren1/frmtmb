# Generate the verification blocks of dev/round-20261028.md from the
# logs of the per-file runs (dev/ciharden-scan.sh, one test file per R
# process). Adapted from dev/rel076-counts.R: the runs are named on the
# command line, each with its label, and a file of a run whose summary
# lists it twice is read once.
#
#   Rscript dev/rel077-counts.R <run> <label> [<run> <label> ...]
# REL077_LIBOK names the library every lib: line must show
# (rellib-r15). Skips are listed by file and test when a run has fewer
# than 40 files with one.
`%||%` <- function(a, b) if (is.null(a)) b else a
L <- "dev/ciharden-log"
a <- commandArgs(TRUE)
runs <- lapply(seq(1, length(a), by = 2), function(i) {
  list(run = a[i], name = a[i + 1])
})
lib_ok <- Sys.getenv("REL077_LIBOK", "rellib-r15")
pk_order <- c("frmtmb", "frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
              "frmtmb.learn", "frmtmb.ode", "frmtmb.sample",
              "frmtmb.spatial", "frmtmb.survival")
read_dir <- function(d) {
  fs <- list.files(file.path(L, paste0(d, "-suite")), "[.]txt$",
                   full.names = TRUE)
  out <- lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    pk <- sub("--.*$", "", basename(f))
    res <- tail(grep("^RESULT ", x, value = TRUE), 1)
    lib <- grep("^lib: ", x, value = TRUE)
    num <- function(k) {
      if (!length(res)) return(NA_integer_)
      as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", res))
    }
    list(pkg = pk, file = sub("^.*--", "", sub("[.]txt$", ".R", basename(f))),
         ok = length(res) == 1L && !grepl("LOADERROR", res),
         lib_ok = length(lib) == 1L &&
           lengths(regmatches(lib, gregexpr(lib_ok, lib))) == 2L,
         pass = num("pass"), fail = num("fail"), err = num("err"),
         skip = num("skip"), warn = num("warn"),
         detail = grep("^DETAIL\t", x, value = TRUE))
  })
  names(out) <- basename(fs)
  out
}
for (r in runs) {
  final <- read_dir(r$run)
  # a run given its own job list writes no .jobs file; its summary's
  # last line, "ran X of Y", has the count
  jf <- file.path(L, paste0(r$run, ".jobs"))
  jobs <- if (file.exists(jf)) readLines(jf) else {
    s <- tail(readLines(file.path(L, paste0(r$run, ".sum"))), 1)
    seq_len(as.integer(sub("^ran [0-9]+ of ([0-9]+)$", "\\1", s)))
  }
  cat("### ", r$name, "\n\n", sep = "")
  cat("From `dev/ciharden-log/", r$run, "-suite/` (", length(jobs),
      " jobs; ", sum(vapply(final, `[[`, NA, "ok")),
      " files with a RESULT line; ",
      sum(vapply(final, `[[`, NA, "lib_ok")), " whose `lib:` line names ",
      lib_ok, " for the package and for frmtmb).\n\n", sep = "")
  cat("| package | files | pass | fail | error | skip | warn |\n")
  cat("|---|---|---|---|---|---|---|\n")
  tot <- numeric(6)
  for (p in pk_order) {
    x <- Filter(function(z) z$pkg == p, final)
    if (!length(x)) next
    v <- c(length(x), sum(vapply(x, `[[`, 0L, "pass")),
           sum(vapply(x, `[[`, 0L, "fail")), sum(vapply(x, `[[`, 0L, "err")),
           sum(vapply(x, `[[`, 0L, "skip")), sum(vapply(x, `[[`, 0L, "warn")))
    tot <- tot + v
    cat("|", p, "|", paste(v, collapse = " | "), "|\n")
  }
  cat("| **total** |", paste0("**", tot, "**", collapse = " | "), "|\n\n")
  sk <- Filter(function(z) isTRUE(z$skip > 0), final)
  cat("Files with a skip: ", length(sk), ".\n\n", sep = "")
  if (length(sk) && length(sk) < 40) {
    for (z in sk) {
      d <- grep("^DETAIL\tSKIP", z$detail, value = TRUE)
      s <- vapply(strsplit(d, "\t"), function(v) paste0("[", v[4], "] ",
                                                        substr(v[5], 1, 90)), "")
      cat("- ", z$pkg, " `", z$file, "` (", z$skip, "): ",
          paste(unique(s), collapse = " | "), "\n", sep = "")
    }
    cat("\n")
  }
  bad <- Filter(function(z) !isTRUE(z$ok) || !isTRUE(z$fail == 0) ||
                  !isTRUE(z$err == 0) || !isTRUE(z$warn == 0), final)
  cat("Files with a failure, an error, an escaped warning or no RESULT",
      "line:", length(bad), "\n\n")
  for (z in bad) {
    cat("- ", z$pkg, " `", z$file, "`: pass ", z$pass, " fail ", z$fail,
        " err ", z$err, " warn ", z$warn, "\n", sep = "")
    for (d in grep("^DETAIL\t(FAIL|ERROR|WARN)", z$detail, value = TRUE)) {
      v <- strsplit(d, "\t")[[1]]
      cat("    - ", v[2], " ", v[3], " [", v[4], "] ",
          substr(v[5], 1, 300), "\n", sep = "")
    }
  }
  cat("\n")
}
