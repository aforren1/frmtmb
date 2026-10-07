# Generate the verification block of dev/round-20261007.md from the logs
# of the full runs (dev/ciharden-scan.sh, one test file per process,
# every gate on: NOT_CRAN, FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS,
# FRMTMB_FUZZ and FRMTMB_SCALE_TESTS). For each configuration, a rerun
# directory (the files run again after the ledger regeneration) replaces
# the first result of the files it holds, and both are listed.
#
#   Rscript dev/rel069-counts.R > dev/rel069-log/counts.md
L <- "dev/ciharden-log"
runs <- list(
  list(name = "reference BLAS and LAPACK", run = "r069-ref",
       rerun = "r069-ref-rr"),
  list(name = "OpenBLAS 0.3.32 as BLAS and LAPACK, 4 threads",
       run = "r069-ob32", rerun = "r069-ob32-rr"))
lib_ok <- "rellib-r7"
pk_order <- c("frmtmb", "frmtmb.eam", "frmtmb.sample", "frmtmb.coupling",
              "frmtmb.learn", "frmtmb.latent", "frmtmb.ode",
              "frmtmb.spline")
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
    det <- grep("^DETAIL\t", x, value = TRUE)
    list(pkg = pk, file = sub("^.*--", "", sub("[.]txt$", ".R", basename(f))),
         ok = length(res) == 1L && !grepl("LOADERROR", res),
         lib_ok = length(lib) == 1L &&
           lengths(regmatches(lib, gregexpr(lib_ok, lib))) == 2L,
         pass = num("pass"), fail = num("fail"), err = num("err"),
         skip = num("skip"), warn = num("warn"), detail = det)
  })
  names(out) <- basename(fs)
  out
}
for (r in runs) {
  a <- read_dir(r$run)
  b <- if (dir.exists(file.path(L, paste0(r$rerun, "-suite")))) {
    read_dir(r$rerun)
  } else list()
  final <- a
  for (k in names(b)) final[[k]] <- b[[k]]
  jobs <- readLines(file.path(L, paste0(r$run, ".jobs")))
  cat("### All eight suites,", r$name, "\n\n")
  cat("From `dev/ciharden-log/", r$run, "-suite/` (", length(jobs),
      " jobs; ", sum(vapply(a, `[[`, NA, "ok")), " files with a RESULT line; ",
      sum(vapply(a, `[[`, NA, "lib_ok")), " whose `lib:` line names ",
      lib_ok, " for the package and for frmtmb).\n\n", sep = "")
  cat("| package | files | pass | fail | error | skip | warn |\n")
  cat("|---|---|---|---|---|---|---|\n")
  tot <- c(0, 0, 0, 0, 0, 0)
  for (p in pk_order) {
    x <- Filter(function(z) z$pkg == p, final)
    v <- c(length(x), sum(vapply(x, `[[`, 0L, "pass")),
           sum(vapply(x, `[[`, 0L, "fail")), sum(vapply(x, `[[`, 0L, "err")),
           sum(vapply(x, `[[`, 0L, "skip")), sum(vapply(x, `[[`, 0L, "warn")))
    tot <- tot + v
    cat("|", p, "|", paste(v, collapse = " | "), "|\n")
  }
  cat("| **total** |", paste0("**", tot, "**", collapse = " | "), "|\n\n")
  sk <- Filter(function(z) isTRUE(z$skip > 0), final)
  if (length(sk)) {
    cat("Every skip, by file and test:\n\n")
    for (z in sk) {
      d <- grep("^DETAIL\tSKIP", z$detail, value = TRUE)
      s <- vapply(strsplit(d, "\t"), function(v) paste0("[", v[4], "] ",
                                                        v[5]), "")
      cat("- ", z$pkg, " `", z$file, "` (", z$skip, "): ",
          paste(s, collapse = " | "), "\n", sep = "")
    }
    cat("\n")
  } else cat("Skips: 0.\n\n")
  if (length(b)) {
    cat("Files run again after the ledger regeneration (`",
        r$rerun, "-suite/`):\n\n", sep = "")
    for (k in names(b)) {
      f <- a[[k]]; g <- b[[k]]
      cat(sprintf("- %s `%s`: first pass %s fail %s err %s, rerun pass %s fail %s err %s warn %s\n",
                  g$pkg, g$file, f$pass %||% NA, f$fail %||% NA,
                  f$err %||% NA, g$pass, g$fail, g$err, g$warn))
    }
    cat("\n")
  }
  bad <- Filter(function(z) !isTRUE(z$ok) || !isTRUE(z$fail == 0) ||
                  !isTRUE(z$err == 0) || !isTRUE(z$warn == 0), final)
  cat("Files with a failure, an error or an escaped warning in the final",
      "results:", length(bad), "\n\n")
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
