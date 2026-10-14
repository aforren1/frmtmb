f <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-findings.md"
x <- paste(readLines(f), collapse = "\n")
old <- "   `final-install-sample2.txt`); no test reads an Rd file, so the suite
   runs above stand for it."
new <- "   `final-install-sample2.txt`). No test reads a package Rd file (grep),
   so the suite runs above stand for it; four core files were rerun on
   the reinstall as a check (`final-core-*.txt`): PASTE-CORE4"
stopifnot(grepl(old, x, fixed = TRUE))
r <- vapply(c("test-families-extreme.R", "test-conditions-census.R",
              "test-compat.R", "test-brms-families.R"), function(t) {
  l <- grep("^RESULT", readLines(file.path(
    "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out",
    paste0("final-core-", t, ".txt"))), value = TRUE)
  sub("^RESULT ", "", l)
}, "")
x <- sub(old, sub("PASTE-CORE4", paste0("\n\n   ```text\n   ",
                                        paste(r, collapse = "\n   "),
                                        "\n   ```\n"), new, fixed = TRUE),
         x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], f)
