# The tests section of dev/splitsurv-findings.md, from the logs.
#   Rscript dev/splitsurv-tests-md.R > dev/splitsurv-out/tests.md
sumlog <- function(path) {
  x <- readLines(path)
  r <- grep("^[^ ]+ RESULT ", x, value = TRUE)
  pkg <- sub(" RESULT.*", "", r)
  num <- function(k) as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
  data.frame(pkg = pkg, file = sub(".* RESULT ([^ ]+) .*", "\\1", r),
             pass = num("pass"), fail = num("fail"), err = num("err"),
             skip = num("skip"), warn = num("warn"))
}
cat("### Every test file, one per process, lane build\n\n")
cat("`dev/splitsurv-par2.sh` over `dev/splitsurv-after.lst` (every test file\n")
cat("of core and of the seven extensions, frmtmb.spline retired),\n")
cat("`NOT_CRAN=true`, library `wt-splitsurv-lib` first, run whole after\n")
cat("fix round 1, which changed paths every suite reaches (frame\n")
cat("assembly, the objective, the standard errors). Every log's `lib:`\n")
cat("line names the lane's frmtmb.\n\n")
a <- sumlog("dev/splitsurv-after-out/summary.log")
lst <- readLines("dev/splitsurv-after.lst")
cat(sprintf("Files: %d listed, %d with a RESULT line.\n\n", length(lst),
            nrow(a)))
agg <- aggregate(cbind(files = 1L, pass, fail, err, skip, warn) ~ pkg, a, sum)
cat("| package | files | pass | fail | err | skip | warn |\n|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(agg))) {
  with(agg[i, ], cat(sprintf("| %s | %d | %d | %d | %d | %d | %d |\n", pkg,
                             files, pass, fail, err, skip, warn)))
}
miss <- setdiff(basename(sub("^[^ ]+ ", "", lst)), a$file)
cat("\nWithout a RESULT line (every block gated off, as on the base):",
    paste(sort(unique(miss)), collapse = ", "), "\n\n")
bad <- a[a$fail > 0 | a$err > 0 | a$warn > 0, ]
cat("Files with a failure, an error or an escaped warning:",
    if (nrow(bad)) paste(bad$pkg, bad$file, collapse = "; ") else "none", "\n\n")
b <- sumlog("dev/splitsurv-before-out/summary.log")
cat("The base run (rellib-r11, base tree, same driver): ",
    nrow(b), " files with a RESULT line; failing: ",
    paste(b$file[b$fail > 0 | b$err > 0], collapse = ", "),
    " (a 0.73.1 test of the `fitted()` clamp, which rellib-r11's 0.73.0 ",
    "lacks). The other ", length(readLines("dev/splitsurv-before.lst")) -
      nrow(b), " listed files are the gated brms-suite files, on which ",
    "that first runner (`dev/splitsurv-run-tests.R`) stopped at an empty ",
    "result; `dev/splitsurv-run-tests2.R` reports them.\n\n",
    sep = "")
cat("Every new or moved test file, with its counts:\n\n")
cat("| package | file | pass | fail | err | skip | warn |\n|---|---|---|---|---|---|---|\n")
mine <- a[a$pkg == "frmtmb.survival" |
            a$file %in% c(grep("^test-curve", a$file, value = TRUE),
                          "test-simplex-pars.R", "test-cox-home.R",
                          "test-extension-home.R",
                          "test-simplex-family-notes.R",
                          "test-famgaps.R", "test-brms-families.R",
                          "test-hurdle-nb-zoib.R", "test-mv-gaps.R",
                          "test-simulate-density.R", "test-update-pool.R",
                          "test-compat.R", "test-compat-extensions.R",
                          "test-arg-refusal.R", "test-cox-draws.R",
                          "test-simulators.R", "test-ci-siblings.R"), ]
mine <- mine[order(mine$pkg, mine$file), ]
for (i in seq_len(nrow(mine))) {
  with(mine[i, ], cat(sprintf("| %s | %s | %d | %d | %d | %d | %d |\n", pkg,
                              file, pass, fail, err, skip, warn)))
}
cat("\n### Tests, accounted\n\n")
cat("`dev/splitsurv-account.R`, matching tests by name between the base\n")
cat("run and the lane run. A change in a core file's count is the moved\n")
cat("assertions; the list after the table says where each went.\n\n")
cat(readLines("dev/splitsurv-out/account.md"), sep = "\n")
cat("\n\nThe 21 assertions core's `test-brms-families.R` lost are cox()'s\n")
cat("rows of three registry loops: 4 links built, 13 roster links\n")
cat("refused (the roster has 17 and cox() takes 4), and 4 of \"link and\n")
cat("link_<dpar> agree with links on every registry family\". The first\n")
cat("record said 4 + 17 and missed the third test (review finding 1).\n")
cat("frmtmb.survival's `test-cox-links.R` runs all three for cox().\n")
cat("\n### One session per suite, three orders\n\n")
cat("`bash dev/release/run-sessions.sh splitsurv 8` with `FRMTMB_RELLIB` the\n")
cat("lane library and `NOT_CRAN=true` (log `dev/release/splitsurv.log`,\n")
cat("copied to `dev/splitsurv-out/sessions-all.log`). LEAK lines without\n")
cat("`[with load of ...]`: none.\n\n```\n")
s <- readLines("dev/splitsurv-out/sessions-all.log")
cat(grep("^SESSION|ran .* of", s, value = TRUE), sep = "\n")
cat("```\n\nTagged LEAK lines that name a package this lane added:\n\n```\n")
cat(grep("frmtmb.survival\\]", s, value = TRUE), sep = "\n")
cat("```\n\n`files=262/272` and `59/64` are the 10 and 5 ported brms-suite\n")
cat("files, gated off without `FRMTMB_BRMS_FIT_TESTS`, which report no\n")
cat("RESULT line.\n")
cat("\n### OpenBLAS\n\n")
cat("`dev/ciharden-openblas.sh 0.3.26 lapack`, then the new numerical\n")
cat("tests and the moved curve and cox() tests through that R\n")
cat("(`dev/splitsurv-par3.sh`, `OPENBLAS_NUM_THREADS=4`):\n\n```\n")
cat(readLines("dev/splitsurv-ob-out/summary.log"), sep = "\n")
cat("```\n")
