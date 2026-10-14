# Generated blocks for section 9 (fix round 1) of
# dev/famsimplex-findings.md. Output: dev/famsimplex-summary-fix1.md
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
out <- character(0)
block <- function(title, lines) {
  out <<- c(out, paste0("<!-- BEGIN generated: ", title, " -->"), "```",
            lines, "```", "<!-- END generated -->", "")
}
rd <- function(f) readLines(file.path(wt, f), warn = FALSE)
strip <- function(x) sub(" [|] lib: .*$", "", x)

p <- rd("famsimplex-pfix-fixed.log")
keep <- grep("^---|^logLik|^ERROR", p, value = TRUE)
block("constants on locations, fixed build (dev/famsimplex-rev-pfix.R lane)",
      substr(keep, 1, 110))

block("pre-fix variant (dev/famsimplex-prefix.R off, dev/famsimplex-run.sh prefix2)",
      c(strip(rd("famsimplex-tests/prefix2.log")),
        substr(sub("^.*DETAIL ", "", grep("DETAIL", readLines(
          list.files(file.path(wt, "famsimplex-tests", "prefix2"),
                     full.names = TRUE)), value = TRUE)), 1, 150)))

block("scaled density error (dev/famsimplex-ulps.R)",
      rd("famsimplex-ulps.log")[-1L])

f1 <- rd("famsimplex-tests/fix1.log")
res <- grep("^RESULT", f1, value = TRUE)
jobs <- readLines(file.path(wt, "famsimplex-tests-fix1.jobs"))
pk <- sub(" .*$", "", jobs)
num <- function(k) {
  as.integer(sub(paste0("^.* ", k, "=([0-9]+).*$"), "\\1", res))
}
tab <- aggregate(cbind(pass = num("pass"), fail = num("fail"),
                       err = num("err"), skip = num("skip"),
                       warn = num("warn")) ~ pk, FUN = sum)
tab$files <- as.integer(table(pk)[tab$pk])
block("parser-reaching files (dev/famsimplex-run.sh fix1)",
      c(utils::capture.output(print(tab, row.names = FALSE)), tail(f1, 1),
        strip(grep("fail=[1-9]|err=[1-9]|skip=[1-9]|warn=[1-9]", res,
                   value = TRUE))))
block("final build, simplex and gated parser tiers (dev/famsimplex-run.sh fix2)",
      strip(rd("famsimplex-tests/fix2.log")))
block("OpenBLAS 0.3.26 LAPACK (dev/famsimplex-run-ob.sh ob326b)",
      strip(rd("famsimplex-tests/ob326b.log")))
for (pkg in c("frmtmb", "frmtmb.sample")) {
  f <- file.path(wt, "famsimplex-check", paste0("check-", pkg, ".log"))
  if (file.exists(f)) {
    l <- readLines(f, warn = FALSE)
    block(paste0("R CMD check ", pkg, " (dev/famsimplex-check.sh)"),
          grep("Status:|NOTE|WARNING|ERROR|Running 'testthat", l,
               value = TRUE))
  }
}
writeLines(out, file.path(wt, "famsimplex-summary-fix1.md"))
cat("wrote", length(out), "lines\n")
