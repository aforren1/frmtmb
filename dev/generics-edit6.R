root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}
p <- "tests/testthat/test-generic-collision.R"
sub1(p,
  "    \"cat('FIXEF:', names(fixef(fit)), '\\\\n')\",",
  "    \"cat('FIXEF:', paste(names(fixef(fit)), collapse = ','), '\\\\n')\",")
sub1(p,
  "  expect_equal(parse_field(out, \"FIXEF\"), \"mu\")",
  paste0("  expect_equal(parse_field(out, \"FIXEF\"), c(\"mu\", \"sigma\"))\n",
         "  expect_equal(parse_field(out, \"VARCORR\"), \"1\")"))
cat("DONE\n")
