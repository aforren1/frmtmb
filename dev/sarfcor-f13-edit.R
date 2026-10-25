p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial/tests/testthat/test-guardrail.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  rds <- tempfile(fileext = \".rds\")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(list(car = fc, spde = fs), rds)",
"  # and a sar() residual term, on a row-standardized path graph
  P <- matrix(0, 10, 10)
  P[cbind(1:9, 2:10)] <- 1
  P <- P + t(P)
  Wp <- P / rowSums(P)
  dp <- data.frame(x = rnorm(10))
  dp$y <- as.vector(solve(diag(10) - 0.4 * Wp, dp$x + rnorm(10)))
  fsar <- frm(y ~ x + sar(Wp), data = dp, data2 = list(Wp = Wp))
  rds <- tempfile(fileext = \".rds\")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(list(car = fc, spde = fs, sar = fsar), rds)")
rep1("    expect_length(res, 22L)
    for (r in res) {
      nm <- strsplit(r, \" \", fixed = TRUE)[[1L]][2L]
      expect_match(r, paste0(\"TRUE The \", nm, \"() structure\"), fixed = TRUE)",
"    expect_length(res, 33L)
    for (r in res) {
      nm <- strsplit(r, \" \", fixed = TRUE)[[1L]][2L]
      # a field is a structure of the model, sar() a residual term of it
      expect_match(r, paste0(\"TRUE The \", nm, \"() \",
                             if (nm == \"sar\") \"term\" else \"structure\"),
                   fixed = TRUE)")
rep1("test_that(\"a saved fit read without frmtmb.spatial is refused by name\", {",
     "test_that(\"a saved fit read without frmtmb.spatial is refused by name\", {")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
