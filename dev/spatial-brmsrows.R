# Ported brms rows (dev/brmsport-verdicts.tsv) that depend on car() or
# spde(): a row id is "<file stem>:<line>" in brms 2.23.0's test suite
# (dev/brms-suite/brms/tests/testthat/tests.<stem>.R). A row depends on
# car() when its line, or the test_that() block holding it, calls
# car( or cor_car(. spde() is not a brms special, so no brms row can
# name it; the count is checked rather than assumed.
for (vf in c("dev/brmsport-verdicts.tsv", "dev/brmsport-bin2-verdicts.tsv")) {
cat("LEDGER", vf, "
")
v <- utils::read.delim(vf, stringsAsFactors = FALSE)
dir <- "dev/brms-suite/brms/tests/testthat"
stem <- sub(":.*$", "", v$id)
line <- suppressWarnings(as.integer(sub("^.*:", "", v$id)))
hit <- function(pat) {
  out <- logical(nrow(v))
  for (s in unique(stem)) {
    f <- file.path(dir, paste0("tests.", s, ".R"))
    if (!file.exists(f)) next
    x <- readLines(f, warn = FALSE)
    ex <- parse(text = x, keep.source = TRUE)
    sr <- attr(ex, "srcref")
    for (i in seq_along(ex)) {
      a <- sr[[i]][1]; b <- sr[[i]][3]
      if (!any(grepl(pat, x[a:b]))) next
      k <- which(stem == s & !is.na(line) & line >= a & line <= b)
      out[k] <- TRUE
    }
  }
  out
}
for (w in c(car = "(^|[^a-zA-Z_.])(cor_)?car\\(", spde = "spde\\(")) {
  h <- hit(w)
  cat("ROWS depending on", w, ":", sum(h), "\n")
  if (any(h)) print(v[h, c("id", "pkg", "verdict")], row.names = FALSE)
}
}
