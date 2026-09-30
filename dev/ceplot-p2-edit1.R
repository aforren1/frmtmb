# Lane ceplot punch 2 (P1-B1): predict() draws the "old_levels" choice
# once per call, over every response.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/predict-brms.R"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(s, old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("expected one match, got ", n, ": ", substr(old, 1, 50))
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, "                                   sample_new_levels = \"gaussian\") {",
"                                   sample_new_levels = \"gaussian\",
                                   pick = list()) {")
s <- rep1(s, "  if (identical(sample_new_levels, \"old_levels\")) {
    pick <- list()
    for (e in out) pick <- old_level_pick_add(pick, object, e$ed)",
"  if (identical(sample_new_levels, \"old_levels\")) {
    # `pick` carries the choices another response of the same call made,
    # so a grouping factor two responses share is chosen once
    for (e in out) pick <- old_level_pick_add(pick, object, e$ed)")
s <- rep1(s, "    nls[[nm]] <- predict_new_level_spec(object, rspec, newdata, re_formula,
                                        allow_new_levels, sample_new_levels)",
"    nls[[nm]] <- predict_new_level_spec(object, rspec, newdata, re_formula,
                                        allow_new_levels, sample_new_levels,
                                        pick = pick_all)
    pick_all <- attr(nls[[nm]], \"old_pick\") %||% pick_all")
s <- rep1(s, "  av <- list()
  nls <- list()",
"  av <- list()
  nls <- list()
  # brms's \"old_levels\" choice is one per call and grouping factor,
  # whichever response a term belongs to
  pick_all <- list()")
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
