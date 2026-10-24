# Port, onto a reconstructed stage tree, the fixes to parts A and B made
# after that stage's snapshot (dev/spatial-out/*.diff), so that the
# stage patches each stand alone at their final state. Run with the
# stage tree as the working directory:
#   Rscript <worktree>/dev/spatial-stage-fix.R <stage: AB|Ccar>
stage <- commandArgs(trailingOnly = TRUE)[1]
sub_file <- function(p, old, new, fixed = TRUE, must = TRUE) {
  if (!file.exists(p)) { if (must) stop("missing ", p); return(invisible()) }
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  hit <- grepl(old, s, fixed = fixed)
  if (!hit && must) stop(p, ": no match: ", substr(old, 1, 60))
  writeLines(gsub(old, new, s, fixed = fixed), p)
}
# the build slot's resolver is named lookup, not data: codetools reads a
# call to a function named data() as utils::data() (R CMD check NOTE)
for (p in c("R/covstruct.R", "tests/testthat/test-covstruct-seam.R",
            "extensions/frmtmb.spatial/R/car.R")) {
  if (!file.exists(p)) next
  sub_file(p, "function(spec, gr, data)", "function(spec, gr, lookup)",
           must = FALSE)
  sub_file(p, "car_term_build(spec, gr, data)",
           "car_term_build(spec, gr, lookup)", must = FALSE)
  sub_file(p, "spde_term_build(spec, gr, data)",
           "spde_term_build(spec, gr, lookup)", must = FALSE)
  sub_file(p, "data(spec[[\"M_expr\"]], \"car()\")",
           "lookup(spec[[\"M_expr\"]], \"car()\")", must = FALSE)
  sub_file(p, "data(spec[[\"fem_expr\"]], \"spde()\")",
           "lookup(spec[[\"fem_expr\"]], \"spde()\")", must = FALSE)
  sub_file(p, "A <- data(spec[[\"A\"]], \"toyfld()\")",
           "A <- lookup(spec[[\"A\"]], \"toyfld()\")", must = FALSE)
}
sub_file("R/covstruct-seam.R", "`function(spec, gr, data)`",
         "`function(spec, gr, lookup)`", must = FALSE)
sub_file("R/covstruct-seam.R", "`data(expr, what)` resolves",
         "`lookup(expr, what)` resolves")
sub_file("R/covstruct-seam.R", "\"spec, gr, data\"", "\"spec, gr, lookup\"")
sub_file("R/frame.R", "#' and `data(expr, what)`, which resolves",
         "#' and `lookup(expr, what)`, which resolves")
sub_file("tests/testthat/test-guardrail1.R",
         "  expect_identical(guard_classify(\"Addition term `zz()` is not supported (supported: weights())\",\n                                  \"zz\"), \"generic\")",
         "  expect_identical(guard_classify(\n    \"Addition term `zz()` is not supported (supported: weights())\", \"zz\"),\n    \"generic\")")
news <- "extensions/frmtmb.latent/NEWS.md"
x <- readLines(news, warn = FALSE)
if (!any(grepl("development version", x))) {
  writeLines(c("# frmtmb.latent (development version)", "",
               "* `hmm()` reads core's one list of field terms (`car()` and `spde()`,",
               "  now parsed through core's covariance-structure registry) when it",
               "  decides whether a distributional parameter is constant, so it needs",
               "  the frmtmb release that carries `frmtmb_register_covstruct()`.",
               "", x), news)
}
cat("fixed", stage, "\n")
