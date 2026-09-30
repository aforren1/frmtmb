# One-off source edit (punch round 2, P1-M3): the reuse key carries
# re_formula.
f <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/R/conditional-effects.R"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(from, to) {
  n <- lengths(regmatches(s, gregexpr(from, s, fixed = TRUE)))
  stopifnot(n == 1)
  s <<- sub(from, to, s, fixed = TRUE)
}
rep1("#' `serialize()` rather than a hash: no new dependency, and an exact",
"#' `re_formula` goes in too: on `(1 + x | g)` at an observed `g`, the\n#' calls under `NULL` and under `~ (1 | g)` read the same grid and hold\n#' the same term in the simulation, and predict different curves.\n#'\n#' `serialize()` rather than a hash: no new dependency, and an exact")
rep1("ce_boot_key <- function(grids, categorical, resp, dpar, lens,
                        plan = NULL, kept = integer(0)) {
  serialize(list(
    resp = resp, dpar = dpar, categorical = categorical, lens = lens,",
"ce_boot_key <- function(grids, categorical, resp, dpar, lens,
                        plan = NULL, kept = integer(0), re_form = NA) {
  serialize(list(
    resp = resp, dpar = dpar, categorical = categorical, lens = lens,
    re_formula = ce_re_key(re_form),")
rep1("  key <- ce_boot_key(grids, categorical, resp, dpar, lens, plan, kept)",
"  key <- ce_boot_key(grids, categorical, resp, dpar, lens, plan, kept,
                     re_form)")
rep1("               } else if (!identical(boot$ce_kept %||% integer(0), kept)) {",
"               } else if (!identical(boot$ce_re %||% NA, ce_re_key(re_form))) {
                 paste0(\"predictions under a different re_formula (\",
                        ce_re_label(boot$ce_re %||% NA), \" there, \",
                        ce_re_label(ce_re_key(re_form)), \" here), which \",
                        \"keep different group-level terms\")
               } else if (!identical(boot$ce_kept %||% integer(0), kept)) {")
rep1("  bs$ce_kept <- kept\n",
"  bs$ce_kept <- kept\n  bs$ce_re <- ce_re_key(re_form)\n")
rep1("#' ONE parametric bootstrap for every grid of the call.",
"#' `re_formula` as a comparable value: a formula's text, without the\n#' environment that would make two equal formulas differ.\n#'\n#' @noRd\nce_re_key <- function(re_form) {\n  if (inherits(re_form, \"formula\")) {\n    paste(deparse(re_form, width.cutoff = 500L), collapse = \" \")\n  } else if (is.null(re_form)) {\n    \"NULL\"\n  } else {\n    \"NA\"\n  }\n}\n\nce_re_label <- function(k) {\n  if (k %in% c(\"NULL\", \"NA\")) paste(\"re_formula =\", k) else k\n}\n\n#' ONE parametric bootstrap for every grid of the call.")
writeLines(s, f)
cat("ok\n")
