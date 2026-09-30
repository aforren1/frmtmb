# Lane ceplot: one-off rewrite of the "New group levels" paragraphs of
# ?conditional_effects for the crossed and mm(by = ) draws.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/conditional-effects.R"
s <- readLines(f)
a <- grep("^#'   With `band = \"boot\"`, a row with one member at an$", s)
b <- grep("^#'   refused by name; the Wald band answers it[.]$", s)
stopifnot(length(a) == 1, length(b) == 1)
s <- c(s[seq_len(a - 1)],
  "#'   With `band = \"boot\"`, a row with one member at an observed level",
  "#'   and one at a new level of the same block is refused by name, as",
  "#'   are rows that mix the two.",
  "#' - `mm(g1, g2, by = cbind(f1, f2))` is one block per by-level. A",
  "#'   member at an observed group reads the block that holds the group,",
  "#'   and a new member reads the block of its own by-value, where its",
  "#'   new level is drawn with that by-level's covariance. Two new",
  "#'   members with the same value and the same by-value read one new",
  "#'   level, and their weights add.",
  s[(b + 1):length(s)])
a <- grep("^#' placeholder exists, the call is refused by$", s)
b <- grep("^#' brms answers it[.] `band = \"wald\"` answers it too[.]$", s)
stopifnot(length(a) == 1, length(b) == 1)
s <- c(s[seq_len(a - 1)],
  "#' placeholder moves the rows without moving such a column, the rows",
  "#' keep their columns and the term's level is renamed instead: in a",
  "#' copy of the fit used for those rows alone, one level of the term",
  "#' takes the label the rows carry, and the draw replaces its effects.",
  "#' That is how crossed terms such as `(1 | g) + (1 | h) + (1 | g:h)`",
  "#' are drawn when `g` and `h` are set to observed levels that the data",
  "#' never has together: `g` and `h` keep their levels, and `g:h` reads a",
  "#' new level. On draws the band equals brms's at the same draws and",
  "#' the same seed, with `sample_new_levels = \"gaussian\"`",
  "#' (`dev/ceplot-crossed-brms.R`).",
  s[(b + 1):length(s)])
writeLines(s, f)
cat("patched\n")
