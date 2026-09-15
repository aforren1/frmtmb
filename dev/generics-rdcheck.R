# Verify the scale-contract page by RENDERING it, not by reading the
# source: `%` starts a comment in Rd even inside verbatim macros, and a
# rendered page is the only proof a table survived.
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
out <- file.path(root, "dev", "generics-out2", "scales-final.txt")
tools::Rd2txt(file.path(root, "man", "frmtmb-scales.Rd"), out = out)
x <- readLines(out)
rd <- readLines(file.path(root, "man", "frmtmb-scales.Rd"))
chk <- c(
  `conditional_effects() row, response` =
    any(grepl("conditional_effects()", x, fixed = TRUE) &
          grepl("response", x, fixed = TRUE)),
  `sigma() row says NA when it varies` =
    any(grepl("NA' if it varies", x, fixed = TRUE)),
  `distributional-sigma condition` =
    any(grepl("sigma ~ x", x, fixed = TRUE)),
  `truncation condition` = any(grepl("TRUNCATED", x, fixed = TRUE)),
  # round 2: the single figure was mislabeled and is gone; the page
  # states the closed form, which depends on where the bound sits
  `truncation closed form, no single figure` =
    any(grepl("pnorm(sigma", x, fixed = TRUE)) &&
      !any(grepl("0.1150", x, fixed = TRUE)),
  `hypothesis() row` = any(grepl("hypothesis()", x, fixed = TRUE)),
  # roxygen's own two header lines are Rd comments and belong there; a
  # % anywhere in the BODY would silently truncate the rendered line
  `no percent sign in the Rd body` =
    !any(grepl("%", rd[!grepl("^% ", rd)], fixed = TRUE)),
  `not keywords internal` =
    !any(grepl("\\keyword{internal}", rd, fixed = TRUE)))
cat(sprintf("rendered lines: %d\n", length(x)))
for (k in names(chk)) cat(sprintf("  %-40s %s\n", k, chk[[k]]))
cat(sprintf("all checks: %s\n", all(chk)))
