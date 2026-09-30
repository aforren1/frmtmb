# Lane ceplot punch 1 (m7): ce_plot_ignored() and the doc lines.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/"
rep_all <- function(f, old, new, n_expect) {
  s <- paste(readLines(f), collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != n_expect) stop(basename(f), ": expected ", n_expect, " got ", n,
                          ": ", substr(old, 1, 60))
  s <- gsub(old, new, s, fixed = TRUE)
  writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
}
fn <- readLines(paste0(wt, "../dev/ceplot-p1-ignored.txt"))
rep_all(paste0(wt, "conditional-effects.R"),
        "#' brms's `theme` is a ggplot2 theme. Base graphics has no use for one,",
        paste(c(fn, "#' brms's `theme` is a ggplot2 theme. Base graphics has no use for one,"),
              collapse = "\n"), 1L)
for (f in paste0(wt, c("conditional-effects.R", "confint.R"))) {
  rep_all(f, "#'   with brms's warning. Base graphical parameters are accepted and\n#'   ignored, as by other `plot()` methods. Any other argument is an\n#'   error that names it.",
          "#'   with brms's warning. Base graphical parameters (`col`, `main`,\n#'   ...) are accepted, as every `plot()` method must accept them, and\n#'   ignored with a warning that names them. Any other argument is an\n#'   error that names it.", 1L)
}
rep_all(paste0(wt, "conditional-effects.R"),
        "#' @param ask If `TRUE` (the default), prompt before each new page after\n#'   the first on an interactive device.",
        "#' @param ask If `TRUE` (the default), prompt before each new page after\n#'   the first on an interactive device. `NULL`, the default before\n#'   this version, is taken as `TRUE`.", 1L)
rep_all(paste0(wt, "confint.R"),
        "#' @param ask If `TRUE` (the default), prompt before each new page\n#'   after the first on an interactive device.",
        "#' @param ask If `TRUE` (the default), prompt before each new page\n#'   after the first on an interactive device. `NULL`, the default\n#'   before this version, is taken as `TRUE`.", 1L)
cat("edited\n")
