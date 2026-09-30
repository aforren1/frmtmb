# Lane ceplot punch 1 (m7): plot() names the graphical parameters it
# ignores, hides do_plot from its argument list, and takes ask = NULL.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/"
rep_all <- function(f, old, new, n_expect) {
  s <- paste(readLines(f), collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != n_expect) stop(basename(f), ": expected ", n_expect, " got ", n,
                          ": ", substr(old, 1, 60))
  s <- gsub(old, new, s, fixed = TRUE)
  writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
}
for (f in paste0(wt, c("conditional-effects.R", "confint.R"))) {
  rep_all(f, "  frm_check_dots(..., .allow = \"do_plot\")\n",
          "  frm_check_dots(..., .hidden = \"do_plot\")\n  ce_plot_ignored(...)\n", 1L)
  rep_all(f, "  check_flag(ask, \"ask\")\n",
          "  # NULL, the default of 0.66.0, asks as brms's TRUE does\n  if (is.null(ask)) ask <- TRUE\n  check_flag(ask, \"ask\")\n", 1L)
}
rep_all(paste0(wt, "conditional-effects.R"),
        "  frm_check_dots(...)\n  ce_draw(x)\n",
        "  frm_check_dots(...)\n  ce_plot_ignored(...)\n  ce_draw(x)\n", 2L)
rep_all(paste0(wt, "confint.R"),
        "  frm_check_dots(...)\n  hyp_draw_page(x)\n",
        "  frm_check_dots(...)\n  ce_plot_ignored(...)\n  hyp_draw_page(x)\n", 2L)
cat("edited\n")
