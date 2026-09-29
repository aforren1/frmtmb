# Reviewer, claim 7: render every changed Rd topic and grep the OUTPUT,
# not the source, because `%` starts a comment in Rd even inside
# \preformatted{}.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-resmooth"
files <- c(file.path(root, "man",
                     c("conditional_effects.Rd", "frm_bootstrap.Rd",
                       "frm_linpred.Rd", "frm_lp_basis.Rd", "pp_check.Rd",
                       "simulate.frmtmb_fit.Rd")),
           file.path(root, "extensions/frmtmb.spline/man",
                     c("frm_curve.Rd", "frm_curve_deriv.Rd",
                       "frm_curve_feature.Rd")))
dir.create("dev/resmooth-rev-rd", showWarnings = FALSE)
for (f in files) {
  out <- file.path("dev/resmooth-rev-rd",
                   sub("\\.Rd$", ".txt", basename(f)))
  r <- tryCatch({tools::Rd2txt(f, out = out); "OK"},
                error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  txt <- if (file.exists(out)) readLines(out, warn = FALSE) else character(0)
  cat(sprintf("%-26s %-10s lines=%d\n", basename(f), r, length(txt)))
  # a stray % would have eaten the rest of its line
  bad <- grep("%", txt, value = TRUE)
  if (length(bad)) {
    cat("   lines with a literal % in the RENDERED text:\n")
    cat(paste0("     ", bad, collapse = "\n"), "\n")
  }
}
cat("\n== phrases the new prose must carry, grepped on the RENDERED text ==\n")
want <- list(
  "frm_linpred.txt" = c("Every SMOOTH stays in", "as brms does",
                        "population curve", "allow_new_levels"),
  "pp_check.txt" = c("holds every SMOOTH at its estimate"),
  "simulate.frmtmb_fit.txt" = c("A SMOOTH is never redrawn"),
  "frm_bootstrap.txt" = c("holds EVERY smooth"),
  "conditional_effects.txt" = c("keeps every SMOOTH"),
  "frm_curve.txt" = c("KEEPS every"),
  "frm_lp_basis.txt" = c("drops the"))
for (nm in names(want)) {
  txt <- paste(readLines(file.path("dev/resmooth-rev-rd", nm),
                         warn = FALSE), collapse = " ")
  txt <- gsub("[[:space:]]+", " ", txt)
  for (w in want[[nm]]) {
    cat(sprintf("  %-26s %-40s %s\n", nm, w,
                if (grepl(w, txt, fixed = TRUE)) "present" else "ABSENT"))
  }
}
cat("\n== style scan on the rendered Rd and on NEWS.md ==\n")
scan_style <- function(path, label) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  em <- grep("\u2014|\u2013", txt)
  sp <- grep(" - ", txt)
  brit <- grep(paste0("behaviour|colour|centre[ds]?\\b|analyse|",
                     "normalis|organis|recognis|modell(ing|ed)"),
               txt, ignore.case = TRUE)
  long <- which(nchar(txt) > 80)
  cat(sprintf("%-30s em/en dash lines: %d | ' - ' lines: %d | ",
              label, length(em), length(sp)))
  cat(sprintf("British spellings: %d | >80 col: %d\n", length(brit),
              length(long)))
  if (length(em)) cat("   dash at:", paste(em, collapse = ","), "\n")
  if (length(brit)) {
    cat("   British at:", paste(brit, collapse = ","), "\n")
    cat(paste0("     ", txt[brit], collapse = "\n"), "\n")
  }
  if (length(long)) cat("   >80 col at:", paste(long, collapse = ","), "\n")
}
scan_style(file.path(root, "NEWS.md"), "NEWS.md")
for (f in list.files("dev/resmooth-rev-rd", full.names = TRUE)) {
  scan_style(f, basename(f))
}
