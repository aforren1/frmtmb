# Part B4: what it would take to put rr()'s factor-space expansion,
# gp()/hsgp() as specials, and the Student-t sharing of us_t/diag_t on
# the published covariance-structure entry shape. Counts, from the
# sources at the lane's base commit (18b536fc), of:
#   - the lines of R/ that special-case each feature by name or by its
#     own helpers (the code a port would move onto entry slots),
#   - the files and the top-level functions those lines sit in (the
#     consumers a port would touch),
#   - the commits of the last five rounds that changed such a line
#     (git log -G over 4f5ea39f..18b536fc, 0.68.1 to 0.74.0).
# Run: Rscript dev/sarfcor-b4-measure.R > dev/sarfcor-out/b4-measure.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
base <- "18b536fc"
files <- system2("git", c("-C", wt, "ls-tree", "--name-only", base, "R/"),
                 stdout = TRUE)
src <- lapply(files, function(f) {
  system2("git", c("-C", wt, "show", paste0(base, ":", f)), stdout = TRUE)
})
names(src) <- files
features <- list(
  rr = paste0("\"rr\"|rr_jacobians|rr_loadings|rr_npar|rr_start|",
              "[[]\"rank\"[]]|expand_b_rr|rr_expand"),
  gp_special = paste0("parse_gp_call|gpterms|gp_by|gp_pred|hsgp_basis|",
                      "hsgp_|gp_krige|gp_new|gp_corr|\"gp\"|\"hsgp\""),
  student_share = paste0("student_share|us_t|diag_t|t_mix|mixing scale|",
                         "dist = \"student\"|dist_cs")
)
top_fun <- function(lines, i) {
  # the nearest top-level definition at or above line i
  j <- i
  while (j > 0 && !grepl("^[A-Za-z_.][A-Za-z0-9_.]* <- ", lines[j])) {
    j <- j - 1L
  }
  if (j == 0) "(file top)" else sub(" <- .*$", "", lines[j])
}
for (nm in names(features)) {
  pat <- features[[nm]]
  hits <- list()
  for (f in files) {
    l <- src[[f]]
    code <- !grepl("^\\s*#", l)
    i <- which(code & grepl(pat, l, perl = TRUE))
    if (!length(i)) next
    hits[[f]] <- data.frame(file = f, line = i,
                            fun = vapply(i, function(k) top_fun(l, k), ""))
  }
  h <- do.call(rbind, hits)
  # compat rules are notes, not consumers
  h_code <- h[h$file != "R/compat.R", ]
  cat("====", nm, "\n")
  cat(sprintf("code lines naming it: %d in %d files (compat.R notes: %d lines)\n",
              nrow(h_code), length(unique(h_code$file)),
              sum(h$file == "R/compat.R")))
  funs <- unique(paste(h_code$file, h_code$fun, sep = "::"))
  cat(sprintf("top-level functions holding them: %d\n", length(funs)))
  tab <- sort(table(h_code$file), decreasing = TRUE)
  cat("by file:", paste(sprintf("%s %d", names(tab), tab), collapse = ", "),
      "\n")
  log <- system2("git", c("-C", wt, "log", "--oneline", "-G",
                          shQuote(pat), "4f5ea39f..18b536fc", "--", "R/"),
                 stdout = TRUE)
  cat(sprintf("commits 4f5ea39f..18b536fc changing such a line: %d\n",
              length(log)))
  if (length(log)) cat(paste0("  ", log), sep = "\n")
}
# the extensions that read these features
ext <- system2("git", c("-C", wt, "grep", "-n", "-E",
                        shQuote("\"rr\"|rr_jacobians|us_t|diag_t|student_share|\"gp\"|\"hsgp\"|gpterms"),
                        base, "--", "extensions/*/R/*.R"), stdout = TRUE)
cat("==== extension R/ lines naming any of them:", length(ext), "\n")
if (length(ext)) cat(paste0("  ", ext), sep = "\n")
