## Blast radius of lane brmsnames: every source line in the monorepo
## that READS one of the surfaces this lane renames or reshapes. Run on
## the unchanged tree first (the record) and again at the end (what is
## left pointing at an old spelling should be zero or explained).
##
##   Rscript dev/brmsnames-blast.R > dev/brmsnames-log/blast-<tag>.txt
##
## Scans R/, tests/, vignettes/ and the same three under every
## extensions/<pkg>/. Roxygen lines (#') count, because an example or a
## documented name that goes stale is part of the radius.
roots <- c(".", list.dirs("extensions", recursive = FALSE))
files <- unlist(lapply(roots, function(r) {
  unlist(lapply(c("R", "tests", "vignettes"), function(d) {
    p <- file.path(r, d)
    if (!dir.exists(p)) return(character())
    list.files(p, pattern = "[.](R|Rmd)$", recursive = TRUE,
               full.names = TRUE)
  }))
}))
files <- sub("^[.]/", "", files)

## One pattern per surface. Written as bracket classes, not escapes, so
## the file survives any shell that eats a backslash.
pats <- list(
  "VarCorr element read" =
    "VarCorr[(][^)]*[)] *([[][[]|[$])",
  "VarCorr names or length" =
    "(names|length|expect_named|expect_length)[(] *VarCorr",
  "as.data.frame(VarCorr())" = "as[.]data[.]frame[(] *VarCorr",
  "variables() call" = "variables[(]",
  "draws matrix by name" = "[$]draws *[[] *, *[\"'c]",
  "draws column names" = "colnames[(][^)]*[$]draws",
  "old draws name b[i]" = "\"b[[][0-9]",
  "old draws name theta_k" = "\"theta_[0-9]",
  "bare coefficient names in a draws read" =
    "(rhat|neff_ratio|posterior_summary|as_draws_[a-z]+|summary)[(][^)]*[)] *[[][^]]*\"(Intercept|x|sigma_Intercept)\"",
  "hypothesis() result column" =
    "[$](estimate|se|lwr|upr|evid_ratio|post_prob)[^a-zA-Z_]",
  "hypothesis() call" = "hypothesis[(]",
  "posterior_summary() call" = "posterior_summary[(]",
  "fixef/ranef/coef/VarCorr summary arg" =
    "(fixef|ranef|coef|VarCorr)[(][^)]*(summary|robust|probs) *="
)

tab <- list()
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  for (nm in names(pats)) {
    n <- sum(grepl(pats[[nm]], txt))
    if (n) tab[[length(tab) + 1L]] <- data.frame(surface = nm, file = f,
                                                  lines = n)
  }
}
tab <- do.call(rbind, tab)
pkg_of <- function(f) {
  if (startsWith(f, "extensions/")) strsplit(f, "/")[[1L]][2L] else "frmtmb"
}
tab$package <- vapply(tab$file, pkg_of, "")
tab$kind <- ifelse(grepl("/tests/", paste0("/", tab$file)), "tests",
            ifelse(grepl("/vignettes/", paste0("/", tab$file)), "vignettes",
                   "R"))

cat("files scanned:", length(files), "\n\n")
cat("== lines per surface, by kind ==\n")
agg <- stats::aggregate(lines ~ surface + kind, tab, sum)
agg <- agg[order(agg$surface, agg$kind), ]
print(agg, row.names = FALSE)
cat("\n== lines per surface, by package ==\n")
agg2 <- stats::aggregate(lines ~ surface + package, tab, sum)
agg2 <- agg2[order(agg2$surface, agg2$package), ]
print(agg2, row.names = FALSE)
cat("\n== files per surface ==\n")
for (nm in names(pats)) {
  s <- tab[tab$surface == nm, ]
  if (!nrow(s)) next
  cat("--", nm, ":", nrow(s), "files,", sum(s$lines), "lines\n")
  cat(paste0("   ", s$file, " (", s$lines, ")"), sep = "\n")
}
cat("DONE\n")
