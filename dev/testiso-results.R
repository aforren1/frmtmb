# Write the "One-session results" section of dev/testiso-findings.md
# from the logs, to dev/testiso-log/results-section.md.
rs <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
counts <- function(tier) {
  system2(rs, c("dev/testiso-counts.R", tier), stdout = TRUE)
}
block <- function(lines, src) {
  c(paste0("<!-- BEGIN GENERATED: ", src, " -->"), "", lines, "",
    "<!-- END GENERATED -->")
}
cl <- function(dir) {
  out <- character(0)
  for (g in sort(list.files(dir, "[.]txt$", full.names = TRUE))) {
    y <- readLines(g, warn = FALSE)
    a <- grep("^test-compat", y, value = TRUE)
    r <- sub("registries unchanged: ", "",
             trimws(grep("^registries unchanged", y, value = TRUE)))
    out <- c(out, paste0("| ", sub("[.]txt$", "", basename(g)), " | ",
                         trimws(sub("^test-compat.R ", "", a[1])), " | ",
                         trimws(sub("^test-compat-register.R ", "", a[2])),
                         " | ", r, " |"))
  }
  c(paste("| preloaded | test-compat.R | test-compat-register.R |",
          "registries unchanged |"), "|---|---|---|---|", out)
}
x <- c(
  "## One-session results, before and after", "",
  "`NOT_CRAN=true`, no gate, every extension installed",
  "(`dev/release/run-sessions.sh`, summarized by `dev/testiso-counts.R`).",
  "\"own LEAK\" counts LEAK lines that no namespace load explains.", "",
  "### Before (rellib-r9, base tests from the main checkout)", "",
  paste("The forward and reverse sessions ran the first version of the",
        "runner. The preload sessions were added in fix round 1 and ran",
        "the final runner on the same base tests."), "",
  block(counts("testiso-before"), "dev/testiso-counts.R testiso-before"),
  "", "### After fix round 1 (this worktree, final tree)", "",
  block(counts("testiso-after2"), "dev/testiso-counts.R testiso-after2"),
  "", paste("The first after run (16 sessions, `testiso-after`, before",
            "fix round 1) also read fail=0 err=0 everywhere with no own",
            "LEAK line."), "",
  "### The two files with each extension preloaded", "",
  paste("`dev/testiso-compatload.sh` (one fresh session per preload set;",
        "logs `dev/testiso-log/compatload-before/` and",
        "`compatload-after5/`):"), "",
  "Before:", "", block(cl("dev/testiso-log/compatload-before"),
                       "dev/testiso-compatload.sh before"),
  "", "After:", "", block(cl("dev/testiso-log/compatload-after5"),
                          "dev/testiso-compatload.sh after5"))
writeLines(x, "dev/testiso-log/results-section.md")
