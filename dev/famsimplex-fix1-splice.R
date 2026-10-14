# Put the fix-round text and its generated blocks into section 9 of the
# findings, in place of its PENDING line.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
f <- file.path(wt, "famsimplex-findings.md")
s <- readLines(f)
i <- which(s == "## 9. Fix round 1")
stopifnot(length(i) == 1L, s[i + 2L] == "PENDING.")
txt <- readLines(file.path(wt, "famsimplex-fix1-text.md"))
blk <- readLines(file.path(wt, "famsimplex-summary-fix1.md"))
cut_at <- function(title) grep(paste0("BEGIN generated: ", title), blk,
                               fixed = TRUE)
b1 <- cut_at("constants on locations")
b2 <- cut_at("pre-fix variant")
b3 <- cut_at("scaled density error")
b4 <- cut_at("parser-reaching files")
seg <- function(a, b) blk[a:(b - 1L)]
new <- c(txt, seg(b1, b2),
  "Seen to fail before the fix: `dev/famsimplex-prefix.R off` switches",
  "the four guards off in the source, the lane library is rebuilt and",
  "`test-simplex-families.R` run, and `dev/famsimplex-prefix.R on`",
  "restores the files (checked byte for byte with `cmp`):",
  "", seg(b2, b3),
  "The scaled error, on the tests' own draws and on 400 wider cases",
  "(\"old\" is the former measure, over `max(1, |ll|)`):", "",
  seg(b3, b4),
  "Reruns. `fix1` ran every file that names `categorical()`,",
  "`mixture()`, `multinomial()`, `mvbf()` or a simplex family, and",
  "`test-parse.R`, in core and the extensions, gated. It ran on the",
  "build before the refusal's last rewording and the dharma test's",
  "split. `fix2` and the OpenBLAS run are on the final build:", "",
  blk[b4:length(blk)])
s <- c(s[seq_len(i + 1L)], new, s[(i + 3L):length(s)])
writeLines(s, f)
cat("spliced", length(new), "lines\n")
