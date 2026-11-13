# Lane nlstop: the bigx grid (dev/nlstop-bigx-grid.R) on two builds, cell
# by cell. Usage: Rscript dev/nlstop-bigxcmp.R <before.tsv> <after.tsv>
f <- commandArgs(TRUE)
a <- read.delim(f[1], stringsAsFactors = FALSE)
b <- read.delim(f[2], stringsAsFactors = FALSE)
k <- c("family", "re", "m", "seed")
stopifnot(identical(a[k], b[k]))
cat("cells", nrow(a),
    "| dll differs by > 1e-8:", sum(abs(a$dll - b$dll) > 1e-8, na.rm = TRUE),
    "| slope_z differs by > 1e-6:",
    sum(abs(a$slope_z - b$slope_z) > 1e-6, na.rm = TRUE),
    "| warning count differs:", sum(a$w_warn != b$w_warn),
    "| code differs:", sum(a$w_code != b$w_code, na.rm = TRUE),
    "| max |dll after - dll before|:",
    format(max(abs(a$dll - b$dll), na.rm = TRUE), digits = 3), "\n")
