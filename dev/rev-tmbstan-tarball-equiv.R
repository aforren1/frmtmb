# Reviewer check: does the checked tarball behave like the tree?
#
# The check that produced PASS 1036 ran on a tarball built before a
# later edit. Only two files differ from the tree: the guard test file
# and DESCRIPTION. This asks what the difference is, rather than taking
# the claim that it is comment-only.
tar <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
              "c--Users-adf44-source-r-frmtmb/",
              "529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/",
              "rev-tar2/frmtmb.sample")
tree <- paste0("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/",
               "extensions/frmtmb.sample")
rel <- "tests/testthat/test-tmbstan-build-guard.R"

a <- parse(file.path(tar, rel), keep.source = FALSE)
b <- parse(file.path(tree, rel), keep.source = FALSE)
cat("expressions: tarball", length(a), " tree", length(b), "\n")
cat("identical(parse, parse):", identical(a, b), "\n")
cat("identical(deparse, deparse):",
    identical(deparse(a), deparse(b)), "\n")

# what actually changed in the bytes
la <- readLines(file.path(tar, rel), warn = FALSE)
lb <- readLines(file.path(tree, rel), warn = FALSE)
cat("lines: tarball", length(la), " tree", length(lb), "\n")
d <- setdiff(lb, la)
cat("lines present only in the tree:", length(d), "\n")
cat(paste0("  + ", d), sep = "\n")
d2 <- setdiff(la, lb)
cat("lines present only in the tarball:", length(d2), "\n")
if (length(d2)) cat(paste0("  - ", d2), sep = "\n")
cat("all tree-only lines are comments or blank:",
    all(grepl("^\\s*(#|$)", d)), "\n")

# and the non-ASCII scan R CMD check runs over R sources reads comments
cat("non-ASCII in the tree file:",
    any(grepl("[^\x01-\x7f]", lb, useBytes = TRUE)), "\n")

# DESCRIPTION: which fields differ
da <- read.dcf(file.path(tar, "DESCRIPTION"))
db <- read.dcf(file.path(tree, "DESCRIPTION"))
only_tar <- setdiff(colnames(da), colnames(db))
shared <- intersect(colnames(da), colnames(db))
changed <- shared[vapply(shared,
                         function(k) !identical(da[, k], db[, k]), NA)]
cat("DESCRIPTION fields only in the tarball:",
    paste(only_tar, collapse = ", "), "\n")
cat("DESCRIPTION shared fields that differ:",
    if (length(changed)) paste(changed, collapse = ", ") else "none",
    "\n")
