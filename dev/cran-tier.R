# Lane cran: build a proposed CRAN test tier from the per-block timings
# of dev/cran-runall.sh (NOT_CRAN unset) and map every export of core to
# the tier blocks and the examples that call it.
#
#   Rscript dev/cran-tier.R <worktree> <threshold s> > dev/cran-log/tier.txt
#
# Rule: keep every test_that() block that ran (not skipped on CRAN) in
# under <threshold> seconds elapsed; then, for each export that no kept
# block, no kept file's top-level code and no example names, add the
# cheapest moved block that names it (greedy). A name counts when it is
# one of all.names() of the block's parsed code, so a call, a bare
# reference (lapply(x, fitted)) and a method name all count; a generic
# such as predict() counts on any class, which over-reports.
local({
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
root <- a[1]
X <- as.numeric(a[2])
logd <- file.path(root, "dev", "cran-log")
helper_cost <- 0.76  # elapsed of a file that skips at once (files-sorted)

b <- utils::read.delim(file.path(logd, "blocks.tsv"), header = FALSE,
                       quote = "", col.names = c("tag", "file", "desc", "nb",
                                                 "skipped", "cpu", "real"))
b$skipped <- as.logical(b$skipped)

# every test_that() call in a file, in source order, with its code
find_tt <- function(e, out) {
  if (is.call(e)) {
    fn <- e[[1L]]
    if (identical(fn, quote(test_that)) ||
        identical(fn, quote(testthat::test_that))) {
      out[[length(out) + 1L]] <- e
      return(out)
    }
    for (k in seq_along(e)[-1L]) {
      if (is.call(e[[k]])) out <- find_tt(e[[k]], out)
    }
  }
  out
}
tdir <- file.path(root, "tests", "testthat")
files <- sort(list.files(tdir, "^test-.*[.]R$"))
parsed <- list()
toplevel <- list()
for (f in files) {
  ex <- parse(file.path(tdir, f), keep.source = FALSE)
  tt <- list()
  top <- character()
  for (e in as.list(ex)) {
    got <- find_tt(e, list())
    if (length(got)) tt <- c(tt, got) else top <- c(top, all.names(e))
  }
  parsed[[f]] <- tt
  toplevel[[f]] <- unique(top)
}

# names each block's code uses; match a block to its call by literal
# description, else to the file's only non-literal call, else the file
desc_of <- function(cl) {
  d <- cl[[2L]]
  if (is.character(d)) d else NA_character_
}
b$names <- vector("list", nrow(b))
b$how <- NA_character_
for (i in seq_len(nrow(b))) {
  tt <- parsed[[b$file[i]]]
  ds <- vapply(tt, desc_of, "")
  j <- which(ds == b$desc[i])
  if (length(j)) {
    b$names[[i]] <- unique(all.names(tt[[j[1L]]]))
    b$how[i] <- "literal"
  } else if (sum(is.na(ds)) == 1L) {
    b$names[[i]] <- unique(all.names(tt[[which(is.na(ds))]]))
    b$how[i] <- "computed"
  } else {
    b$names[[i]] <- unique(unlist(lapply(tt, all.names)))
    b$how[i] <- "file"
  }
}
cat("== block matching\n")
print(table(b$how))

nsf <- parseNamespaceFile("frmtmb", dirname(root))
exports <- sort(unique(nsf$exports))
cat("\n== exports in NAMESPACE export():", length(exports), "\n")

# names the examples use (examples run on CRAN as well)
rd <- tools::Rd_db(dir = root)
ex_names <- unique(unlist(lapply(rd, function(r) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f))
  tools::Rd2ex(r, f, commentDontrun = TRUE, commentDonttest = FALSE)
  if (!file.exists(f)) return(character())
  tryCatch(all.names(parse(f, keep.source = FALSE)),
           error = function(e) character())
})))

ran <- !b$skipped
keep <- ran & b$real < X
cat("\n== threshold", X, "s\n")
cat("blocks run on CRAN:", sum(ran), " kept:", sum(keep),
    " moved:", sum(ran & !keep), "\n")

kept_files <- unique(b$file[keep])
top_names <- unique(unlist(toplevel[kept_files]))
covered_by <- function(k) {
  vapply(exports, function(e) {
    any(vapply(b$names[k], function(n) e %in% n, TRUE))
  }, TRUE)
}
cov_all <- covered_by(which(ran))
cov_any <- covered_by(seq_len(nrow(b)))
cov_keep <- covered_by(which(keep)) | exports %in% top_names
cov_ex <- exports %in% ex_names
cat("exports named by any block (incl. skipped on CRAN):", sum(cov_any), "\n")
cat("exports named by a block that runs on CRAN:", sum(cov_all), "\n")
cat("exports named by the threshold tier (blocks + top level):",
    sum(cov_keep), "\n")
cat("exports named by an example:", sum(cov_ex), "\n")

# greedy: cheapest moved block for each export the tier and examples miss
add <- integer()
need <- exports[!cov_keep & !cov_ex & cov_all]
cand <- which(ran & !keep)
while (length(need)) {
  gain <- vapply(cand, function(k) sum(need %in% b$names[[k]]), 0)
  if (!any(gain > 0)) break
  score <- b$real[cand] / pmax(gain, 1e-9)
  score[gain == 0] <- Inf
  k <- cand[which.min(score)]
  add <- c(add, k)
  need <- setdiff(need, b$names[[k]])
  cand <- setdiff(cand, k)
}
keep2 <- keep
keep2[add] <- TRUE
cov_final <- covered_by(which(keep2)) |
  exports %in% unique(unlist(toplevel[unique(b$file[keep2])]))
cat("\n== greedy additions for exports no tier block or example names:",
    length(add), "blocks,", sprintf("%.1f", sum(b$real[add])), "s\n")
for (k in add) {
  cat("ADD\t", b$file[k], "\t", b$desc[k], "\t", sprintf("%.2f", b$real[k]),
      "\t", paste(intersect(exports, b$names[[k]])[
        !cov_keep[intersect(exports, b$names[[k]])]], collapse = ","),
      "\n", sep = "")
}
nf <- length(unique(b$file[keep2]))
cat("\n== proposed tier: blocks", sum(keep2), " files", nf,
    " expectations", sum(b$nb[keep2]),
    " block time", sprintf("%.1f", sum(b$real[keep2])), "s\n")
cat("moved behind NOT_CRAN: blocks", sum(ran & !keep2), " expectations",
    sum(b$nb[ran & !keep2]), " block time",
    sprintf("%.1f", sum(b$real[ran & !keep2])), "s, in",
    length(unique(b$file[ran & !keep2])), "files\n")
whole <- setdiff(unique(b$file[ran]), unique(b$file[keep2]))
cat("files moved whole:", length(whole), "\n")

cat("\n== exports with no tier block and no example\n")
miss <- exports[!cov_final & !cov_ex]
for (e in miss) {
  why <- if (!cov_any[exports == e]) "no test names it at all" else
    if (!cov_all[exports == e]) "only blocks already skipped on CRAN" else
      "?"
  cat("MISS\t", e, "\t", why, "\n", sep = "")
}
cat("\n== exports reached only by examples (no tier block)\n")
cat(exports[!cov_final & cov_ex], sep = ", ")
cat("\n")

cat("\n== per-export tier blocks (first three files)\n")
for (e in exports) {
  k <- which(keep2 & vapply(b$names, function(n) e %in% n, TRUE))
  fl <- unique(b$file[k])
  cat("MAP\t", e, "\tblocks=", length(k), "\tfiles=", length(fl), "\t",
      paste(utils::head(sub("[.]R$", "", fl), 3), collapse = " "),
      if (e %in% ex_names) "\t+example" else "", "\n", sep = "")
}

# the moved blocks, slowest first, for the findings
mv <- b[ran & !keep2, c("file", "desc", "nb", "real")]
mv <- mv[order(-mv$real), ]
utils::write.table(mv, file.path(logd, "tier-moved.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
kp <- b[keep2, c("file", "desc", "nb", "real")]
utils::write.table(kp, file.path(logd, "tier-kept.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
pf <- do.call(rbind, lapply(split(seq_len(nrow(b)), b$file), function(ix) {
  data.frame(file = b$file[ix[1L]], run = sum(ran[ix]),
             kept = sum(keep2[ix]), moved = sum(ran[ix] & !keep2[ix]),
             kept_s = sum(b$real[ix][keep2[ix]]),
             moved_s = sum(b$real[ix][ran[ix] & !keep2[ix]]))
}))
pf <- pf[order(-pf$moved_s), ]
utils::write.table(pf, file.path(logd, "tier-perfile.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
})
