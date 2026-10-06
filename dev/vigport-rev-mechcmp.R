# Reviewer: compare the reviewer's mechanical reruns (seeded with the
# lane's runner, unseeded with the reviewer's copy) with the lane's
# results, expression by expression, and re-derive "transform
# unchanged" against the tracked 0.34.0 baseline.
#
#   Rscript dev/vigport-rev-mechcmp.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
L <- readRDS(file.path(root, "vigport-port-out/r5/results-merged.rds"))
S <- readRDS(file.path(root, "vigport-rev-out/seeded/results-merged.rds"))
U <- readRDS(file.path(root, "vigport-rev-out/noseed/results-merged.rds"))
B <- readRDS(file.path(root, "brms-port/results-merged.rds"))
cmp <- function(a, b, tag) {
  k <- intersect(a$id, b$id)
  ca <- a$class[match(k, a$id)]; cb <- b$class[match(k, b$id)]
  sa <- a$status_raw[match(k, a$id)]; sb <- b$status_raw[match(k, b$id)]
  cat(sprintf(paste0("%-22s ids %d/%d common %d; class differs %d;",
                     " raw status differs %d\n"),
              tag, nrow(a), nrow(b), length(k), sum(ca != cb), sum(sa != sb)))
  for (i in which(ca != cb | sa != sb)) cat("   ", k[i], sa[i], ca[i], "|",
                                           sb[i], cb[i], "\n")
}
cmp(L, S, "lane vs rev seeded")
cmp(L, U, "lane vs rev unseeded")
cat("\n## transform unchanged\n")
for (nm in c("L", "S", "U")) {
  X <- get(nm)
  k <- intersect(X$id, B$id)
  same <- X$src[match(k, X$id)] == B$src[match(k, B$id)]
  cat(sprintf("%s: ids %d, baseline ids %d, common %d, same src %d\n", nm,
              nrow(X), nrow(B), length(k), sum(same)))
}
# Independently of any run: transform the vignettes now with the
# lane's port-lib.R and with the 0.34.0 port-lib.R from git, and compare.
HERE <- file.path(root, "brms-port")
e_new <- new.env(); sys.source(file.path(HERE, "port-lib.R"), e_new)
old <- system2("git", c("-C", root, "show",
                        "9e902909:dev/brms-port/port-lib.R"), stdout = TRUE)
tf <- tempfile(fileext = ".R"); writeLines(old, tf)
e_old <- new.env(); sys.source(tf, e_old)
VIGS <- unique(B$vignette)
n <- 0; same <- 0
for (v in VIGS) {
  cn <- e_new$extract_vignette(v); co <- e_old$extract_vignette(v)
  stopifnot(length(cn) == length(co))
  for (i in seq_along(cn)) {
    tn <- e_new$transform_code(cn[[i]]$code)
    to <- e_old$transform_code(co[[i]]$code)
    for (j in seq_along(tn)) {
      n <- n + 1
      same <- same + identical(tn[[j]]$src, to[[j]]$src)
      id <- sprintf("%s.%d.%d", v, cn[[i]]$idx, j)
      bsrc <- B$src[B$id == id]
      # summarize.R stores src with each newline and its indent as a space
      if (length(bsrc) && !identical(bsrc, gsub("\n *", " ", tn[[j]]$src)))
        cat("  differs from baseline:", id, "\n")
    }
  }
}
cat(sprintf("old vs new port-lib transform: %d of %d expressions identical\n",
            same, n))
