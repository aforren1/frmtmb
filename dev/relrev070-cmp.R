# Reviewer of the 0.70.0 consolidation: the cox() and mo() studies rerun
# on rellib-r8 against the records they must reproduce.
#   Rscript dev/relrev070-cmp.R
R <- "C:/Users/adf44/source/r"
L <- file.path(R, "frmtmb-wt-release/dev/relrev070-log")
cmp_cox <- function(a, b, tag) {
  stopifnot(identical(a$seed, b$seed))
  cols <- setdiff(names(a), character(0))
  same <- vapply(cols, function(k) identical(a[[k]], b[[k]]), NA)
  cat(tag, ": seeds", nrow(a), "| columns identical:",
      paste(cols[same], collapse = ","), "| differ:",
      paste(cols[!same], collapse = ","), "\n")
  d <- a$ll - b$ll
  cat("   logLik identical", sum(a$ll == b$ll), "| max |d|",
      format(max(abs(d)), digits = 3), "| evals", sum(a$evals), sum(b$evals),
      "| gap > 1e-6", sum(a$gap > 1e-6), sum(b$gap > 1e-6), "| max gap",
      format(max(a$gap), digits = 4), format(max(b$gap), digits = 4), "\n")
}
r8 <- readRDS(file.path(L, "cox-study-r8.rds"))
r8i <- readRDS(file.path(L, "cox-study-r8-instr.rds"))
lane <- readRDS(file.path(R, "frmtmb-wt-moparity/dev/moparity-log/cox-study-lane.rds"))
first <- readRDS(file.path(R, "frmtmb-wt-release/dev/rel070-log/studies/moparity/dev/moparity-log/cox-study-lane.rds"))
cmp_cox(r8, lane, "cox r8 vs lane moparity (its rds)")
cmp_cox(r8, first, "cox r8 vs the first release build")
cmp_cox(r8, r8i, "cox r8 vs r8 instrumented")
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE)
mine <- do.call(rbind, lapply(c("1-100", "101-200", "201-300", "301-400"),
                              function(r) rd(file.path(L, "mo", paste0("r8-", r, ".tsv")))))
mi <- do.call(rbind, lapply(c("1-200", "201-400"), function(r)
  rd(file.path(L, "mo", paste0("r8i-", r, ".tsv")))))
c3 <- do.call(rbind, lapply(c("1-100", "101-200", "201-300", "301-400"),
                            function(r) rd(file.path(R, "frmtmb-wt-release/dev/rel070-log/studies/optrap/c3", paste0("mo-", r, ".tsv")))))
ofiles <- list.files(file.path(R, "frmtmb-wt-optrap/dev/optrap-log/mo-final"),
                     "[.]tsv$", full.names = TRUE)
ol <- do.call(rbind, lapply(ofiles, rd))
ol <- ol[order(ol$seed), ]
cmp_mo <- function(a, b, tag) {
  b <- b[match(a$seed, b$seed), ]
  ok <- !is.na(b$seed)
  a <- a[ok, ]; b <- b[ok, ]
  cols <- setdiff(names(a), "sec")
  rowsame <- Reduce(`&`, lapply(cols, function(k) {
    x <- a[[k]]; y <- b[[k]]
    (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  }))
  cat(tag, ": rows", nrow(a), "| identical in every non-timing column",
      sum(rowsame), "\n")
  if (any(!rowsame)) {
    bad <- which(!rowsame)
    for (k in cols) {
      dif <- !((is.na(a[[k]]) & is.na(b[[k]])) | (a[[k]] == b[[k]]) %in% TRUE)
      if (any(dif)) cat("   column", k, "differs on", sum(dif), "rows\n")
    }
  }
}
cat("mo r8: rows", nrow(mine), "| gap > 1e-6", sum(mine$gap > 1e-6),
    "| max gap", format(max(mine$gap), digits = 3), "| codes",
    paste(names(table(mine$code)), table(mine$code), sep = "="),
    "| warned", sum(mine$nwarn > 0), "| NaN SE", sum(mine$nse_nan > 0), "\n")
cmp_mo(mine, c3, "mo r8 vs the integration's c3 rerun")
cmp_mo(mine, mi, "mo r8 vs r8 instrumented")
cmp_mo(mine, ol, paste0("mo r8 vs lane optrap mo-final (", nrow(ol), " rows)"))
p2 <- do.call(rbind, lapply(Sys.glob(file.path(R, "frmtmb-wt-optrap/dev/optrap-log/p2/mo/p2-*.tsv")), rd))
cmp_mo(mine, p2, paste0("mo r8 vs lane optrap p2 final (", nrow(p2), " rows)"))
