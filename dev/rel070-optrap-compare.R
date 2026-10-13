# Release 0.70.0: compare the merged build's reruns of lane optrap's
# studies (dev/rel070-log/studies/optrap/) with the lane's final logs
# (build 27, frmtmb-wt-optrap/dev/optrap-log/), row for row, because a
# summary line can agree while rows move.
#   Rscript dev/rel070-optrap-compare.R
lane <- "C:/Users/adf44/source/r/frmtmb-wt-optrap/dev"
mine <- "dev/rel070-log/studies/optrap"
ulps <- function(a, b) {
  # distance in units of the larger operand's last place
  d <- abs(a - b)
  e <- pmax(abs(a), abs(b)) * .Machine$double.eps
  ifelse(d == 0, 0, d / e)
}
cmp <- function(A, B, key, label) {
  ka <- do.call(paste, A[key]); kb <- do.call(paste, B[key])
  cat(sprintf("\n== %s: lane rows %d, merged rows %d, keys matched %d\n",
              label, nrow(A), nrow(B), sum(ka %in% kb)))
  B <- B[match(ka, kb), , drop = FALSE]
  same_row <- rep(TRUE, nrow(A))
  for (v in setdiff(intersect(names(A), names(B)), key)) {
    a <- A[[v]]; b <- B[[v]]
    if (is.numeric(a) && is.numeric(b)) {
      eq <- (a == b) | (is.na(a) & is.na(b))
      eq[is.na(eq)] <- FALSE
      d <- abs(a - b)
      cat(sprintf(paste0("%-12s identical %3d of %3d | max |diff| %.3g |",
                         " max ulps %.3g\n"), v, sum(eq), length(eq),
                  max(d, na.rm = TRUE), max(ulps(a, b), na.rm = TRUE)))
    } else {
      eq <- as.character(a) == as.character(b)
      eq[is.na(eq)] <- is.na(a)[is.na(eq)] & is.na(b)[is.na(eq)]
      cat(sprintf("%-12s identical %3d of %3d (text)\n", v, sum(eq),
                  length(eq)))
    }
    if (!v %in% c("sec", "secs")) same_row <- same_row & eq
  }
  cat(sprintf("rows identical in every non-timing column: %d of %d\n",
              sum(same_row), length(same_row)))
  invisible(list(A = A, B = B, same = same_row))
}

# 1. The scale study (dev/optrap-scale.R, lane tag p2, seeds 1:40)
A <- read.delim(file.path(lane, "optrap-log/scale-p2.tsv"))
B <- read.delim(file.path(mine, "optrap-log/scale-r8.tsv"))
r <- cmp(A, B, c("design", "seed", "scale"), "scale study")
if (any(!r$same)) {
  cat("rows that differ (lane | merged dll4, dll1, code):\n")
  i <- which(!r$same)
  print(data.frame(r$A[i, c("design", "seed", "scale")],
                   dll4_lane = r$A$dll4[i], dll4_r8 = r$B$dll4[i],
                   dll1_lane = r$A$dll1[i], dll1_r8 = r$B$dll1[i],
                   code_lane = r$A$code[i], code_r8 = r$B$code[i]),
        digits = 4, row.names = FALSE)
}

# 2. The mo() study (dev/optrap-mo-study.R, seeds 1:400)
rd <- function(f) do.call(rbind, lapply(f, read.delim))
A <- rd(Sys.glob(file.path(lane, "optrap-log/p2/mo/p2-*.tsv")))
B <- rd(Sys.glob(file.path(mine, "mo/r8-*.tsv")))
A <- A[order(A$seed), ]; B <- B[order(B$seed), ]
for (rng in list(1:200, 201:400)) for (arm in c("p2", "r8")) {
  X <- if (arm == "p2") A else B
  X <- X[X$seed %in% rng, ]
  cat(sprintf(paste0("%-6s seeds %d-%d n %d | at the exact max (gap < 1e-6)",
                     " %d | max gap %.4g | gap > 1e-6: %s | code != 0 %d ",
                     "| warned %d | evals %d\n"),
              arm, min(rng), max(rng), nrow(X), sum(X$gap < 1e-6),
              max(X$gap), paste(X$seed[X$gap >= 1e-6], collapse = ","),
              sum(X$code != 0), sum(X$nwarn > 0), sum(X$nfev)))
}
r <- cmp(A, B, "seed", "mo() study")
i <- which(r$A$nfev != r$B$nfev)
cat(sprintf("seeds whose evaluation count differs: %d\n", length(i)))
if (length(i)) {
  print(data.frame(seed = r$A$seed[i], nfev_lane = r$A$nfev[i],
                   nfev_r8 = r$B$nfev[i], gap_lane = r$A$gap[i],
                   gap_r8 = r$B$gap[i],
                   dll = r$B$ll_fit[i] - r$A$ll_fit[i]),
        digits = 4, row.names = FALSE)
}
cat(sprintf("sec total: lane %.1f | merged %.1f\n", sum(A$sec), sum(B$sec)))

# 3. The wide fit (dev/optrap-rev3-b1.R wide, seeds 1:3)
A <- read.delim(file.path(lane, "optrap-rev3-log/b1-wide-lane.tsv"))
B <- read.delim(file.path(mine, "rev3-log/b1-wide-r8.tsv"))
r <- cmp(A, B, c("case", "seed"), "wide fit")
print(data.frame(seed = A$seed, dll_lane = A$dll, dll_r8 = r$B$dll,
                 code_r8 = r$B$code, nwarn_r8 = r$B$nwarn,
                 unit_r8 = r$B$unit, secs_lane = A$secs,
                 secs_r8 = r$B$secs), digits = 6, row.names = FALSE)
