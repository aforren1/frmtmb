# Summarize grid arms, generated into dev/impsamp-out/grid-summary-*.txt:
#   Rscript dev/impsamp-grid-summary.R base after
args <- commandArgs(TRUE)
rd <- function(arm) {
  fs <- Sys.glob(sprintf("dev/impsamp-out/grid-%s-*.tsv", arm))
  d <- do.call(rbind, lapply(fs, utils::read.delim, stringsAsFactors = FALSE))
  d$arm <- arm
  d$cell <- paste(d$family, d$shape, d$size, d$seed, sep = "/")
  d
}
arms <- lapply(args, rd)
names(arms) <- args
tab <- function(d, by) {
  sp <- split(d, d[by], drop = TRUE)
  out <- do.call(rbind, lapply(sp, function(x) {
    ok <- x$status == "ok"
    z <- (x$ll[ok] - x$ll_exact_at[ok]) / x$mcse[ok]
    gap <- abs(x$ll[ok] - x$ref[ok])
    data.frame(
      x[1L, by, drop = FALSE], fits = nrow(x), ok = sum(ok),
      rose = sum(x$status == "rose"), other_err = sum(x$status == "error"),
      capped = sum(x$capped %in% TRUE, na.rm = TRUE),
      med_rounds = stats::median(x$rounds, na.rm = TRUE),
      max_rounds = suppressWarnings(max(x$rounds, na.rm = TRUE)),
      med_secs = stats::median(x$secs),
      med_gap_ref = if (any(ok)) stats::median(gap, na.rm = TRUE) else NA,
      max_gap_ref = if (any(ok)) suppressWarnings(max(gap, na.rm = TRUE)) else NA,
      lap_gap = stats::median(abs(x$lap - x$ref), na.rm = TRUE),
      med_mcse = if (any(ok)) stats::median(x$mcse[ok]) else NA,
      med_absz = if (any(ok)) stats::median(abs(z), na.rm = TRUE) else NA,
      ess_min = if (any(ok)) min(x$ess_min[ok]) else NA, row.names = NULL)
  }))
  out[do.call(order, out[by]), ]
}
op <- options(width = 200, digits = 4)
for (a in names(arms)) {
  d <- arms[[a]]
  cat("\n=====", a, "=====\n")
  cat("fits", nrow(d), "cells", length(unique(d$cell)), "\n")
  cat("status: "); print(table(d$status))
  print(tab(d, "N"), row.names = FALSE)
  print(tab(d, "family"), row.names = FALSE)
  print(tab(d, "shape"), row.names = FALSE)
  print(tab(d, "size"), row.names = FALSE)
  cat("reference check against quadrature (int cells): max |ref - quad| ",
      format(max(abs(d$ref - d$quad), na.rm = TRUE), digits = 3), "\n")
  e <- d[d$status == "error", ]
  if (nrow(e)) {
    cat("other errors:\n")
    print(utils::head(as.data.frame(table(substr(e$msg, 1, 90))), 12),
          row.names = FALSE)
    print(utils::head(e[, c("family", "shape", "size", "seed", "N", "msg")],
                      12), row.names = FALSE)
  }
  w <- d[d$status == "ok" & d$nwarn > 0, ]
  if (nrow(w)) {
    cat("warnings on", nrow(w), "fits:\n")
    print(utils::head(sort(table(substr(w$msg, 1, 70)), decreasing = TRUE), 8))
  }
}
if (length(arms) == 2L) {
  a <- arms[[1L]]
  b <- arms[[2L]]
  m <- merge(a, b, by = c("family", "shape", "size", "seed", "N"),
             suffixes = c(".a", ".b"))
  cat("\n===== paired", names(arms)[1L], "->", names(arms)[2L], "=====\n")
  print(table(a = m$status.a, b = m$status.b))
  both <- m[m$status.a == "ok" & m$status.b == "ok", ]
  cat("paired ok fits:", nrow(both), "\n")
  cat("median secs:", stats::median(m$secs.a), "->", stats::median(m$secs.b),
      "; total secs:", sum(m$secs.a), "->", sum(m$secs.b), "\n")
  cat("median |ll - ref|:",
      stats::median(abs(both$ll.a - both$ref.a), na.rm = TRUE), "->",
      stats::median(abs(both$ll.b - both$ref.b), na.rm = TRUE), "\n")
  cat("median max|beta - ref beta|:", stats::median(both$d_beta0.a, na.rm = TRUE), "->",
      stats::median(both$d_beta0.b, na.rm = TRUE), "; median max|log sd - ref|:",
      stats::median(both$d_sd.a, na.rm = TRUE), "->",
      stats::median(both$d_sd.b, na.rm = TRUE), "
")
  cat("rounds, median:", stats::median(both$rounds.a), "->",
      stats::median(both$rounds.b), " capped:",
      sum(both$capped.a %in% TRUE), "->", sum(both$capped.b %in% TRUE), "\n")
  regress <- m[m$status.a == "ok" & m$status.b != "ok", ]
  cat("ok -> not ok:", nrow(regress), "\n")
  if (nrow(regress)) {
    print(utils::head(regress[, c("family", "shape", "size", "seed", "N",
                                  "status.b", "msg.b")], 20), row.names = FALSE)
  }
}
options(op)
