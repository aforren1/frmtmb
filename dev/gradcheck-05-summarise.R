# Counts from the design-set replicates. Emits marked blocks that go into
# dev/gradcheck-findings.md verbatim, so no count in that file is typed
# by hand.
#
#   Rscript dev/gradcheck-05-summarise.R
#
# Reads both builds when both are present, so `warned` is what each build
# ACTUALLY did rather than a criterion re-derived here.

load_build <- function(which_lib) {
  fs <- Sys.glob(file.path("dev",
                           paste0("gradcheck-03-", which_lib, "-g*.rds")))
  if (!length(fs)) return(NULL)
  rows <- unlist(lapply(fs, readRDS), recursive = FALSE)
  col <- function(nm, mode = "numeric") {
    vapply(rows, function(r) {
      v <- r[[nm]]
      if (is.null(v) || length(v) != 1L) {
        if (identical(mode, "character")) NA_character_ else NA_real_
      } else if (identical(mode, "character")) {
        as.character(v)
      } else as.numeric(v)
    }, if (identical(mode, "character")) "" else 0)
  }
  d <- data.frame(design = col("design", "character"), n = col("n"),
                  rep = col("rep"), seed = col("seed"),
                  err = col("err", "character"),
                  warned = col("warned"), np = col("np"),
                  gmax = col("gmax"), gproj = col("gproj"),
                  nactive = col("nactive"), decr = col("decr"),
                  gscaled = col("gscaled"), grel = col("grel"),
                  obj = col("obj"), conv = col("conv"),
                  stringsAsFactors = FALSE)
  d[order(d$design, d$n, d$rep), ]
}

B <- load_build("base")
L <- load_build("lane")
stopifnot(!is.null(B), !is.null(L), nrow(B) == nrow(L),
          identical(B$seed, L$seed), !anyNA(B$seed),
          anyDuplicated(B$seed) == 0L)
cat("rows per build:", nrow(B), " distinct seeds:", length(unique(B$seed)),
    " fit errors:", sum(!is.na(B$err)), sum(!is.na(L$err)), "\n")

gt <- 1e-3
cat("\n<<<BLOCK design-identical>>>\n")
cat("Every fit in the design set, both builds, same seed:\n")
cat("  objectives bitwise identical :",
    identical(B$obj, L$obj), "\n")
cat("  max|grad| bitwise identical  :",
    identical(B$gmax, L$gmax), "\n")
cat("  optimizer codes identical    :", identical(B$conv, L$conv), "\n")
cat("  largest objective difference :",
    format(max(abs(B$obj - L$obj)), digits = 3), "\n")
cat("<<<END design-identical>>>\n")

cat("\n<<<BLOCK design-falsealarm>>>\n")
cat("Warnings raised on CORRECT fits. `base` and `lane` are what the\n")
cat("build did; `proj` and `head` are the two stages, scored here.\n\n")
cat(sprintf("%-11s %6s %5s %5s %5s %5s %5s %11s %11s %11s\n",
            "design", "n", "reps", "base", "lane", "proj", "head",
            "max gmax", "max gproj", "max head"))
key <- paste(B$design, B$n)
tot <- c(base = 0L, lane = 0L, proj = 0L, head = 0L)
for (k in unique(key)) {
  i <- key == k
  a <- B[i, ]
  b <- L[i, ]
  pj <- sum(is.finite(a$gproj) & a$gproj > gt)
  hd <- sum(is.finite(a$gproj) & a$gproj > gt &
              (!is.finite(a$decr) | a$decr > gt))
  tot["base"] <- tot["base"] + sum(a$warned > 0, na.rm = TRUE)
  tot["lane"] <- tot["lane"] + sum(b$warned > 0, na.rm = TRUE)
  tot["proj"] <- tot["proj"] + pj
  tot["head"] <- tot["head"] + hd
  cat(sprintf("%-11s %6d %5d %5d %5d %5d %5d %11.3e %11.3e %11.3e\n",
              a$design[1L], a$n[1L], nrow(a),
              sum(a$warned > 0, na.rm = TRUE),
              sum(b$warned > 0, na.rm = TRUE), pj, hd,
              max(a$gmax, na.rm = TRUE), max(a$gproj, na.rm = TRUE),
              max(a$decr, na.rm = TRUE)))
}
cat(sprintf("\n%-11s %6s %5d %5d %5d %5d %5d\n", "TOTAL", "", nrow(B),
            tot[["base"]], tot[["lane"]], tot[["proj"]], tot[["head"]]))
cat("\nfalse-alarm rate on the reference build:", tot[["base"]], "/",
    nrow(B), "=", format(tot[["base"]] / nrow(B), digits = 3), "\n")
cat("false-alarm rate on the lane build     :", tot[["lane"]], "/",
    nrow(L), "=", format(tot[["lane"]] / nrow(L), digits = 3), "\n")
cat("the two stages scored here, for the mechanism:",
    tot[["proj"]], "after the bound projection,", tot[["head"]],
    "after the headroom\n")
cat("largest headroom over all", nrow(B), "correct fits:",
    format(max(B$decr, na.rm = TRUE), digits = 4),
    "at", B$design[which.max(B$decr)], "n =", B$n[which.max(B$decr)],
    "seed", B$seed[which.max(B$decr)], "\n")
cat("  tolerance", gt, ", margin",
    format(gt / max(B$decr, na.rm = TRUE), digits = 4), "x\n")
cat("headrooms that could not be computed:", sum(!is.finite(B$decr)),
    "\n")
cat("<<<END design-falsealarm>>>\n")

# the n-dependence, which is the mechanism: nlminb stops on a RELATIVE
# objective change, so max|grad| grows with the objective while the
# headroom does not
cat("\n<<<BLOCK design-ndependence>>>\n")
cat(sprintf("%6s %6s %11s %11s %11s %11s\n", "n", "fits", "med gmax",
            "med |obj|", "med head", "warn rate"))
for (nn in sort(unique(B$n))) {
  a <- B[B$n == nn & is.na(B$err), ]
  cat(sprintf("%6d %6d %11.3e %11.3e %11.3e %11.3f\n", nn, nrow(a),
              stats::median(a$gmax, na.rm = TRUE),
              stats::median(abs(a$obj), na.rm = TRUE),
              stats::median(a$decr, na.rm = TRUE),
              mean(a$warned > 0, na.rm = TRUE)))
}
cat("<<<END design-ndependence>>>\n")
