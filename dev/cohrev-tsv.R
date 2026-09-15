## Re-derive the lane's headline table from its own TSV, independently
## of dev/coh-summarize.R, and answer the four questions the review
## asked of it: the drops, the bias, the too-wide rung, and Wilson.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

read_kv <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  rows <- lapply(ln, function(l) {
    kv <- strsplit(strsplit(l, "\t", fixed = TRUE)[[1L]], "=",
                   fixed = TRUE)
    setNames(trimws(vapply(kv, function(x) paste(x[-1L], collapse = "="),
                           character(1))),
             trimws(vapply(kv, `[`, character(1), 1L)))
  })
  nms <- unique(unlist(lapply(rows, names)))
  d <- as.data.frame(do.call(rbind, lapply(rows, function(r) r[nms])),
                     stringsAsFactors = FALSE)
  names(d) <- nms
  for (cl in c("rep", "seed", "est", "se", "lo", "hi", "width", "sd_id",
               "sd_idcond", "loglik", "conv", "maxgrad", "nbadse",
               "secs")) {
    if (cl %in% names(d)) d[[cl]] <- as.numeric(d[[cl]])
  }
  for (cl in c("ok", "covers", "pdhess")) {
    if (cl %in% names(d)) d[[cl]] <- d[[cl]] == "TRUE"
  }
  d
}
wilson <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  cen <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(cen - hw, cen + hw)
}
truth <- 0.5
rungs <- c("cond", "smooth", "id", "idcond", "full")

m <- read_kv("dev/coh-recovery-main.tsv")
n <- read_kv("dev/coh-recovery-null.tsv")

## --- 1. the drops -------------------------------------------------
cat("==== completeness, main arm ====\n")
tb <- table(m$seed)
part <- names(tb)[tb != 5L]
cat("seeds with fewer than 5 rungs:", length(part), "\n")
for (s in part) {
  r <- m[as.character(m$seed) == s, ]
  cat(sprintf("  %s: %s | all ok: %s | conv: %s\n", s,
              paste(r$rung, collapse = ","), all(r$ok),
              paste(r$conv, collapse = ",")))
}
cat("is every partial a PREFIX of the fit order?",
    all(vapply(part, function(s) {
      r <- m[as.character(m$seed) == s, ]
      identical(r$rung, rungs[seq_len(nrow(r))])
    }, logical(1))), "\n")
cat("rows with ok = FALSE anywhere:", sum(!m$ok), "\n")
cat("does 'full' appear in ANY partial replicate?",
    any(m$rung[as.character(m$seed) %in% part] == "full"), "\n")

## --- 2. the headline table, complete replicates only ---------------
keep <- names(tb)[tb == 5L]
d <- m[as.character(m$seed) %in% keep, ]
cat("\n==== complete-replicate table (n =", length(keep), ") ====\n")
cat(sprintf("%-8s %4s %8s %-16s %9s %9s %8s %8s %9s %8s\n", "rung",
            "n", "cover", "wilson95", "width_rel", "mean_est",
            "sd_est", "mean_se", "bias_abs", "sd_z"))
fw <- mean(d$width[d$rung == "full"])
for (rg in rungs) {
  s <- d[d$rung == rg, ]
  ci <- wilson(sum(s$covers), nrow(s))
  z <- (s$est - truth) / s$se
  cat(sprintf(paste0("%-8s %4d %8.3f (%.3f, %.3f) %9.3f %9.4f %8.4f",
                     " %8.4f %9.4f %8.3f\n"),
              rg, nrow(s), mean(s$covers), ci[1], ci[2],
              mean(s$width) / fw, mean(s$est), sd(s$est), mean(s$se),
              mean(s$est) - truth, sd(z)))
}
cat("\nbias in units of the estimator's own spread:\n")
for (rg in rungs) {
  s <- d[d$rung == rg, ]
  cat(sprintf("  %-8s %+.4f absolute, %+.3f of sd(est), t = %.2f\n",
              rg, mean(s$est) - truth,
              (mean(s$est) - truth) / sd(s$est),
              (mean(s$est) - truth) / (sd(s$est) / sqrt(nrow(s)))))
}

## --- 3. what the drops would do if handled the OTHER way -----------
cat("\n==== every available row, partial replicates INCLUDED ====\n")
for (rg in rungs) {
  s <- m[m$rung == rg, ]
  t2 <- d[d$rung == rg, ]
  ci <- wilson(sum(s$covers), nrow(s))
  cat(sprintf(paste0("%-8s all rows %3d cover %.4f (%.3f, %.3f) |",
                     " complete %3d cover %.4f\n"),
              rg, nrow(s), mean(s$covers), ci[1], ci[2], nrow(t2),
              mean(t2$covers)))
}

## --- 4. the too-wide rung: is the id variance absorbed? ------------
cat("\n==== variance components by rung, main arm ====\n")
for (rg in rungs) {
  s <- d[d$rung == rg, ]
  cat(sprintf("%-8s sd_id %s  sd_idcond %s\n", rg,
              if (all(is.na(s$sd_id))) "     -" else
                sprintf("%.4f +/- %.4f", mean(s$sd_id), sd(s$sd_id)),
              if (all(is.na(s$sd_idcond))) "     -" else
                sprintf("%.4f +/- %.4f", mean(s$sd_idcond),
                        sd(s$sd_idcond))))
}
tot <- sqrt(0.35^2 + 0.20^2)
ic <- d[d$rung == "idcond", ]
fl <- d[d$rung == "full", ]
cat(sprintf("\nsqrt(sd_id^2 + sd_idcond^2) at the truth = %.4f\n", tot))
cat(sprintf("idcond rung's sd_idcond           = %.4f\n",
            mean(ic$sd_idcond)))
cat(sprintf("full rung's sqrt(sd_id^2+sd_ic^2) = %.4f\n",
            mean(sqrt(fl$sd_id^2 + fl$sd_idcond^2))))
cat(sprintf("predicted width ratio idcond/full = %.3f\n",
            mean(ic$sd_idcond) / mean(fl$sd_idcond)))
mm <- merge(ic[, c("seed", "width", "se", "sd_idcond")],
            fl[, c("seed", "width", "se", "sd_idcond")], by = "seed",
            suffixes = c("_ic", "_full"))
r <- mm$width_ic / mm$width_full
cat(sprintf(paste0("observed paired width ratio idcond/full:",
                   " min %.3f median %.3f max %.3f mean %.3f\n"),
            min(r), median(r), max(r), mean(r)))
cat(sprintf(paste0("per-seed width ratio against per-seed component",
                   " ratio: r = %.3f, slope %.3f\n"),
            cor(r, mm$sd_idcond_ic / mm$sd_idcond_full),
            coef(lm(r ~ I(mm$sd_idcond_ic / mm$sd_idcond_full)))[2L]))

## --- 5. does "(1 | id) cancels out of the contrast" hold? ----------
cat("\n==== does adding (1 | id) move the contrast? ====\n")
sm <- d[d$rung == "smooth", ]
idr <- d[d$rung == "id", ]
p <- merge(sm[, c("seed", "est", "se", "width")],
           idr[, c("seed", "est", "se", "width")], by = "seed",
           suffixes = c("_sm", "_id"))
cat(sprintf(paste0("smooth -> id  (adds (1|id) only): paired se ratio",
                   " median %.4f, range %.4f to %.4f\n"),
            median(p$se_id / p$se_sm), min(p$se_id / p$se_sm),
            max(p$se_id / p$se_sm)))
p2 <- merge(ic[, c("seed", "se")], fl[, c("seed", "se")], by = "seed",
            suffixes = c("_ic", "_full"))
cat(sprintf(paste0("idcond -> full (adds (1|id) only): paired se ratio",
                   " median %.4f, range %.4f to %.4f\n"),
            median(p2$se_full / p2$se_ic), min(p2$se_full / p2$se_ic),
            max(p2$se_full / p2$se_ic)))

## --- 6. null arm ---------------------------------------------------
cat("\n==== null arm ====\n")
tbn <- table(n$seed)
cat("replicates:", sum(tbn == 5L), " partial:", sum(tbn != 5L), "\n")
dn <- n[as.character(n$seed) %in% names(tbn)[tbn == 5L], ]
fwn <- mean(dn$width[dn$rung == "full"])
for (rg in rungs) {
  s <- dn[dn$rung == rg, ]
  z <- (s$est - truth) / s$se
  cat(sprintf(paste0("%-8s cover %2d/%2d  width_rel %.3f  sd_z %.3f",
                     "  mean_w %.4f  bias %+.4f\n"),
              rg, sum(s$covers), nrow(s), mean(s$width) / fwn, sd(z),
              mean(s$width), mean(s$est) - truth))
}
an <- dn[dn$rung == "id", ]
bn <- dn[dn$rung == "full", ]
mn <- merge(an[, c("seed", "width")], bn[, c("seed", "width")],
            by = "seed", suffixes = c("_id", "_full"))
rn <- mn$width_full / mn$width_id
cat(sprintf(paste0("null arm paired width ratio full/id: min %.3f",
                   " median %.3f max %.3f\n"),
            min(rn), median(rn), max(rn)))
cat(sprintf(paste0("null arm full rung sd_idcond: %.4f +/- %.4f;",
                   " below 0.01 on %d of %d\n"),
            mean(bn$sd_idcond), sd(bn$sd_idcond),
            sum(bn$sd_idcond < 0.01), nrow(bn)))
cat(sprintf(paste0("main arm survey mean width %.4f vs null arm full",
                   " mean width %.4f, ratio %.4f\n"),
            mean(d$width[d$rung == "id"]), mean(bn$width),
            mean(bn$width) / mean(d$width[d$rung == "id"])))

## --- 7. main arm paired full vs id, and the 53 percent -------------
a <- d[d$rung == "id", ]
mm4 <- merge(a[, c("seed", "width", "covers")],
             fl[, c("seed", "width", "covers")], by = "seed",
             suffixes = c("_id", "_full"))
rr <- mm4$width_full / mm4$width_id
cat(sprintf(paste0("\nmain arm paired width ratio full/id: min %.3f",
                   " median %.3f max %.3f\n"),
            min(rr), median(rr), max(rr)))
cat(sprintf(paste0("mean width ratio id/full = %.4f, that is %.1f",
                   " percent narrower\n"),
            mean(a$width) / mean(fl$width),
            100 * (1 - mean(a$width) / mean(fl$width))))
cat(sprintf(paste0("mean of the PAIRED ratio id/full = %.4f (%.1f",
                   " percent narrower)\n"),
            mean(mm4$width_id / mm4$width_full),
            100 * (1 - mean(mm4$width_id / mm4$width_full))))
cat(sprintf(paste0("id missed while full covered: %d; full missed",
                   " while id covered: %d\n"),
            sum(!mm4$covers_id & mm4$covers_full),
            sum(mm4$covers_id & !mm4$covers_full)))

## --- 8. convergence nuisance --------------------------------------
cat("\n==== convergence, main arm, complete replicates ====\n")
bad <- d[d$conv != 0 | !d$pdhess, ]
cat("rows with conv != 0 or no pd Hessian:", nrow(bad), "\n")
if (nrow(bad)) {
  print(bad[, c("rung", "seed", "conv", "maxgrad", "pdhess", "covers")])
}
fb <- fl[fl$conv != 0, ]
cat(sprintf(paste0("full rung: %d of %d with conv != 0; dropping",
                   " them: %d/%d = %.4f (vs %d/%d = %.4f)\n"),
            nrow(fb), nrow(fl), sum(fl$covers[fl$conv == 0]),
            sum(fl$conv == 0), mean(fl$covers[fl$conv == 0]),
            sum(fl$covers), nrow(fl), mean(fl$covers)))

## --- 9. Wilson half-width resolution ------------------------------
cat("\n==== Wilson resolution at n = 148 ====\n")
for (pp in c(0.95, 0.5)) {
  ci <- wilson(round(pp * 148), 148)
  cat(sprintf("p = %.2f: (%.4f, %.4f), half-width %.4f\n", pp, ci[1],
              ci[2], diff(ci) / 2))
}
ci <- wilson(141, 148)
cat(sprintf("141/148 = %.4f, Wilson (%.4f, %.4f)\n", 141 / 148, ci[1],
            ci[2]))
ci <- wilson(78, 148)
cat(sprintf(" 78/148 = %.4f, Wilson (%.4f, %.4f)\n", 78 / 148, ci[1],
            ci[2]))
