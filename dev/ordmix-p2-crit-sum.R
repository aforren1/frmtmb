# Punch round 2, B3: read dev/ordmix-p2-crit-log/ and report, per
# population and label, the spread of each candidate statistic, and the
# false alarms and misses of candidate rules. Per fit, a statistic is
# taken at its most degenerate component: the largest reach, pm and
# near, and the largest (least negative) dll2.
# Usage: Rscript dev/ordmix-p2-crit-sum.R
dir <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-p2-crit-log"
rows <- NULL
for (f in list.files(dir, "[.]txt$", full.names = TRUE)) {
  for (l in grep("^REP .* label=", readLines(f), value = TRUE)) {
    kv <- strsplit(sub("^REP ", "", l), " ")[[1]]
    v <- setNames(sub("^[^=]*=", "", kv), sub("=.*$", "", kv))
    num <- function(k) as.numeric(strsplit(v[[k]], ",")[[1]])
    rows <- rbind(rows, data.frame(
      cfg = v[["cfg"]], seed = as.integer(v[["seed"]]), label = v[["label"]],
      reach = max(num("reach")),
      collapsed = any(strsplit(v[["collapsed"]], ",")[[1]] == "TRUE"),
      pm9 = max(num("pm9")), dll2 = max(num("dll2")),
      near = max(num("near")), ll1err = as.numeric(v[["ll1err"]]),
      maxse = if ("maxse" %in% names(v)) as.numeric(v[["maxse"]]) else NA,
      maxabs = as.numeric(v[["maxabs"]]),
      degwarn = if ("degwarn" %in% names(v)) v[["degwarn"]] == "TRUE" else NA))
  }
}
# the cs() designs are read apart (end of this script): the rule
# computed here has no per-threshold probe, which the build has
cs_rows <- rows[startsWith(rows$cfg, "cs_"), ]
rows <- rows[!startsWith(rows$cfg, "cs_"), ]
rows$pop <- ifelse(startsWith(rows$cfg, "margin_"), "margin",
                   ifelse(startsWith(rows$cfg, "rev_"), "rev160", "id140"))
cat("largest |ll1err| (the evaluation against logLik()):",
    signif(max(abs(rows$ll1err), na.rm = TRUE), 3), "\n")
cat("fits per population and label:\n")
print(table(rows$pop, rows$label))
for (st in c("reach", "pm9", "dll2", "near")) {
  cat("\n==", st, "quantiles (0, .1, .5, .9, 1) by label, collapsed excluded\n")
  for (lb in c("sound", "poor", "degenerate")) {
    v <- rows[[st]][rows$label == lb & !rows$collapsed]
    cat(sprintf("  %-10s n=%3d  %s\n", lb, length(v),
                paste(signif(quantile(v, c(0, .1, .5, .9, 1), na.rm = TRUE), 3),
                      collapse = "  ")))
  }
}
cat("\nsound fits by config: max reach, pm9, dll2, near\n")
s <- rows[rows$label == "sound", ]
print(aggregate(cbind(reach, pm9, dll2, near) ~ cfg, s, max), digits = 3)
cat("\ndegenerate fits, not collapsed:\n")
dgn <- rows[rows$label == "degenerate" & !rows$collapsed, ]
print(dgn[order(dgn$dll2), c("cfg", "seed", "reach", "pm9", "dll2", "near")],
      digits = 3, row.names = FALSE)
cat("\npoor fits (finite SEs, a slope off by more than 50%):\n")
pr <- rows[rows$label == "poor", ]
print(pr[order(pr$dll2), c("cfg", "seed", "reach", "pm9", "dll2", "near")],
      digits = 3, row.names = FALSE)
rule <- function(name, fire) {
  cat(sprintf("\nrule %s\n", name))
  for (p in c("margin", "rev160", "id140")) {
    i <- rows$pop == p
    cat(sprintf(paste0("  %-7s sound %3d, warns on %2d; poor %2d, warns on ",
                       "%2d; degenerate %3d, warns on %2d (misses %2d)\n"),
                p, sum(i & rows$label == "sound"),
                sum(i & rows$label == "sound" & fire),
                sum(i & rows$label == "poor"),
                sum(i & rows$label == "poor" & fire),
                sum(i & rows$label == "degenerate"),
                sum(i & rows$label == "degenerate" & fire),
                sum(i & rows$label == "degenerate" & !fire)))
  }
}
rule("round 1: reach > 50 or collapsed", rows$reach > 50 | rows$collapsed)
for (tol in c(1e-6, 1e-3, 1e-2, 0.1)) {
  rule(sprintf("dll2 > -%g or collapsed", tol),
       (rows$dll2 > -tol) %in% TRUE | rows$collapsed)
}
for (nr in c(10, 20, 50)) {
  rule(sprintf("dll2 > -1e-3 or near > %g or collapsed", nr),
       (rows$dll2 > -1e-3) %in% TRUE | rows$near > nr | rows$collapsed)
}
# the rule taken: a step function (sharpening costs under 0.1), a
# threshold no row is near (over 50), or a collapsed gap
take <- (rows$dll2 > -0.1) %in% TRUE | rows$near > 50 | rows$collapsed
rule("TAKEN: dll2 > -0.1 or near > 50 or collapsed", take)
# the probit fits whose standard errors are NaN at a latent distance of
# about 38, where the probit's log-odds form underflows
# (dev/test-backlog.md): degenerate by the label, not by a component
ps <- grepl("^margin_probit", rows$cfg) & rows$label == "degenerate"
cat(sprintf("\nprobit saturation (label degenerate, probit margin designs): %d fits, max reach %.1f, max near %.3g, dll2 NaN in %d; the rule warns on %d\n",
            sum(ps), max(rows$reach[ps]), max(rows$near[ps]),
            sum(is.na(rows$dll2[ps])), sum(take[ps])))
rule("TAKEN, probit saturation left out", take)
keep <- !ps
for (p in c("margin", "rev160", "id140")) {
  i <- rows$pop == p & keep
  cat(sprintf("  %-7s degenerate %3d, warns on %2d (misses %2d)\n", p,
              sum(i & rows$label == "degenerate"),
              sum(i & rows$label == "degenerate" & take),
              sum(i & rows$label == "degenerate" & !take)))
}
cat("\nmisses of the rule taken (probit saturation left out):\n")
m <- rows[rows$label == "degenerate" & !take & !ps, ]
print(m[, c("cfg", "seed", "reach", "pm9", "dll2", "near")], digits = 3,
      row.names = FALSE)
cat(sprintf("\nmargins: sound and poor fits, largest dll2 %.3f, largest near %.3g; collapsed among sound %d\n",
            max(rows$dll2[rows$label != "degenerate"], na.rm = TRUE),
            max(rows$near[rows$label != "degenerate"]),
            sum(rows$collapsed[rows$label == "sound"])))
# the installed build's own warning, where the log has it: it must
# agree with the rule computed here
dw <- NULL
for (f in list.files(dir, "[.]txt$", full.names = TRUE)) {
  for (l in grep("^REP .* degwarn=", readLines(f), value = TRUE)) {
    dw <- c(dw, setNames(grepl("degwarn=TRUE", l),
                         sub("^REP cfg=([^ ]*) seed=([0-9]*).*$", "\\1#\\2", l)))
  }
}
if (length(dw)) {
  key <- paste0(rows$cfg, "#", rows$seed)
  fire <- unname(dw[key])
  cat(sprintf("\ninstalled build's warning: logged on %d fits; agrees with the rule on %d\n",
              sum(!is.na(fire)), sum(fire == take, na.rm = TRUE)))
  rule("INSTALLED BUILD's degenerate warning", fire %in% TRUE)
}
# the cs() designs, by the installed build's own warning, under the
# review's label and under one that also reads the standard errors (an
# SE above 30 means nothing, as an estimate above 30 does)
if (nrow(cs_rows)) {
  cat("\ncs() designs, the installed build's warning:\n")
  se_lab <- ifelse(cs_rows$label == "degenerate" | !is.finite(cs_rows$maxse) |
                     cs_rows$maxse > 30, "degenerate", "sound")
  for (cf in unique(cs_rows$cfg)) {
    i <- cs_rows$cfg == cf
    for (lab in list(list("review label", cs_rows$label),
                     list("label with SE > 30", se_lab))) {
      L <- lab[[2]]
      cat(sprintf("  %-10s %-18s sound %2d, warns on %d; degenerate %2d, warns on %2d (misses %d)\n",
                  cf, lab[[1]], sum(i & L == "sound"),
                  sum(i & L == "sound" & cs_rows$degwarn),
                  sum(i & L == "degenerate"),
                  sum(i & L == "degenerate" & cs_rows$degwarn),
                  sum(i & L == "degenerate" & !cs_rows$degwarn)))
    }
  }
  cat("  fits where the two disagree with the warning:\n")
  print(cs_rows[(cs_rows$label == "sound") == cs_rows$degwarn |
                  (se_lab == "sound") == cs_rows$degwarn,
                c("cfg", "seed", "label", "maxabs", "maxse", "degwarn")],
        row.names = FALSE, digits = 3)
}
