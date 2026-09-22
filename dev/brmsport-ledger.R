# Turn the recorded runs of the brms suite tier into the ledger.
#
#   sh dev/brmsport-record.sh          # record every generated file
#   Rscript dev/brmsport-ledger.R [--partial]
#   Rscript dev/brmsport-gen.R         # regenerate with the verdicts
#
# Inputs: dev/brmsport-log/rec-<pkg>-<topic>.tsv, one row per assertion
# per package, written by the harness in record mode; and
# dev/brmsport-verdicts-manual.tsv, the verdict and reason of every
# assertion that does NOT hold.
#
# A PASS is not typed: it is an assertion observed to hold, genuinely
# (not on a missing function, not on a stale object, not reaching brms's
# namespace). The script STOPS when
#   - an assertion that did not hold has no manual verdict,
#   - a manual verdict names an assertion that held (the verdict is
#     stale: frmtmb changed, or the verdict was wrong),
#   - a manual verdict is "pass",
#   - a both-tier assertion holds in one package and not the other
#     while its verdict is shared,
#   - an assertion of bin 1 was never recorded (unless --partial).
# Outputs: dev/brmsport-verdicts.tsv (read by the generator),
# dev/brmsport-ledger.tsv (one row per assertion), and the generated
# summary dev/brmsport-log/ledger-summary.md.
source("dev/brmsport-blocks.R")
partial <- "--partial" %in% commandArgs(trailingOnly = TRUE)

outcomes <- c("pass", "defect", "divergence", "pending 2.6d",
              "pending 2.6e", "cannot transfer")

recs <- list.files("dev/brmsport-log", pattern = "^rec-.*[.]tsv$",
                   full.names = TRUE)
cols <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
          "caught", "raw_held")
rec <- do.call(rbind, lapply(recs, function(f) {
  x <- utils::read.delim(f, header = FALSE, quote = "", col.names = cols,
                         colClasses = "character", na.strings = NULL)
  x[x$kind == "assert", ]
}))
rec$held <- rec$held == "TRUE"
rec$vacuous <- rec$vacuous == "TRUE"
rec$raw_held <- rec$raw_held == "TRUE"
stopifnot(!anyDuplicated(rec[c("pkg", "id")]))

man <- do.call(rbind, lapply(
  c("dev/brmsport-verdicts-manual.tsv",
    "dev/brmsport-verdicts-manual-fit.tsv"),
  utils::read.delim, quote = "", colClasses = "character",
  na.strings = NULL))
own <- utils::read.delim("dev/brmsport-verdicts-own.tsv", quote = "",
                         colClasses = "character", na.strings = NULL)
# A pass that is WEAKER than the assertion brms wrote. It still holds,
# so it is not a defect and no verdict can be given for it, but a
# reader has to be told: `brmsfit-methods:340` compares
# `c(19.122, NA, NA, NA)` with itself now, 3 of its 4 cells NA on both
# sides, where before the shapes item it compared two numbers. Without
# this file such a row moves from defect to pass and reads like a fix.
notes <- utils::read.delim("dev/brmsport-notes.tsv", quote = "",
                           colClasses = "character", na.strings = NULL)
stopifnot(!anyDuplicated(notes[c("id", "pkg")]), all(nzchar(notes$note)))
stopifnot(!anyDuplicated(own[c("id", "pkg")]),
          !any(paste(own$id) %in% paste(man$id)))
stopifnot(!anyDuplicated(man[c("id", "pkg")]),
          all(man$verdict %in% outcomes), !any(man$verdict == "pass"),
          all(nzchar(man$reason)))

# every bin-1 assertion, from brms's source
blocks <- brms_bin1_blocks()
topic_of <- function(f) sub("^tests[.](.*)[.]R$", "\\1", f)
exp_rows <- list()
for (b in blocks) {
  for (j in seq_along(b$stmts)) {
    if (n_expect(b$stmts[[j]]) == 0L) next
    txt <- gsub("[\r\n\t]+", " ", b$texts[j])
    txt <- gsub(" +", " ", txt)
    exp_rows[[length(exp_rows) + 1L]] <- data.frame(
      file = b$file, block = gsub("[\r\n\t]+", " ", b$label),
      line = b$lines[j], tier = b$tier,
      id = sprintf("%s:%d", topic_of(b$file), b$lines[j]),
      assertion = substr(txt, 1, 160), stringsAsFactors = FALSE)
  }
}
led <- do.call(rbind, exp_rows)
stopifnot(nrow(led) == 494L, !anyDuplicated(led$id))

pkgs_of <- function(tier) switch(tier, core = "frmtmb",
                                 sample = "frmtmb.sample",
                                 both = c("frmtmb", "frmtmb.sample"))
problems <- character()
out_v <- list()
led$outcome <- NA_character_
led$class <- ""
led$reason <- ""
led$held <- ""
led$message <- ""
for (i in seq_len(nrow(led))) {
  id <- led$id[i]
  runs <- rec[rec$id == id, ]
  want <- pkgs_of(led$tier[i])
  missing <- setdiff(want, runs$pkg)
  if (length(missing)) {
    problems <- c(problems, sprintf("NOT RECORDED %s in %s", id,
                                    paste(missing, collapse = ",")))
    next
  }
  runs <- runs[match(want, runs$pkg), ]
  led$held[i] <- paste(sprintf("%s=%s", sub("frmtmb[.]?", "", runs$pkg),
                               runs$held), collapse = " ")
  led$held[i] <- sub("^=", "core=", led$held[i])
  led$message[i] <- paste(unique(runs$msg[nzchar(runs$msg)]),
                          collapse = " || ")
  for (k in seq_len(nrow(runs))) {
    p <- runs$pkg[k]
    m <- man[man$id == id & man$pkg %in% c(p, "*"), ]
    if (nrow(m) > 1L) m <- m[m$pkg == p, ]
    o <- own[own$id == id & own$pkg %in% c(p, "*"), ]
    if (nrow(o)) {
      if (!runs$held[k]) {
        problems <- c(problems, sprintf(
          "OWN-WORDS row does not hold: %s in %s: %s", id, p,
          substr(runs$msg[k], 1, 150)))
        next
      }
      v <- data.frame(id = id, pkg = p, verdict = "pass",
                      reason = paste("frmtmb's own words:", o$pattern,
                                     "|", o$note))
      out_v[[length(out_v) + 1L]] <- v
      cls <- "own-words"
    } else if (runs$held[k] && nrow(m) && identical(m$class, "hollow")) {
      # held, but on something that is not frmtmb's answer (a NULL slot
      # read as an empty attribute); the harness asserts it keeps holding
      v <- data.frame(id = id, pkg = p, verdict = "hollow",
                      reason = m$reason)
      out_v[[length(out_v) + 1L]] <- v
      v$verdict <- m$verdict
      cls <- "hollow"
    } else if (runs$held[k]) {
      if (nrow(m)) {
        problems <- c(problems, sprintf(
          "HOLDS but has manual verdict '%s': %s in %s", m$verdict, id, p))
      }
      nt <- notes[notes$id == id & notes$pkg %in% c(p, "*"), ]
      v <- data.frame(id = id, pkg = p, verdict = "pass",
                      reason = if (nrow(nt)) nt$note[1L] else "")
      out_v[[length(out_v) + 1L]] <- v
      cls <- if (nrow(nt)) "weak-pass" else ""
    } else {
      if (!nrow(m)) {
        problems <- c(problems, sprintf("UNCLASSIFIED %s in %s: %s", id, p,
                                        substr(runs$msg[k], 1, 150)))
        next
      }
      if (identical(m$class, "hollow")) {
        problems <- c(problems, sprintf(
          "marked hollow but does NOT hold: %s in %s", id, p))
      }
      v <- data.frame(id = id, pkg = p, verdict = m$verdict,
                      reason = m$reason)
      out_v[[length(out_v) + 1L]] <- v
      cls <- m$class
    }
    if (is.na(led$outcome[i])) {
      led$outcome[i] <- v$verdict
      led$reason[i] <- v$reason
      led$class[i] <- cls
    } else if (!identical(led$outcome[i], v$verdict)) {
      problems <- c(problems, sprintf(
        "TIERS DISAGREE on %s: %s against %s", id, led$outcome[i],
        v$verdict))
    }
  }
}

if (length(problems)) {
  cat(problems, sep = "\n")
  n_unrec <- sum(grepl("^NOT RECORDED", problems))
  if (!partial || length(problems) > n_unrec) {
    stop(length(problems), " problem(s); the ledger was not written")
  }
}

verd <- do.call(rbind, out_v)
write.table(verd, "dev/brmsport-verdicts.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
ledger <- led[c("file", "block", "line", "tier", "outcome", "class",
                "reason", "id", "held", "assertion", "message")]
ledger$outcome[is.na(ledger$outcome)] <- "not recorded"
write.table(ledger, "dev/brmsport-ledger.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)

# ---- the generated summary ----
ledger$outcome <- factor(ledger$outcome,
                         levels = c(outcomes, "not recorded"))
tab <- table(ledger$file, ledger$outcome)
tab <- tab[, colSums(tab) > 0, drop = FALSE]
lines <- c("<!-- generated by dev/brmsport-ledger.R; do not edit by hand -->",
           "", "### Per file", "",
           paste0("| file | assertions | ",
                  paste(colnames(tab), collapse = " | "), " |"),
           paste0("|---|---|", strrep("---|", ncol(tab))))
for (f in rownames(tab)) {
  lines <- c(lines, sprintf("| `%s` | %d | %s |", f, sum(tab[f, ]),
                            paste(tab[f, ], collapse = " | ")))
}
lines <- c(lines, sprintf("| **total** | **%d** | %s |", sum(tab),
                          paste0("**", colSums(tab), "**",
                                 collapse = " | ")))
ctab <- table(ledger$outcome, ledger$class)
lines <- c(lines, "", "### Outcome by class", "",
           "| outcome | class | assertions |", "|---|---|---|")
for (a in seq_len(nrow(ctab))) for (b in seq_len(ncol(ctab))) {
  cl <- colnames(ctab)[b]
  if (ctab[a, b] > 0) {
    lines <- c(lines, sprintf("| %s | %s | %d |", rownames(ctab)[a],
                              if (nzchar(cl)) cl else "-", ctab[a, b]))
  }
}
np <- sum(ledger$outcome == "pass")
lines <- c(lines, "",
           sprintf("Bin 1 passes: %d of 494 (%.1f%%).", np, 100 * np / 494),
           sprintf("Against bins 1 and 2: %d of 823 (%.1f%%); bin 2 was not ported.",
                   np, 100 * np / 823),
           sprintf(paste("Runs testthat alone would count as a pass and the",
                         "harness does not (missing function or object,",
                         "stale object, argument-name refusal, a NULL",
                         "read through a partial $ match), over both",
                         "packages: %d; hollow passes marked by hand: %d."),
                   sum(rec$raw_held & !rec$held),
                   sum(ledger$class == "hollow")))
# the sample half, its own totals: every run in frmtmb.sample
vs <- verd[verd$pkg == "frmtmb.sample", ]
vs$verdict[vs$verdict == "hollow"] <- "pass (hollow)"
vs$file <- led$file[match(vs$id, led$id)]
vs$tier <- led$tier[match(vs$id, led$id)]
stab <- table(paste0(vs$file, " (", vs$tier, ")"), vs$verdict)
lines <- c(lines, "", "### The frmtmb.sample half", "",
           paste0("| file (tier) | assertions | ",
                  paste(colnames(stab), collapse = " | "), " |"),
           paste0("|---|---|", strrep("---|", ncol(stab))))
for (f in rownames(stab)) {
  lines <- c(lines, sprintf("| `%s` | %d | %s |", f, sum(stab[f, ]),
                            paste(stab[f, ], collapse = " | ")))
}
lines <- c(lines, sprintf("| **total** | **%d** | %s |", sum(stab),
                          paste0("**", colSums(stab), "**",
                                 collapse = " | ")),
           "", sprintf(paste("frmtmb.sample passes %d of its %d runs (%.1f%%):",
                             "%d of the 61 sample-tier and %d of the 104",
                             "both-tier assertions."),
                       sum(vs$verdict == "pass"), nrow(vs),
                       100 * mean(vs$verdict == "pass"),
                       sum(vs$verdict == "pass" & vs$tier == "sample"),
                       sum(vs$verdict == "pass" & vs$tier == "both")))
writeLines(lines, "dev/brmsport-log/ledger-summary.md")
cat(lines, sep = "\n")
