# Turn the recorded bin-2 runs and the 2P survey into the bin-2 ledger,
# and compose its totals with bin 1's.
#
#   sh dev/brmsport-bin2-record.sh       # record every generated file
#   Rscript dev/brmsport-bin2-ledger.R
#   Rscript dev/brmsport-bin2-gen.R      # regenerate with the verdicts
#
# The rules are bin 1's (dev/brmsport-ledger.R): a pass is an assertion
# OBSERVED to hold, genuinely; every run that does not hold needs a
# manual verdict; a manual verdict on a run that holds is stale, unless
# it marks the run hollow. The script STOPS on any of these. Bin 2 adds:
#
# - a TARGET per ported row (dev/brmsport-bin2-targets.tsv, or the one
#   package a 2P adapter block runs in): a 2T row runs on the ML fit and
#   on draws, and the row's outcome is the target's verdict. The other
#   run keeps its own verdict, shown in the `other` column.
# - every 2P row that is not ported takes its group's class from
#   dev/brmsport-bin2-survey.tsv (dev/brmsport-bin2-survey.R).
#
# A ledger row is one STATEMENT of brms's block; `n` is the number of
# assertions it holds (a loop or SW(expect_...) holds several), and the
# totals count assertions, as the audit and bin 1 do.
#
# Outputs: dev/brmsport-bin2-verdicts.tsv (read by the generator),
# dev/brmsport-bin2-ledger.tsv, dev/brmsport-bin2-summary.md. For a dry
# run, FRMTMB_BIN2_REC_DIR names another directory of recorded runs and
# FRMTMB_BIN2_LEDGER_OUT another directory for the outputs; the bin-1
# ledger is read where it is.
#
# A row that holds only because its adapter substituted brms's input
# (dev/brmsport-bin2-adapters-spec.R, `bin2_adapted_input`) is a pass of
# class "adapted-input", with the substitution as its reason, so the
# summary counts it apart from the literal passes.
source("dev/brmsport-bin2-blocks.R")
source("dev/brmsport-bin2-adapters-spec.R")
topic_of <- function(f) sub("^tests[.](.*)[.]R$", "\\1", f)
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
outcomes <- c("pass", "defect", "divergence", "cannot transfer")
# Classes of a run that HOLDS on something other than frmtmb's answer:
# "hollow" (a NULL or empty read, bin 1's meaning) and "adapter" (the
# assertion tests the bin-2P adapter's own arithmetic). The tier asserts
# such a row keeps holding, and the ledger gives it the manual verdict.
held_classes <- c("hollow", "adapter")

rec_dir <- Sys.getenv("FRMTMB_BIN2_REC_DIR", "dev/brmsport-bin2-log")
out_dir <- Sys.getenv("FRMTMB_BIN2_LEDGER_OUT", "dev")
out_path <- function(f) file.path(out_dir, f)
recs <- list.files(rec_dir, pattern = "^rec-.*[.]tsv$",
                   full.names = TRUE)
cols <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
          "caught", "raw_held")
rec_all <- do.call(rbind, lapply(recs, function(f) {
  utils::read.delim(f, header = FALSE, quote = "", col.names = cols,
                    colClasses = "character", na.strings = NULL)
}))
rec <- rec_all[rec_all$kind == "assert", ]
rec$held <- rec$held == "TRUE"
rec$raw_held <- rec$raw_held == "TRUE"
stopifnot(!anyDuplicated(rec[c("pkg", "id")]))

man_files <- list.files("dev", "^brmsport-bin2-verdicts-manual.*[.]tsv$",
                        full.names = TRUE)
man <- do.call(rbind, lapply(man_files, rd))
stopifnot(!anyDuplicated(man[c("id", "pkg")]), all(man$verdict %in% outcomes),
          !any(man$verdict == "pass"), all(nzchar(man$reason)))
tg_files <- list.files("dev", "^brmsport-bin2-targets.*[.]tsv$",
                       full.names = TRUE)
tg <- do.call(rbind, lapply(tg_files, rd))
stopifnot(!anyDuplicated(tg$id))
sv <- rd("dev/brmsport-bin2-survey.tsv")
sv$n <- as.integer(sv$n)
groups <- rd("dev/brmsport-bin2-groups.tsv")

rows <- list()
for (bin in c("2T", "2P")) {
  for (b in brms_bin_blocks(bin)) {
    ported <- bin == "2T" || (b$file %in% names(bin2_adapted) &&
                                b$label %in% bin2_adapted[[b$file]]$blocks)
    pk <- switch(b$tier, core = "frmtmb", sample = "frmtmb.sample",
                 both = c("frmtmb", "frmtmb.sample"))
    if (ported && b$file %in% names(bin2_adapted) &&
        !is.null(bin2_adapted[[b$file]]$pkgs)) {
      pk <- intersect(pk, bin2_adapted[[b$file]]$pkgs)
    }
    for (j in seq_along(b$stmts)) {
      k <- n_expect(b$stmts[[j]])
      if (k == 0L) next
      txt <- gsub(" +", " ", gsub("[\r\n\t]+", " ", b$texts[j]))
      rows[[length(rows) + 1L]] <- data.frame(
        file = b$file, block = gsub("[\r\n\t]+", " ", b$label),
        line = b$lines[j], bin = bin, tier = b$tier, n = k,
        id = sprintf("%s:%d", topic_of(b$file), b$lines[j]),
        ported = ported, pkgs = paste(pk, collapse = " "),
        assertion = substr(txt, 1, 160), stringsAsFactors = FALSE)
    }
  }
}
led <- do.call(rbind, rows)
stopifnot(sum(led$n[led$bin == "2T"]) == 59L,
          sum(led$n[led$bin == "2P"]) == 270L, !anyDuplicated(led$id))

problems <- character()
out_v <- list()
led$target <- ""; led$outcome <- NA_character_; led$class <- ""
led$reason <- ""; led$held <- ""; led$other <- ""; led$message <- ""
for (i in seq_len(nrow(led))) {
  id <- led$id[i]
  if (!led$ported[i]) {
    s <- sv[sv$id == id, ]
    if (nrow(s) != 1L) {
      problems <- c(problems, paste("NOT IN THE SURVEY", id))
      next
    }
    g <- groups[groups$group == s$group, ]
    led$outcome[i] <- "cannot transfer"
    # a survey row may carry its own class (a block that differs from
    # its group); otherwise the group's
    own <- !is.null(s$class) && nzchar(s$class)
    led$class[i] <- if (own) s$class else g$class
    led$reason[i] <- paste0(s$group, ": ", if (own) s$reason else g$reason)
    next
  }
  want <- strsplit(led$pkgs[i], " ", fixed = TRUE)[[1]]
  t <- tg$target[tg$id == id]
  if (!length(t)) t <- if (length(want) == 1L) want else NA
  if (is.na(t) || !t %in% want) {
    problems <- c(problems, paste("NO TARGET", id))
    next
  }
  led$target[i] <- t
  runs <- rec[rec$id == id, ]
  missing <- setdiff(want, runs$pkg)
  if (length(missing)) {
    problems <- c(problems, sprintf("NOT RECORDED %s in %s", id,
                                    paste(missing, collapse = ",")))
    next
  }
  runs <- runs[match(want, runs$pkg), ]
  led$held[i] <- paste(sprintf("%s=%s", ifelse(runs$pkg == "frmtmb",
                                                "fit", "draws"),
                               runs$held), collapse = " ")
  for (r in seq_len(nrow(runs))) {
    p <- runs$pkg[r]
    m <- man[man$id == id & man$pkg == p, ]
    if (runs$held[r] && nrow(m) && m$class %in% held_classes) {
      out_v[[length(out_v) + 1L]] <- data.frame(
        id = id, pkg = p, verdict = "hollow", reason = m$reason)
      v <- m$verdict; cls <- m$class; why <- m$reason
    } else if (runs$held[r]) {
      if (nrow(m)) {
        problems <- c(problems, sprintf(
          "HOLDS but has manual verdict '%s': %s in %s", m$verdict, id, p))
      }
      why <- unname(bin2_adapted_input[id])
      cls <- if (is.na(why)) "" else "adapted-input"
      if (is.na(why)) why <- ""
      out_v[[length(out_v) + 1L]] <- data.frame(
        id = id, pkg = p, verdict = "pass", reason = why)
      v <- "pass"
    } else {
      if (!nrow(m)) {
        problems <- c(problems, sprintf("UNCLASSIFIED %s in %s: %s", id, p,
                                        substr(runs$msg[r], 1, 150)))
        next
      }
      if (m$class %in% held_classes) {
        problems <- c(problems, sprintf(
          "marked hollow but does NOT hold: %s in %s", id, p))
      }
      out_v[[length(out_v) + 1L]] <- data.frame(
        id = id, pkg = p, verdict = m$verdict, reason = m$reason)
      v <- m$verdict; cls <- m$class; why <- m$reason
    }
    if (p == t) {
      led$outcome[i] <- v; led$class[i] <- cls; led$reason[i] <- why
      led$message[i] <- runs$msg[r]
    } else {
      led$other[i] <- sprintf("%s: %s%s", ifelse(p == "frmtmb", "fit",
                                                 "draws"), v,
                              if (nzchar(cls)) paste0(" (", cls, ")") else "")
    }
  }
}
# an adapted input whose row no longer holds, or is not ported, is a
# stale entry: the reason would describe a pass that is not there
held_ids <- rec$id[rec$held]
for (id in setdiff(names(bin2_adapted_input), held_ids)) {
  problems <- c(problems, paste("ADAPTED INPUT FOR A ROW THAT DOES NOT",
                                "HOLD:", id))
}
stale <- setdiff(paste(man$id, man$pkg), paste(rec$id, rec$pkg))
if (length(stale)) {
  problems <- c(problems, paste("MANUAL VERDICT FOR A RUN NOT RECORDED:",
                                stale))
}

# The value checks (lane bin2port, dev/brmsport-bin2-adapters.R): the
# value a line computed through the adapter, compared with brms's own
# function on the same arguments, after a brms assertion that checks
# only its shape. Each is recorded under the brms row's id and is judged
# beside that row, never in its verdict. On a row that passes on a
# shape, a check that holds needs no manual verdict and one that does
# not hold needs one (dev/brmsport-bin2-valuechecks-manual.tsv, written
# by dev/brmsport-bin2-verdicts-2p.R); the rules are the rows' own.
val <- rec_all[rec_all$kind == "value", ]
val$held <- val$held == "TRUE"
stopifnot(!anyDuplicated(val[c("pkg", "id")]))
val$outcome <- ifelse(val$held, "holds",
                      ifelse(grepl("^DIFFERS", val$msg), "differs", "error"))
vmf <- "dev/brmsport-bin2-valuechecks-manual.tsv"
vman <- if (file.exists(vmf)) rd(vmf) else
  data.frame(id = character(), pkg = character(), verdict = character(),
             class = character(), msg_has = character(),
             reason = character())
stopifnot(!anyDuplicated(vman[c("id", "pkg")]),
          all(vman$verdict %in% outcomes), all(nzchar(vman$reason)))
led$value_check <- ""; led$value_verdict <- ""; led$value_class <- ""
led$value_reason <- ""; led$value_msg <- ""
shape_pass <- led$outcome %in% "pass" &
  grepl(bin2_shape_re, led$assertion)
used_vman <- character()
for (i in which(shape_pass)) {
  r <- val[val$id == led$id[i] & val$pkg == led$target[i], ]
  if (!nrow(r)) next
  m <- vman[vman$id == led$id[i] & vman$pkg == led$target[i], ]
  led$value_check[i] <- r$outcome
  led$value_msg[i] <- r$msg
  if (r$held) {
    if (nrow(m)) {
      problems <- c(problems, sprintf(
        "VALUE CHECK HOLDS but has manual verdict '%s': %s", m$verdict,
        led$id[i]))
    }
    led$value_verdict[i] <- "pass"
  } else if (!nrow(m)) {
    problems <- c(problems, sprintf("VALUE CHECK UNCLASSIFIED %s: %s",
                                    led$id[i], substr(r$msg, 1, 150)))
  } else {
    # a manual verdict names the message it was written for; the same
    # label with another message is a stale verdict
    if (!nzchar(m$msg_has) || !grepl(m$msg_has, r$msg, fixed = TRUE)) {
      problems <- c(problems, sprintf(
        "VALUE CHECK '%s' no longer matches its verdict's message '%s': %s",
        led$id[i], m$msg_has, substr(r$msg, 1, 150)))
    }
    led$value_verdict[i] <- m$verdict
    led$value_class[i] <- m$class
    led$value_reason[i] <- m$reason
  }
  used_vman <- c(used_vman, paste(led$id[i], led$target[i]))
}
vstale <- setdiff(paste(vman$id, vman$pkg), used_vman)
if (length(vstale)) {
  problems <- c(problems, paste(
    "VALUE-CHECK VERDICT FOR NO CHECK ON A SHAPE PASS:", vstale))
}
if (length(problems)) {
  cat(problems, sep = "\n")
  stop(length(problems), " problem(s); the ledger was not written")
}
verd <- do.call(rbind, out_v)
utils::write.table(verd, out_path("brmsport-bin2-verdicts.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
# every recorded check's outcome, which the generator gives the tier to
# assert (a check on a row that is not a shape pass included)
val$msg_has <- vman$msg_has[match(paste(val$id, val$pkg),
                                  paste(vman$id, vman$pkg))]
val$msg_has[is.na(val$msg_has)] <- ""
utils::write.table(val[c("id", "pkg", "outcome", "msg_has", "msg")],
                   out_path("brmsport-bin2-valuechecks.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)

# A statement that holds several assertions (a loop) is run whole by the
# harness, so its record is one outcome. Where its assertions differ,
# dev/brmsport-bin2-split.tsv (written by dev/brmsport-bin2-verdicts-2p.R
# from a probe that judges each assertion line over every iteration)
# gives each its own verdict, class "per-line". The statement's recorded
# run must NOT hold (else the split is stale), and the parts must carry
# its whole assertion count.
splf <- "dev/brmsport-bin2-split.tsv"
if (file.exists(splf)) {
  sp <- rd(splf)
  sp$n <- as.integer(sp$n)
  for (sid in unique(sp$id)) {
    i <- which(led$id == sid)
    p <- sp[sp$id == sid, ]
    stopifnot(length(i) == 1L, sum(p$n) == led$n[i],
              !identical(led$outcome[i], "pass"))
    new <- led[rep(i, nrow(p)), ]
    new$id <- paste0(sid, "/", p$part)
    new$n <- p$n
    new$outcome <- p$verdict
    new$class <- p$class
    new$reason <- p$reason
    new$assertion <- p$assertion
    led <- rbind(led[seq_len(i - 1L), ], new,
                 led[seq_len(nrow(led))[-seq_len(i)], ])
  }
  stopifnot(sum(led$n) == 329L)
}

# What a pass asserts. brms wrote most of these blocks as "runs without
# errors": a check of a shape, a length, names or a type, which says the
# call answered and not that the answer is right. "value" is an
# assertion that compares a computed number with a reference or with
# another computation. Read from brms's assertion text, with the
# pattern dev/brmsport-bin2-adapters-spec.R shares with the generator.
led$kind <- ifelse(led$outcome %in% "pass",
                   ifelse(grepl(bin2_shape_re, led$assertion), "shape",
                          "value"),
                   "")
# The classes of "cannot transfer", split by what they mean. A decision
# is a question frmtmb does not answer on purpose or cannot be asked,
# each with its reason in the row: a Stan or MCMC object, a draws-only
# quantity of a fit, a brms internal with no counterpart (parsed
# formulas, link-array maps, insert_refcat()), no multivariate density
# outside a fit and no quantile function (SPEC.md:153-154 puts numeric
# d/p/q/r functions in gamlss.dist or extraDistr), a comparison of
# brms's own two backends. A gap is a function, family or CDF frmtmb
# lacks. "cost" is a row frmtmb could answer that this port did not
# reach. An artifact is the port's own limit (a fixture term, a
# reference that cannot run here, a row that tests the adapter, an
# upstream change that breaks brms's own assertion). The rate below
# excludes the decisions only.
ct_kind <- c(stan = "decision", mcmc = "decision", "no-draws" = "decision",
             "brms-internal" = "decision", "no-density" = "decision",
             "no-quantile" = "decision", "brms-backends" = "decision",
             absent = "gap", "no-cdf" = "gap", "not-ported" = "cost",
             fixture = "artifact", upstream = "artifact",
             reference = "artifact", adapter = "artifact")
led$ct_kind <- ""
ct <- led$outcome %in% "cannot transfer"
stopifnot(all(led$class[ct] %in% names(ct_kind)))
led$ct_kind[ct] <- ct_kind[led$class[ct]]
ledger <- led[c("file", "block", "line", "bin", "tier", "n", "id", "ported",
                "target", "outcome", "class", "ct_kind", "kind",
                "value_check", "value_verdict", "value_class", "reason",
                "held", "other", "assertion", "message")]
utils::write.table(ledger, out_path("brmsport-bin2-ledger.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)

## the generated summary --------------------------------------------------
L <- ledger
L$outcome <- factor(L$outcome, levels = outcomes)
md <- c("<!-- generated by dev/brmsport-bin2-ledger.R; do not edit by hand -->",
        "")
sumn <- function(d, ...) sum(d$n[Reduce(`&`, list(...))])
md <- c(md, "### Bin 2 per file and bin", "",
        paste("| file (bin) | assertions | pass (value / shape) | defect |",
              "divergence | cannot transfer: decision / gap / cost /",
              "artifact |"),
        "|---|---|---|---|---|---|")
keys <- unique(L[c("file", "bin")])
keys <- keys[order(keys$file, keys$bin), ]
row_of <- function(d, label) {
  sprintf("| %s | %d | %d (%d / %d) | %d | %d | %d / %d / %d / %d |",
          label, sum(d$n),
          sumn(d, d$outcome == "pass"),
          sumn(d, d$outcome == "pass", d$kind == "value"),
          sumn(d, d$outcome == "pass", d$kind == "shape"),
          sumn(d, d$outcome == "defect"), sumn(d, d$outcome == "divergence"),
          sumn(d, d$ct_kind == "decision"), sumn(d, d$ct_kind == "gap"),
          sumn(d, d$ct_kind == "cost"), sumn(d, d$ct_kind == "artifact"))
}
for (r in seq_len(nrow(keys))) {
  d <- L[L$file == keys$file[r] & L$bin == keys$bin[r], ]
  md <- c(md, row_of(d, sprintf("`%s` (%s)", keys$file[r], keys$bin[r])))
}
md <- c(md, sub("^[|] ([^|]*) [|]", "| **total** |", row_of(L, "total")))
md <- c(md, "", "### Bin 2 outcome by class", "",
        "| bin | outcome | class | kind | assertions |",
        "|---|---|---|---|---|")
for (bn in c("2T", "2P")) {
  d <- L[L$bin == bn, ]
  key <- unique(d[c("outcome", "class", "ct_kind")])
  key <- key[order(match(key$outcome, outcomes), key$class), ]
  for (r in seq_len(nrow(key))) {
    s <- sum(d$n[d$outcome == key$outcome[r] & d$class == key$class[r]])
    md <- c(md, sprintf("| %s | %s | %s | %s | %d |", bn, key$outcome[r],
                        if (nzchar(key$class[r])) key$class[r] else "-",
                        if (nzchar(key$ct_kind[r])) key$ct_kind[r] else "-",
                        s))
  }
}
p2t <- sum(L$n[L$bin == "2T" & L$outcome == "pass"])
p2p <- sum(L$n[L$bin == "2P" & L$outcome == "pass"])
ported2p <- sum(L$n[L$bin == "2P" & L$ported])
# a guard on the guard: this sentence once read a column the ledger did
# not carry and printed 0 ported
stopifnot(ported2p == sum(led$n[led$bin == "2P" & led$ported]))
# the 2T rows on the ML fit: the restatement, not counted
fit2t <- L[L$bin == "2T" & L$target == "frmtmb.sample", ]
fit_v <- sub(" [(].*$", "", sub("^fit: ", "", fit2t$other))
fit_v[grepl("fit=TRUE", fit2t$held)] <- "pass"
pv <- sum(L$n[L$outcome %in% "pass" & L$kind == "value"])
pa <- sum(L$n[L$outcome %in% "pass" & L$class == "adapted-input"])
md <- c(md, "",
        sprintf(paste("Bin 2T passes: %d of 59 (%.1f%%), judged on the",
                      "target object."), p2t, 100 * p2t / 59),
        sprintf(paste("Bin 2P: %d of 270 assertions run through the",
                      "adapter, %d of them pass (the rest it refused or",
                      "they failed); the other %d are classified from the",
                      "survey."), ported2p, p2p, 270L - ported2p),
        sprintf(paste("Of the %d bin-2 passes, %d compare a value and %d",
                      "check only a shape, length, names or type."),
                p2t + p2p, pv, p2t + p2p - pv),
        sprintf(paste("%d of the bin-2 passes hold on an adapted input",
                      "(class adapted-input: the port replaced brms's",
                      "input, which frmtmb refuses by design); %d are",
                      "literal."), pa, p2t + p2p - pa),
        sprintf(paste("The 58 methods rows of 2T on the maximum-likelihood",
                      "fit (the restatement, not counted): %s."),
                paste(sprintf("%s %d", names(table(fit_v)), table(fit_v)),
                      collapse = ", ")))

# the value checks beside the shape passes (lane bin2port)
sp <- L[L$outcome %in% "pass" & L$kind == "shape", ]
vrow <- function(d, label) {
  sprintf("| %s | %d | %d | %d | %d | %d | %d |", label, sum(d$n),
          sum(d$n[nzchar(d$value_check)]),
          sum(d$n[d$value_verdict == "pass"]),
          sum(d$n[d$value_verdict == "defect"]),
          sum(d$n[d$value_verdict == "divergence"]),
          sum(d$n[d$value_verdict == "cannot transfer"]))
}
md <- c(md, "", "### Value checks beside the shape passes", "",
        paste("A value check compares the value the line before a shape",
              "assertion computed through the adapter with brms's own",
              "function on the same arguments (brms:::log_lik_<family>(),",
              "or brms's d or p function). It is recorded beside the brms",
              "row and does not change that row's verdict. Simulators are",
              "not checked: two of them agree in distribution, not draw by",
              "draw."), "",
        paste("| file | shape passes | with a value check | holds |",
              "defect | divergence | reference cannot run |"),
        "|---|---|---|---|---|---|---|")
for (f in unique(sp$file)) {
  md <- c(md, vrow(sp[sp$file == f, ], sprintf("`%s`", f)))
}
md <- c(md, sub("^[|] ([^|]*) [|]", "| **total** |", vrow(sp, "total")))
vh <- sum(sp$n[sp$value_verdict == "pass"])
vbad <- sp[nzchar(sp$value_verdict) & sp$value_verdict != "pass", ]
md <- c(md, "",
        sprintf(paste("Of the %d bin-2 passes, %d compare a value as brms",
                      "wrote them, %d more check a shape and carry a value",
                      "check that holds against brms's own function, and %d",
                      "check only a shape."),
                p2t + p2p, pv, vh, p2t + p2p - pv - vh),
        if (nrow(vbad)) {
          paste0("Value checks that do not hold: ",
                 paste(sprintf("`%s` (%s, %s)", vbad$id, vbad$value_verdict,
                               vbad$value_class), collapse = "; "), ".")
        })

# compose with bin 1, whose cannot-transfer classes map the same way
b1 <- utils::read.delim("dev/brmsport-ledger.tsv", quote = "",
                        stringsAsFactors = FALSE)
stopifnot(nrow(b1) == 494L)
b1$n <- 1L
b1$ct_kind <- ""
c1 <- b1$outcome == "cannot transfer"
stopifnot(all(b1$class[c1] %in% names(ct_kind)))
b1$ct_kind[c1] <- ct_kind[b1$class[c1]]
cell <- function(d, o, k = NULL) {
  if (is.null(k)) sum(d$n[d$outcome %in% o]) else
    sum(d$n[d$outcome %in% o & d$ct_kind == k])
}
parts <- list(b1, L[L$bin == "2T", ], L[L$bin == "2P", ])
rows <- list(c("pass", "pass", NA), c("defect", "defect", NA),
             c("divergence", "divergence", NA),
             c("cannot transfer: no counterpart (decision)",
               "cannot transfer", "decision"),
             c("cannot transfer: absent (gap)", "cannot transfer", "gap"),
             c("cannot transfer: not ported (cost)", "cannot transfer",
               "cost"),
             c("cannot transfer: port artifact", "cannot transfer",
               "artifact"))
md <- c(md, "", "### Bins 1 and 2 together", "",
        "| outcome | bin 1 | bin 2T | bin 2P | bins 1 and 2 |",
        "|---|---|---|---|---|")
tot <- list()
for (r in rows) {
  k <- if (is.na(r[3])) NULL else r[3]
  v <- vapply(parts, cell, 0L, o = r[2], k = k)
  tot[[r[1]]] <- sum(v)
  md <- c(md, sprintf("| %s | %d | %d | %d | %d |", r[1], v[1], v[2], v[3],
                      sum(v)))
}
md <- c(md, sprintf("| **total** | **%d** | **%d** | **%d** | **%d** |",
                    nrow(b1), 59L, 270L, nrow(b1) + 329L))
np <- tot[["pass"]]
nd <- tot[["cannot transfer: no counterpart (decision)"]]
ng <- tot[["cannot transfer: absent (gap)"]]
nc <- tot[["cannot transfer: not ported (cost)"]]
na <- tot[["cannot transfer: port artifact"]]
stopifnot(sum(unlist(tot)) == 823L)
# the classes the rate leaves out, read from the rows, so the sentence
# cannot name a class that holds none
all_ct <- rbind(b1[c("class", "ct_kind", "n")],
                as.data.frame(L)[c("class", "ct_kind", "n")])
ex <- tapply(all_ct$n[all_ct$ct_kind == "decision"],
             all_ct$class[all_ct$ct_kind == "decision"], sum)
md <- c(md, "",
        sprintf(paste("Bins 1 and 2 pass: %d of 823 (%.1f%%); bin 1 alone",
                      "%d of 494, as dev/brmsport-ledger.tsv records it",
                      "(regenerated once at consolidation)."),
                np, 100 * np / 823, cell(b1, "pass")),
        sprintf(paste("Over every assertion but the %d that frmtmb has no",
                      "counterpart for by decision, %d of %d pass",
                      "(%.1f%%). The %d gaps, %d rows not ported (cost) and",
                      "%d port artifacts stay in that denominator."),
                nd, np, 823L - nd, 100 * np / (823 - nd), ng, nc, na),
        sprintf(paste("The rate's denominator excludes exactly the",
                      "\"cannot transfer\" rows of these classes: %s."),
                paste(sprintf("%s %d", names(ex), ex), collapse = ", ")))
writeLines(md, out_path("brmsport-bin2-summary.md"))
cat(md, sep = "\n")
