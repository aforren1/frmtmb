# The survey's rates, generated from the per-pair results. Nothing in
# the summary is typed by hand.
#
#   Rscript dev/equiv-summary.R [outdir] [--record]
#
# Reads dev/equiv-data/sample.tsv, strata.tsv, weights.tsv,
# weights-drawn.tsv, corpus.tsv, nl-sweep.tsv and adjudication.tsv, and
# every <outdir>/<id>.tsv. Writes results.tsv, compat-proposals.tsv and
# summary.md into <outdir>. Only with --record does it write them into
# dev/equiv-data/ instead and refresh the generated blocks of
# dev/equiv-findings.md: a rerun into a scratch outdir leaves the record
# alone (review 2026-10-28, finding 6).
args <- commandArgs(trailingOnly = TRUE)
record <- "--record" %in% args
args <- setdiff(args, "--record")
outdir <- if (length(args) >= 1) args[1] else "dev/equiv-out"
dest <- if (record) "dev/equiv-data" else outdir
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
rd <- function(p) read.delim(p, stringsAsFactors = FALSE,
                             na.strings = c("NA"))
smp <- rd("dev/equiv-data/sample.tsv")
fr <- rd("dev/equiv-data/strata.tsv")
w_now <- rd("dev/equiv-data/weights.tsv")
w_drawn <- rd("dev/equiv-data/weights-drawn.tsv")
corpus <- rd("dev/equiv-data/corpus.tsv")

rows <- lapply(seq_len(nrow(smp)), function(i) {
  p <- file.path(outdir, paste0(smp$id[i], ".tsv"))
  if (file.exists(p)) return(rd(p))
  lg <- file.path(outdir, "logs", paste0(smp$id[i], ".log"))
  tail_txt <- if (file.exists(lg)) {
    paste(utils::tail(readLines(lg, warn = FALSE), 3), collapse = " | ")
  } else "no log"
  data.frame(smp[i, c("id", "feature_a", "kind_a", "feature_b", "kind_b",
                      "status", "draws", "weight")],
             verdict = "cannot measure", level = "",
             agree = NA, reason = paste("no result:", substr(tail_txt, 1,
                                                            200)),
             stringsAsFactors = FALSE)
})
allc <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(r) {
  for (m in setdiff(allc, names(r))) r[[m]] <- NA
  r[, allc]
})
res <- do.call(rbind, rows)
res$verdict_measured <- res$verdict

adj_p <- "dev/equiv-data/adjudication.tsv"
res$adjudicated <- ""
if (file.exists(adj_p)) {
  adj <- rd(adj_p)
  for (k in seq_len(nrow(adj))) {
    i <- match(adj$id[k], res$id)
    if (is.na(i)) stop("adjudication names an unsampled id ", adj$id[k])
    res$verdict[i] <- adj$verdict[k]
    res$agree[i] <- as.logical(adj$agree[k])
    res$reason[i] <- adj$reason[k]
    if (!is.na(adj$proposal[k]) && nzchar(adj$proposal[k])) {
      res$proposal[i] <- adj$proposal[k]
    }
    if (!is.na(adj$level[k]) && nzchar(adj$level[k])) {
      res$level[i] <- adj$level[k]
    }
    res$adjudicated[i] <- "yes"
  }
}
verdicts <- c("equivalent", "refusal correct", "compat row wrong",
              "divergence by decision", "defect", "cannot measure")
bad <- setdiff(res$verdict, verdicts)
if (length(bad)) stop("unknown verdict: ", paste(bad, collapse = ", "))
res$agree[is.na(res$agree)] <- res$verdict[is.na(res$agree)] %in%
  c("equivalent", "refusal correct")
res$agree <- as.logical(res$agree)
res$nodefect <- res$agree | res$verdict == "divergence by decision"
res$measured <- res$verdict != "cannot measure"

lv <- res$level
lv[is.na(lv)] <- ""
top <- ifelse(grepl("L1", lv) & grepl("L2", lv), "L1+L2",
         ifelse(grepl("L1", lv), "L1",
           ifelse(grepl("L2", lv), "L2",
             ifelse(grepl("L3", lv), "L3",
               ifelse(grepl("L4", lv), "L4", "none")))))
res$top_level <- top
# Weak evidence (review 2026-10-28, finding 4): an agreement read only
# off the flat-prior posterior's fitted means (the route that scored a
# known identity at 0.851 sd), or a REML fit, which no brms quantity
# can check beyond the model itself.
route <- if (is.null(res$l4_route)) rep("", nrow(res)) else res$l4_route
route[is.na(route)] <- ""
res$weak_evidence <- res$agree & (
  (top == "L4" & route == "posterior") |
    res$feature_a == "REML" | res$feature_b == "REML")

wilson <- function(x, n, z = stats::qnorm(0.975)) {
  if (n == 0) return(c(NA, NA, NA))
  p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(p, ctr - half, ctr + half)
}
pct <- function(v) sprintf("%.1f", 100 * v)
strata <- c("untested", "conditional", "works", "refused")

stratum_tab <- function(ok) {
  do.call(rbind, lapply(strata, function(s) {
    i <- res$status == s
    n <- sum(res$draws[i])
    a <- sum(res$draws[i & ok])
    w <- wilson(a, n)
    data.frame(stratum = s, draws = n, x = a, rate = w[1], lo = w[2],
               hi = w[3])
  }))
}
W <- fr$W_share[match(strata, fr$status)]
Nsh <- fr$N_share[match(strata, fr$status)]

# The population estimate is stratified: each stratum's rate weighted
# by its share of the reach weight (W) of the whole pair table. Beside
# the Wald interval, a parametric bootstrap at the point estimate and a
# nonparametric stratified bootstrap that resamples the draws within
# each stratum (both 10,000 replicates, fixed seeds).
pop <- function(ok, wts = W) {
  st <- stratum_tab(ok)
  p <- sum(wts * st$rate)
  v <- sum(wts^2 * st$rate * (1 - st$rate) / st$draws)
  set.seed(20261028)
  pb <- replicate(10000, sum(wts * stats::rbinom(4, st$draws, st$rate) /
                               st$draws))
  set.seed(20261029)
  dr <- lapply(strata, function(s) {
    i <- res$status == s
    rep(ok[i], res$draws[i])
  })
  nb <- replicate(10000, sum(wts * vapply(dr, function(v) {
    mean(sample(v, length(v), replace = TRUE))
  }, 0)))
  q <- function(b) unname(stats::quantile(b, c(0.025, 0.975)))
  c(est = p, lo = p - 1.96 * sqrt(v), hi = p + 1.96 * sqrt(v),
    q(pb), q(nb))
}
p_agree <- pop(res$agree)
p_nd <- pop(res$nodefect)
p_weak <- pop(res$agree & !res$weak_evidence)
st_a <- stratum_tab(res$agree)
st_nd <- stratum_tab(res$nodefect)
st_w <- stratum_tab(res$agree & !res$weak_evidence)
# measured only: each stratum over its measured draws
st_m <- do.call(rbind, lapply(strata, function(s) {
  i <- res$status == s & res$measured
  w <- wilson(sum(res$draws[i & res$agree]), sum(res$draws[i]))
  data.frame(stratum = s, draws = sum(res$draws[i]), rate = w[1],
             lo = w[2], hi = w[3])
}))
p_m <- sum(W * st_m$rate)
p_count <- sum(Nsh * st_a$rate)

# The bar_crossing pattern was corrected after the draw (weights.tsv
# against weights-drawn.tsv). Its effect, by reweighting each drawn pair
# by its new over its old weight within the strata and recomputing the
# strata shares on the pair table (needs dev/equiv-frame-out/compat.tsv).
wd <- stats::setNames(w_drawn$weight, w_drawn$feature)
wn <- stats::setNames(w_now$weight, w_now$feature)
reweight_line <- NULL
if (!identical(wd[order(names(wd))], wn[order(names(wn))])) {
  ratio <- (wn[res$feature_a] * wn[res$feature_b]) /
    (wd[res$feature_a] * wd[res$feature_b])
  rr <- vapply(strata, function(s) {
    i <- res$status == s
    sum((res$draws * ratio * res$agree)[i]) / sum((res$draws * ratio)[i])
  }, 0)
  fp_frame <- "dev/equiv-frame-out/compat.tsv"
  W2 <- W
  if (file.exists(fp_frame)) {
    cf <- rd(fp_frame)
    cf <- cf[cf$in_core_only, ]
    pw <- wn[cf$feature_a] * wn[cf$feature_b]
    W2 <- vapply(strata, function(s) sum(pw[cf$status == s]), 0)
    W2 <- W2 / sum(W2)
  }
  changed <- names(wn)[wn != wd[names(wn)]]
  reweight_line <- sprintf(paste0(
    "- agree, with the corrected reach weights of %s (the sample stays ",
    "the one drawn; each draw reweighted by new over old weight%s): %s"),
    paste(changed, collapse = ", "),
    if (file.exists(fp_frame)) ", strata shares recomputed" else "",
    pct(sum(W2 * rr)))
}

utils::write.table(res, file.path(dest, "results.tsv"), sep = "\t",
                   row.names = FALSE, quote = TRUE, qmethod = "double")

# Proposals: the sampled pairs', then the nl x family sweep's rows that
# no sampled pair already covers.
prop <- res[!is.na(res$proposal) & nzchar(res$proposal),
            c("id", "feature_a", "feature_b", "status", "proposal",
              "verdict", "reason")]
names(prop)[names(prop) == "status"] <- "status_now"
prop$source <- "sampled pair"
sw_p <- "dev/equiv-data/nl-sweep.tsv"
if (file.exists(sw_p)) {
  sw <- rd(sw_p)
  sw <- sw[!is.na(sw$proposal) & nzchar(sw$proposal), ]
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  have <- key(prop$feature_a, prop$feature_b)
  sw <- sw[!(key("nl", sw$family) %in% have), ]
  if (nrow(sw)) {
    prop <- rbind(prop, data.frame(
      id = "", feature_a = "nl", feature_b = sw$family,
      status_now = sw$status_now, proposal = sw$proposal,
      verdict = ifelse(sw$outcome == "nl body DROPPED", "defect (D1)",
                       "compat row wrong"),
      # the sweep compares frm() maxima of a LINEAR body with the
      # linear fit; it is not an L1 identity with brms
      reason = paste0(sprintf(paste0(
        "nl x family sweep (dev/equiv-nl-sweep.R): %s; logLik linear %s, ",
        "nl %s, intercept-only %s; brms nl spelling %s. A frm() maximum ",
        "of a linear body, not an L1 identity with brms"),
        sw$outcome, sw$frm_linear, sw$frm_nl, sw$frm_intercept,
        substr(sw$brms_nl, 1, 40)),
        ifelse(nzchar(sw$default_outcome %||% ""), paste0(
          "; under the log link. With the default link: nl from frm()'s ",
          "default starts fails (", sw$default_nl_default_start_err,
          "); from starts that keep mu positive, logLik ",
          sw$default_nl_positive_start, " against the linear ",
          sw$default_linear, ". The default starts are a defect candidate"),
          "")),
      source = "nl x family sweep"))
  }
}
utils::write.table(prop, file.path(dest, "compat-proposals.tsv"),
                   sep = "\t", row.names = FALSE, quote = TRUE,
                   qmethod = "double")

out <- character()
say <- function(...) out <<- c(out, paste0(...))
say("<!-- generated by dev/equiv-summary.R from ", outdir,
    "; do not edit -->")
say("")
say("frmtmb measured: ", paste(unique(stats::na.omit(res$frmtmb[
  nzchar(res$frmtmb)])), collapse = ", "), "; ", nrow(res),
  " unique pairs, ", sum(res$draws), " draws; ",
  sum(res$adjudicated == "yes"), " verdicts adjudicated by hand ",
  "(dev/equiv-data/adjudication.tsv).")
say("")
say("Verdicts by stratum (draws, with unique pairs in brackets):")
say("")
say("| stratum | ", paste(verdicts, collapse = " | "), " |")
say("|---|", paste(rep("---", length(verdicts)), collapse = "|"), "|")
for (s in strata) {
  r <- res[res$status == s, ]
  say("| ", s, " | ", paste(vapply(verdicts, function(v) {
    sprintf("%d (%d)", sum(r$draws[r$verdict == v]), sum(r$verdict == v))
  }, ""), collapse = " | "), " |")
}
say("")
say("Rates per stratum, percent, Wilson 95% interval, over draws:")
say("")
say("| stratum | draws | agree with brms | no defect | ",
    "agree, measured only | agree, weak evidence out |")
say("|---|---|---|---|---|---|")
ci <- function(t, i) sprintf("%s [%s, %s]", pct(t$rate[i]), pct(t$lo[i]),
                             pct(t$hi[i]))
for (i in seq_along(strata)) {
  say("| ", strata[i], " | ", st_a$draws[i], " | ", ci(st_a, i), " | ",
      ci(st_nd, i), " | ", ci(st_m, i), " (", st_m$draws[i], ") | ",
      ci(st_w, i), " |")
}
say("")
say("Population estimate, weighted by reach (stratum shares of the ",
    "reach weight: ", paste(sprintf("%s %.3f", strata, W),
                            collapse = ", "), "). Wald 95% interval, then ",
    "the parametric and the nonparametric stratified bootstrap:")
say("")
f3 <- function(p) sprintf(paste0("%s, Wald [%s, %s]; bootstrap parametric ",
                                 "[%s, %s], nonparametric [%s, %s]"),
                          pct(p[1]), pct(p[2]), pct(p[3]), pct(p[4]),
                          pct(p[5]), pct(p[6]), pct(p[7]))
say("- **agree with brms (the headline): ", f3(p_agree), "**")
say("- no defect (agree, or a divergence by decision): ", f3(p_nd))
say("- sensitivity, weak evidence out: ", f3(p_weak), ". The ",
    sum(res$weak_evidence), " pairs (", sum(res$draws[res$weak_evidence]),
    " draws) scored agree only by L4 fitted means against a flat-prior ",
    "posterior, or under REML, count as not agree: ",
    paste(res$id[res$weak_evidence], collapse = ", "), ".")
say("- agree, among measured pairs only: ", pct(p_m))
say("- agree, strata weighted by pair COUNT instead of reach: ",
    pct(p_count))
if (!is.null(reweight_line)) say(reweight_line)
say("")
say("'agree' is equivalent, refusal correct, or a compat row whose ",
    "behavior agrees with brms. Cannot-measure pairs stay in the ",
    "denominator except in the 'measured only' column.")
say("")
say("Strongest level reached, unique pairs (draws):")
say("")
tl <- table(factor(top, c("L1+L2", "L1", "L2", "L3", "L4", "none")))
for (k in names(tl)) {
  say("- ", k, ": ", tl[[k]], " (", sum(res$draws[top == k]), ")")
}
l4r <- table(route[top == "L4"])
say("- of the L4-only pairs, by route: ",
    paste(sprintf("%s %d", names(l4r), as.integer(l4r)), collapse = ", "))
say("")
say("Compat-row proposals: ", nrow(prop), " (", sum(prop$source ==
  "sampled pair"), " from sampled pairs, ", sum(prop$source ==
  "nl x family sweep"), " from the nl x family sweep; ",
  file.path(dest, "compat-proposals.tsv"), ").")
if (file.exists(sw_p)) {
  sw_all <- rd(sw_p)
  tt <- table(sw_all$outcome[sw_all$status_now == "refused"])
  say("nl x family sweep, rows the table says refused: ",
      paste(sprintf("%s %d", names(tt), as.integer(tt)), collapse = ", "),
      "; of the honored ones, brms builds the nl spelling for ",
      sum(sw_all$status_now == "refused" & sw_all$outcome == "nl honored" &
            grepl("^builds", sw_all$brms_nl)), " (", sw_p, ").")
}
say("Pairs touching other lanes' areas: ",
    paste(sprintf("%s %s x %s [%s]", res$id, res$feature_a, res$feature_b,
                  res$overlap)[!is.na(res$overlap) & nzchar(res$overlap)],
          collapse = "; "))
noan <- grepl("no brms analogue", paste(res$reason, res$brms_error))
say("")
say("Pairs with no brms spelling (frmtmb-only structures and terms): ",
    sum(noan), " unique (", sum(res$draws[noan]), " draws): ",
    paste(sprintf("%s %s x %s", res$id, res$feature_a,
                  res$feature_b)[noan], collapse = "; "), ".")
sec <- suppressWarnings(as.numeric(res$seconds))
if (any(is.finite(sec))) {
  k <- which.max(sec)
  say("")
  say(sprintf(paste0("Cost: %d pairs timed; per-pair wall seconds sum ",
                     "%.0f, median %.1f, max %.0f (%s %s x %s); with 8 at ",
                     "once the run is bounded below by the longest pair."),
              sum(is.finite(sec)), sum(sec, na.rm = TRUE),
              stats::median(sec, na.rm = TRUE), sec[k], res$id[k],
              res$feature_a[k], res$feature_b[k]))
}
writeLines(out, file.path(dest, "summary.md"))
cat(out, sep = "\n")

# The frame of section 2, generated from the files it describes.
fo <- character()
fsay <- function(...) fo <<- c(fo, paste0(...))
fsay("<!-- generated by dev/equiv-summary.R; do not edit -->")
fsay("")
fp_frame <- "dev/equiv-frame-out/compat.tsv"
if (file.exists(fp_frame)) {
  cf <- rd(fp_frame)
  fsay("- pair table: ", sum(cf$in_core_only), " core pairs (",
       nrow(cf), " with all eight extensions loaded); core statuses ",
       paste(sprintf("%s %d", names(table(cf$status[cf$in_core_only])),
                     as.integer(table(cf$status[cf$in_core_only]))),
             collapse = ", "), ".")
} else {
  fsay("- pair table: ", sum(fr$N), " core pairs (strata.tsv; rerun ",
       "dev/equiv-dump-compat.R for the full table).")
}
fsay("- corpus: ", nrow(corpus), " documents (",
     paste(sprintf("%s %d", names(table(corpus$source)),
                   as.integer(table(corpus$source))), collapse = ", "),
     ", corpus.tsv).")
fsay("- features: ", nrow(w_now), " core features; ", sum(w_now$d == 0),
     " appear in no document (d = 0), ", sum(w_now$d_brms == 0),
     " in no brms document.")
top10 <- utils::head(w_now[order(-w_now$weight), ], 11)
fsay("- top weights (weights.tsv): ", paste(sprintf("%s %d", top10$feature,
                                                     top10$weight),
                                             collapse = ", "), ".")
fsay("")
fsay("| stratum | pairs N | reach share W | draws | unique pairs |")
fsay("|---|---|---|---|---|")
for (s in strata) {
  j <- match(s, fr$status)
  fsay(sprintf("| %s | %d | %.3f | %d | %d |", s, fr$N[j], fr$W_share[j],
               fr$draws[j], fr$unique[j]))
}
fsay("")
fsay("Sample: ", nrow(smp), " unique pairs, ", sum(smp$draws), " draws, ",
     length(unique(c(smp$feature_a, smp$feature_b))), " features.")

if (record) {
  # The findings carry each block verbatim between two markers, so the
  # record's numbers are this script's and are never retyped.
  fp <- "dev/equiv-findings.md"
  put <- function(tag, lines) {
    fx <- readLines(fp, warn = FALSE)
    b <- grep(paste0("^<!-- ", tag, " begin -->$"), fx)
    e <- grep(paste0("^<!-- ", tag, " end -->$"), fx)
    if (length(b) == 1L && length(e) == 1L && e > b) {
      writeLines(c(fx[seq_len(b)], lines, fx[e:length(fx)]), fp)
    }
  }
  esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
  put("equiv-summary", out)
  put("equiv-frame", fo)
  put("equiv-proposals", c(
    "<!-- generated by dev/equiv-summary.R; do not edit -->", "",
    "| id | pair | status now | proposed | verdict | source |",
    "|---|---|---|---|---|---|",
    sprintf("| %s | %s x %s | %s | %s | %s | %s |", prop$id,
            esc(prop$feature_a), esc(prop$feature_b), prop$status_now,
            prop$proposal, prop$verdict, prop$source)))
  cat("record updated: dev/equiv-data/ and the findings blocks\n")
}
