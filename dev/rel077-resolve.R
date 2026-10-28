# Consolidation edits of the 0.77.0 release tree, applied after the
# textual merge (dev/rel077-merge.sh) and before the bin-2 generator and
# the NEWS merge (dev/rel077-news-merge.R). Each edit checks that its old
# text is there exactly once (or the count it names), so a rerun on a
# fresh merge either applies all of them or stops on the first that no
# longer fits. A file no lane changed is restored from its HEAD blob
# (435e7ce6) before it is edited, so every edit here is applied once
# however often this runs. Every hand resolution of the round is here,
# each with its reason. Adapted from dev/rel076-resolve.R.
#
#   Rscript dev/rel077-resolve.R
rel <- Sys.getenv("REL077_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel077")
git_dir <- Sys.getenv("REL077_GIT", rel)
setwd(rel)
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)),
               stdout = TRUE)
  wr(sub("\r$", "", x), f)
}
edit_file <- function(f, old, new, n = 1L) {
  s <- paste(rd(f), collapse = "\n")
  hits <- gregexpr(old, s, fixed = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ",
                   substr(old, 1, 70))
  s <- gsub(old, new, s, fixed = TRUE)
  wr(strsplit(s, "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}
merged <- readLines("dev/rel077-log/merge-work/files")
copied <- unlist(lapply(list.files("dev/rel077-log", "^copied-.*[.]lst$",
                                   full.names = TRUE), readLines))
restore <- function(f) {
  if (f %in% c(merged, copied)) {
    stop(f, " is a merged or copied file; edit it, do not restore it")
  }
  from_head(f)
}
# Replace the one conflict block of `f` whose two sides contain `ours`
# and `theirs` (fixed strings) with `new` (a character vector of lines,
# or a function of the two sides' lines). Stops unless exactly one
# block matches.
resolve_block <- function(f, ours, theirs, new) {
  x <- rd(f)
  a <- which(startsWith(x, "<<<<<<< "))
  m <- which(x == "=======")
  z <- which(startsWith(x, ">>>>>>> "))
  stopifnot(length(a) == length(z))
  hit <- integer(0)
  for (i in seq_along(a)) {
    mid <- m[m > a[i] & m < z[i]]
    stopifnot(length(mid) == 1L)
    o <- paste(x[(a[i] + 1L):(mid - 1L)], collapse = "\n")
    t <- if (z[i] > mid + 1L) paste(x[(mid + 1L):(z[i] - 1L)],
                                    collapse = "\n") else ""
    if (grepl(ours, o, fixed = TRUE) && grepl(theirs, t, fixed = TRUE)) {
      hit <- c(hit, i)
    }
  }
  if (length(hit) != 1L) stop(f, ": expected 1 conflict block, found ",
                              length(hit), ": ", ours)
  i <- hit
  mid <- m[m > a[i] & m < z[i]]
  if (is.function(new)) {
    new <- new(x[seq_len(mid - a[i] - 1L) + a[i]],
               if (z[i] > mid + 1L) x[(mid + 1L):(z[i] - 1L)] else
                 character(0))
  }
  wr(c(x[seq_len(a[i] - 1L)], new, x[-seq_len(z[i])]), f)
  cat("resolved a conflict in", f, "\n")
}
both <- function(ours, theirs) c(ours, theirs)

# ---- 1. R/priors.R: mointerop x reslot, the simplex branch ------------
# Both lanes added a skip at the top of the loop over a simo or family
# simplex target: mointerop for a simplex of one weight (a two-value
# mo(), brms's simplex[1], whose Dirichlet is a constant), reslot for a
# specification with no density (a bounds-only or (flat) row, which a
# simplex reads as no statement). The two are independent; both stay,
# reslot's first since it does not read the target.
resolve_block("R/priors.R", "if (tm$D < 2L) next", "is.null(s$dist)) next",
              function(ours, theirs) c(theirs, ours))
# The overlap the reviews named, reslot's hunk at base line 3839 (the
# ordinal Intercept bound message) next to sigmainn's removal of the
# innov block in the loop over target_coefs(s), merged clean; the loop
# body is sigmainn's (checked: the merged loop reads no `innov`, and
# R/priors.R has no `innov_logc`, `lc =` or `$innov`).
pr <- rd("R/priors.R")
stopifnot(!any(grepl("innov_logc|[$]innov\\b|lc = ", pr)))

# ---- 2. dev/brmsport-bin2-adapters-spec.R: sigmainn x mointerop -------
# sigmainn added the ARMA blocks (gaussian_time, student_time) and
# mointerop the multivariate ones (gaussian_mv, student_mv) at the same
# anchors: the header's list of what goes through an adapter, the
# function vectors, both block lists and both substitution lists. Every
# block keeps both. The tier files are regenerated from this spec by
# dev/brmsport-bin2-gen.R (dev/rel077-merge.sh), never merged by hand.
sf <- "dev/brmsport-bin2-adapters-spec.R"
resolve_block(sf, "The ARMA blocks", "The multivariate", function(o, t) {
  n <- length(o)
  stopifnot(n == 12L, length(t) == 7L,
            o[n - 1L] == paste("# the value checks read brms's ARMA",
                               "covariance. Left to the survey: the"),
            o[n] == "# multivariate blocks and custom families.",
            startsWith(t[6], "# joint draw of predict(). Left to the"))
  c(o[seq_len(n - 2L)],
    "# the value checks read brms's ARMA covariance. The multivariate",
    t[2:5],
    "# joint draw of predict(). Left to the survey: custom families.")
})
resolve_block(sf, "\"gaussian_time\", \"student_time\")",
              "brms_adapt_mv_spec_funs <-", function(o, t) {
  c(o, t[-1L])
})
resolve_block(sf, "\"log_lik for ARMA models runs without errors\"),",
              "\"log_lik for multivariate linear models runs without",
              c("      # the covariance form of a first-order ARMA",
                "      \"log_lik for ARMA models runs without errors\",",
                paste0("      # lane mointerop: multivariate gaussian and ",
                       "student-t with rescor"),
                paste0("      \"log_lik for multivariate linear models ",
                       "runs without errors\"),")))
# the posterior_predict block first: its release side also holds
# "brms_adapt_ac_spec_funs)),", the log_lik block's anchor
resolve_block(sf, "posterior_predict for ARMA covariance models",
              "posterior_predict for multivariate linear models", c(
  "      # the covariance form of a first-order ARMA",
  paste0("      \"posterior_predict for ARMA covariance models runs ",
         "without errors\","),
  "      # lane mointerop: multivariate gaussian and student-t with rescor",
  "      paste(\"posterior_predict for multivariate linear models runs\",",
  "            \"without errors\")),",
  "    pre = brms_adapt_pre_links,",
  "    subst = c(brms_adapt_subst(\"posterior_predict\",",
  "                               c(brms_adapt_funs, brms_adapt_ac_spec_funs,",
  "                                 brms_adapt_mv_spec_funs)),",
  "              brms_adapt_sar_subst, brms_adapt_time_subst)),"))
resolve_block(sf, "brms_adapt_ac_spec_funs)),",
              "brms_adapt_mv_spec_funs)),", c(
  "                                            brms_adapt_ac_spec_funs,",
  "                                            brms_adapt_mv_spec_funs)),",
  "              brms_adapt_sar_subst, brms_adapt_time_subst)),"))
stopifnot(!any(grepl("^(<<<<<<<|=======$|>>>>>>>)", rd(sf))))

# ---- 3. frmtmb.sample R/sample.R: reslot x priordraws, ?frm_sample -----
# Both lanes rewrote the lines between `seed` and `prior`: priordraws the
# `seed` sentence (it fixes the prior draws too), `save_warmup` and the
# new @param sample_prior; reslot the @param prior text (brms's (flat)
# for an empty specification; the named list is gone). Each keeps its
# own part: priordraws' lines up to @param prior, reslot's from there.
sa <- "extensions/frmtmb.sample/R/sample.R"
resolve_block(sa, "#' @param prior Priors: a [frmtmb::set_prior()] specification, a",
              "#' @param sample_prior brms's `sample_prior`", function(o, t) {
  ko <- grep("^#' @param prior Priors", o)
  kt <- grep("^#' @param prior Priors", t)
  stopifnot(length(ko) == 1L, length(kt) == 1L)
  c(t[seq_len(kt - 1L)], o[ko:length(o)])
})

# ---- 4. the overlaps the reviews named, edited -------------------------
# sigmainn x reslot: flat_spec() cleared the `innov_sd` field that the
# 0.76.0 prior layer set; sigmainn removed that field, so the line is a
# dead store
edit_file(sa, "  s$prior <- \"\"\n  s$innov_sd <- NULL\n  p[[1L]] <- s",
          "  s$prior <- \"\"\n  p[[1L]] <- s")
# reslot x priordraws: nothing sets the `overrides` attribute since the
# named-list prior went (reslot), and priordraws rewrote the evidence
# ratio without er_entries(); the two checks left read a dead attribute
md <- "extensions/frmtmb.sample/R/methods-draws.R"
edit_file(md,
"  if (is.null(pl) || (!length(unclass(pl)) &&
                        !length(attr(pl, \"overrides\")))) {",
"  if (is.null(pl) || !length(unclass(pl))) {")
edit_file(sa,
"#' program writes it (`prior_plan()`). The named-list prior spelling
#' takes over internal parameters that no `set_prior()` row describes,
#' so brms's draws of the rows would not be the prior that was sampled
#' under; it is refused.
#'
#' @noRd
sample_prior_plan <- function(fit, rp, sample_prior) {
  eff <- rp$effective
  if (length(attr(eff, \"overrides\") %||% list())) {
    frm_stop(\"frm_sample(sample_prior = \\\"\", sample_prior, \"\\\") reads the \",
             \"prior as set_prior() rows, as brms does, and the named-list \",
             \"prior spelling takes over internal parameters no row \",
             \"describes. Write the prior with set_prior()\", call. = FALSE)
  }
  pl <-",
"#' program writes it (`prior_plan()`).
#'
#' @noRd
sample_prior_plan <- function(fit, rp, sample_prior) {
  eff <- rp$effective
  pl <-")
for (f in c("extensions/frmtmb.sample/R/evidence-ratio.R", md, sa)) {
  stopifnot(!any(grepl("overrides", rd(f), fixed = TRUE)))
}
# sigmainn x reslot: reslot's test asserted the `innov` field that
# sigmainn removed. (flat) on class sigma under ar(cov = TRUE) is flat on
# the innovation sd itself, brms's real<lower=0> sigma with no prior: no
# innov field, and a -log prior that is the log-Jacobian of exp alone,
# so the same at any ar
edit_file("extensions/frmtmb.sample/tests/testthat/test-flat-prior.R",
"  expect_identical(sg[[1L]]$dist$kind, \"flat\")
  expect_false(is.null(sg[[1L]]$innov))
})",
"  expect_identical(sg[[1L]]$dist$kind, \"flat\")
  # (flat) on the innovation sd itself, brms's real<lower=0> sigma with
  # no prior: the entry reads no other parameter, and its -log prior at
  # a fixed log sigma is minus that log sigma (the log-Jacobian of exp)
  # at any ar
  expect_null(sg[[1L]]$innov)
  nlp <- vapply(c(-0.8, 0.9), function(a) {
    fp_nlp(fit, sg, list(betad = 0.3, thetaac = a))
  }, 0)
  expect_identical(nlp[1], nlp[2])
  expect_lt(abs(nlp[1] + 0.3), 1e-12 * 0.3)
})")

# ---- 5. dev/upstream-bugs.md: brms's lb on class cor --------------------
# reslot found it (dev/reslot-brms-corlb.R); not posted. No lane changed
# the file; the numbers stay unique
ub <- "dev/upstream-bugs.md"
restore(ub)
edit_file(ub, "\n## RTMB\n", paste0("\n", paste(rd(
  "dev/rel077-upstream-brms42.md"), collapse = "\n"), "\n\n## RTMB\n"))
nums <- sub("^### (brms-[0-9]+)[.].*", "\\1",
            grep("^### brms-[0-9]+[.]", rd(ub), value = TRUE))
if (anyDuplicated(nums)) stop("dev/upstream-bugs.md: a number twice: ",
                              nums[duplicated(nums)][1])
cat("dev/upstream-bugs.md:", length(nums), "entries, last",
    nums[length(nums)], ", no number twice\n")

# ---- 6. versions and floors ---------------------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
dcore <- "DESCRIPTION"
restore(dcore)
ver(dcore, "0.76.0", "0.77.0")
cm <- "codemeta.json"
restore(cm)
# no lane changed a dependency of core, so only the version moves
edit_file(cm, '"version": "0.76.0",', '"version": "0.77.0",')
# frmtmb.sample: reads build_objective()'s prior_only (on 0.76.0 core,
# sample_prior = "only" samples the posterior: priordraws' review,
# finding 3), writes the (flat) density kind and the flat LKJ that the
# 0.77.0 resolver reads, and reads the route mark; minor, breaking
dsa <- "extensions/frmtmb.sample/DESCRIPTION"
restore(dsa)
ver(dsa, "0.23.0", "0.24.0")
edit_file(dsa, "frmtmb (>= 0.76.0)", "frmtmb (>= 0.77.0)")
# frmtmb.spatial and the other six: see dev/rel077-versions.R, which the
# measurement of their floors decides (dev/round-20261028.md)
vf <- "dev/rel077-versions.R"
if (file.exists(vf)) source(vf, local = TRUE)

# ---- 7. the release harness reads rellib-r15 -----------------------------
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R", "dev/brmsport-bin2-all.sh",
            "dev/brmsport-bin2-record.sh", "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r14", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("rellib-r14", "rellib-r15", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}

# ---- 8. the ported brms suite: the verdicts the round moves -------------
# Bin 1. Rows that now hold lose their manual verdict (the ledger stops
# on a manual verdict for a row that holds): brmsfit-methods:794, 797
# and 798 (priordraws: prior_draws()), priors:5 and priors:123 (reslot:
# one slot per parameter, no theta rows), brmsfit-methods:924
# (mointerop: mo() on one value refits, as brms's update() builds it).
# Rows 107 and 802 keep "cannot transfer" for another reason now: each
# reads fit1$fit, which is brms's stanfit and here the frmtmb fit.
vm <- "dev/brmsport-verdicts-manual.tsv"
vmf <- "dev/brmsport-verdicts-manual-fit.tsv"
restore(vm)
restore(vmf)
drop_rows <- function(f, ids) {
  x <- rd(f)
  hit <- vapply(ids, function(id) sum(startsWith(x, paste0(id, "\t"))), 0L)
  if (any(hit != 1L)) stop(f, ": not one row for ",
                           paste(ids[hit != 1L], collapse = ", "))
  wr(x[!sub("\t.*", "", x) %in% ids], f)
  cat("dropped", length(ids), "verdicts from", f, "\n")
}
drop_rows(vm, c("priors:5", "priors:123", "brmsfit-methods:924"))
drop_rows(vmf, c("brmsfit-methods:794", "brmsfit-methods:797",
                 "brmsfit-methods:798"))
set_row <- function(f, id, new) {
  x <- rd(f)
  i <- which(startsWith(x, paste0(id, "\t")))
  stopifnot(length(i) == 1L, !grepl("\t", new$reason))
  x[i] <- paste(id, new$pkg, new$verdict, new$class, new$reason, sep = "\t")
  wr(x, f)
  cat("set the verdict of", id, "in", f, "\n")
}
set_row(vmf, "brmsfit-methods:107", list(
  pkg = "frmtmb.sample", verdict = "cannot transfer", class = "stan",
  reason = paste(
    "reads fit1$fit@sim$iter: in brms fit1$fit is the stanfit, here it is",
    "the frmtmb fit, which has no @sim. The warmup itself is stored now:",
    "dim(as.mcmc(fit1, inc_warmup = TRUE)[[1]]) is 75 x the variable",
    "count, brms's shape for the fixture's 1 chain of 75 with 50 warmup",
    "(dev/priordraws-findings.md, section 7)")))
set_row(vmf, "brmsfit-methods:802", list(
  pkg = "frmtmb.sample", verdict = "cannot transfer", class = "stan",
  reason = paste(
    "prior_draws() on fit1$fit: in brms the stanfit, here the frmtmb fit,",
    "which carries no draws. prior_draws(fit1) itself answers with brms's",
    "names and values (brmsfit-methods:794, 797 and 798 pass;",
    "dev/priordraws-findings.md, section 7)")))

# Bin 2. The 18 rows that test brms's own numeric helpers (classes
# no-density and no-quantile, "cannot transfer / decision" until now)
# are divergences by policy (the coordinator's decision at this
# consolidation): six run through the adapter and take their verdict in
# dev/brmsport-bin2-verdicts-2p.R; twelve are classified from the survey,
# and dev/brmsport-bin2-ledger.R now reads a policy verdict for an
# unported row from the spec (bin2_policy_unported). The two rescor
# posterior_predict rows get a note on what their pass says.
v2p <- "dev/brmsport-bin2-verdicts-2p.R"
restore(v2p)
edit_file(v2p,
"stopifnot(!anyDuplicated(V$id), all(nzchar(V$reason)),
          !any(grepl(\"\\t\", V$reason)))",
"# the 0.77.0 consolidation: brms's own numeric helpers are a divergence
# by policy, not a decision of \"cannot transfer\" (bin2_numeric_helper_reason
# in dev/brmsport-bin2-adapters-spec.R, which also lists the unported
# rows of the multivariate block)
for (k in c(8, 68, 71, 97, 100, 139)) {
  V <- bp_set(V, paste0(\"distributions:\", k), \"divergence\", \"policy\",
              bin2_numeric_helper_reason)
}
stopifnot(!anyDuplicated(V$id), all(nzchar(V$reason)),
          !any(grepl(\"\\t\", V$reason)))")
edit_file(v2p, "v <- list()\nadd <- function(id, verdict, class, reason) {",
          "source(\"dev/brmsport-bin2-adapters-spec.R\")\nv <- list()\nadd <- function(id, verdict, class, reason) {")
anchor <- paste("# brms's exported d, r and p functions that",
                "dev/brmsport-bin2-adapters.R")
edit_file(sf, anchor, paste(c(rd("dev/rel077-bin2-policy-spec.R"), "",
                              anchor), collapse = "\n"))
l2 <- "dev/brmsport-bin2-ledger.R"
restore(l2)
edit_file(l2,
"    g <- groups[groups$group == s$group, ]
    led$outcome[i] <- \"cannot transfer\"",
"    g <- groups[groups$group == s$group, ]
    if (id %in% bin2_policy_unported) {
      # a divergence by policy on a row the port does not run (brms's
      # own numeric helpers, the 0.77.0 consolidation)
      led$outcome[i] <- \"divergence\"
      led$class[i] <- \"policy\"
      led$reason[i] <- paste0(s$group, \": \", bin2_numeric_helper_reason)
      next
    }
    led$outcome[i] <- \"cannot transfer\"")
edit_file(l2,
"      why <- unname(bin2_adapted_input[id])
      cls <- if (is.na(why)) \"\" else \"adapted-input\"
      if (is.na(why)) why <- \"\"",
"      why <- unname(bin2_adapted_input[id])
      cls <- if (is.na(why)) \"\" else \"adapted-input\"
      # a literal pass may carry a note on what it does not show
      if (is.na(why)) why <- unname(bin2_pass_note[id])
      if (is.na(why)) why <- \"\"")
edit_file(l2,
"# an adapted input whose row no longer holds, or is not ported, is a
# stale entry: the reason would describe a pass that is not there",
"# a policy verdict names an unported row, and a pass note a row that
# holds
for (id in setdiff(bin2_policy_unported, led$id[!led$ported])) {
  problems <- c(problems, paste(\"POLICY VERDICT FOR A ROW THAT IS PORTED\",
                                \"OR ABSENT:\", id))
}
# an adapted input whose row no longer holds, or is not ported, is a
# stale entry: the reason would describe a pass that is not there")
edit_file(l2,
"for (id in setdiff(names(bin2_adapted_input), held_ids)) {",
"for (id in setdiff(names(bin2_pass_note), held_ids)) {
  problems <- c(problems, paste(\"PASS NOTE FOR A ROW THAT DOES NOT HOLD:\",
                                id))
}
for (id in setdiff(names(bin2_adapted_input), held_ids)) {")

# ---- 9. release fixes -----------------------------------------------------
# the lanes' reviews and the coordinator's instructions; the fixes and
# their reasons are dev/rel077-release-fix.R, sourced here so a merge
# rerun applies them
source("dev/rel077-release-fix.R", local = TRUE)

# ---- 10. dev/test-backlog.md and dev/round-handoff.md ---------------------
bl <- "dev/test-backlog.md"
stopifnot(bl %in% merged)
bsec <- "dev/rel077-backlog-section.md"
if (file.exists(bsec)) {
  edit_file(bl, "\n## Reference\n",
            paste0("\n", paste(rd(bsec), collapse = "\n"),
                   "\n\n## Reference\n"))
} else cat("NOTE: no", bsec, "yet\n")
hf <- "dev/round-handoff.md"
restore(hf)
hsec <- "dev/rel077-handoff-section.md"
if (file.exists(hsec)) {
  edit_file(hf,
"0.74.0 release, 2026-10-24 at the 0.75.0 release and 2026-10-27 at the
0.76.0 release. Read this, then",
"0.74.0 release, 2026-10-24 at the 0.75.0 release, 2026-10-27 at the
0.76.0 release and 2026-10-28 at the 0.77.0 release. Read this, then")
  edit_file(hf, "\n## The 0.76.0 round, ",
            paste0("\n", paste(rd(hsec), collapse = "\n"),
                   "\n\n## The 0.76.0 round, "))
} else cat("NOTE: no", hsec, "yet\n")
cat("RESOLVE DONE\n")

# ---- 11. the repro of the rescor posterior_predict() defect ---------------
# The backlog's MAJOR item cites the mointerop reviewer's script, which
# the merge does not copy (it skips dev/review-*); it and the helper it
# sources are copied as written (both read the lane tree's paths)
for (f in c("dev/review-mointerop-rescorpp.R",
            "dev/review-mointerop-common.R")) {
  src <- file.path("C:/Users/adf44/source/r/frmtmb-wt-mointerop", f)
  wr(rd(src), f)
  cat("copied", f, "from the mointerop tree\n")
}
cat("RESOLVE DONE\n")
