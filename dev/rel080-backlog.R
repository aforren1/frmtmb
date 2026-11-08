# dev/test-backlog.md: the sections the lanes of the round of 2026-11-06
# added, gathered into one section "Filed at the 0.80.0 release", with the
# items the reviews left for later (brief, "Survey and records", item 5).
#
#   Rscript dev/rel080-backlog.R
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
p <- file.path(REL, "dev/test-backlog.md")
x <- readLines(p, warn = FALSE)
lane_heads <- c(
  ci = "## Filed by lane ci, round of 2026-11-06",
  latperf = "## Filed in the round of 2026-11-06",
  fitedge = "## Filed in the round opened 2026-11-06",
  nlstart = "## Filed by lane nlstart, round of 2026-11-06",
  bigx = "## Filed in the round opened 2026-11-06")
h2 <- grep("^## ", x)
st <- which(x %in% unique(lane_heads))
stopifnot(length(st) == 5L)
take <- list()
drop <- integer()
for (s in st) {
  e <- min(c(h2[h2 > s], length(x) + 1L)) - 1L
  body <- x[(s + 1L):e]
  # bigx's last item was joined to the one before it (no newline at the
  # end of its file)
  body <- unlist(strsplit(sub("this round\\.- \\*\\*", "this round.\n- **",
                              body), "\n", fixed = TRUE))
  body <- body[!grepl("^### (fitedge|Filed by lane bigx)", body)]
  k <- which(nzchar(trimws(body)))
  take[[length(take) + 1L]] <- body[min(k):max(k)]
  drop <- c(drop, s:e)
}
names(take) <- c("ci", "latperf", "fitedge", "nlstart", "bigx")
x <- x[-drop]
rel_items <- c(
"### From the reviews, filed at the release",
"",
"- **MINOR. A binding `uniform(0, 0.2)` on a slope with no `lb`/`ub` fails",
"  with \"NA/NaN gradient evaluation\"**, on 0.79.1 and on 0.80.0 alike.",
"  With `lb = 0, ub = 0.2` the fit holds the slope at 0.2. brms would want",
"  the bounds; frmtmb could set them from the uniform's support",
"  (dev/reviews/2026-11-06-loglik.md, short check, item 2;",
"  `dev/review-loglik-unifbase.R`).",
"- **MINOR. A switch in `frmtmb_control()` for the nonlinear sign probe.**",
"  0.80.0 has `options(frmtmb.nl_sign_probe = FALSE)` only; a control entry",
"  would travel with the fit (nlstart decision, brief of the 0.80.0",
"  release).",
"- **MINOR, cost. `pbeta()` may be cheaper than the summed CDF** for the",
"  binomial and negative binomial OSA residuals on the tape. Default OSA on",
"  a 2000-row poisson `(1 | g)` fit took 41.9 s (dev/reviews/2026-11-06-osa.md,",
"  short check, item 2).",
"- **MINOR. Fractional weights on a continuous response give `NA`.** The",
"  one-step CDF of a fresh copy of the row is defined there (`U = 0`) and",
"  could be returned; for a count `NA` stays right (osa review, item 8).",
"- **MINOR. Nested quadrature: an unsettled fit can carry a",
"  log-likelihood error of a few 1e-2** (1.7e-2 near separation, 6.2e-2 on",
"  70 of 72 ones); the warning now gives the rule's own error. The memory",
"  fallback factor 42.9 over-predicts binomial and lognormal by 4.7 to 5.9",
"  times, so the default `quad_max_gb = 4` refuses them at 2000 to 4000",
"  rows where they need under 1 GB: calibrate them. `profile(fit)` with no",
"  `parm` stops on the missing argument before the nested refusal",
"  (dev/reviews/2026-11-06-nestquad.md, remaining items).",
"- **MINOR. `fitted()` under a draw on a periodic body can drop the draw:**",
"  `y ~ sin(a)^2 + c` at `a = 0` with a draw sd of `10 pi` reads 0.509",
"  against 1.009, silent, because the steps 0.5, 1 and 1.7 sds land on zeros",
"  of the sine. Steps at irrational multiples would close it",
"  (dev/reviews/2026-11-06-readers.md, confirmation). Q1 of that review",
"  stays for the user: the 11-node value under `frmtmb_mean_not_finite`",
"  reads 222685 on a lognormal mixture at link sd 0.555, where",
"  brms-style draws read about 2.5.",
"- **MINOR. A probe win can return a quieter lower maximum:** bernoulli",
"  Michaelis-Menten seed 55, 1.84 short at code 0, silent, where 0.79.1",
"  was 2.29 short with 2 warnings; NEWS says so",
"  (dev/reviews/2026-11-06-nlstart.md, third short check, MINOR 1).",
"- **QUESTION. A user `vcov =` over the coefficients alone holds the",
"  group-level effects fixed, with a warning** (decided for 0.80.0). The",
"  alternative is to carry it through `J` with the model's own",
"  variance-parameter block (dev/reviews/2026-11-06-mfxse.md, m3).",
"- **QUESTION. Pooled D1 under a prior** runs as a Wald test at the MAP; its",
"  nominal size is not guaranteed (statistic 1.125 against 2.084 without",
"  the prior, seed 702). Kept with a caveat for 0.80.0",
"  (dev/reviews/2026-11-06-loglik.md, short check, remaining).",
"- **MINOR, CI. stan-next has a 150-minute job limit with a 120-minute",
"  step limit**, so a wedged identity tier after a slow setup cancels the",
"  job and skips the gp and methods tiers (dev/reviews/2026-11-06-ci.md,",
"  fix round, MINOR 1). A cold identity tier that hits its limit saves no",
"  compiles (same review, NIT).")
sec <- c("## Filed at the 0.80.0 release (2026-11-06)", "",
         "Gathered from the lanes' own sections of the round of 2026-11-06,",
         "then the reviews' items left for later. Lane records:",
         "`dev/<lane>-findings.md`, reviews `dev/reviews/2026-11-06-<lane>.md`.",
         "")
for (l in names(take)) {
  sec <- c(sec, paste0("### From lane ", l), "", take[[l]], "")
}
sec <- c(sec, rel_items, "")
at <- grep("^## Filed at the 0\\.79\\.1 release", x)
stopifnot(length(at) == 1L)
x <- c(x[seq_len(at - 1L)], sec, x[at:length(x)])
# the file keeps one blank line between sections
x <- x[!(c(FALSE, !nzchar(x[-1L]) & !nzchar(x[-length(x)])))]
x <- x[seq_len(max(which(nzchar(x))))]
con <- file(p, open = "wb")
writeLines(x, con, sep = "\r\n")
close(con)
cat("BACKLOG section:", length(sec), "lines\n")
