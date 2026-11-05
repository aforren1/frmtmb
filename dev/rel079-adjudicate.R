# The adjudication rows at the 0.79.0 release, sourced by dev/rel079-resolve.R.
# Rows the round's lanes have answered go, so that the measurement reads. The
# exact gp() pairs are divergences by decision on one line of the program, the
# nugget, with their reason in NEWS; the quote of the withdrawn gap goes.
# sample.tsv stays as drawn: its status is the stratum a pair was drawn from,
# which the estimator reads, so p010, p011, p059, p146, p169, p158, p166, p171,
# p014 and p167 read "compat row wrong" against the frozen status, and they count
# as agreement (dev/round-20261104.md, "The survey").
f <- "dev/equiv-data/adjudication.tsv"
if (!f %in% merged) {
  h <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)), stdout = TRUE)
  wr(sub("\r$", "", h), f)
}
x <- rd(f)
drop <- c("p013", "p124", "p069", "p043", "p169", "p167", "p049", "p158", "p166")
ids <- sub("\t.*$", "", x)
stopifnot(all(drop %in% ids))
x <- x[!ids %in% drop]
reason <- function(what, const) {
  paste0(what, ": L1 reads brms's own program, whose nugget line is an absolute ",
         "1e-12 where frmtmb's is 1e-12 * sdgp^2, relative because brms's breaks on ",
         "y in other units (NEWS 0.79.0, ?\"frmtmb-gp\"). The constant is ", const,
         " at 1.1 to 8.5 times the survey's tolerance, and the gradients 0.013 to ",
         "0.064 against 1e-3; nothing else in the two models differs ",
         "(dev/ncgp-findings.md, section 11, block l1eq).")
}
rows <- list(
  p022 = reason("exact gp() under mvbf", "-6.3e-4"),
  p023 = reason("exact gp() with a prior on sdgp",
                "0.671, the free prior's constant, with -6.3e-4 to -2.6e-3 on the others"),
  p094 = reason("exact gp() under lognormal", "-2.5e-4"),
  p109 = reason("exact gp() under student", "-2.6e-3"))
for (id in names(rows)) {
  k <- which(sub("\t.*$", "", x) == id)
  stopifnot(length(k) == 1L)
  x[k] <- paste(id, "divergence by decision", "FALSE", "", "", rows[[id]], sep = "\t")
}
# p069 keeps a row with the review's ruling, since the flow reads L4 alone (it
# measures defect L1+L4, at z 0.627, as the control does)
x <- c(x, "p069\tequivalent\tTRUE\tL4\t\tus x importance: agree, weak evidence, the review's ruling. L4 reads max z 0.627 on this data seed, with the correlation at its boundary, and the control reads 0.627 too: on 11 seeds the correction's reading is the control's within 0.022, so the route has no power here in either direction. Against an exact marginal at seed 69, fits at 200, 512 and 2000 draws lie within 0.37 MCSE on all nine, 1 round, no warning, and the L1 identity reads -1.9e-10 (dev/reviews/2026-11-02-impsamp.md, p069). The compat row stays conditional.")
wr(x, f)
cat("adjudication rows left:", length(x) - 1L, "\n")
# the survey keeps its own L1: gpgap's rewrite of it is not taken
edit_file("dev/equiv-findings.md",
          paste0("An edited copy of\n  brms's program never scores as L1: the exact gp() pairs also read their\n",
                 "  program with its one nugget line rewritten to frmtmb's, and that reading\n",
                 "  attributes the gap (`dev/gpgap-findings.md`). Under a prior,"),
          "Under a prior,")
# the release's section of the backlog
bl <- "dev/test-backlog.md"
x <- rd(bl)
a <- which(x == "## Reference")
stopifnot(length(a) == 1L)
wr(c(x[seq_len(a - 1L)], rd("dev/rel079-backlog-section.md"), "", x[a:length(x)]), bl)
cat("edited", bl, "\n")

# The findings records the reviews correct (impsamp, record only; latres RF3). They
# are the lanes' files as copied, and the lanes' own trees stay as reviewed.
edit_file("dev/impsamp-findings.md",
          paste0("2, 0.867 to 18.27 (review 0.867 to 18.63), with no warning, against that\n",
                 "design's reference on the uncentered column. The year table's centered reading\n",
                 "is the one to trust there."),
          paste0("2, 0.867 to 18.27 (review 0.867 to 18.63), with no warning. That reference was\n",
                 "the reviewer's own wrong difference, 1e-3 on a raw column near 2000: against the\n",
                 "exact Hessian on the centered column, by the delta rule, seed 2 reads SE 1.0011\n",
                 "(dev/review-impsamp-fr2-year.R), and no warning is due."))
edit_file("dev/impsamp-findings.md",
          "So only the test's wrapper reaches the refusal now, and no design I tried does.\n",
          paste0("So only the test's wrapper reached the refusal at this reading. Corrected by the\n",
                 "review's short check: a wider search finds 8 natural firings in 1260 fits, 5 at\n",
                 "the default rounds (`dev/review-impsamp-fr2-natural.R`), and they are chaotic.\n"))
edit_file("dev/latres-findings.md",
          paste0("where it does not bind. In the L1 sweep, `unstr()` on `Beta()` and on\n",
                 "`Gamma(link = \"log\")` now stop at a correlation matrix with no Cholesky factor,\n",
                 "so they have no point in brms's coordinates. Those two are the 2 cells with an\n",
                 "error in the report, and `test-latres-brms.R` allows at most 2 such cells. Each\n",
                 "fit is reported as it stands."),
          paste0("where it does not bind. In the L1 sweep, `unstr()` on `Beta()` and on `Gamma()`\n",
                 "(this section first misnamed it `Gamma(link = \"log\")`, which reads 1.7e-13) stop\n",
                 "at a correlation matrix with no Cholesky factor. The test now reads the factor\n",
                 "by construction (`us_chol_L()`), so `Beta()` has a number, const -3.34e-08.\n",
                 "`Gamma()` has no number in that sweep: with L read directly, and with no box,\n",
                 "the check is off by 301 nats, as Stan's `cholesky_factor_corr` does not hold a\n",
                 "matrix this close to singular through unconstraining. Each fit is reported as\n",
                 "it stands."))
