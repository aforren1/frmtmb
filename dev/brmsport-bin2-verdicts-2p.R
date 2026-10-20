# The verdicts of the bin-2P rows ported through the thin adapter
# (dev/brmsport-bin2-adapters.R): groups G3, tests.log_lik.R, and G4,
# tests.posterior_predict.R. Each block runs in frmtmb alone, so a row's
# target is frmtmb and this script writes only the manual verdicts,
# dev/brmsport-bin2-verdicts-manual-2p.tsv, which
# dev/brmsport-bin2-ledger.R reads beside the 2T file. Written as a
# script so the reasons are reviewable R strings.
#
#   Rscript dev/brmsport-bin2-verdicts-2p.R
#
# Every row that holds is a pass by observation and has no entry here.
# Constructions: dev/remeasure-bin2-probe2.R, log
# dev/remeasure-out/probe/bin2-probe2.txt (seeds are the blocks' own).
v <- list()
add <- function(id, verdict, class, reason) {
  v[[length(v) + 1L]] <<- data.frame(id = id, pkg = "frmtmb",
                                     verdict = verdict, class = class,
                                     reason = reason)
}

## families frmtmb does not have ----------------------------------------
absent_family <- function(fam, why) paste0(
  fam, " is not a frmtmb family: the adapter's brmsfamily(\"", fam,
  "\") is refused ('", fam, " is not a supported family'), so the row ",
  "reads the value its last successful line left. ", why)
for (id in c("log_lik:271", "posterior_predict:159")) {
  add(id, "cannot transfer", "absent", absent_family("discrete_weibull",
    "A gap, not a decision"))
}
for (id in c("log_lik:288", "posterior_predict:169")) {
  add(id, "cannot transfer", "absent", absent_family("negbinomial2",
    paste("brms's sigma = 1 / shape spelling of the negative binomial;",
          "frmtmb's negbinomial() is the same density in shape, and the",
          "block's negbinomial row (log_lik:285) holds")))
}
for (id in c("log_lik:324", "posterior_predict:184")) {
  add(id, "cannot transfer", "absent", absent_family("frechet",
    "A gap, not a decision"))
}
for (id in c("log_lik:349", "posterior_predict:190")) {
  add(id, "cannot transfer", "absent", absent_family("gen_extreme_value",
    "A gap, not a decision"))
}

## the reference side cannot run ---------------------------------------
add("log_lik:267", "cannot transfer", "reference", paste(
  "brms's exported reference dbeta_binomial() is extraDistr::dbbinom(),",
  "and extraDistr is not installed here ('Please install the",
  "'extraDistr' package'), so the expected value never exists. frmtmb's",
  "side runs: against the closed form lchoose + lbeta at the block's",
  "seed it agrees to 1.7e-14 absolute (dev/remeasure-bin2-probe2.R,",
  "section 2). Installing extraDistr makes this row run"))
add("log_lik:302", "divergence", "parameterization", paste(
  "a known parameterization divergence (findings R6, fix round 2): the",
  "reference brms:::dcom_poisson() takes mu = lambda^(1 / shape), while",
  "frmtmb's compois mu is the mean (RTMB::dcompois2), so the two log",
  "densities differ at the same numbers (max 94 log units at the",
  "block's seed; at mu = 3, shape = 2 brms's mean is 2.737,",
  "dev/remeasure-bin2-probe2.R section 6). The row's text reaches",
  "brms:::, so the harness records it and does not run it"))

## a feature frmtmb's family lacks --------------------------------------
add("log_lik:392", "cannot transfer", "absent", paste(
  "censored von_mises: frmtmb's von_mises() has no CDF and accepts only",
  "weights, so cens() is refused (the adapter by the family's",
  "accepts_aterms; frm() with 'cens()/trunc() need a family with a",
  "CDF'); brms censors it through pvon_mises(). The row's cens[i] is -1",
  "at the block's seed. The uncensored row (log_lik:389) holds, and the",
  "density matches brms::dvon_mises() to 8.9e-16"))

## ordinal densities outside the ordinal model's domain (fix round 1) --
# The statement is one loop holding 8 assertions, so the harness records
# it whole (it does not hold); dev/brmsport-bin2-split.tsv below gives
# each assertion line its own verdict
d238_why <- paste(
  "brms's block draws its thresholds unsorted (rnorm), and at a",
  "decreasing pair brms's dcumulative() returns a NEGATIVE category",
  "probability (4 of 20 cells at seed 238), which its closed-form",
  "helper reproduces; frmtmb's cumulative builds the categories from",
  "the threshold gaps and returns NaN with a warning, in categories on",
  "either side of the decreasing pair (at draw 5 of bin2-probe3.txt",
  "categories 1 and 2 are NaN where brms's category 1 is 0.805). With",
  "sorted thresholds every cell agrees, and sratio, cratio and acat",
  "agree everywhere, to 2.2e-16 or better (dev/remeasure-bin2-probe3.R,",
  "log dev/remeasure-out/probe/bin2-probe3.txt). A fit never reaches",
  "decreasing thresholds: frmtmb's threshold map orders them")
add("distributions:238", "divergence", "domain", d238_why)
split <- data.frame(
  id = "distributions:238",
  part = c("cumulative-equivalent", "cumulative-dim", "sratio-equivalent",
           "sratio-dim", "cratio-equivalent", "cratio-dim",
           "acat-equivalent", "acat-dim"),
  n = 1L,
  verdict = c("divergence", rep("pass", 7)),
  class = c("domain", rep("per-line", 7)),
  reason = c(d238_why, rep(paste(
    "holds in 32 of 32 iterations of brms's loops, each assertion line",
    "judged alone (dev/remeasure-rev2-d238.R, the review's script, rerun",
    "by the lane: dev/remeasure-out/probe/d238.txt); the cumulative",
    "identity holds in 16 of 32"), 7)),
  assertion = c("expect_equivalent(d_cumul, d_cumul_ch)",
                "expect_equal(dim(d_cumul), c(ndraws, ncat))",
                "expect_equivalent(d_sratio, d_sratio_ch)",
                "expect_equal(dim(d_sratio), c(ndraws, ncat))",
                "expect_equivalent(d_cratio, d_cratio_ch)",
                "expect_equal(dim(d_cratio), c(ndraws, ncat))",
                "expect_equivalent(d_acat, d_acat_ch)",
                "expect_equal(dim(d_acat), c(ndraws, ncat))"))
# fix 2b: brms's inv_link_<ordinal>() array maps (distributions:300),
# one loop of 8 assertion lines, through brms_adapt_ilord(); the
# review's dev/remeasure-rev3-g2.R (seed 300), rerun by the lane:
# dev/remeasure-out/probe/g2-inv-link.txt
d300_why <- paste(
  "brms's block draws its array x = thres - eta unsorted (rnorm), so",
  "for cumulative a decreasing pair gives brms's inv_link_cumulative() a",
  "NEGATIVE probability, which frmtmb's cumulative returns as NaN: the",
  "identity holds in 16 of 32 iterations, the ordered-threshold domain",
  "limit of distributions:238. A fit never reaches decreasing",
  "thresholds")
split <- rbind(split, data.frame(
  id = "distributions:300",
  part = c("cumulative-equivalent", "cumulative-dim", "sratio-equivalent",
           "sratio-dim", "cratio-equivalent", "cratio-dim",
           "acat-equivalent", "acat-dim"),
  n = 1L,
  verdict = c("divergence", rep("pass", 7)),
  class = c("domain", rep("per-line", 7)),
  reason = c(d300_why, rep(paste(
    "holds in 32 of 32 iterations of brms's loops, each assertion line",
    "judged alone (dev/remeasure-rev3-g2.R, the review's script, rerun",
    "by the lane: dev/remeasure-out/probe/g2-inv-link.txt)"), 7)),
  assertion = c("expect_equivalent(il_cumul, il_cumul_ch)",
                "expect_equal(dim(il_cumul), c(ndraws, nobsv, ncat))",
                "expect_equivalent(il_sratio, il_sratio_ch)",
                "expect_equal(dim(il_sratio), c(ndraws, nobsv, ncat))",
                "expect_equivalent(il_cratio, il_cratio_ch)",
                "expect_equal(dim(il_cratio), c(ndraws, nobsv, ncat))",
                "expect_equivalent(il_acat, il_acat_ch)",
                "expect_equal(dim(il_acat), c(ndraws, nobsv, ncat))")))
add("distributions:300", "divergence", "domain", paste(
  "the statement is one loop of 8 assertions; recorded whole it does",
  "not hold, and dev/brmsport-bin2-split.tsv gives each line its",
  "verdict:", d300_why))
utils::write.table(split, "dev/brmsport-bin2-split.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)

## G1 (fix round 2): brms's p and q functions, and two absent families --
no_cdf <- c(6, 52, 64, 83, 125, 150, 164, 169, 174, 179, 185, 195, 200,
            205, 210)
for (k in no_cdf) add(paste0("distributions:", k), "cannot transfer",
  "no-cdf", paste(
  "brms's p function of a family frmtmb has, and the frmtmb family has",
  "no CDF (no $lcdf), so there is nothing to read: a gap.",
  "inverse.gaussian, the one family here with a CDF, is shimmed and its",
  "row holds (distributions:111)"))
for (k in c(8, 68, 71)) add(paste0("distributions:", k), "cannot transfer",
  "no-quantile", paste(
  "decision: brms's q function; no frmtmb family object carries a",
  "quantile function, and SPEC.md:153-154 puts numeric d/p/q/r",
  "functions for post-processing outside frmtmb, in gamlss.dist or",
  "extraDistr"))
for (k in c(92, 95, 97, 100, 102, 134, 137, 139, 141)) {
  add(paste0("distributions:", k), "cannot transfer", "absent", paste(
    "frechet and gen_extreme_value are not frmtmb families (lane famuni",
    "is adding them; the consolidation rerun may flip these rows)"))
}

## the categorical block (fix round 2: tried with no new adapter code) --
for (id in c("log_lik:495", "posterior_predict:344")) {
  add(id, "cannot transfer", "not-ported", paste(
    "categorical: brms's prep holds one linear predictor per non-reference",
    "category (mu1, mu2), which the thin adapter does not map onto",
    "frmtmb's categorical family ('prep$dpars has no mu'); frmtmb has the",
    "family, so the row is transferable at a cost"))
}
for (id in c("log_lik:504", "posterior_predict:349")) {
  add(id, "cannot transfer", "absent", paste(
    "multinomial: the block builds the family as brms spells it,",
    "multinomial(), which frmtmb refuses ('needs the number of",
    "categories'); frmtmb needs K, which brms reads from the data (gap",
    "rank 26)"))
}
for (id in c("log_lik:510", "log_lik:516", "log_lik:523", "log_lik:534",
             "posterior_predict:355", "posterior_predict:360",
             "posterior_predict:361", "posterior_predict:368",
             "posterior_predict:369", "posterior_predict:380",
             "posterior_predict:381")) {
  add(id, "cannot transfer", "absent", paste(
    "dirichlet, dirichlet_multinomial, dirichlet2 or logistic_normal is",
    "not a frmtmb family (lane famsimplex is adding dirichlet and",
    "logistic_normal; the consolidation rerun may flip these rows)"))
}

## a row that tests the adapter, not frmtmb (fix round 1) -------------
add("log_lik:32", "cannot transfer", "adapter", paste(
  "holds, but tests the adapter: the assertion is brms's",
  "log_lik_weight(), ll * weights[i], and the adapter computes that",
  "product itself (dev/brmsport-bin2-adapters.R), because frmtmb weights",
  "a row outside the family, in the objective and in log_lik() of a",
  "fit. The weighted value is not frmtmb's answer here, so it is not a",
  "pass. Removing the adapter's weights line breaks it",
  "(dev/remeasure-rev-noweights.txt), which shows exactly that"))

V <- do.call(rbind, v)
# the verdict moves of the 0.71.0 consolidation
V <- source("dev/rel071-bin2-moves.R")$value(V, "2p")

## lane bin2port: the categorical-like adapter and the extreme-value
## and categorical shims (dev/brmsport-bin2-adapters.R) ----------------
# Rows that hold now carry no manual verdict: the simplex rows and
# categorical through brms_adapt_cat_*() (values against brms's own
# log_lik_<family>() in dev/bin2port-probe1.R, and as value checks in
# the ledger), the frechet and GEV p, d and r rows through brms_g1_shims
# (pfrechet, dgen_extreme_value, pgen_extreme_value and
# rgen_extreme_value).
held <- c("log_lik:495", "log_lik:516", "log_lik:523", "log_lik:534",
          "posterior_predict:344", "posterior_predict:360",
          "posterior_predict:361", "posterior_predict:368",
          "posterior_predict:369", "posterior_predict:380",
          "posterior_predict:381", "distributions:95", "distributions:134",
          "distributions:137", "distributions:141")
stopifnot(all(held %in% V$id))
V <- V[!V$id %in% held, , drop = FALSE]
bp_set <- function(V, id, verdict, class, reason) {
  i <- which(V$id == id)
  stopifnot(length(i) == 1L)
  V$verdict[i] <- verdict
  V$class[i] <- class
  V$reason[i] <- reason
  V
}
for (id in c("log_lik:504", "posterior_predict:349")) {
  V <- bp_set(V, id, "cannot transfer", "absent", paste(
    "multinomial: the block builds the family as brms spells it,",
    "multinomial(), which frmtmb 0.71.0 refuses ('needs the number of",
    "categories'), and the adapter builds frmtmb's family from that same",
    "bare constructor, resolved on brms's ncat categories",
    "(brms_adapt_cat_family()). A gap (rank 26), blocked on lane",
    "catmulti, which makes multinomial() read its category count from",
    "the data. The density the row would then reach agrees with brms's",
    "log_lik_multinomial() on this block's prep to 4.3e-14",
    "(multinomial(K = 3), dev/bin2port-probe1.R)"))
}
for (k in c(92, 102)) {
  # "domain", as log_lik:349 (the GEV support): review finding 8
  V <- bp_set(V, paste0("distributions:", k), "divergence",
              "domain", paste(
    sprintf("brms's %s takes the Frechet scale and the row passes %s;",
            if (k == 92) "dfrechet()" else "rfrechet()",
            if (k == 92) "shape = 1" else "shape = 1:10"),
    "frmtmb's frechet(), like brms's frechet family, takes the mean",
    "mu = scale * gamma(1 - 1 / nu), which is infinite at shape <= 1,",
    "so there is no frmtmb parameter for it. The shim refuses the call",
    "(brms_g1_frechet_dp()) rather than hand the density mu = Inf, on",
    "which brms's length-only assertion would hold over NaN values. At",
    "shape 3 the same shims agree with brms's pfrechet() to 1.1e-16",
    "(distributions:95 and its value check)"))
}
for (k in c(97, 100, 139)) {
  V <- bp_set(V, paste0("distributions:", k), "cannot transfer",
              "no-quantile", paste(
    sprintf("decision: brms's %s; no frmtmb family object carries a",
            if (k == 139) "qgen_extreme_value()" else "qfrechet()"),
    "quantile",
    "function, and SPEC.md:153-154 puts numeric d/p/q/r functions for",
    "post-processing outside frmtmb, in gamlss.dist or extraDistr (as",
    "distributions:8, :68 and :71)"))
}
# posterior_predict:420 (review finding 7, fix round 2) holds and has
# no manual verdict: the adapter returns brms's
# posterior_predict_wiener() container, ndraws x 1
# (dev/brmsport-bin2-adapters.R), and the row's nrow() then counts
# frmtmb.eam's draws, a shape check as for every other family. A
# difference the row does not test: frmtmb.eam draws the time given the
# boundary the row ended at (brms's dec, c(1, 0, 1)), and brms draws
# the boundary too; frmtmb.eam documents it (?frmtmb.eam::wiener,
# 'Post-fit quantities against brms'). brms 2.23.0 skips its own block
# ('skip as long as RWiener fails on R-devel'), so the port runs it in
# place of the skip.
stopifnot(!anyDuplicated(V$id), all(nzchar(V$reason)),
          !any(grepl("\t", V$reason)))

## lane bin2port: the value checks that do not hold -------------------
# A value check compares the adapter's value with brms's own function on
# the same arguments (dev/brmsport-bin2-adapters.R, brms_value_check()).
# One that holds needs no entry; the ledger stops on one that does not
# hold without an entry here, and on an entry for one that holds.
# Constructions: dev/bin2port-probe2.R, log dev/bin2port-out/probe2.txt
# (the blocks' own seeds).
#
# msg_has is the text of the recorded message the verdict was written
# for; the ledger stops, and the tier fails, when the check keeps its
# label but its message no longer contains it (review finding 3)
W <- list()
addv <- function(id, verdict, class, msg_has, reason) {
  W[[length(W) + 1L]] <<- data.frame(id = id, pkg = "frmtmb",
                                     verdict = verdict, class = class,
                                     msg_has = msg_has, reason = reason)
}
addv("log_lik:271", "divergence", "precision",
     "DIFFERS: 25 values (0 NA, 17 infinite)", paste(
  "discrete_weibull at y = 12: brms's log_lik_discrete_weibull() forms",
  "log(mu^(y^shape) - mu^((y + 1)^shape)), which underflows to -Inf in",
  "17 of the 25 draws (brms's block says so: 'they will be -Inf for",
  "this data'); frmtmb's density works on the log scale and is finite",
  "there (down to -3.6e9), matching the same closed form taken on the",
  "log scale to 5.7e-6 absolute, about 1e-15 relative. On the 8 draws",
  "where brms's value is finite the two agree to 2.4e-9 relative,",
  "brms's cancellation. frmtmb is the more accurate side"))
addv("log_lik:462", "defect", "tied-thresholds",
     "NA in 200 values from the adapter, 200 from brms", paste(
  "cumulative with all three thresholds 0 (brms's prep): frmtmb's log",
  "probability of category 1 and of category 4 is NaN, where it is",
  "finite (log F(-mu) and log(1 - F(-mu)), brms's values); the two",
  "zero-width middle categories are -Inf in frmtmb, the exact value,",
  "and NaN in brms. The cause is ord_cumulative_logpmf() (R/families.R):",
  "an end category's row still evaluates the clamped interior pair,",
  "which is log(0) = -Inf at tied thresholds, and masks it as",
  "0 * -Inf. hurdle_cumulative calls the same function on pmax(y, 1)",
  "and masks it with (1 - i0), so its y = 0 rows are NaN too at a tie",
  "of the lowest pair (review finding 4; at tau (0, 0, 1), mu 0, hu 0.2,",
  "y = 0, 1, 4 gives NaN, NaN, -1.536405 where the exact values are",
  "-1.609438, -0.916291, -1.536405; dev/bin2port-defect-cumulative.R).",
  "With",
  "distinct thresholds (-1, 0, 1) every category agrees with brms to",
  "6.2e-14. Not reached by the optimizer (the tape folds the masked",
  "term); reachable by log_lik() on doubles at such draws. Not fixed in",
  "this lane: dev/bin2port-findings.md has the reproduction and a",
  "one-line fix. After the fix this check stays 'differs' with a new",
  "message (frmtmb 0 NA, brms 200), the ledger stops on it, and it is",
  "re-judged: a divergence (domain), as log_lik:478"))
addv("log_lik:478", "divergence", "domain",
     "NA in 0 values from the adapter, 50 from brms", paste(
  "hurdle_cumulative at y = 3 with all thresholds 0: the category has",
  "probability 0, which frmtmb returns as -Inf, the exact log, and",
  "brms's log_diff_exp() as NaN. With distinct thresholds the two",
  "agree to 6.0e-14 (dev/bin2port-probe2.R)"))
for (id in c("distributions:123", "distributions:177", "log_lik:436")) {
  addv(id, "cannot transfer", "reference",
       "Please install the 'extraDistr' package", paste(
    "the reference cannot run here: brms's",
    if (id == "log_lik:436") "log_lik_zero_inflated_beta_binomial()" else
      if (id == "distributions:123") "dbeta_binomial()" else
        "dzero_inflated_beta_binomial()",
    "needs extraDistr, which is not installed ('Please install the",
    "'extraDistr' package'), as for log_lik:267"))
}
addv("log_lik:580", "cannot transfer", "reference",
     "Please install the 'betareg' package", paste(
  "the reference cannot run here: brms's log_lik_xbeta() needs betareg,",
  "which is not installed ('Please install the 'betareg' package')"))
W <- do.call(rbind, W)
# the verdict moves of the 0.72.0 consolidation
mv072 <- source("dev/rel072-bin2-moves.R")$value
W <- mv072(W, "values")
V <- mv072(V, "2p")
stopifnot(!anyDuplicated(W$id), !any(grepl("\t", W$reason)),
          all(nzchar(W$msg_has)))
utils::write.table(W, "dev/brmsport-bin2-valuechecks-manual.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)

utils::write.table(V, "dev/brmsport-bin2-verdicts-manual-2p.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("bin-2P adapter verdicts:", nrow(V), "(",
    paste(sprintf("%s %d", names(table(V$verdict)), table(V$verdict)),
          collapse = ", "), ")\n")
