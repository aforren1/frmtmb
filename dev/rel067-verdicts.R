# The verdict moves of the 0.67.0 consolidation to the ported brms tier
# (dev/brmsport-ledger.tsv), applied idempotently. The three lanes of
# 2026-09-30 each moved rows with their own code; a row is recorded
# only once every lane's code is in one build, so the moves are made
# here, against the record of the merged tree.
#
#   FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r5 sh dev/defects-record.sh
#   Rscript dev/rel067-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# dev/round-20260930.md has the before and after counts;
# dev/rel067-log/ledger-stale.txt is the builder's list of the 22 rows
# that hold, and dev/rel067-log/msgdiff.txt the rows whose message
# moved, which is how the reasons below were checked.
man_f <- "dev/brmsport-verdicts-manual.tsv"
fit_f <- "dev/brmsport-verdicts-manual-fit.tsv"
own_f <- "dev/brmsport-verdicts-own.tsv"

rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
man <- rd(man_f)
fit <- rd(fit_f)
own <- rd(own_f)

# Rows that now HOLD on brms's own assertion, as the lanes reported
to_pass <- c(
  # ordinal: disc, threshold structures, acat links, the standata view's
  # constant dpars and Z columns (families:37 is one row for both
  # packages)
  "priors:14", paste0("standata:", c(724, 725, 738, 739, 740, 745)),
  "families:37",
  # ceplot: effect validity, parnames() on a fit, nsamples(incl_warmup),
  # posterior_samples(pars = ) order
  paste0("brmsfit-methods:", c(205, 995, 595, 635)),
  # formrobust: ar() and friends as functions, autocor(), the arma fill,
  # bernoulli coding, 0 + intercept
  "brm:112", paste0("brmsfit-methods:", c(112, 113, 747)),
  paste0("standata:", c(75, 970, 974)))
# own-words rows whose brms assertion now holds as brms wrote it
own_pass <- c("brmsfit-methods:203", "standata:112")

ggplot_reason <- function(what) paste0(
  what, " returns frmtmb's base-graphics plot objects now (lane ceplot, ",
  "classes frmtmb_ce_plot and frmtmb_hyp_plot), with brms's arguments; ",
  "the row asserts is(object, 'ggplot'), and frmtmb does not depend on ",
  "ggplot2, so it cannot return one")

recl <- list(
  list(id = "brmsfit-methods:154", verdict = "divergence", class = "ggplot",
       reason = ggplot_reason("plot(ce, rug = TRUE, plot = FALSE)")),
  list(id = "brmsfit-methods:162", verdict = "divergence", class = "ggplot",
       reason = ggplot_reason("plot() of a surface display with plot = FALSE")),
  list(id = "brmsfit-methods:164", verdict = "divergence", class = "ggplot",
       reason = ggplot_reason("plot(stype = , plot = FALSE)")),
  list(id = "brmsfit-methods:169", verdict = "divergence", class = "ggplot",
       reason = ggplot_reason(paste0(
         "plot() of the spaghetti display of brmsfit-methods:167 ",
         "with plot = FALSE"))),
  list(id = "brmsfit-methods:391", verdict = "divergence", class = "ggplot",
       reason = ggplot_reason("plot() of a hypothesis with plot = FALSE")),
  list(id = "brmsfit-methods:396", verdict = "divergence", class = "ggplot",
       reason = ggplot_reason(paste0(
         "plot() of a hypothesis with ignore_prior and plot = FALSE ",
         "(ignore_prior changes nothing: a fit has no prior draws)"))),
  list(id = "brmsfit-methods:217", verdict = "divergence", class = "policy",
       reason = paste0(
         "the ordinal default of conditional_effects() stays the ",
         "per-category display, and brms's 'Predictions are treated as ",
         "continuous variables' warning comes with categorical = FALSE ",
         "(user decision, 2026-09-30; lane ceplot 1.5): brms's own warning ",
         "calls its default display likely invalid for ordinal families ",
         "and asks for categorical = TRUE, frmtmb's default")),
  list(id = "standata:1133", verdict = "divergence", class = "policy",
       reason = paste0(
         "frm(drop_unused_levels = FALSE) is accepted now (lane ",
         "formrobust); the unused level's all-zero column xc is dropped ",
         "with frmtmb's rank-deficiency message, as lm() drops it, where ",
         "brms keeps a column that only its prior places (as ",
         "standata:928)")),
  list(id = "brmsfit-methods:314", verdict = "defect", class = "argument",
       reason = paste0(
         "fitted() takes sample_new_levels = 'old_levels' now (lane ",
         "ceplot); the call also asks ndraws, which a maximum-likelihood ",
         "fit has no draws to thin, and the setup dies on it (D3's ",
         "reason), so the assertion reads a stale object")),
  # the lane recorded this row as a class-name divergence like :953; the
  # merged build, and the lane's own library, stop in the refit instead
  # (dev/rel067-update955.R, dev/rel067-log/update955.txt)
  list(id = "brmsfit-methods:955", verdict = "defect", class = "start-values",
       reason = paste0(
         "update(fit2, formula. = bf(count ~ a + b, nl = TRUE)) is no ",
         "longer refused and keeps fit2's parameter formulas with brms's ",
         "message (lane formrobust), but the refit stops from the default ",
         "starting values ('NA/NaN gradient evaluation'); brms's testmode ",
         "does not fit, so brms never meets the start"))
)

man <- man[!(man$id %in% to_pass), , drop = FALSE]
fit <- fit[!(fit$id %in% to_pass), , drop = FALSE]
own <- own[!(own$id %in% own_pass), , drop = FALSE]
for (r in recl) {
  hit <- 0L
  for (nm in c("man", "fit")) {
    x <- get(nm)
    i <- which(x$id == r$id)
    if (!length(i)) next
    x$verdict[i] <- r$verdict
    x$class[i] <- r$class
    x$reason[i] <- r$reason
    assign(nm, x)
    hit <- hit + length(i)
  }
  if (!hit) stop("no verdict row for ", r$id)
}
wr(man, man_f)
wr(fit, fit_f)
wr(own, own_f)
cat("removed", length(to_pass), "manual ids and", length(own_pass),
    "own-words ids; reclassified", length(recl), "\n")
