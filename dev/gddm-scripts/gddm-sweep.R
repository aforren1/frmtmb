# lane gddm: the false-alarm and true-positive sweep for the
# within-condition constancy refusal.
#
# Two paired arms over the SAME designs, so that neither can be read on
# its own:
#
#   ok   the condition index is built from every variable in the model,
#        which is what ?gddm asks for. The check must NOT fire. A count
#        of zero here is only worth reading beside the other arm.
#   drop one variable is removed from the index and left in the model,
#        which is the defect. The check MUST fire and must name that
#        variable and the parameter carrying it.
#
# A design with no variable in any formula cannot produce a `drop` case
# and is counted separately, because a refusal that cannot fire is not
# evidence that it does not false-alarm.
#
# The models are the shapes a response-time paper writes: drift by
# coherence, boundary by block, start point by cue, non-decision time by
# subject, subject random effects, a monotonic ordered predictor, a
# spline in coherence, an offset, and the coherence drift with its
# vreal() covariate.
#
# Seed 202609. Arm chosen by GDDM_LIB. Everything runs at
# dry_run = "frame", so no Fokker-Planck solve is taped: the check runs
# at frame assembly and that is all this measures.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
gd_varying_groups_ <- function(...) frmtmb.eam:::gd_varying_groups(...)
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

seed <- 202609
set.seed(seed)
ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
n <- 480L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
# The three crossed factors are DRAWN rather than cycled. A cycled
# subject of period 6 beside a cycled block of period 2 makes the block
# a function of the subject, and the first version of this script was
# cycled: dropping `blk` from the index then changed nothing, which the
# drop arm reported as a miss when it was an aliased design.
d$subj <- factor(sample.int(6L, n, replace = TRUE))
d$cohn <- rep(c(0, 0.128, 0.256, 0.512), each = n / 4L)
d$coh <- factor(d$cohn)
d$cohord <- factor(d$cohn, ordered = TRUE)
d$blk <- factor(sample.int(2L, n, replace = TRUE))
d$cue <- factor(sample(c("left", "right"), n, replace = TRUE))
d$o <- as.numeric(d$blk) / 10
d$trial <- seq_len(n)
# the aliasing the first version had, measured rather than asserted
cat("block within subject, distinct values per subject:",
    paste(tapply(d$blk, d$subj, function(z) length(unique(z))),
          collapse = " "), "\n")

# The models. `vars` names every variable in the model, which is what
# gddm_conditions() has to be given; `fam` is the family; `dp` is the
# formula list after the response.
M <- list()
add <- function(label, vars, dp, mu = NULL,
                fam = gddm(control = ctl),
                aterm = "vint(upper, cond)") {
  M[[length(M) + 1L]] <<- list(label = label, vars = vars, dp = dp,
                               mu = mu, fam = fam, aterm = aterm)
}

add("intercept only", character(0),
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5))
add("drift by coherence", "coh",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5), mu = quote(~coh))
add("drift by coherence, boundary by block", c("coh", "blk"),
    list(quote(bs ~ blk), quote(ndt ~ 1), bias = 0.5), mu = quote(~coh))
add("drift by coherence times block", c("coh", "blk"),
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~coh * blk))
add("non-decision time by subject", "subj",
    list(quote(bs ~ 1), quote(ndt ~ subj), bias = 0.5))
add("start point by cue", "cue",
    list(quote(bs ~ 1), quote(ndt ~ 1), quote(bias ~ cue)))
add("all four parameters", c("coh", "blk", "cue", "subj"),
    list(quote(bs ~ blk), quote(ndt ~ subj), quote(bias ~ cue)),
    mu = quote(~coh))
add("subject random intercept on the drift", "subj",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~1 + (1 | subj)))
add("subject random slope in coherence", c("coh", "subj"),
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~coh + (coh | subj)))
add("subject random intercept on the boundary", "subj",
    list(quote(bs ~ 1 + (1 | subj)), quote(ndt ~ 1), bias = 0.5))
add("drift polynomial in coherence", "cohn",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~poly(cohn, 2)))
add("drift spline in coherence", "cohn",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~s(cohn, k = 4)))
add("drift monotonic in ordered coherence", "cohord",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~mo(cohord)))
add("drift with an offset", "o",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~offset(o)))
add("cell means, no intercept", "coh",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5),
    mu = quote(~0 + coh))
add("drift by numeric coherence", "cohn",
    list(quote(bs ~ 1), quote(ndt ~ 1), bias = 0.5), mu = quote(~cohn))
add("collapsing boundary, rate by block", "blk",
    list(quote(bs ~ 1), quote(tau ~ blk), quote(ndt ~ 1), bias = 0.5),
    fam = gddm(bound = gddm_bound_exponential(), control = ctl))
add("uniform start width by block", "blk",
    list(quote(bs ~ 1), quote(sz ~ blk), quote(ndt ~ 1), bias = 0.5),
    fam = gddm(start = gddm_start_uniform(), control = ctl))
add("lapse rate by block", "blk",
    list(quote(bs ~ 1), quote(ndt ~ 1), quote(lapse ~ blk), bias = 0.5),
    fam = gddm(lapse = "uniform", control = ctl))
add("coherence drift nonlinearity", "blk",
    list(quote(bs ~ blk), quote(ndt ~ 1), quote(alpha ~ 1), bias = 0.5),
    fam = gddm(drift = gddm_drift_coherence(cmax = 0.512),
               control = ctl),
    aterm = "vint(upper, cond) + vreal(cohn)")
add("leaky integration, leak by block", "blk",
    list(quote(bs ~ 1), quote(leak ~ blk), quote(ndt ~ 1), bias = 0.5),
    fam = gddm(drift = list(gddm_drift_constant(), gddm_drift_leak()),
               control = ctl))

build <- function(m, idx_vars) {
  dd <- d
  dd$cond <- if (length(idx_vars)) {
    do.call(gddm_conditions, c(list(dd), lapply(idx_vars, as.name)))
  } else rep(1L, nrow(dd))
  # the coherence drift reads the covariate through vreal(), which must
  # also be constant within a condition: that is the family's own older
  # check and it is part of what a correct index has to satisfy
  if (grepl("vreal", m[["aterm"]], fixed = TRUE)) {
    dd$cond <- as.integer(factor(paste(dd$cond, dd$cohn, sep = "\r")))
  }
  mu <- m[["mu"]] %||% quote(~1)
  lhs <- str2lang(paste0("rt | ", m[["aterm"]]))
  f0 <- as.formula(call("~", lhs, mu[[2L]]))
  args <- c(list(f0), m[["dp"]])
  list(data = dd, spec = do.call(bf, args))
}
`%||%` <- function(x, y) if (is.null(x)) y else x

run <- function(m, idx_vars) {
  b <- build(m, idx_vars)
  r <- tryCatch({
    frm(b[["spec"]], family = m[["fam"]], data = b[["data"]],
        dry_run = "frame")
    ""
  }, error = function(e) conditionMessage(e))
  list(msg = r, ncond = length(unique(b[["data"]]$cond)))
}

ok_fail <- 0L; ok_n <- 0L; ok_live <- 0L
drop_n <- 0L; drop_miss <- 0L; drop_wrongname <- 0L
rows <- list()
for (m in M) {
  vars <- m[["vars"]]
  # is `vars` really every variable in the model?  taken from the
  # formulas rather than trusted, so a mislabelled row cannot make the
  # ok arm look clean
  allv <- unique(unlist(lapply(c(list(m[["mu"]]), m[["dp"]]),
                               function(z) all.vars(z))))
  stopifnot(setequal(intersect(allv, names(d)), vars))

  a <- run(m, vars)
  ok_n <- ok_n + 1L
  if (nzchar(a[["msg"]])) {
    ok_fail <- ok_fail + 1L
    cat("FALSE ALARM [", m[["label"]], "]: ", a[["msg"]], "\n", sep = "")
  }
  # could the check have fired on this design at all?
  live <- length(vars) > 0L && a[["ncond"]] < nrow(d)
  if (live) ok_live <- ok_live + 1L

  for (v in vars) {
    drop_n <- drop_n + 1L
    b <- run(m, setdiff(vars, v))
    if (!nzchar(b[["msg"]])) {
      drop_miss <- drop_miss + 1L
      cat("MISSED [", m[["label"]], "] dropping `", v, "`\n", sep = "")
    } else if (!grepl(paste0("`", v, "`"), b[["msg"]], fixed = TRUE)) {
      drop_wrongname <- drop_wrongname + 1L
      cat("NAMED WRONG [", m[["label"]], "] dropping `", v, "`: ",
          b[["msg"]], "\n", sep = "")
    }
  }
  rows[[length(rows) + 1L]] <-
    data.frame(model = m[["label"]], nvar = length(vars),
               ncond_ok = a[["ncond"]], live = live,
               refused_ok = nzchar(a[["msg"]]))
}

cat("\n== designs\n")
print(do.call(rbind, rows), row.names = FALSE)
cat("\n== ok arm: index built from every variable in the model\n")
cat("  designs:", ok_n, "  of which the check could fire on:", ok_live,
    "\n")
cat("  FALSE ALARMS:", ok_fail, "\n")
cat("\n== drop arm: one variable removed from the index\n")
cat("  cases:", drop_n, "  missed:", drop_miss,
    "  refused without naming the variable:", drop_wrongname, "\n")

# A condition of one row cannot distinguish a constant parameter from a
# varying one, and nothing is refused there.
d1 <- d
d1$cond <- seq_len(nrow(d1))
r1 <- tryCatch({
  frm(bf(rt | vint(upper, cond) ~ coh + cue, bs ~ blk, ndt ~ subj,
         bias = 0.5),
      family = gddm(control = ctl), data = d1, dry_run = "frame")
  "accepted"
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("\none row per condition, every parameter varying between rows:",
    r1, "with", nrow(d1), "conditions\n")

# The check must be silent for every other family in the package and
# for a core family, whatever the data look like.
other <- list(
  gaussian = list(f = quote(gaussian()), s = quote(bf(rt ~ coh + blk))),
  wiener = list(f = quote(wiener()),
                s = quote(bf(rt | dec(upper) ~ coh + blk))),
  lba = list(f = quote(lba(2)),
             s = quote(bf(rt | vint(win) ~ coh + blk))),
  rdm = list(f = quote(rdm(2)),
             s = quote(bf(rt | vint(win) ~ coh + blk))),
  wiener_gng = list(f = quote(wiener_gng(deadline = 2.5)),
                    s = quote(bf(rt | dec(upper) ~ coh + blk))))
d2 <- d
d2$win <- d2$upper + 1L
cat("\n== other families, a covariate varying across every grouping\n")
for (nm in names(other)) {
  r <- tryCatch({
    frm(eval(other[[nm]][["s"]]), family = eval(other[[nm]][["f"]]),
        data = d2, dry_run = "frame")
    "accepted"
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat("  ", nm, ": ", substr(r, 1, 90), "\n", sep = "")
}

# ---------------------------------------------------------------------
# The tolerance: what it absorbs, and what it lets through.
#
# An exact comparison is what gd_check_response() makes on the vreal()
# covariates, and it is the wrong shape here, because a design column is
# COMPUTED. poly() orthogonalizes over the whole column, so two rows
# built from bitwise identical inputs come back different.
P <- poly(d$cohn, 2)
gi <- match(d$cohn, sort(unique(d$cohn)))
ref <- match(gi, gi)
dev <- abs(P - P[ref, , drop = FALSE])
noise <- max(dev)
cat("\n== the tolerance\n")
cat(sprintf(paste0("  poly(cohn, 2) rows from identical inputs",
                   " differ by %.4g, on a column of max %.4g\n"),
            noise, max(abs(P))))
cat(sprintf(paste0("  that is %.4g relative to the column maximum,",
                   " and %.1f ulp of it\n"),
            noise / max(abs(P)),
            noise / (.Machine$double.eps * max(abs(P)))))
# The tolerance the guard actually applies to each pair, taken from
# the package rather than retyped: the largest ratio of a real
# deviation to the tolerance it is compared against. Anything below 1
# is accepted, so this is the headroom.
tolm <- 1e-8 * pmax(abs(P), abs(P[ref, , drop = FALSE])) +
  1e-11 * max(abs(P))
cat(sprintf("  worst deviation / its own tolerance: %.4g\n",
            max(dev / tolm)))
cat(sprintf("  headroom, the reciprocal: %.4g x\n",
            1 / max(dev / tolm)))

# And the smallest within-condition difference it still refuses. The
# probe range goes to 1e-16 because the constant moved: under the
# first shipped tolerance nothing below 1e-9 was refused, and a probe
# that stopped at 1e-12 would report the same floor for both.
smallest <- NA_real_
for (e in 10^-(4:16)) {
  dd <- d
  dd$cond <- gddm_conditions(dd, cohn)
  dd$w <- dd$cohn
  dd$w[seq(2L, nrow(dd), by = 2L)] <-
    dd$w[seq(2L, nrow(dd), by = 2L)] + e
  r <- tryCatch({
    frm(bf(rt | vint(upper, cond) ~ w, bs ~ 1, ndt ~ 1, bias = 0.5),
        family = gddm(control = ctl), data = dd, dry_run = "frame")
    FALSE
  }, error = function(z) TRUE)
  cat(sprintf("  a step of %.0e on a column of max %.3g: %s\n",
              e, max(dd$w), if (r) "refused" else "ACCEPTED"))
  if (r) smallest <- e
}
cat(sprintf("  smallest step refused: %.0e, %.1e relative to the",
            smallest, smallest / max(d$cohn)))
cat(" column maximum\n")

# The wide-range column the review built, checked here too so the
# sweep carries it: the guard's answer must not depend on centering.
wide <- c(rep(c(1, 2, 3), each = 10L), rep(3.15e8, 30L))
wg <- rep(1:2, each = 30L)
wf <- c(1L, 31L)
cat("  wide-range column (1, 2, 3 s beside ten years):",
    length(gd_varying_groups_(wide, wg, wf)), "conditions flagged,",
    length(gd_varying_groups_(wide - mean(wide), wg, wf)),
    "when centered\n")
