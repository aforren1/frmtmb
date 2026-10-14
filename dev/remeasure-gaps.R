# Rank what still blocks vignette calls, by how many calls each cause
# blocks, with its status (deliberate or open), the lane of the round
# of 2026-10-13 that works on it, and the row ids.
#
#   Rscript dev/remeasure-gaps.R <tag>
#
# The successor of dev/vigport-gaps.R, reading dev/remeasure-out/<tag>.
# Every failing row of the mechanical port's spell pass (class FAIL or
# CASCADE) and every open row of the hand translation (STILL-OPEN, or
# FAILS-NOW that errors) gets ONE root cause from its error text. A
# cascade ("object 'fit4' not found", or a hand row dispatching on the
# NULL an earlier failure left) is charged to the cause of the row that
# should have made the object, so a gap is charged for all it blocks.
# The script stops if a row matches no cause.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
tag <- commandArgs(trailingOnly = TRUE)[1]
d <- file.path(root, "remeasure-out", tag)
M <- readRDS(file.path(d, "port", "results-merged.rds"))
D <- utils::read.csv(file.path(d, "bv", "drift.csv"), stringsAsFactors = FALSE)
ver <- sub("^# frmtmb ([^ ]+) from.*$", "\\1", readLines(
  file.path(d, "port", "results", "brms_overview.log"), n = 1))

# name, pattern, status, lane. "deliberate" means frmtmb refuses on
# purpose for a stated reason (a maximum-likelihood fit has no
# posterior; there is no Stan program; brms refuses it too); "open"
# means a brms call that frmtmb could serve and does not yet.
gap <- function(cause, pattern, status, lane) {
  data.frame(cause = cause, pattern = pattern, status = status,
             lane = lane, stringsAsFactors = FALSE)
}
C <- rbind(
  gap("loo(), LOO(), waic() on a maximum-likelihood fit",
      paste0(
        "^loo[(][)] is a posterior|^LOO[(][)] is the deprecated ",
        "brms spelling, and on a maximum|^WAIC[(][)] is .* ",
        "maximum"),
      "deliberate: no posterior",
      "no lane"),
  gap("add_criterion(fit, 'loo') on a maximum-likelihood fit",
      "add_criterion[(]criterion = .loo.[)] has no meaning",
      "deliberate: no posterior",
      "no lane"),
  gap("bayes_R2() on a maximum-likelihood fit",
      "^bayes_R2[(][)] is computed per posterior draw",
      "deliberate: no posterior",
      "no lane"),
  gap(paste0(
        "plot() of a fit refuses plot.brmsfit()'s N, variable, ",
        "regex (posterior panels)"),
      "^plot[(][)] cannot honor",
      "deliberate: no posterior",
      "no lane"),
  gap("frm_multiple() result has no draws (plot, as_draws_array, nchains)",
      paste0(
        "has no pooled display for a frm_multiple|needs draws, ",
        "and a frm_multiple[(][)] result has none"),
      "deliberate: no posterior",
      "no lane"),
  gap("brms's Stan-side custom_family()/stanvar()/expose_functions()",
      paste0(
        "unused arguments [(]lb|stanvar|no Stan ",
        "source|expose_functions[(][)] has ",
        "no|expose_functions[(][)] has nothing"),
      "deliberate: Stan-side",
      "no lane"),
  gap("stancode(), standata() on frmtmb objects",
      "^stancode[(][)] has no meaning|^standata[(][)] has no meaning",
      "deliberate: Stan-side",
      "no lane"),
  gap(paste0(
        "brms-only names with no analog (threading, backend, ",
        "shiny, marginal_effects)"),
      paste0(
        "threading|`threads`|`backend`|grainsize|reduce_sum|",
        "launch_shiny|marginal_effects"),
      "deliberate: Stan-side or retired in brms",
      "no lane"),
  gap("LOO(), WAIC() on draws (spellings brms deprecated)",
      paste0(
        "^LOO[(][)] is the deprecated brms spelling and frmtmb ",
        "never had it|^WAIC[(][)] is the deprecated"),
      "deliberate: deprecated in brms",
      "no lane"),
  gap("hypothesis() with no relation (brms refuses it too)",
      "states no relation",
      "deliberate: as brms",
      "no lane"),
  gap("class = sd without nlpar on a nonlinear model (brms refuses it too)",
      "No random-effect SDs match class=sd",
      "deliberate: as brms",
      "no lane"),
  gap(paste0(
        "class Intercept, dpar sigma without a sigma formula ",
        "(brms refuses it too)"),
      "dpar = .sigma. is a density on the log-scale",
      "deliberate: as brms",
      "no lane"),
  gap("a simo prior other than dirichlet (brms refuses it too)",
      "Currently 'dirichlet' is the only valid prior for simplex",
      "deliberate: as brms",
      "no lane"),
  gap("brms attached after frmtmb masks bf()",
      "built by brms::bf",
      "deliberate: coexistence, the refusal says how",
      "no lane"),
  gap(paste0(
        "rescor_matrix() on draws (a frmtmb function; the ",
        "refusal names the fit)"),
      "rescor_matrix[(][)] reads the fitted point estimate",
      "deliberate: frmtmb API",
      "no lane"),
  gap("summary(waic = TRUE), which brms's summary() swallows",
      "summary[(][)] cannot honor .waic|summary[(][)] has no argument .waic",
      "open",
      "postfix"),
  gap("pp_check() on a multivariate model (resp =)",
      "pp_check[(][)] is not supported yet for multivariate",
      "open",
      "postfix"),
  gap("threshold = as an frm() argument (brms 2.23.0 fails the same call)",
      "unused argument [(]threshold",
      "deliberate: as brms",
      "postfix (at most a refusal naming the family argument)"),
  gap("conditional_effects(method = 'predict') on draws",
      "on draws has no method = .predict",
      "open",
      "postfix"),
  gap("frm_sample() of a frm_multiple() result (pooled draws)",
      paste0(
        "frm_sample[(][)] takes a frmtmb fit or a formula; got ",
        "an object of class frmtmb_multiple"),
      "open",
      "postfix"),
  gap(paste0(
        "families brms has and core lacks: discrete_weibull, ",
        "frechet, gen_extreme_value"),
      paste0(
        "^(discrete_weibull|frechet|gen_extreme_value) is not a ",
        "supported family"),
      "open",
      "famuni"),
  gap("families brms has and core lacks: dirichlet, logistic_normal",
      "^(dirichlet|logistic_normal) is not a supported family",
      "open",
      "famsimplex"),
  gap(paste0(
        "wiener scored missing: the hand row probes core, ",
        "frmtmb.eam exports wiener()"),
      "^wiener is not a supported family",
      "instrument: covered by frmtmb.eam",
      "famuni (parity)"),
  gap("mixture: brms's mu1 ~ formula beside the main formula",
      "dpar[(]s[)] not available for family 'mixture",
      "open",
      "no lane"),
  gap("multinomial() needs K, which brms reads from the data",
      "multinomial[(][)] needs the number of categories",
      "open",
      "no lane"),
  gap(paste("loo() of several models on draws (brms compares them in one",
            "call; dev/brms-api-diff.md (c) lists it as Small open work)"),
      "loo[(][)] takes one model here",
      "open",
      "postfix"),
  gap("update() on draws",
      "update[(][)] has no method for draws",
      "open",
      "no lane"),
  gap("log_lik() on an mi() model",
      "log_lik[(][)] is not defined for a model with in-model",
      "open",
      "no lane"),
  gap("horseshoe prior",
      "horseshoe",
      "open",
      "no lane"),
  gap("custom family without a simulator (pp_check)",
      "has no simulator",
      "open",
      "no lane"),
  gap("fixef() returns brms's matrix since 0.61.0 (hand translation stale)",
      "values must be length 5",
      "instrument: stale translation",
      "no lane")
)
cause_of <- function(msg) {
  hit <- C$cause[vapply(C$pattern, function(p) grepl(p, msg), NA)]
  if (length(hit)) hit[1] else "UNCLASSIFIED"
}

bad <- M[M$bucket %in% c("FAIL", "CASCADE"), ]
creator <- function(obj, id) {
  i <- match(id, M$id)
  same <- which(M$vignette == M$vignette[i] & seq_len(nrow(M)) < i &
                  grepl(paste0("^", obj, " *<-"), M$src_run))
  if (length(same)) M$id[max(same)] else NA_character_
}
root_cause <- function(i, depth = 0) {
  msg <- bad$msg_spell[i]
  m <- regmatches(msg, regexec("^object '([^']+)' not found", msg))[[1]]
  if (length(m) == 2 && depth < 10) {
    j <- match(creator(m[2], bad$id[i]), bad$id)
    if (!is.na(j)) return(root_cause(j, depth + 1))
    return(paste("cascade from", m[2]))
  }
  cause_of(msg)
}
bad$cause <- vapply(seq_len(nrow(bad)), root_cause, "")
bad$cascade <- grepl("^object '", bad$msg_spell)

op <- D[D$drift %in% c("STILL-OPEN", "FAILS-NOW") & !D$ok, ]
D$path <- ifelse(grepl("^SAMPLE:", D$label), "SAMPLE", "ML")
op$path <- ifelse(grepl("^SAMPLE:", op$label), "SAMPLE", "ML")
op$cause <- vapply(op$msg, cause_of, "")
op$cascade <- grepl("applied to an object of class \"NULL\"", op$msg)
for (i in which(op$cascade)) {
  j <- which(D$vignette == op$vignette[i] & D$path == op$path[i] &
               D$kind == "model" & !D$ok &
               seq_len(nrow(D)) < match(rownames(op)[i], rownames(D)))
  op$cause[i] <- if (length(j)) cause_of(D$msg[max(j)]) else
    "cascade with no failing model row before it"
}
un <- c(bad$id[bad$cause == "UNCLASSIFIED"],
        paste(op$vignette, op$label)[op$cause == "UNCLASSIFIED"])
if (length(un)) {
  cat("UNCLASSIFIED rows:\n", paste(" ", un, collapse = "\n"), "\n")
  for (i in which(bad$cause == "UNCLASSIFIED")) cat(bad$msg_spell[i], "\n")
  for (i in which(op$cause == "UNCLASSIFIED")) cat(op$msg[i], "\n")
  stop("a failing row matches no cause")
}

# Hand rows are named by vignette and their order among the vignette's
# non-data rows, the key dev/remeasure-history.R uses.
D2 <- D[D$kind != "data", ]
D2$n <- stats::ave(seq_len(nrow(D2)), D2$vignette, FUN = seq_along)
op$n <- D2$n[match(rownames(op), rownames(D2))]
short <- function(v) sub("^brms[_-]", "", v)
tabc <- do.call(rbind, lapply(C$cause, function(k) {
  b <- bad[bad$cause == k, ]; h <- op[op$cause == k, ]
  data.frame(cause = k,
             status = C$status[C$cause == k], lane = C$lane[C$cause == k],
             model = sum(b$kind == "model"), post = sum(b$kind == "post"),
             other = sum(b$kind == "other"), cascade = sum(b$cascade),
             hand_ml = sum(h$path == "ML"),
             hand_sample = sum(h$path == "SAMPLE"),
             hand_cascade = sum(h$cascade),
             mech_ids = paste(sub("^brms_", "", b$id), collapse = ", "),
             hand_ids = paste(sprintf("%s:%d (%s)", short(h$vignette), h$n,
                                      substr(h$label, 1, 60)),
                              collapse = "; "),
             stringsAsFactors = FALSE)
}))
tabc$calls <- tabc$model + tabc$post + tabc$other + tabc$hand_ml +
  tabc$hand_sample
tabc <- tabc[tabc$calls > 0, ]
tabc <- tabc[order(-(tabc$model + tabc$post + tabc$other), -tabc$calls), ]

cat(sprintf("### frmtmb %s, generated by dev/remeasure-gaps.R %s\n\n",
            ver, tag))
cat(sprintf(paste0("Mechanical port, spell pass: %d rows fail (%d model, %d ",
                   "post, %d other; %d of them cascades). Hand translation: ",
                   "%d open rows (%d ML, %d SAMPLE).\n\n"),
            nrow(bad), sum(bad$kind == "model"), sum(bad$kind == "post"),
            sum(bad$kind == "other"), sum(bad$cascade), nrow(op),
            sum(op$path == "ML"), sum(op$path == "SAMPLE")))
cat("Ranked by mechanical calls blocked (model + post + other), then by",
    "all calls.\n\n")
cat("| rank | cause | status | lane | mech model/post/other (cascade) |",
    "hand ML/SAMPLE (cascade) |\n|---|---|---|---|---|---|\n")
for (r in seq_len(nrow(tabc))) {
  x <- tabc[r, ]
  cat(sprintf("| %d | %s | %s | %s | %d/%d/%d (%d) | %d/%d (%d) |\n", r,
              x$cause, x$status, x$lane, x$model, x$post, x$other, x$cascade,
              x$hand_ml, x$hand_sample, x$hand_cascade))
}
cat("\nRow ids per cause (mechanical: vignette.chunk.expr without `brms_`;",
    "hand: vignette:n, the n-th non-data bv() row of the script):\n\n")
for (r in seq_len(nrow(tabc))) {
  x <- tabc[r, ]
  cat(sprintf("%d. %s\n", r, x$cause))
  if (nzchar(x$mech_ids)) cat("   mechanical:", x$mech_ids, "\n")
  if (nzchar(x$hand_ids)) cat("   hand:", x$hand_ids, "\n")
}
cat("\n### totals by status\n\n")
st <- sub(":.*$", "", tabc$status)
for (s in unique(st)) {
  y <- tabc[st == s, ]
  cat(sprintf("%-11s causes %2d  mechanical rows %3d  hand rows %3d\n", s,
              nrow(y), sum(y$model + y$post + y$other),
              sum(y$hand_ml + y$hand_sample)))
}
cat("\n### open calls per lane\n\n")
o <- tabc[st == "open", ]
for (l in unique(o$lane)) {
  y <- o[o$lane == l, ]
  cat(sprintf("%-18s mechanical %2d  hand %2d  (%s)\n", l,
              sum(y$model + y$post + y$other), sum(y$hand_ml + y$hand_sample),
              paste(y$cause, collapse = "; ")))
}

# Defects that block no vignette call but that this lane's estimate
# check and bin-2 port found, ranked by severity, with their owners
# (assigned at fix round 1). The bin-2 rows come from the bin-2 ledger,
# so a row that stops being a defect drops out.
cat("\n### Defects outside the vignette calls, ranked by severity\n\n")
cat("| rank | defect | owner | evidence |\n|---|---|---|---|\n")
cat(paste("| 1 | R1: two `mi()` fits of brms_missings lose their standard",
          "errors at an unchanged optimum (0.70.0, commit ea711734) |",
          "optstart | `dev/remeasure-repro-missings.R`, secmp |\n"))
bl <- file.path(root, "brmsport-bin2-ledger.tsv")
if (file.exists(bl)) {
  L2 <- utils::read.delim(bl, quote = "", stringsAsFactors = FALSE)
  d <- L2[L2$outcome == "defect", ]
  sev <- c(hollow = 1L, method = 2L, argument = 3L)
  d <- d[order(sev[d$class], d$id), ]
  lab <- c("brmsfit-methods:469" =
             "R2: loo() on draws silently ignores newdata",
           "brmsfit-methods:434" = "R3: logLik() of draws has no method",
           "brmsfit-methods:121" = paste("R4: bayes_R2() on draws refuses",
                                         "newdata and re_formula"))
  grp <- ifelse(d$id %in% names(lab), lab[d$id],
                "loo() and waic() of several models on draws")
  r <- 1L
  for (g in unique(grp)) {
    r <- r + 1L
    cat(sprintf("| %d | %s | postfix | %s (bin 2T, %s) |\n", r, g,
                paste(sub("^brmsfit-methods", "methods", d$id[grp == g]),
                      collapse = ", "),
                paste(unique(d$class[grp == g]), collapse = ", ")))
  }
}
