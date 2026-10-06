# Rank the remaining gaps by how many vignette calls each one blocks.
#
#   Rscript dev/vigport-gaps.R
#
# Why a script: the ranking feeds the next round, so it must come from
# the records and not from a reading of them. Every failing row of the
# mechanical port's spell pass is assigned ONE root cause by its error
# text. A cascade row ("object 'fit4' not found") is charged to the
# cause of the expression that should have made the object, so a gap is
# charged for everything it takes down. The hand translation's open rows
# (STILL-OPEN, FAILS-NOW) are counted against the same causes as a
# second column, ML and SAMPLE paths apart.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
M <- readRDS(file.path(root, "vigport-port-out/r5/results-merged.rds"))
D <- utils::read.csv(file.path(root, "vigport-bv-out/r5/drift.csv"),
                     stringsAsFactors = FALSE)

causes <- c(
  "loo(), LOO(), waic() on a maximum-likelihood fit" =
    paste0("^loo\\(\\) is a posterior|^LOO\\(\\) is the deprecated|",
           "^WAIC\\(\\) is|^loo\\(\\) takes one model"),
  "add_criterion() does not exist" = "add_criterion",
  "plot() of a fit lacks plot.brmsfit()'s N, variable, regex" =
    "^plot\\(\\) has no argument",
  "summary() refuses brms's swallowed waic =" =
    "^summary\\(\\) has no argument .waic",
  "frm_multiple() result has no pooled post-processing" =
    paste0("frm_multiple\\(\\) result|takes a frmtmb fit or a formula; ",
           "got an object of class frmtmb_multiple"),
  "pp_check() refuses a multivariate fit" =
    "pp_check\\(\\) is not supported yet for multivariate",
  "bayes_R2() on a maximum-likelihood fit" =
    "^bayes_R2\\(\\) is computed per posterior draw",
  "group-level mo() is refused" = "mo\\(\\) is only supported",
  "threshold = as an frm() argument (brms 2.23 drops it silently)" =
    "unused argument \\(threshold",
  "brms's Stan-side custom_family()/stanvar() spelling" =
    paste0("unused arguments \\(lb|stanvar|Stan source|no Stan program|",
           "expose_functions\\(\\) has no"),
  "nonlinear start values once the priors are dropped" =
    "default starting values",
  "dirichlet / simo prior on a mo() simplex" =
    "dirichlet|class = \"simo\"",
  "mixture: brms's mu1 ~ formula beside the main formula" =
    "dpar\\(s\\) not available for family 'mixture",
  "pp_mixture(), stancode(), standata() on a fit: no refusal" =
    "pp_mixture|Data must be specified using the 'data' argument",
  "conditional_effects(method = 'predict') on draws" =
    "on draws has no method",
  "log_lik() on an mi() model" =
    "log_lik\\(\\) is not defined for a model with in-model",
  "families brms has and frmtmb lacks" = "is not a supported family",
  "update() on draws" = "update\\(\\) has no method for draws",
  "plot() of draws" = "plot\\(\\) has no display for frmtmb draws",
  "brms-only names with no frmtmb analog (threading, shinystan, ...)" =
    paste0("threading|`threads`|`backend`|grainsize|reduce_sum|",
           "launch_shiny|marginal_effects"),
  "horseshoe prior" = "horseshoe",
  "custom family without a simulator" = "has no simulator",
  "brms attached after frmtmb (bf() masked)" = "built by brms::bf",
  "hypothesis() with no relation (brms refuses it too)" = "states no relation",
  "class = sd without nlpar on a nonlinear model (brms refuses it too)" =
    "No random-effect SDs match class=sd",
  "class Intercept, dpar sigma without a sigma formula (brms refuses it too)" =
    "dpar = \"sigma\" is a density on the log-scale",
  "fixef() returns brms's matrix since 0.61.0 (translation stale)" =
    "values must be length 5",
  "multinomial() needs K, which brms reads from the data" =
    "multinomial\\(\\) needs the number of categories",
  "rescor_matrix() on draws (refused, names the fit)" =
    "rescor_matrix\\(\\) reads the fitted point estimate"
)
cause_of <- function(msg) {
  hit <- names(causes)[vapply(causes, function(p) grepl(p, msg), NA)]
  if (length(hit)) hit[1] else "UNCLASSIFIED"
}

# The spell pass, model, post and other rows that did not run.
# Scored as the 0.34.0 record scored them: a row is open when its class
# is FAIL or CASCADE. A row CLEAN in the raw pass is not open even where
# a stale patch makes it fail in the spell pass.
bad <- M[M$bucket %in% c("FAIL", "CASCADE"), ]
bad$cause <- NA_character_
# the last assignment to `obj` before row `id`, in the same vignette
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
    src <- creator(m[2], bad$id[i])
    j <- match(src, bad$id)
    if (!is.na(j)) return(root_cause(j, depth + 1))
    return(paste("cascade from", m[2]))
  }
  cause_of(msg)
}
for (i in seq_len(nrow(bad))) bad$cause[i] <- root_cause(i)
bad$cascade <- grepl("^object '", bad$msg_spell)

cat("## mechanical port (spell pass, 0.67.0): calls blocked per cause\n\n")
tb <- table(bad$cause, factor(bad$kind, levels = c("model", "post", "other")))
ord <- order(-(tb[, "model"] + tb[, "post"]), -tb[, "other"])
cat("| rank | cause | model | post | other | of which cascade |\n")
cat("|---|---|---|---|---|---|\n")
for (r in seq_along(ord)) {
  k <- rownames(tb)[ord[r]]
  cat(sprintf("| %d | %s | %d | %d | %d | %d |\n", r, k, tb[k, "model"],
              tb[k, "post"], tb[k, "other"],
              sum(bad$cascade & bad$cause == k)))
}
cat(sprintf(paste0("\nrows: %d failing (%d model, %d post, %d other),",
                   " %d unclassified\n"),
            nrow(bad), sum(bad$kind == "model"), sum(bad$kind == "post"),
            sum(bad$kind == "other"), sum(bad$cause == "UNCLASSIFIED")))
for (i in which(bad$cause == "UNCLASSIFIED")) {
  cat("  UNCLASSIFIED", bad$id[i], substr(bad$msg_spell[i], 1, 120), "\n")
}

cat("\n## hand translation (0.67.0): open rows per cause\n\n")
op <- D[D$drift %in% c("STILL-OPEN", "FAILS-NOW") & !D$ok, ]
op$path <- ifelse(grepl("^SAMPLE:", op$label), "SAMPLE", "ML")
op$cause <- vapply(op$msg, cause_of, "")
# A NULL fit downstream of an earlier failure is a cascade in this
# harness: bv() returns NULL and the next call dispatches on it. It is
# charged, as on the mechanical side, to the cause of the nearest
# earlier failing model row of the same vignette and path.
D$path <- ifelse(grepl("^SAMPLE:", D$label), "SAMPLE", "ML")
op$cascade <- grepl("applied to an object of class \"NULL\"", op$msg)
for (i in which(op$cascade)) {
  j <- which(D$vignette == op$vignette[i] & D$path == op$path[i] &
               D$kind == "model" & !D$ok &
               seq_len(nrow(D)) < match(rownames(op)[i], rownames(D)))
  op$cause[i] <- if (length(j)) cause_of(D$msg[max(j)]) else
    "cascade with no failing model row before it"
}
th <- table(op$cause, factor(op$path, levels = c("ML", "SAMPLE")))
oh <- order(-(th[, "ML"] + th[, "SAMPLE"]))
cat("| cause | ML | SAMPLE | of which cascade |\n|---|---|---|---|\n")
for (k in rownames(th)[oh]) {
  cat(sprintf("| %s | %d | %d | %d |\n", k, th[k, "ML"], th[k, "SAMPLE"],
              sum(op$cascade & op$cause == k)))
}
cat(sprintf("\nrows: %d open, %d unclassified\n", nrow(op),
            sum(op$cause == "UNCLASSIFIED")))
for (i in which(op$cause == "UNCLASSIFIED")) {
  cat("  UNCLASSIFIED", op$vignette[i], "/", op$label[i], ":",
      substr(op$msg[i], 1, 120), "\n")
}
