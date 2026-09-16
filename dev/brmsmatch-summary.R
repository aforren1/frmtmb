## Emits the before-and-after block that dev/brmsmatch-findings.md
## pastes verbatim. Counts and numbers come out of the run, not out of
## the author's memory (dev/lane-rules.md).
##
##   Rscript dev/brmsmatch-summary.R
##
## Both arms in ONE process: the BASE build is loaded in a child so the
## two never share a namespace.
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(posterior))
q(requireNamespace("brms", quietly = TRUE))
q(requireNamespace("bayesplot", quietly = TRUE))

cache <- "dev/stan-cache/brmsmatch-draws.rds"
ds <- readRDS(cache)
arr <- posterior::as_draws_array(ds)
nd <- posterior::ndraws(arr)

## brms's two definitions, read off brms in dev/brmsmatch-formals.R
sdr <- posterior::summarise_draws(arr, rhat = posterior::rhat)
brms_rhat <- stats::setNames(sdr$rhat, sdr$variable)
sde <- posterior::summarise_draws(arr, ess_bulk = posterior::ess_bulk,
                                  ess_tail = posterior::ess_tail)
brms_neff <- stats::setNames(pmin(sde$ess_bulk, sde$ess_tail) / nd,
                             sde$variable)
## and the bulk-only figure the first review recorded, so the two can
## be told apart
brms_neff_bulk <- stats::setNames(sde$ess_bulk / nd, sde$variable)

## the BASE answers, from the round's shared reference build
base <- local({
  f <- tempfile(fileext = ".rds")
  src <- c(
    '.libPaths(c("C:/Users/adf44/source/r/rellib-r3",',
    '            "C:/Users/adf44/source/r/pinlib",',
    '            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))',
    'q <- function(e) suppressWarnings(suppressMessages(e))',
    'q(library(frmtmb)); q(library(frmtmb.sample))',
    'ds <- readRDS("dev/stan-cache/brmsmatch-draws.rds")',
    sprintf('saveRDS(list(rhat = rhat(ds), neff = neff_ratio(ds),
      vars = variables(ds), summ = summary(ds)), %s)', deparse(f)))
  sf <- tempfile(fileext = ".R")
  writeLines(src, sf)
  out <- system2(file.path(R.home("bin"), "Rscript"), shQuote(sf),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(f)) {
    cat(paste(out, collapse = "\n"), "\n")
    stop("the base child wrote nothing")
  }
  readRDS(f)
})

now <- list(rhat = rhat(ds), neff = neff_ratio(ds),
            vars = variables(ds), summ = summary(ds))

rel <- function(got, want) {
  cm <- intersect(names(got), names(want))
  if (!length(cm)) return(NA_real_)
  max(abs(got[cm] - want[cm]) / abs(want[cm]))
}

cat("<!-- BEGIN GENERATED: dev/brmsmatch-summary.R -->\n")
cat("```\n")
cat("== the draws ==\n")
cat("dev/stan-cache/brmsmatch-draws.rds, built by\n")
cat("dev/brmsmatch-measure.R: y ~ x + (1 | g), gaussian, n = 120,\n")
cat("data seed 9, frm_sample(chains = 4, iter = 1000, seed = 20260915)\n")
cat(sprintf("draws %d x %d, chains %d\n", nrow(ds$draws),
            ncol(ds$draws), nchains(ds)))
cat(sprintf("brms %s, posterior %s, rstan %s, bayesplot %s\n",
            packageVersion("brms"), packageVersion("posterior"),
            packageVersion("rstan"), packageVersion("bayesplot")))

cat("\n== 1. rhat() against brms's rhat.brmsfit ==\n")
cat(sprintf("%-10s %-30s %12s\n", "arm", "quantity", "value"))
cat(sprintf("%-10s %-30s %12.8f\n", "BEFORE",
            "max relative difference", rel(base$rhat, brms_rhat)))
cat(sprintf("%-10s %-30s %12.8f\n", "AFTER",
            "max relative difference", rel(now$rhat, brms_rhat)))
cmn <- intersect(names(now$rhat), names(brms_rhat))
cat(sprintf("%-10s %-30s %12s\n", "AFTER", "identical() to brms's",
            identical(unname(now$rhat[cmn]), unname(brms_rhat[cmn]))))
cat(sprintf("%-10s %-30s %12.8f\n", "both",
            "max |brms rhat - 1|", max(abs(brms_rhat - 1))))
cmb <- intersect(names(base$rhat), names(brms_rhat))
cat(sprintf("%-10s %-30s %12.8f\n", "BEFORE",
            "diff / (rhat - 1)",
            max(abs(base$rhat[cmb] - brms_rhat[cmb]) /
                  abs(brms_rhat[cmb] - 1))))
cat(sprintf("%-10s %-30s %12d\n", "BEFORE", "names shared with brms's",
            length(cmb)))
cat(sprintf("%-10s %-30s %12d\n", "AFTER", "names shared with brms's",
            length(cmn)))
cat(sprintf("%-10s %-30s %12d\n", "both", "variables in the draws",
            length(brms_rhat)))

cat("\n== 2. neff_ratio() against brms's neff_ratio.brmsfit ==\n")
cat(sprintf("%-10s %-30s %12.8f\n", "BEFORE",
            "max relative difference", rel(base$neff, brms_neff)))
cat(sprintf("%-10s %-30s %12.8f\n", "AFTER",
            "max relative difference", rel(now$neff, brms_neff)))
cmn <- intersect(names(now$neff), names(brms_neff))
cat(sprintf("%-10s %-30s %12s\n", "AFTER", "identical() to brms's",
            identical(unname(now$neff[cmn]), unname(brms_neff[cmn]))))
cat(sprintf("%-10s %-30s %12.8f\n", "BEFORE",
            "vs bulk-only (the R1 figure)",
            rel(base$neff, brms_neff_bulk)))
cat(sprintf("%-10s %-30s %12.8f\n", "ratio",
            "neff gap / rhat gap", rel(base$neff, brms_neff) /
              rel(base$rhat, brms_rhat)))

cat("\n== 3. the names ==\n")
## wrapped: this block is pasted verbatim into the findings, where
## house style is 80 columns
show <- function(lbl, v) {
  txt <- paste(v, collapse = " ")
  if (!nzchar(txt)) txt <- "(none)"
  lines <- strwrap(txt, width = 48)
  cat(sprintf("%-28s %s\n", lbl, lines[[1L]]))
  for (l in lines[-1L]) cat(sprintf("%-28s %s\n", "", l))
}
show("variables(ds)", now$vars)
show("BEFORE names(rhat(ds))", names(base$rhat))
show("AFTER  names(rhat(ds))", names(now$rhat))
show("BEFORE not in variables()", setdiff(names(base$rhat), now$vars))
show("AFTER  not in variables()", setdiff(names(now$rhat), now$vars))
cat(sprintf("%-28s BEFORE %s   AFTER %s\n", "rhat(ds)[\"x\"]",
            format(unname(base$rhat["x"])),
            format(unname(now$rhat["x"]), digits = 7)))
cat(sprintf("%-28s BEFORE %s   AFTER %s\n", "neff_ratio(ds)[\"x\"]",
            format(unname(base$neff["x"])),
            format(unname(now$neff["x"]), digits = 7)))

## Whether the package agrees with ITSELF. brms:::summary.brmsfit adds
## Rhat = posterior::rhat, Bulk_ESS = ess_bulk, Tail_ESS = ess_tail, so
## in brms summary(fit)[, "Rhat"] IS rhat(fit).
cat("\n== 4. summary(ds) against brms's, and against rhat() ==\n")
line <- function(a, q, v) cat(sprintf("%-8s %-28s %s\n", a, q, v))
for (a in c("BEFORE", "AFTER")) {
  s <- if (a == "BEFORE") base$summ else now$summ
  rn <- rownames(s)
  line(a, "columns", paste(colnames(s), collapse = " "))
  ## against brms's definition, row by row, on the frmtmb names the
  ## summary table has always used
  ok <- rn[rn %in% names(brms_rhat)]
  line(a, "rows comparable to brms's", length(ok))
  d <- abs(s[ok, "Rhat"] - brms_rhat[ok])
  line(a, "max |Rhat - brms's rhat|", sprintf("%12.8f", max(d)))
  line(a, "identical() to brms's", identical(unname(s[ok, "Rhat"]),
                                             unname(brms_rhat[ok])))
  ## the sign question a reader actually asks of an R-hat column
  cross <- ok[(s[ok, "Rhat"] < 1) != (brms_rhat[ok] < 1)]
  line(a, "opposite sides of 1",
       if (length(cross)) paste(cross, collapse = " ") else "none")
  ## The effective size the table offers, against brms's neff_ratio,
  ## which is min(bulk, tail). BEFORE there is one column, rstan's
  ## n_eff; AFTER there are two and their minimum is the quantity.
  got <- if (all(c("Bulk_ESS", "Tail_ESS") %in% colnames(s))) {
    pmin(s[ok, "Bulk_ESS"], s[ok, "Tail_ESS"])
  } else {
    s[ok, "n_eff"]
  }
  want <- brms_neff[ok] * nd
  line(a, "max |table ESS - brms's|", sprintf("%12.4f",
                                              max(abs(got - want))))
  line(a, "max relative overstatement",
       sprintf("%12.6f", max((got - want) / want)))
  ## whether the package agrees with ITSELF
  own <- if (a == "BEFORE") base$rhat else now$rhat
  k <- rn[rn %in% names(own)]
  line(a, "rows rhat() can address", paste0(length(k), " of ",
                                            length(rn)))
  if (length(k)) {
    line(a, "identical() to this rhat()",
         identical(unname(s[k, "Rhat"]), unname(own[k])))
  }
}

cat("\n== 5. the ten positional signatures ==\n")
pos_args <- function(f) {
  a <- names(formals(f))
  a[seq_len(match("...", a, nomatch = length(a) + 1L) - 1L)]
}
first_div <- function(b, o) {
  k <- min(length(b), length(o))
  if (!k) return(NA_integer_)
  d <- which(b[seq_len(k)] != o[seq_len(k)])
  if (length(d)) d[[1L]] else NA_integer_
}
ten <- c("as.mcmc", "log_lik", "mcmc_plot", "posterior_epred",
         "posterior_interval", "posterior_linpred", "posterior_predict",
         "pp_mixture", "predictive_error", "psis")
tb <- get(".__S3MethodsTable__.", envir = asNamespace("frmtmb.sample"),
          inherits = FALSE)
## the BEFORE formals, read off the shared reference build
base_fm <- local({
  f <- tempfile(fileext = ".rds")
  sf <- tempfile(fileext = ".R")
  writeLines(c(
    '.libPaths(c("C:/Users/adf44/source/r/rellib-r3",',
    '            "C:/Users/adf44/source/r/pinlib",',
    '            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))',
    'suppressMessages(library(frmtmb.sample))',
    'tb <- get(".__S3MethodsTable__.",',
    '          envir = asNamespace("frmtmb.sample"), inherits = FALSE)',
    sprintf('ten <- %s', paste(deparse(ten), collapse = "")),
    'out <- lapply(ten, function(n)',
    '  names(formals(get(paste0(n, ".frmtmb_draws"), envir = tb))))',
    'names(out) <- ten',
    sprintf('saveRDS(out, %s)', deparse(f))), sf)
  system2(file.path(R.home("bin"), "Rscript"), shQuote(sf),
          stdout = TRUE, stderr = TRUE)
  readRDS(f)
})
cut_dots <- function(a) a[seq_len(match("...", a,
                                        nomatch = length(a) + 1L) - 1L)]
cat(sprintf("%-20s %-5s %-16s %-16s %-5s\n", "generic", "before",
            "brms's", "before ours", "after"))
nbad_before <- 0L
nbad_after <- 0L
for (nm in ten) {
  b <- pos_args(get(paste0(nm, ".brmsfit"), envir = asNamespace("brms")))
  o0 <- cut_dots(base_fm[[nm]])
  o1 <- pos_args(get(paste0(nm, ".frmtmb_draws"), envir = tb))
  p0 <- first_div(b, o0)
  p1 <- first_div(b, o1)
  if (!is.na(p0)) nbad_before <- nbad_before + 1L
  if (!is.na(p1)) nbad_after <- nbad_after + 1L
  cat(sprintf("%-20s %-5s %-16s %-16s %-5s\n", nm,
              if (is.na(p0)) "ok" else as.character(p0),
              if (is.na(p0)) "" else b[[p0]],
              if (is.na(p0)) "" else o0[[p0]],
              if (is.na(p1)) "ok" else as.character(p1)))
}
cat(sprintf("\ndiverging at a positional slot: BEFORE %d of %d, ",
            nbad_before, length(ten)))
cat(sprintf("AFTER %d of %d\n", nbad_after, length(ten)))
cat("```\n")
cat("<!-- END GENERATED -->\n")
cat("DONE\n")
