## Reviewer's re-check of the two BLOCKERs, on the reviewer's OWN draws
## (dev/stan-cache/bmrev-draws.rds: data seed 17, frm_sample(chains = 4,
## iter = 1000, seed = 31337)).
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
q(library(posterior)); q(library(brms))
cat("frmtmb.sample from ", dirname(system.file(package = "frmtmb.sample")),
    "\n")
ds <- readRDS("dev/stan-cache/bmrev-draws.rds")
cat("draws ", nrow(ds$draws), " x ", ncol(ds$draws), "\n\n")

## the shim that lets brms's OWN bodies run on OUR array
arr_all <- posterior::as_draws_array(ds)
ARR <- arr_all
shim <- structure(list(), class = "bmrevshim")
registerS3method("as_draws_array", "bmrevshim",
                 function(x, variable = NULL, ...) {
                   if (is.null(variable)) ARR else
                     posterior::subset_draws(ARR, variable = variable, ...)
                 }, envir = asNamespace("posterior"))
utils::assignInNamespace("contains_draws",
                         function(x, ...) invisible(TRUE), ns = "brms")

## ================= BLOCKER 1: summary(ds) =========================
cat("================ BLOCKER 1: summary(ds) ================\n")
s <- summary(ds)
cat("colnames: ", paste(colnames(s), collapse = " "), "\n")
cat("rownames: ", paste(rownames(s), collapse = " "), "\n\n")

keep <- rownames(s)
## 1a. the three columns against brms's OWN summary measures, run by
## brms's own summarise_draws call shape on the same subset array
ARR <- posterior::subset_draws(arr_all, variable = keep)
bd <- posterior::summarise_draws(ARR, Rhat = posterior::rhat,
                                 Bulk_ESS = posterior::ess_bulk,
                                 Tail_ESS = posterior::ess_tail)
i <- match(keep, bd$variable)
for (cn in c("Rhat", "Bulk_ESS", "Tail_ESS")) {
  cat(sprintf("%-9s identical() to brms's measure: %-6s max|diff| %s\n",
              cn, identical(unname(s[, cn]), unname(bd[[cn]][i])),
              format(max(abs(s[, cn] - bd[[cn]][i])), digits = 17)))
}
## 1b. and against brms's OWN rhat.brmsfit / neff_ratio.brmsfit bodies
ARR <- arr_all
br <- brms:::rhat.brmsfit(shim)
bn <- brms:::neff_ratio.brmsfit(shim)
cat("\nsummary Rhat vs brms:::rhat.brmsfit on the same draws\n")
cat("  identical(): ", identical(unname(s[, "Rhat"]), unname(br[keep])),
    "  max|diff|: ",
    format(max(abs(s[, "Rhat"] - br[keep])), digits = 17), "\n")
cat("  rows summary has of rhat(ds)'s: ", length(keep), " of ",
    length(br), " (the ", length(br) - length(keep),
    " dropped are lp__ and the b[] modes)\n", sep = "")
cat("  variables on opposite sides of 1: ",
    sum((s[, "Rhat"] > 1) != (br[keep] > 1)), "\n")
cat("  addressable by rhat(ds)[rowname]: ",
    sum(!is.na(br[keep])), " of ", length(keep), "\n")
cat("  n_eff column still present:      ",
    "n_eff" %in% colnames(s), "\n")

## 1c. the pmin identity: construction or measurement?
cat("\n-- the pmin(Bulk_ESS, Tail_ESS) identity --\n")
nr <- neff_ratio(ds)
lhs <- pmin(s[, "Bulk_ESS"], s[, "Tail_ESS"])
rhs <- unname(nr[keep]) * posterior::ndraws(ds)
cat("identical(pmin(B,T), neff_ratio(ds)*ndraws(ds)): ",
    identical(unname(lhs), unname(rhs)), "\n")
cat("max |lhs - rhs|: ", format(max(abs(lhs - rhs)), digits = 17), "\n")
cat("max ulp gap:     ",
    max(abs(lhs - rhs) / pmax(.Machine$double.eps * abs(rhs),
                              .Machine$double.xmin)), "\n")
cat("ndraws(ds) == posterior::ndraws(as_draws_array(ds)): ",
    identical(posterior::ndraws(ds),
              posterior::ndraws(arr_all)), "\n")
cat("NOTE the two sides call posterior::ess_bulk/ess_tail on the same\n")
cat("per-variable draws, so this is an IDENTITY up to the x/N*N\n")
cat("round trip, not an independent measurement.\n")
## the round trip, isolated
rt <- (unname(lhs) / posterior::ndraws(ds)) * posterior::ndraws(ds)
cat("x/N*N round trip on these values is exact: ",
    identical(rt, unname(lhs)), "\n")

## 1d. NOTHING ELSE MOVED: mean, sd and the two quantiles
cat("\n-- did any other summary column move? --\n")
m <- ds$draws
ref <- t(vapply(keep, function(nm)
  c(mean = mean(m[, nm]), sd = stats::sd(m[, nm]),
    `2.5%` = unname(stats::quantile(m[, nm], 0.025)),
    `97.5%` = unname(stats::quantile(m[, nm], 0.975))), numeric(4)))
for (cn in c("mean", "sd", "2.5%", "97.5%")) {
  cat(sprintf("%-7s identical() to the raw draws statistic: %s\n",
              cn, identical(unname(s[, cn]), unname(ref[, cn]))))
}
cat("column count: ", ncol(s), " (was 6: mean sd 2.5% 97.5% n_eff Rhat)\n")
## and against the BASE build's first four columns, in a child process
cat("\n-- base build's first four columns, child process --\n")
tf <- tempfile(fileext = ".rds")
scr <- tempfile(fileext = ".R")
writeLines(c(
  '.libPaths(c("C:/Users/adf44/source/r/rellib-r3",',
  '            "C:/Users/adf44/source/r/pinlib",',
  '            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))',
  'suppressMessages(library(frmtmb)); suppressMessages(library(frmtmb.sample))',
  'ds <- readRDS("dev/stan-cache/bmrev-draws.rds")',
  sprintf('saveRDS(summary(ds), %s)', deparse(tf))), scr)
system2(file.path(R.home("bin"), "Rscript"), scr,
        stdout = NULL, stderr = NULL)
if (file.exists(tf)) {
  s0 <- readRDS(tf)
  cat("base colnames: ", paste(colnames(s0), collapse = " "), "\n")
  cat("base rownames identical: ",
      identical(rownames(s0), rownames(s)), "\n")
  cm <- intersect(colnames(s0), colnames(s))
  for (cn in cm) {
    cat(sprintf("  %-7s identical() base vs now: %-6s max|diff| %s\n",
                cn, identical(unname(s0[, cn]), unname(s[, cn])),
                format(max(abs(s0[, cn] - s[, cn])), digits = 17)))
  }
  cat("columns only in base: ",
      paste(setdiff(colnames(s0), colnames(s)), collapse = " "), "\n")
  cat("columns only in now:  ",
      paste(setdiff(colnames(s), colnames(s0)), collapse = " "), "\n")
} else cat("base summary NOT captured\n")

## ================= BLOCKER 2: pars ================================
cat("\n================ BLOCKER 2: pars ================\n")
cat("formals rhat:       ",
    paste(names(formals(getS3method("rhat", "frmtmb_draws"))),
          collapse = " "), "\n")
cat("formals neff_ratio: ",
    paste(names(formals(getS3method("neff_ratio", "frmtmb_draws"))),
          collapse = " "), "\n")
cat("brms rhat:          ",
    paste(names(formals(brms:::rhat.brmsfit)), collapse = " "), "\n")
cat("brms neff_ratio:    ",
    paste(names(formals(brms:::neff_ratio.brmsfit)), collapse = " "),
    "\n\n")

ARR <- arr_all
res <- function(e) tryCatch({
  v <- q(e)
  paste0("OK n=", length(v), " [",
         paste(utils::head(names(v), 3), collapse = ","), "]")
}, error = function(c) paste0("ERROR: ", conditionMessage(c)))
probe <- list(
  list("missing", quote(rhat(ds)), quote(brms:::rhat.brmsfit(shim))),
  list("\"x\"",   quote(rhat(ds, "x")),
       quote(brms:::rhat.brmsfit(shim, "x"))),
  list("\"^x$\"", quote(rhat(ds, "^x$")),
       quote(brms:::rhat.brmsfit(shim, "^x$"))),
  list("NULL",    quote(rhat(ds, NULL)),
       quote(brms:::rhat.brmsfit(shim, NULL))),
  list("NA",      quote(rhat(ds, NA)),
       quote(brms:::rhat.brmsfit(shim, NA))),
  list("TRUE",    quote(rhat(ds, TRUE)),
       quote(brms:::rhat.brmsfit(shim, TRUE))),
  list("0.9",     quote(rhat(ds, 0.9)),
       quote(brms:::rhat.brmsfit(shim, 0.9))),
  list("regex=T", quote(rhat(ds, "^b", regex = TRUE)),
       quote(brms:::rhat.brmsfit(shim, "^b", regex = TRUE))),
  list("c(x,lp__)", quote(rhat(ds, c("x", "lp__"))),
       quote(brms:::rhat.brmsfit(shim, c("x", "lp__"))))
)
cat(sprintf("%-11s %-46s %s\n", "pars", "brms", "here"))
nsame <- 0L
for (p in probe) {
  b <- res(eval(p[[3]])); o <- res(eval(p[[2]]))
  same <- identical(b, o)
  nsame <- nsame + same
  cat(sprintf("%-11s %-46s %s%s\n", p[[1]], substr(b, 1, 46),
              substr(o, 1, 46), if (same) "" else "   <-- DIFFERS"))
}
cat("\nmatching: ", nsame, " of ", length(probe), "\n")

## neff_ratio on the same five
cat("\nneff_ratio, the same probe:\n")
nsame2 <- 0L
for (p in list(list("missing", quote(neff_ratio(ds)),
                    quote(brms:::neff_ratio.brmsfit(shim))),
               list("\"x\"", quote(neff_ratio(ds, "x")),
                    quote(brms:::neff_ratio.brmsfit(shim, "x"))),
               list("\"^x$\"", quote(neff_ratio(ds, "^x$")),
                    quote(brms:::neff_ratio.brmsfit(shim, "^x$"))),
               list("NULL", quote(neff_ratio(ds, NULL)),
                    quote(brms:::neff_ratio.brmsfit(shim, NULL))),
               list("NA", quote(neff_ratio(ds, NA)),
                    quote(brms:::neff_ratio.brmsfit(shim, NA))))) {
  b <- res(eval(p[[3]])); o <- res(eval(p[[2]]))
  nsame2 <- nsame2 + identical(b, o)
  cat(sprintf("%-11s %-46s %s%s\n", p[[1]], substr(b, 1, 46),
              substr(o, 1, 46),
              if (identical(b, o)) "" else "   <-- DIFFERS"))
}
cat("matching: ", nsame2, " of 5\n")

## where does the refusal come from?
cat("\n-- provenance of the two refusals --\n")
for (v in list(TRUE, 0.9, "nosuchvar")) {
  e <- tryCatch(q(rhat(ds, v)), error = function(c) c)
  cat(sprintf("rhat(ds, %-11s) msg:  %s\n",
              paste(deparse(v), collapse = ""),
              sub("\n.*$", "", conditionMessage(e))))
  cat("                          call: ",
      paste(deparse(conditionCall(e)), collapse = " "), "\n")
}
cat("\nposterior::subset_draws() directly, for comparison:\n")
for (v in list(TRUE, 0.9, "nosuchvar")) {
  e <- tryCatch(posterior::subset_draws(arr_all, variable = v),
                error = function(c) c)
  cat(sprintf("  variable = %-11s %s\n",
              paste(deparse(v), collapse = ""),
              sub("\n.*$", "", conditionMessage(e))))
}
cat("\ndoes the package still hand-roll a pars check for these two?\n")
bodies <- paste(c(
  paste(deparse(getS3method("rhat", "frmtmb_draws")), collapse = " "),
  paste(deparse(getS3method("neff_ratio", "frmtmb_draws")),
        collapse = " "),
  paste(deparse(frmtmb.sample:::draws_diag_array), collapse = " ")),
  collapse = " ")
cat("  mentions draws_extract_pars / draws_select_variables: ",
    grepl("draws_extract_pars|draws_select_variables", bodies), "\n")
cat("  mentions 'must be NA or a character vector':          ",
    grepl("must be NA or a character vector", bodies), "\n")
cat("  mentions posterior::subset_draws:                     ",
    grepl("subset_draws", bodies), "\n")

## the removal: does dropping variable= and fixed= break a call brms
## ANSWERS?
cat("\n-- the removal of variable= and fixed= --\n")
cat("brms rhat(x, variable = \"x\"):\n   ")
cat(res(brms:::rhat.brmsfit(shim, variable = "x")), "\n")
cat("here rhat(ds, variable = \"x\"):\n   ")
cat(res(rhat(ds, variable = "x")), "\n")
cat("brms rhat(x, fixed = TRUE):\n   ")
cat(res(brms:::rhat.brmsfit(shim, "x", fixed = TRUE)), "\n")
cat("here rhat(ds, \"x\", fixed = TRUE):\n   ")
cat(res(rhat(ds, "x", fixed = TRUE)), "\n")
cat("brms rhat(x, \"^b\", regex = TRUE):\n   ")
cat(res(brms:::rhat.brmsfit(shim, "^b", regex = TRUE)), "\n")
cat("here rhat(ds, \"^b\", regex = TRUE):\n   ")
cat(res(rhat(ds, "^b", regex = TRUE)), "\n")
cat("brms rhat(x, \"x\", inc_warmup = TRUE):\n   ")
cat(res(brms:::rhat.brmsfit(shim, "x", inc_warmup = FALSE)), "\n")
cat("here rhat(ds, \"x\", inc_warmup = FALSE):\n   ")
cat(res(rhat(ds, "x", inc_warmup = FALSE)), "\n")
cat("DONE\n")

## The identity in both directions is in dev/bmrev-identity.R.
