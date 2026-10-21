# Reviewer: does each new test file see a defect put back? One mutation
# of the lane build per process, then one test file.
#   Rscript dev/nb2dm-rev-mut.R <mutation> <test file>
# Mutations:
#   none       the lane build as installed (control)
#   mnfam      multinomial_families back to "multinomial" only, so every
#              fam_is_multinomial() switch the lane added is undone
#   noeta      negbinomial2's .eta_shape translation removed
#   nopost     fam_reparam() leaves post$* untranslated
#   plain      lgamma_shift_diff() replaced by the plain difference
#   nostart    start_ls_rules$negbinomial2 removed
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
a <- commandArgs(trailingOnly = TRUE)
mut <- a[1]
f <- a[2]
ns <- asNamespace("frmtmb")
put <- function(name, value) {
  utils::assignInNamespace(name, value, "frmtmb")
}
switch(mut,
  none = NULL,
  mnfam = put("multinomial_families", "multinomial"),
  noeta = {
    g <- get("fam_negbinomial2", ns)
    body(g) <- do.call(substitute, list(body(g), list(
      robust_logmu = quote(function(...) NULL))))
    put("fam_negbinomial2", g)
  },
  nopost = {
    # the top-level slots translated, post$* handed over as the inner
    # family has them (a first spelling that edited the deparsed source
    # matched nothing, because deparse() splits that line)
    orig <- get("fam_reparam", ns)
    g <- function(inner, family, dpars, links, to_inner, init_dpars) {
      out <- orig(inner, family, dpars, links, to_inner, init_dpars)
      out[["post"]] <- inner[["post"]]
      out
    }
    put("fam_reparam", g)
  },
  plain = put("lgamma_shift_diff", function(a, s) lgamma(a + s) - lgamma(a)),
  # the probability-scale Dirichlet draw of the first build
  olddraw = put("simplex_rdirichlet", function(alpha, cats) {
    G <- matrix(stats::rgamma(length(alpha), shape = alpha), nrow(alpha))
    out <- G / rowSums(G)
    colnames(out) <- cats
    out
  }),
  # negbinomial2's response check back to negbinomial's name
  oldvalid = {
    # brmsfamily() reads the constructor captured in family_registry,
    # so the registry entry is what has to change
    g <- get("fam_negbinomial2", ns)
    reg <- get("family_registry", ns)
    reg[["negbinomial2"]] <- function(...) {
      f <- g(...)
      f[["valid_y"]] <- get("count_y", ns)("negbinomial")
      f
    }
    put("family_registry", reg)
  },
  # the fitted() count clamp for dirichlet_multinomial undone
  noclamp = put("fam_is_multinomial", function(fam) {
    identical(fam[["family"]], "multinomial")
  }),
  nostart = {
    r <- get("start_ls_rules", ns)
    r$negbinomial2 <- NULL
    put("start_ls_rules", r)
  },
  stop("unknown mutation ", mut))
p <- if (grepl("frmtmb.sample", f)) "frmtmb.sample" else "frmtmb"
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
r <- as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT mut=", mut, " ", basename(f), " pass=", sum(r$passed),
    " fail=", sum(r$failed), " err=", sum(r$error), " skip=",
    sum(r$skipped), " warn=", sum(r$warning), " blocks=", nrow(r), "\n",
    sep = "")
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("FAILING BLOCK:", bad, sep = "\n  ")
cat("\n")
