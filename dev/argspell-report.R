## The generated half of dev/argspell-findings.md. Everything numeric in
## that document comes out of here and is pasted verbatim, because a
## count that is TYPED is a count that can be wrong without any file
## being wrong.
##
## Run AFTER installing the change:  Rscript dev/argspell-report.R
LIB <- "C:/Users/adf44/source/r/argspell-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkgs <- c("frmtmb", "frmtmb.sample")
for (p in pkgs) suppressMessages(library(p, character.only = TRUE))
suppressMessages(library(brms))

cat("frmtmb       :", dirname(system.file(package = "frmtmb")), "\n")
cat("frmtmb.sample:", dirname(system.file(package = "frmtmb.sample")), "\n")
cat("brms         :", as.character(utils::packageVersion("brms")), "\n\n")

# The PARSE TREE, not the deparsed text: three dots inside a string
# literal are not a use of the dots, and reading the text scored
# print.frmtmb_par_template() as a dots user on the strength of
# cat("... ", length(v) - n, " more\n"). Mirrors ar_dots_names in
# tests/testthat/test-arg-refusal.R.
dots_names <- c("...", "...length", "...names", "...elt")
is_dots_name <- function(x) {
  x %in% dots_names | grepl("^[.][.][0-9]+$", x)
}

# A method whose whole body is an unconditional refusal cannot swallow
# anything, so it is exempt by SHAPE rather than by name. Mirrors
# ar_refuses_always() in tests/testthat/test-arg-refusal.R.
refusers <- c("stop", "fit_no_draws", "multiple_no_draws")
refuses_always <- function(fn) {
  b <- body(fn)
  head_is <- function(e) {
    is.call(e) && as.character(e[[1L]])[1L] %in% refusers
  }
  if (!is.call(b) || !identical(as.character(b[[1L]]), "{")) return(head_is(b))
  length(b) == 2L && head_is(b[[2L]])
}

methods_of <- function(p) {
  reg <- parseNamespaceFile(p, dirname(system.file(package = p)))$S3methods
  nm <- ifelse(is.na(reg[, 3]), paste(reg[, 1], reg[, 2], sep = "."),
               reg[, 3])
  gen <- reg[, 1]
  u <- !duplicated(nm)
  data.frame(pkg = p, fun = nm[u], generic = gen[u],
             stringsAsFactors = FALSE)
}
tab <- do.call(rbind, lapply(pkgs, methods_of))
tab$formals <- NA_character_
tab$has_dots <- FALSE
tab$reads_dots <- FALSE
for (i in seq_len(nrow(tab))) {
  o <- tryCatch(get(tab$fun[i], envir = asNamespace(tab$pkg[i])),
                error = function(e) NULL)
  if (is.null(o)) next
  fo <- names(formals(o))
  tab$formals[i] <- paste(fo, collapse = ", ")
  tab$has_dots[i] <- "..." %in% fo
  tab$reads_dots[i] <- tab$has_dots[i] &&
    (any(is_dots_name(all.names(body(o)))) || refuses_always(o))
}
# Measured, not assumed: see dev/argspell-exempt2.R. Seven, not the
# thirteen the first pass carried.
exempt <- c(
  "emm_basis.frmtmb_fit", "get_predict.frmtmb_fit",
  "get_vcov.frmtmb_fit", "get_varcov.frmtmb_fit",
  "get_parameters.frmtmb_fit", "find_formula.frmtmb_fit",
  "find_random.frmtmb_fit",
  "get_coef.frmtmb_fit", "set_coef.frmtmb_fit")
tab$exempt <- tab$fun %in% exempt

brms_method <- function(gen) {
  for (bc in c("brmsfit", "brmsformula", "brmsprior", "mvbrmsformula")) {
    o <- tryCatch(get(paste0(gen, ".", bc), envir = asNamespace("brms")),
                  error = function(e) NULL)
    if (!is.null(o)) return(o)
  }
  NULL
}
tab$brms_formals <- vapply(tab$generic, function(g) {
  o <- brms_method(g)
  if (is.null(o)) NA_character_ else paste(names(formals(o)), collapse = ", ")
}, "")

cat("---- GENERATED: method counts ----\n")
for (p in pkgs) {
  t <- tab[tab$pkg == p, ]
  cat(sprintf("%s\n", p))
  cat(sprintf("  registered S3 methods (deduped) : %d\n", nrow(t)))
  cat(sprintf("  taking `...`                    : %d\n", sum(t$has_dots)))
  cat(sprintf("  with `...` never read           : %d\n",
              sum(t$has_dots & !t$reads_dots)))
  cat(sprintf("  of those, exempt by contract    : %d\n",
              sum(t$has_dots & !t$reads_dots & t$exempt)))
  cat(sprintf("  with a brms counterpart         : %d\n",
              sum(!is.na(t$brms_formals))))
}
bad <- tab[tab$has_dots & !tab$reads_dots & !tab$exempt, ]
cat(sprintf("\nMETHODS STILL SWALLOWING THEIR DOTS: %d\n", nrow(bad)))
if (nrow(bad)) cat(paste0("  ", bad$pkg, " ", bad$fun, collapse = "\n"), "\n")
cat(sprintf("frm_check_dots(...) call sites: frmtmb %d, frmtmb.sample %d\n",
            length(grep("frm_check_dots[(][.][.][.]",
                        unlist(lapply(list.files("R", pattern = "[.]R$",
                                                 full.names = TRUE),
                                      readLines, warn = FALSE)))),
            length(grep("frm_check_dots[(][.][.][.]",
                        unlist(lapply(list.files("extensions/frmtmb.sample/R",
                                                 pattern = "[.]R$",
                                                 full.names = TRUE),
                                      readLines, warn = FALSE))))))
cat("---- END GENERATED: method counts ----\n\n")

cat("---- GENERATED: re.form / re_formula surface ----\n")
surf <- c("predict.frmtmb_fit", "fitted.frmtmb_fit", "simulate.frmtmb_fit",
          "residuals.frmtmb_fit", "conditional_effects.frmtmb_fit",
          "pp_check.frmtmb_fit")
for (nm in surf) {
  fo <- names(formals(getFromNamespace(nm, "frmtmb")))
  cat(sprintf("%-32s re_formula %-3s re.form %s\n", nm,
              "re_formula" %in% fo, "re.form" %in% fo))
}
for (nm in c("frm_bootstrap", "dharma_residuals")) {
  fo <- names(formals(getFromNamespace(nm, "frmtmb")))
  cat(sprintf("%-32s re_formula %-3s re.form %s\n", nm,
              "re_formula" %in% fo, "re.form" %in% fo))
}
dr <- c("posterior_epred", "posterior_linpred", "posterior_predict",
        "predictive_error", "predictive_interval", "pp_check")
# DECLARES is not ACCEPTS. `predictive_interval.brmsfit` declares
# neither spelling and its whole body is `posterior_predict(object,
# ...)`, so brms honors the alias one frame down. Both columns are
# printed, because reading only the first one got this wrong.
brms_accepts_alias <- function(gen) {
  bm <- brms_method(gen)
  if (is.null(bm)) return(NA)
  if ("re.form" %in% names(formals(bm))) return(TRUE)
  bd <- deparse(body(bm))
  if (!any(grepl("...", bd, fixed = TRUE))) return(FALSE)
  targets <- c("posterior_predict", "posterior_epred", "posterior_linpred")
  hit <- targets[vapply(targets,
                        function(t) any(grepl(t, bd, fixed = TRUE)), TRUE)]
  any(vapply(hit, function(t) {
    m <- brms_method(t)
    !is.null(m) && "re.form" %in% names(formals(m))
  }, TRUE))
}
for (g in dr) {
  fo <- names(formals(getFromNamespace(paste0(g, ".frmtmb_draws"),
                                       "frmtmb.sample")))
  bm <- brms_method(g)
  cat(sprintf("%-32s ours %-5s brms declares %-5s accepts %s\n",
              paste0(g, ".frmtmb_draws"), "re.form" %in% fo,
              if (is.null(bm)) NA else "re.form" %in% names(formals(bm)),
              brms_accepts_alias(g)))
}
cat("---- END GENERATED: re.form / re_formula surface ----\n\n")

cat("---- GENERATED: formals diff against brms, prediction family ----\n")
key <- c("fitted.frmtmb_fit", "predict.frmtmb_fit", "residuals.frmtmb_fit",
         "coef.frmtmb_fit", "fixef.frmtmb_fit", "ranef.frmtmb_fit",
         "VarCorr.frmtmb_fit", "hypothesis.frmtmb_fit",
         "conditional_effects.frmtmb_fit",
         "posterior_epred.frmtmb_draws", "posterior_linpred.frmtmb_draws",
         "posterior_predict.frmtmb_draws", "log_lik.frmtmb_draws",
         "predictive_error.frmtmb_draws", "pp_mixture.frmtmb_draws")
k <- tab[tab$fun %in% key, ]
for (i in seq_len(nrow(k))) {
  a <- strsplit(k$formals[i], ", ")[[1]]
  b <- if (is.na(k$brms_formals[i])) character() else
    strsplit(k$brms_formals[i], ", ")[[1]]
  cat(sprintf("%s\n  ours: %s\n  brms: %s\n  brms-only: %s\n",
              k$fun[i], k$formals[i],
              if (is.na(k$brms_formals[i])) "<no brmsfit method>"
              else k$brms_formals[i],
              paste(setdiff(b, c(a, "object", "x", "...")), collapse = ", ")))
}
cat("---- END GENERATED: formals diff ----\n")
