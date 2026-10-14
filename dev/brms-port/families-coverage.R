# brms_families.Rmd is a parameterization reference with no runnable
# chunks, so it is audited by family: every family the document names,
# looked up by name and then fitted by frm() on simulated data.
#
#   PORT_LIB=<lib> Rscript families-coverage.R
#
# Two counts. "core" is the 0.34.0 probe unchanged, so it compares
# across builds. "core plus extensions" also counts a family that an
# installed frmtmb extension exports as a constructor: up to 0.70.0 the
# probe read core only and scored `wiener` missing although frmtmb.eam
# exports wiener(). The column `covered_by` names the package. Both are
# NAME-level; every named family is then FITTED (below), and the
# FITTED counts are the coverage this instrument reports.
HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
source(file.path(HERE, "env.R"))
suppressMessages(library(frmtmb))
cat("#", port_build(), "\n")
src <- readLines(file.path(system.file("doc", package = "brms"),
                           "brms_families.Rmd"), warn = FALSE)
bold <- unlist(regmatches(src, gregexpr("\\*\\*[A-Za-z_0-9.]+\\*\\*", src)))
bold <- unique(gsub("\\*", "", bold))

catalogue <- tryCatch({
  f <- getFromNamespace("family_names", "brms")
  sort(unique(f()))
}, error = function(e) character())

# The probe the 0.34.0 audit used, kept so the counts compare.
accepted <- function(nm) {
  tryCatch({
    frmtmb:::as_frmtmb_family(nm)
    TRUE
  }, error = function(e) FALSE)
}

# An extension covers a family when it exports a function of that name
# whose default call frm() accepts as a family. Every installed frmtmb.*
# package is probed, so a new extension is counted without editing this.
ext_pkgs <- unique(sort(grep("^frmtmb[.]", rownames(utils::installed.packages(
  lib.loc = .libPaths())), value = TRUE)))
ext_covers <- function(nm) {
  for (p in ext_pkgs) {
    ns <- tryCatch(suppressMessages(loadNamespace(p)),
                   error = function(e) NULL)
    if (is.null(ns) || !nm %in% getNamespaceExports(p)) next
    ok <- tryCatch({
      f <- getExportedValue(p, nm)
      is.function(f) && inherits(frmtmb:::as_frmtmb_family(f()),
                                 "frmtmb_family")
    }, error = function(e) FALSE)
    if (ok) return(p)
  }
  NA_character_
}

# Keep only tokens that look like family names (drop bolded prose)
cand <- bold[nchar(bold) > 2]
cand <- setdiff(cand, c("mu", "sigma", "not", "Note", "brms", "hurdle"))
reg <- names(frmtmb:::family_registry)
# multinomial is registered but needs its category count, so the
# bare-name probe fails where frm(family = multinomial(K)) works
core <- vapply(cand, function(n) accepted(n) || n %in% reg, logical(1))
by <- ifelse(core, "frmtmb", vapply(cand, ext_covers, character(1)))
res <- data.frame(family = cand, core = unname(core),
                  covered_by = ifelse(is.na(by), "-", unname(by)),
                  row.names = NULL)
res <- res[order(!res$core, res$covered_by == "-", res$family), ]

# A name is not a fit. Each covered family is FITTED, y ~ x on data
# simulated from a model of its kind (seed 20261013, n = 300; the
# generators are the review's, dev/remeasure-rev-famfit.R), the family
# passed by brms's name through frmtmb::brmsfamily(). A family counts as
# covered only when frm() returns a fit with a finite log-likelihood.
# multinomial is fitted through multinomial(K = 3): frmtmb needs the
# category count, which brms reads from the response.
set.seed(20261013)
n <- 300
x <- rnorm(n)
pos <- function() rgamma(n, 4, 4 / exp(0.5 + 0.3 * x))
cnt <- function() rpois(n, exp(0.7 + 0.3 * x))
p01 <- function() plogis(0.2 + 0.5 * x)
ord <- function() cut(0.5 * x + rlogis(n), c(-Inf, -1, 0, 1, Inf),
                      labels = FALSE)
zi <- function(v, p = 0.2) ifelse(runif(n) < p, 0, v)
bet <- function() rbeta(n, 2 * p01() * 5, 2 * (1 - p01()) * 5)
gen <- list(
  acat = ord, cratio = ord, sratio = ord, cumulative = ord,
  categorical = function() sample(1:3, n, TRUE),
  asym_laplace = function() 1 + x + rnorm(n),
  bernoulli = function() rbinom(n, 1, p01()),
  Beta = bet,
  binomial = function() rbinom(n, 10, p01()),
  exgaussian = function() 1 + x + rnorm(n) + rexp(n, 2),
  exponential = function() rexp(n, 1 / exp(0.5 + 0.3 * x)),
  Gamma = pos, gaussian = function() 1 + x + rnorm(n),
  geometric = cnt, hurdle_gamma = function() zi(pos()),
  hurdle_lognormal = function() zi(exp(rnorm(n, 0.3 * x))),
  hurdle_negbinomial = function() zi(pmax(1, cnt())),
  hurdle_poisson = function() zi(pmax(1, cnt())),
  inverse.gaussian = pos, lognormal = function() exp(rnorm(n, 0.3 * x)),
  negbinomial = cnt, poisson = cnt,
  shifted_lognormal = function() 0.2 + exp(rnorm(n, 0.3 * x)),
  skew_normal = function() 1 + x + rnorm(n),
  student = function() 1 + x + rt(n, 5),
  von_mises = function() {
    v <- 0.3 * x + rnorm(n, 0, 0.5)
    atan2(sin(v), cos(v))
  },
  weibull = function() rweibull(n, 2, exp(0.5 + 0.3 * x)),
  zero_inflated_beta = function() zi(bet()),
  zero_inflated_binomial = function() zi(rbinom(n, 10, p01())),
  zero_inflated_negbinomial = function() zi(cnt()),
  zero_inflated_poisson = function() zi(cnt()),
  zero_one_inflated_beta = function() {
    v <- bet()
    u <- runif(n)
    v[u < 0.1] <- 0
    v[u > 0.9] <- 1
    v
  },
  cox = function() rexp(n, exp(0.3 * x)))
fit_one <- function(fam_name, pkg) {
  out <- tryCatch({
    if (identical(fam_name, "multinomial")) {
      d <- data.frame(x = x, tr = 10L)
      d$Y <- t(sapply(seq_len(n), function(i) {
        stats::rmultinom(1, 10, c(0.2, 0.3, 0.5))
      }))
      f <- suppressWarnings(frm(bf(Y | trials(tr) ~ x),
                                family = multinomial(K = 3), data = d))
    } else if (identical(fam_name, "wiener")) {
      d <- getExportedValue(pkg, "ddm_simulate")(n, mu = 0.8, bs = 1.4,
                                                 ndt = 0.3, bias = 0.5)
      f <- suppressWarnings(frm(bf(rt | dec(upper) ~ 1),
                                family = getExportedValue(pkg, fam_name)(),
                                data = d))
    } else {
      if (is.null(gen[[fam_name]])) stop("no data generator")
      d <- data.frame(x = x, y = gen[[fam_name]](), tr = 10L)
      fo <- if (fam_name %in% c("binomial", "zero_inflated_binomial")) {
        bf(y | trials(tr) ~ x)
      } else {
        bf(y ~ x)
      }
      f <- suppressWarnings(frm(fo, family = frmtmb::brmsfamily(fam_name),
                                data = d))
    }
    ll <- as.numeric(logLik(f))
    if (!is.finite(ll)) stop("non-finite log-likelihood")
    sprintf("fit, code %d, logLik %.3f", f$opt$convergence, ll)
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
  out
}
res$fit <- vapply(seq_len(nrow(res)), function(i) {
  if (res$covered_by[i] == "-") return("-")
  fit_one(res$family[i], res$covered_by[i])
}, "")
res$fits <- startsWith(res$fit, "fit,")
print(res, row.names = FALSE)
cat("\nextensions probed:",
    if (length(ext_pkgs)) paste(sprintf("%s %s", ext_pkgs,
      vapply(ext_pkgs, function(p) format(utils::packageVersion(p)), "")),
      collapse = ", ") else "none", "\n")
cat("named by core (the 0.34.0 name probe):", sum(res$core), "of",
    nrow(res), "\n")
cat("named by core plus extensions:", sum(res$covered_by != "-"), "of",
    nrow(res), "\n")
cat("FITTED by core:", sum(res$fits & res$core), "of", nrow(res), "\n")
cat("FITTED by core plus extensions:", sum(res$fits), "of", nrow(res), "\n")
cat("named but not fitted:",
    paste(res$family[res$covered_by != "-" & !res$fits], collapse = ", "),
    "\n")
cat("missing from core:", paste(res$family[!res$core], collapse = ", "), "\n")
cat("missing from core plus extensions:",
    paste(res$family[res$covered_by == "-"], collapse = ", "), "\n")
cat("\nfrmtmb registry:\n")
cat(paste(sort(unique(reg)), collapse = ", "), "\n")
if (length(catalogue)) {
  cat("\nbrms catalogue (", length(catalogue), " names) not in the frmtmb",
      " registry:\n", sep = "")
  cat(paste(setdiff(catalogue, reg), collapse = ", "), "\n")
}
