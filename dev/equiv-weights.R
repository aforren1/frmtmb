# Reach weights for the equivalence survey.
#
# "Common" is read off the code users copy: the examples and vignettes
# of brms 2.23.0 (the tarball under dev/brms-suite/) and frmtmb's own
# vignettes. One document is one vignette's code, or one Rd file's
# \examples section. A feature's count d is the number of documents
# whose code uses it; its weight is 1 + d, so a feature no document
# shows still has a nonzero chance to be drawn. A pair's weight is the
# product of its two features' weights, which assumes a user reaches
# the two features independently.
#
# Usage: Rscript dev/equiv-weights.R [outdir]
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "dev/equiv-data"
brms_dir <- "dev/brms-suite/brms"

rmd_code <- function(path) {
  x <- readLines(path, warn = FALSE)
  on <- FALSE
  keep <- character()
  for (ln in x) {
    if (!on && grepl("^\\s*```\\s*[{]r", ln)) {
      on <- TRUE
      next
    }
    if (on && grepl("^\\s*```\\s*$", ln)) {
      on <- FALSE
      next
    }
    if (on) keep <- c(keep, ln)
  }
  keep
}
ltx_code <- function(path) {
  x <- readLines(path, warn = FALSE)
  on <- FALSE
  keep <- character()
  for (ln in x) {
    if (grepl("\\\\begin[{](Sinput|Code|CodeInput)[}]", ln)) {
      on <- TRUE
      next
    }
    if (grepl("\\\\end[{](Sinput|Code|CodeInput)[}]", ln)) {
      on <- FALSE
      next
    }
    if (on) keep <- c(keep, ln)
  }
  # the overview's comparison table writes calls in \code{} cells
  cells <- regmatches(x, gregexpr("\\\\code[{][^}]*[}]", x))
  keep <- c(keep, gsub("\\$\\\\sim\\$", "~", unlist(cells)))
  keep
}
rd_code <- function(path) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f))
  ok <- try(tools::Rd2ex(path, f, commentDontrun = FALSE,
                            commentDonttest = FALSE), silent = TRUE)
  if (inherits(ok, "try-error") || !file.exists(f)) return(character())
  readLines(f, warn = FALSE)
}
strip_comments <- function(x) sub("#.*$", "", x)

docs <- list()
for (p in list.files(file.path(brms_dir, "vignettes"), "[.]Rmd$",
                     full.names = TRUE)) {
  docs[[paste0("brms:vignette:", basename(p))]] <- rmd_code(p)
}
for (p in list.files(file.path(brms_dir, "vignettes"), "[.]ltx$",
                     full.names = TRUE)) {
  docs[[paste0("brms:vignette:", basename(p))]] <- ltx_code(p)
}
for (p in list.files(file.path(brms_dir, "man"), "[.]Rd$",
                     full.names = TRUE)) {
  cd <- rd_code(p)
  if (length(cd)) docs[[paste0("brms:Rd:", basename(p))]] <- cd
}
for (p in list.files("vignettes", "[.]Rmd$", full.names = TRUE)) {
  docs[[paste0("frmtmb:vignette:", basename(p))]] <- rmd_code(p)
}
docs <- lapply(docs, function(x) paste(strip_comments(x), collapse = "\n"))

# The model calls inside a document, with balanced parentheses, so that
# a call that names no family can count for gaussian, brms's default.
model_calls <- function(txt) {
  st <- gregexpr("\\b(brm|frm|brm_multiple)\\(", txt, perl = TRUE)[[1]]
  if (st[1] < 0) return(character())
  ch <- strsplit(txt, "")[[1]]
  vapply(st, function(s) {
    i <- s + attr(st, "match.length")[match(s, st)] - 1L
    depth <- 0L
    j <- i
    while (j <= length(ch)) {
      if (ch[j] == "(") depth <- depth + 1L
      if (ch[j] == ")") {
        depth <- depth - 1L
        if (depth == 0L) break
      }
      j <- j + 1L
    }
    paste(ch[s:min(j, length(ch))], collapse = "")
  }, "")
}

L <- "(?<![A-Za-z0-9_.])"
pat <- c(
  # families
  gaussian = paste0(L, "gaussian(?![A-Za-z_])"),
  student = paste0(L, "student(?![A-Za-z_])"),
  lognormal = paste0(L, "lognormal(?![A-Za-z_])"),
  shifted_lognormal = "shifted_lognormal",
  skew_normal = "skew_normal",
  exgaussian = "exgaussian",
  asym_laplace = "asym_laplace",
  Gamma = paste0(L, "Gamma(?![A-Za-z_])"),
  weibull = paste0(L, "weibull(?![A-Za-z_])"),
  exponential = "family\\s*=\\s*(brmsfamily\\(\")?exponential",
  inverse.gaussian = "inverse[.]gaussian",
  beta = paste0(L, "(Beta|beta)\\(\\s*(link|\\))"),
  tweedie = "tweedie",
  poisson = paste0(L, "poisson(?![A-Za-z_])"),
  negbinomial = paste0(L, "negbinomial(?![A-Za-z_2])"),
  negbinomial2 = paste0(L, "negbinomial2"),
  nbinom1 = "nbinom1",
  geometric = paste0(L, "geometric(?![A-Za-z_])"),
  compois = "com_poisson|compois",
  binomial = paste0(L, "binomial(?![A-Za-z_])"),
  bernoulli = paste0(L, "bernoulli(?![A-Za-z_])"),
  beta_binomial = paste0(L, "beta_binomial"),
  zero_inflated_beta_binomial = "zero_inflated_beta_binomial",
  multinomial = paste0(L, "multinomial(?![A-Za-z_])"),
  dirichlet_multinomial = "dirichlet_multinomial",
  zero_inflated_poisson = "zero_inflated_poisson",
  zero_inflated_negbinomial = "zero_inflated_negbinomial",
  zero_inflated_binomial = "zero_inflated_binomial",
  zero_inflated_beta = "zero_inflated_beta(?![A-Za-z_])",
  zero_one_inflated_beta = "zero_one_inflated_beta",
  xbeta = paste0(L, "xbeta"),
  hurdle_poisson = "hurdle_poisson",
  hurdle_negbinomial = "hurdle_negbinomial",
  hurdle_gamma = "hurdle_gamma",
  hurdle_lognormal = "hurdle_lognormal",
  cumulative = paste0(L, "cumulative(?![A-Za-z_])"),
  hurdle_cumulative = "hurdle_cumulative",
  sratio = paste0(L, "sratio"),
  cratio = paste0(L, "cratio"),
  acat = paste0(L, "acat"),
  categorical = paste0(L, "categorical\\("),
  von_mises = "von_mises",
  dirichlet = "family\\s*=\\s*(brmsfamily\\(\")?dirichlet(?![A-Za-z_2])",
  dirichlet2 = "dirichlet2",
  logistic_normal = "logistic_normal",
  frechet = "frechet",
  gen_extreme_value = "gen_extreme_value",
  discrete_weibull = "discrete_weibull",
  # covariance structures; us is any single-bar group-level term
  us = "\\([^()|]+[^|]\\|[^|][^()]*\\)",
  diag = "\\|\\||cor\\s*=\\s*FALSE|(?<![A-Za-z_])diag\\(",
  homdiag = "homdiag\\(",
  cs = paste0(L, "cs\\([^)]*\\|"),
  ar1 = paste0(L, "ar1\\("),
  hetar1 = "hetar1\\(",
  ou = paste0(L, "ou\\("),
  toep = paste0(L, "toep\\("),
  homtoep = "homtoep\\(",
  homcs = "homcs\\(",
  exp = paste0(L, "exp\\([^)]*\\|"),
  gau = paste0(L, "gau\\("),
  mat = paste0(L, "mat\\([^)]*\\|"),
  rr = paste0(L, "rr\\("),
  equalto = "equalto\\(",
  gr_cov = "gr\\([^)]*cov\\s*=",
  gr_prec = "gr\\([^)]*prec\\s*=",
  smooth = paste0(L, "(s|t2)\\("),
  gp = paste0(L, "gp\\("),
  hsgp = "gp\\([^)]*\\bk\\s*=",
  us_t = "dist\\s*=\\s*\"student\"",
  diag_t = paste0("dist\\s*=\\s*\"student\"[^)]*cor\\s*=\\s*FALSE|",
                  "cor\\s*=\\s*FALSE[^)]*dist\\s*=\\s*\"student\""),
  # addition terms
  `weights()` = "\\|\\s*weights\\(",
  `trials()` = "trials\\(",
  `cens()` = "cens\\(",
  `trunc()` = "\\|[^~]*trunc\\(",
  `se()` = "\\|[^~]*(?<![A-Za-z_])se\\(",
  `mi()` = "\\|[^~]*(?<![A-Za-z_])mi\\(",
  `vint()` = "vint\\(",
  `vreal()` = "vreal\\(",
  `thres()` = "thres\\(",
  `rate()` = "\\|[^~]*rate\\(",
  `subset()` = "\\|[^~]*subset\\(",
  `index()` = "\\|[^~]*index\\(",
  # predictor specials
  `s()` = paste0(L, "s\\("),
  `t2()` = paste0(L, "t2\\("),
  `mo()` = paste0(L, "mo\\("),
  `mi_pred()` = "~[^\\n]*(?<![A-Za-z_])mi\\(",
  `gp_pred()` = paste0(L, "gp\\("),
  `cs_pred()` = paste0(L, "cs\\([^|)]*\\)"),
  `ps()` = paste0(L, "ps\\("),
  `me()` = paste0(L, "me\\("),
  # residual correlation
  `ar()` = paste0(L, "ar\\("),
  `ma()` = paste0(L, "ma\\("),
  `arma()` = paste0(L, "arma\\("),
  `cosy()` = paste0(L, "cosy\\("),
  `unstr()` = paste0(L, "unstr\\("),
  `fcor()` = paste0(L, "fcor\\("),
  # estimation modes and options
  REML = "REML",
  quadrature = "quadrature",
  importance = "importance\\s*=",
  profile = "profile\\s*=\\s*TRUE",
  autoscale = "autoscale",
  sparse_x = "sparse_x|sparse\\s*=\\s*TRUE",
  prior = "prior\\s*=|set_prior\\(|prior\\(",
  bounds = "\\b(lb|ub)\\s*=|bounds\\s*=",
  verbose = "verbose",
  # structures
  mvbf = "mvbf\\(|mvbind\\(|bf\\([^\\n]*\\)\\s*\\+\\s*bf\\(",
  rescor = "rescor",
  `|ID|` = "\\|\\s*[A-Za-z0-9_.]+\\s*\\|",
  nl = "nl\\s*=\\s*TRUE|nlf\\(",
  mixture = "mixture\\(",
  disc = paste0(L, "disc\\s*~"),
  equidistant = "equidistant",
  sum_to_zero = "sum_to_zero",
  # methods
  fitted = paste0(L, "(fitted|posterior_epred|posterior_linpred)\\("),
  predict = paste0(L, "predict\\("),
  simulate = paste0(L, "(simulate|posterior_predict|pp_check)\\("),
  residuals = paste0(L, "(residuals|resid)\\("),
  residuals_osa = "osa",
  emmeans = "emmeans|emtrends",
  confint_profile = paste0(L, "confint\\("),
  hypothesis_profile = paste0(L, "hypothesis\\("),
  frm_lp_basis = "frm_lp_basis",
  frm_curve = "frm_curve\\(|conditional_effects\\(",
  frm_curve_contrast = "frm_curve_contrast",
  # grammar
  # the refused crossing of a whole bar term, x * (1 | g) or (1 | g):x;
  # nested grouping inside the bar, (1 | g1/g2), is not this feature
  bar_crossing = paste0("[*:]\\s*\\([^()|]*\\|[^()]*\\)|",
                        "\\([^()|]*\\|[^()]*\\)\\s*[*:]"),
  call_group = "\\|\\s*(factor|interaction)\\(",
  double_bar = "\\|\\|",
  `mm()` = paste0(L, "mm\\("),
  `mmc()` = "mmc\\(",
  `0 + Intercept` = "0\\s*\\+\\s*Intercept",
  gr_by = "gr\\([^)]*by\\s*="
)

cp <- read.delim("dev/equiv-frame-out/compat.tsv")
core <- cp[cp$in_core_only, ]
feats <- unique(c(core$feature_a, core$feature_b))
miss <- setdiff(feats, names(pat))
if (length(miss)) stop("no pattern for: ", paste(miss, collapse = ", "))

hit <- function(p, txt) grepl(p, txt, perl = TRUE)
src <- sub(":.*$", "", names(docs))
res <- data.frame(feature = feats, d_brms = 0L, d_frmtmb = 0L)
for (i in seq_along(feats)) {
  h <- vapply(docs, function(t) hit(pat[[feats[i]]], t), NA)
  if (feats[i] == "gaussian") {
    # brms's default family: a model call that names no family
    h <- h | vapply(docs, function(t) {
      mc <- model_calls(t)
      any(!grepl("family", mc))
    }, NA)
  }
  res$d_brms[i] <- sum(h & src == "brms")
  res$d_frmtmb[i] <- sum(h & src == "frmtmb")
}
res$d <- res$d_brms + res$d_frmtmb
res$weight <- 1 + res$d
res$kind <- cp$kind_a[match(res$feature, cp$feature_a)]
res$kind[is.na(res$kind)] <- cp$kind_b[match(res$feature[is.na(res$kind)],
                                             cp$feature_b)]
res$pattern <- unname(pat[res$feature])
cat("documents:", length(docs), " brms:", sum(src == "brms"),
    " frmtmb:", sum(src == "frmtmb"), "\n")
res <- res[order(-res$weight), ]
print(res[, c("feature", "kind", "d_brms", "d_frmtmb", "weight")],
      row.names = FALSE)
# the corpus itself, so its size is a count of rows, not a typed number
utils::write.table(data.frame(document = names(docs), source = src),
                   file.path(out, "corpus.tsv"), sep = "\t",
                   row.names = FALSE, quote = TRUE, qmethod = "double")
utils::write.table(res, file.path(out, "weights.tsv"), sep = "\t",
                   row.names = FALSE, quote = TRUE, qmethod = "double")
