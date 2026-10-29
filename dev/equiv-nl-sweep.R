# The nl x family rows of frm_compat(), measured family by family.
#
#   Rscript dev/equiv-nl-sweep.R [outdir] [--record]
#
# The table refuses nl = TRUE for every core family but gaussian (rule
# nl x kind:family). This fits each core family three ways on one data
# set: the predictor written linearly (y ~ x), the same predictor as a
# nonlinear body (y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), and the
# intercept-only model. A body frm() honors has the linear fit's
# maximum; a body it drops has the intercept-only one. brms 2.23.0's
# make_stancode() is run on the same nl spelling, so a row says whether
# brms takes it. Seed 20261029, n = 300.
#
# Gamma and inverse.gaussian are fitted under a log link, and ALSO under
# their default links (inverse, 1/mu^2): from frm()'s default starts
# (0 for a nonlinear coefficient) and from starts that keep mu positive.
#
# Writes <outdir>/nl-sweep.tsv (default dev/equiv-out), or
# dev/equiv-data/nl-sweep.tsv with --record; equiv-summary.R reads the
# latter into the proposals.
args <- commandArgs(trailingOnly = TRUE)
record <- "--record" %in% args
args <- setdiff(args, "--record")
outdir <- if (length(args) >= 1) args[1] else "dev/equiv-out"
if (!record) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- if (record) "dev/equiv-data/nl-sweep.tsv" else
  file.path(outdir, "nl-sweep.tsv")
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r15")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

ft <- frm_compat_features()
fams <- unique(ft$name[ft$kind == "family"])
cp <- frm_compat("nl")
cp <- cp[cp$kind_b == "family", ]
fams <- intersect(fams, cp$feature_b)

set.seed(20261029)
n <- 300
x <- round(rnorm(n), 3)
eta <- 0.3 + 0.5 * x
pos <- exp(0.4 + 0.3 * x + rnorm(n, 0, 0.4))
cnt <- rpois(n, exp(eta))
prob <- plogis(eta)
ord <- cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = 1:4)
simplex <- function() {
  a <- cbind(2, 2 * exp(0.4 * x), 2 * exp(-0.3 + 0.2 * x))
  s <- t(apply(a, 1, function(r) {
    v <- rgamma(3, r)
    v / sum(v)
  }))
  colnames(s) <- c("a", "b", "c")
  s
}
counts3 <- function() {
  m <- t(vapply(seq_len(n), function(i) {
    stats::rmultinom(1, 10, c(1, exp(0.4 * x[i]), 1))[, 1]
  }, numeric(3)))
  colnames(m) <- c("a", "b", "c")
  m
}
zero <- runif(n) < 0.2
ydat <- function(f) {
  switch(f,
    gaussian = , student = , skew_normal = , asym_laplace = ,
    gen_extreme_value = eta + rnorm(n),
    exgaussian = eta + rnorm(n, 0, 0.5) + rexp(n, 2),
    lognormal = , Gamma = , weibull = , exponential = ,
    inverse.gaussian = , frechet = pos,
    shifted_lognormal = 0.1 + pos,
    beta = , xbeta = plogis(eta + rnorm(n, 0, 0.5)),
    zero_inflated_beta = ifelse(zero, 0, plogis(eta + rnorm(n, 0, 0.5))),
    zero_one_inflated_beta = ifelse(zero, rbinom(n, 1, 0.5),
                                    plogis(eta + rnorm(n, 0, 0.5))),
    tweedie = ifelse(zero, 0, pos),
    hurdle_gamma = , hurdle_lognormal = ifelse(zero, 0, pos),
    poisson = , negbinomial = , negbinomial2 = , nbinom1 = , geometric = ,
    compois = , discrete_weibull = cnt,
    zero_inflated_poisson = , zero_inflated_negbinomial = ifelse(zero, 0L,
                                                                 cnt),
    hurdle_poisson = , hurdle_negbinomial = ifelse(zero, 0L, cnt + 1L),
    binomial = , beta_binomial = rbinom(n, 10, prob),
    zero_inflated_binomial = ,
    zero_inflated_beta_binomial = ifelse(zero, 0L, rbinom(n, 10, prob)),
    bernoulli = rbinom(n, 1, prob),
    cumulative = , sratio = , cratio = , acat = ordered(ord),
    hurdle_cumulative = ordered(ifelse(zero, 0L, as.integer(ord)),
                                levels = 0:4),
    categorical = factor(sample(c("a", "b", "c"), n, TRUE)),
    multinomial = , dirichlet_multinomial = counts3(),
    dirichlet = , dirichlet2 = , logistic_normal = simplex(),
    von_mises = (eta + rnorm(n, 0, 0.5) + pi) %% (2 * pi) - pi,
    NULL)
}
lhs <- function(f) {
  if (f %in% c("binomial", "beta_binomial", "zero_inflated_binomial",
               "zero_inflated_beta_binomial", "multinomial",
               "dirichlet_multinomial")) "y | trials(ntr)" else "y"
}
brms_name <- c(Gamma = "gamma", compois = "com_poisson")
log_link <- c("Gamma", "inverse.gaussian")
frm_family <- function(f) {
  g <- tryCatch(get(f, envir = asNamespace("frmtmb")), error = function(e) {
    tryCatch(get(f, envir = asNamespace("stats")), error = function(e) NULL)
  })
  if (identical(f, "beta")) g <- frmtmb::Beta
  # a log link where the default (inverse, 1/mu^2) is undefined at the
  # nonlinear coefficients' default start of 0, on both sides
  if (f %in% log_link) return(g(link = "log"))
  if (is.function(g)) g() else frmtmb::brmsfamily(f)
}
ll <- function(expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    structure(NA_real_, err = conditionMessage(r))
  } else {
    as.numeric(stats::logLik(r))
  }
}

rows <- list()
for (f in fams) {
  y <- ydat(f)
  if (is.null(y)) {
    rows[[f]] <- data.frame(family = f, status_now = NA, frm_linear = NA,
                            frm_nl = NA, frm_intercept = NA, outcome =
                              "no generator", brms_nl = NA, note = "")
    next
  }
  d <- data.frame(x = x, ntr = 10L)
  d$y <- y
  fam <- tryCatch(frm_family(f), error = function(e) NULL)
  if (is.null(fam)) {
    rows[[f]] <- data.frame(family = f, status_now = NA, frm_linear = NA,
                            frm_nl = NA, frm_intercept = NA, outcome =
                              "no constructor", brms_nl = NA, note = "")
    next
  }
  L <- lhs(f)
  f_lin <- stats::as.formula(paste(L, "~ x"))
  f_nl <- stats::as.formula(paste(L, "~ b0 + b1 * x"))
  f_int <- stats::as.formula(paste(L, "~ 1"))
  a <- ll(frm(frmtmb::bf(f_lin), family = fam, data = d))
  b <- ll(frm(frmtmb::bf(f_nl, b0 + b1 ~ 1, nl = TRUE), family = fam,
              data = d))
  c0 <- ll(frm(frmtmb::bf(f_int), family = fam, data = d))
  tol <- function(v) 1e-4 * max(1, abs(v))
  outcome <- if (is.na(b)) {
    "frm refuses nl"
  } else if (is.na(a)) {
    "linear fit fails"
  } else if (abs(b - a) < tol(a)) {
    "nl honored"
  } else if (!is.na(c0) && abs(b - c0) < tol(c0)) {
    "nl body DROPPED"
  } else {
    "nl differs from both"
  }
  bn <- if (f %in% names(brms_name)) brms_name[[f]] else f
  bb <- tryCatch({
    bfam <- if (f %in% log_link) brms::brmsfamily(bn, link = "log") else
      brms::brmsfamily(bn)
    suppressMessages(brms::make_stancode(
      brms::bf(f_nl, b0 + b1 ~ 1, nl = TRUE), data = d, family = bfam))
    "builds"
  }, error = function(e) paste("refuses:", substr(conditionMessage(e), 1,
                                                  80)))
  note <- paste(stats::na.omit(c(attr(a, "err"), attr(b, "err"))),
                collapse = " | ")
  # the default link: the inverse link (Gamma) and 1/mu^2
  # (inverse.gaussian) are undefined at a linear predictor of 0, where
  # frm() starts a nonlinear coefficient; a start of eta = 1/mean(y)
  # (or 1/mean(y)^2) keeps mu positive at the first evaluation
  dl <- list(lin = NA_real_, nl0 = NA_real_, nl0_err = "", nls = NA_real_,
             outcome = "")
  if (f %in% log_link) {
    g <- get(f, envir = asNamespace("stats"))
    e0 <- if (identical(f, "Gamma")) 1 / mean(d$y) else 1 / mean(d$y)^2
    st <- list(beta = c("b0_(Intercept)" = e0, "b1_(Intercept)" = 0))
    dl$lin <- ll(frm(frmtmb::bf(f_lin), family = g(), data = d))
    n0 <- ll(frm(frmtmb::bf(f_nl, b0 + b1 ~ 1, nl = TRUE), family = g(),
                 data = d))
    dl$nl0 <- n0
    dl$nl0_err <- substr(attr(n0, "err") %||% "", 1, 120)
    dl$nls <- ll(frm(frmtmb::bf(f_nl, b0 + b1 ~ 1, nl = TRUE), family = g(),
                     data = d, start = st))
    dl$outcome <- if (!is.na(dl$nls) && !is.na(dl$lin) &&
                        abs(dl$nls - dl$lin) < tol(dl$lin)) {
      "default link: nl honored from positive starts"
    } else {
      "default link: nl not honored from positive starts"
    }
  }
  rows[[f]] <- data.frame(
    family = f, status_now = cp$status[cp$feature_b == f],
    frm_linear = round(a, 4), frm_nl = round(b, 4),
    frm_intercept = round(c0, 4), outcome = outcome, brms_nl = bb,
    note = gsub("[\t\r\n]+", " ", substr(note, 1, 200)),
    default_linear = round(dl$lin, 4), default_nl_default_start =
      round(dl$nl0, 4), default_nl_default_start_err =
      gsub("[\t\r\n]+", " ", dl$nl0_err), default_nl_positive_start =
      round(dl$nls, 4), default_outcome = dl$outcome)
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
res$proposal <- ifelse(
  res$status_now == "refused" & res$outcome == "nl honored",
  ifelse(grepl("^builds", res$brms_nl), "works", "works, no brms reference"),
  ifelse(res$outcome == "nl body DROPPED", "broken",
         ifelse(res$status_now == "refused" &
                  res$outcome == "nl differs from both", "untested", "")))
# Gamma and inverse.gaussian: the honored body above is the log link's;
# the default link decides the condition
dlk <- res$family %in% log_link & res$status_now == "refused"
res$proposal[dlk] <- ifelse(
  res$default_outcome[dlk] == "default link: nl honored from positive starts",
  "conditional (default link needs starting values that keep mu positive)",
  "conditional (log link)")
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = TRUE,
                   qmethod = "double")
print(res[, c("family", "status_now", "frm_linear", "frm_nl",
              "frm_intercept", "outcome", "brms_nl", "proposal")],
      row.names = FALSE)
cat("\n")
print(table(res$outcome, res$status_now, useNA = "ifany"))
print(res[res$family %in% log_link, c("family", "default_linear",
                                      "default_nl_default_start",
                                      "default_nl_positive_start",
                                      "default_outcome", "proposal")],
      row.names = FALSE)
cat("wrote", out, "\n")
