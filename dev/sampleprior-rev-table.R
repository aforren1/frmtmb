# Reviewer of lane sampleprior, fix round 1: the lane's table script
# (dev/sampleprior-table.R, copied unchanged below except that its
# apply_check() is renamed apply_check_lane()) with an apply check that
# calls frm_sample() itself, tmbstan captured, so the objective is the
# one the sampler gets (non-centered blocks, gev and mo retapes included).
#   Rscript dev/sampleprior-rev-table.R <lib> <out.tsv>
CAPX <- new.env()
cap_sweep <- function(bform, d, fam) {
  real <- tmbstan::tmbstan
  utils::assignInNamespace("tmbstan", function(obj, ...) {
    CAPX$obj <- obj
    stop("REV_CAPTURED")
  }, "tmbstan")
  on.exit(utils::assignInNamespace("tmbstan", real, "tmbstan"))
  CAPX$obj <- NULL
  r <- tryCatch(suppressMessages(suppressWarnings(frmtmb.sample::frm_sample(
    bform, data = d, family = fam, chains = 1, iter = 10, refresh = 0))),
    error = function(e) conditionMessage(e))
  if (is.null(CAPX$obj)) stop(r)
  CAPX$obj
}
apply_check <- function(bform, d, fam, route = "formula") {
  mine <- tryCatch({
    obj <- cap_sweep(bform, d, fam)
    v <- obj$fn(obj$par)
    if (is.finite(v)) "SAMPLER OK" else paste("SAMPLER NONFINITE", v)
  }, error = function(e) paste("SAMPLER ERROR",
                               substr(conditionMessage(e), 1, 100)))
  paste(mine, "|", apply_check_lane(bform, d, fam, route))
}
# Lane sampleprior: frm_sample()'s default priors against brms 2.23.0's
# default_prior(), family by family and dpar by dpar, MEASURED.
#
#   Rscript dev/sampleprior-table.R <lib or ""> <out.tsv>
#
# For every brms family frmtmb fits (and the mixture), three kinds of
# model on one small data set (seed 20261014, n = 60, 6 groups):
#   base:  y ~ x + (1 | g), every other dpar without a formula
#   <dp>:  the same, with `<dp> ~ x` (the prior moves to its Intercept)
# Each side's table is read with its own default_prior(); frmtmb's with
# route = "sample", which lists what frm_sample() applies. Rows are
# matched on class, coef, group, resp, dpar and nlpar. brms's empty
# prior is its flat one and is written "(flat)", as frmtmb writes it.
a <- commandArgs(TRUE)
lib <- a[1]
out <- a[2]
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(frmtmb.eam)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
cat("brms", format(packageVersion("brms")), "\n")
# whether this frmtmb.sample puts brms's dirichlet(1) on mixture weights
HAS_MIX <- exists("mix_theta_nlp", asNamespace("frmtmb.sample"))
cat("mixture theta dirichlet applied:", HAS_MIX, "\n")

set.seed(20261014)
n <- 60
x <- rnorm(n)
g <- factor(rep(1:6, length.out = n))
pos <- function() rgamma(n, 4, 4 / exp(0.5 + 0.3 * x))
cnt <- function() rpois(n, exp(0.7 + 0.3 * x))
p01 <- function() plogis(0.2 + 0.5 * x)
ord <- function() cut(0.5 * x + rlogis(n), c(-Inf, -1, 0, 1, Inf),
                      labels = FALSE)
zi <- function(v, p = 0.2) ifelse(runif(n) < p, 0, v)
bet <- function() rbeta(n, 2 * p01() * 5, 2 * (1 - p01()) * 5)
simplex3 <- function() {
  e <- cbind(0, 0.3 * x, -0.2 + 0.2 * x) + matrix(rnorm(3 * n, 0, 0.3), n)
  y <- exp(e) / rowSums(exp(e))
  colnames(y) <- c("a", "b", "c")
  y
}
gen <- list(
  acat = ord, cratio = ord, sratio = ord, cumulative = ord,
  hurdle_cumulative = function() zi(ord()),
  categorical = function() sample(1:3, n, TRUE),
  asym_laplace = function() 1 + x + rnorm(n),
  zero_inflated_asym_laplace = function() zi(1 + x + rnorm(n)),
  bernoulli = function() rbinom(n, 1, p01()),
  beta = bet,
  binomial = function() rbinom(n, 10, p01()),
  beta_binomial = function() rbinom(n, 10, p01()),
  zero_inflated_beta_binomial = function() zi(rbinom(n, 10, p01())),
  exgaussian = function() 1 + x + rnorm(n) + rexp(n, 2),
  exponential = function() rexp(n, 1 / exp(0.5 + 0.3 * x)),
  gamma = pos, gaussian = function() 1 + x + rnorm(n),
  geometric = cnt, com_poisson = cnt, hurdle_gamma = function() zi(pos()),
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
  xbeta = function() {
    v <- bet()
    v[1:3] <- c(0, 1, 0)
    v
  },
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
  cox = function() rexp(n, exp(0.3 * x)),
  frechet = function() exp(0.5 + 0.3 * x) / gamma(1 - 1 / 5) *
    (-log(runif(n)))^(-1 / 5),
  gen_extreme_value = function() 1 + 0.3 * x +
    ((-log(runif(n)))^(-0.1) - 1) / 0.1,
  discrete_weibull = function() {
    q <- plogis(1 + 0.3 * x)
    ceiling((log(runif(n)) / log(q))^(1 / 1.5)) - 1
  })

# every family's data drawn up front, in a fixed order, so adding a
# family to the list does not move another's response
D <- list()
for (nm in names(gen)) {
  D[[nm]] <- data.frame(x = x, g = g, y = gen[[nm]](), tr = 10L)
}
dm <- data.frame(x = x, g = g, tr = 10L)
dm$Y <- t(sapply(seq_len(n), function(i) stats::rmultinom(1, 10,
                                                          c(0.2, 0.3, 0.5))))
D$multinomial <- dm
ds <- data.frame(x = x, g = g)
ds$Y <- simplex3()
D$dirichlet <- D$dirichlet2 <- D$logistic_normal <- ds
dw <- frmtmb.eam::ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3,
                               bias = 0.5)
dw$x <- x
dw$g <- g
D$wiener <- dw
D$mixture <- data.frame(x = x, g = g, y = c(rnorm(n / 2, -1),
                                            rnorm(n / 2, 2)))

lhs <- function(nm) {
  switch(nm,
    binomial = , zero_inflated_binomial = , beta_binomial = ,
    zero_inflated_beta_binomial = "y | trials(tr)",
    multinomial = "Y | trials(tr)",
    dirichlet = , dirichlet2 = , logistic_normal = "Y",
    wiener = "rt | dec(upper)",
    "y")
}
brms_fam <- function(nm) {
  if (identical(nm, "mixture")) {
    return(brms::mixture("gaussian", "gaussian"))
  }
  brms::brmsfamily(nm)
}
frm_fam <- function(nm) {
  if (identical(nm, "mixture")) {
    return(frmtmb::mixture(stats::gaussian(), stats::gaussian()))
  }
  if (identical(nm, "wiener")) return(frmtmb.eam::wiener())
  if (identical(nm, "multinomial")) return(frmtmb::multinomial(K = 3))
  frmtmb::brmsfamily(nm)
}

rows_of <- function(tab) {
  tab <- as.data.frame(tab)
  pr <- ifelse(is.na(tab$prior) | !nzchar(tab$prior), "(flat)", tab$prior)
  lb <- as.character(tab$lb)
  ub <- as.character(tab$ub)
  lb[is.na(lb)] <- ""
  ub[is.na(ub)] <- ""
  data.frame(key = paste(tab$class, tab$coef, tab$group, tab$resp,
                         tab$dpar, tab$nlpar, sep = "|"),
             prior = pr, lb = lb, ub = ub)
}

# variant "base", or a dpar's formula such as "shape ~ x" / "shape ~ 1"
one <- function(nm, variant) {
  fo_mu <- paste(lhs(nm), "~ x + (1 | g)")
  extra <- if (identical(variant, "base")) character(0) else variant
  # frmtmb names com_poisson's dispersion nu where brms names it shape
  extra_f <- if (identical(nm, "com_poisson")) sub("^shape", "nu", extra)
             else extra
  compare(nm, variant, c(fo_mu, extra), c(fo_mu, extra_f), brms_fam(nm),
          frm_fam(nm), D[[nm]])
}

min_resp <- function(d) {
  y <- if (!is.null(d$rt)) d$rt else d$y
  if (is.numeric(y) && !is.matrix(y)) min(y) else NA_real_
}

compare <- function(nm, variant, fb, ff_, fam_b, fam_f, d) {
  bb <- tryCatch({
    bf <- do.call(brms::bf, lapply(fb, stats::as.formula))
    rows_of(suppressWarnings(suppressMessages(
      brms::default_prior(bf, data = d, family = fam_b))))
  }, error = function(e) conditionMessage(e))
  ff <- tryCatch({
    bf <- do.call(frmtmb::bf, lapply(ff_, stats::as.formula))
    rows_of(suppressWarnings(suppressMessages(
      frmtmb::default_prior(bf, data = d, family = fam_f,
                            route = "sample"))))
  }, error = function(e) conditionMessage(e))
  r <- compare_rows(nm, variant, bb, ff, min_resp(d))
  r$applied <- apply_check(do.call(frmtmb::bf,
                                   lapply(ff_, stats::as.formula)),
                           d, fam_f)
  r
}

# Whether the sampler APPLIES the defaults default_prior() lists (fix
# round 1, after a listed default the sampler refused): the formula
# route's assembly, the defaults resolved, and the objective
# frm_sample() would hand tmbstan built from them; then each default
# on its own resolved, so that one that resolves to nothing is counted
# rather than passing in silence.
apply_check_lane <- function(bform, d, fam, route = "formula") {
  ns <- asNamespace("frmtmb.sample")
  out <- tryCatch({
    fit <- suppressWarnings(suppressMessages(if (route == "formula") {
      ns$sample_assemble(
        bform, d, fam, data2 = list(), start = NULL,
        fit_control = frmtmb::frmtmb_control(), na.action = stats::na.omit,
        REML = FALSE)
    } else {
      frmtmb::frm(bform, data = d, family = fam)
    }))
    rp <- suppressMessages(ns$sample_resolve_priors(fit, NULL))
    fit$prior <- rp$effective
    obj <- ns$prior_augmented_obj(fit, rp$ri$entries %||% list())
    v <- obj$fn(obj$par)
    defs <- unclass(ns$default_priors_for(fit))
    none <- 0L
    for (s in defs) {
      one <- structure(list(s), class = "frmtmb_priorlist")
      ri <- frmtmb::resolve_prior_input(fit, one)
      if (!length(ri$entries) && is.na(s$lb) && is.na(s$ub)) {
        none <- none + 1L
      }
    }
    sprintf("%s: %d defaults, %d resolve to nothing, objective %s",
            if (none == 0L && is.finite(v)) "applied" else "NOT APPLIED",
            length(defs), none, format(v, digits = 6))
  }, error = function(e) paste("NOT APPLIED: ERROR",
                               substr(conditionMessage(e), 1, 100)))
  # a scalar GEV xi cannot be sampled from a formula at all, on
  # rellib-r9 too (dev/sampleprior-gevform.R): its chart's start reads
  # frm_linpred() of a fit. Its defaults are checked on the fit route
  if (route == "formula" && grepl("needs a fitted model", out)) {
    out <- paste(apply_check_lane(bform, d, fam, route = "fit"),
                 "(fit route; the formula route fails on rellib-r9 too)")
  }
  out
}

compare_rows <- function(nm, variant, bb, ff, miny = NA_real_) {
  if (is.character(bb) || is.character(ff)) {
    return(data.frame(family = nm, variant = variant, key = "(error)",
                      brms = if (is.character(bb)) substr(bb, 1, 120) else "ok",
                      brms_lb = "", brms_ub = "",
                      frmtmb = if (is.character(ff)) substr(ff, 1, 120)
                               else "ok",
                      frmtmb_lb = "", frmtmb_ub = "", status = "error"))
  }
  # frmtmb names com_poisson's dispersion nu where brms names it shape
  if (identical(nm, "com_poisson")) {
    ff$key <- vapply(strsplit(ff$key, "|", fixed = TRUE), function(k) {
      k <- c(k, rep("", 6L - length(k)))
      k[k == "nu"] <- "shape"
      paste(k, collapse = "|")
    }, "")
  }
  keys <- union(bb$key, ff$key)
  bm <- match(keys, bb$key)
  fm <- match(keys, ff$key)
  st <- ifelse(is.na(fm), "frmtmb lacks row",
         ifelse(is.na(bm), "brms lacks row",
          ifelse(bb$prior[bm] == ff$prior[fm], "same", "differs")))
  # brms leaves a group's or a coefficient's row empty where the class
  # row above it speaks for it (its "vectorized" rows); frmtmb's table
  # repeats the class row's density there. Same density applied
  kp <- lapply(strsplit(keys, "|", fixed = TRUE), function(k) {
    c(k, rep("", 6L - length(k)))
  })
  cls_key <- vapply(kp, function(k) {
    paste(k[1], "", "", k[4], k[5], k[6], sep = "|")
  }, "")
  sub_row <- vapply(kp, function(k) nzchar(k[2]) || nzchar(k[3]), NA)
  inh <- st == "differs" & sub_row & bb$prior[bm] %in% "(flat)" &
    ff$prior[fm] %in% ff$prior[match(cls_key, ff$key)]
  st[inh] <- "same (brms inherits)"
  # the same density written two ways: brms keeps the digits as typed
  # (normal(5.0, 0.8)) and its symbol min_Y, which is the smallest
  # response; compared as numbers, min_Y as measured on the data
  same_num <- mapply(function(a, b) {
    if (is.na(a) || is.na(b)) return(FALSE)
    a <- gsub("min_Y", sprintf("%.17g", miny), a, fixed = TRUE)
    pa <- regmatches(a, regexec("^([a-z_]+)[(](.*)[)]$", a))[[1]]
    pb <- regmatches(b, regexec("^([a-z_]+)[(](.*)[)]$", b))[[1]]
    if (length(pa) != 3L || length(pb) != 3L || pa[2] != pb[2]) {
      return(FALSE)
    }
    na <- suppressWarnings(as.numeric(strsplit(pa[3], ",")[[1]]))
    nb <- suppressWarnings(as.numeric(strsplit(pb[3], ",")[[1]]))
    length(na) == length(nb) && !anyNA(c(na, nb)) &&
      all(abs(na - nb) <= 1e-14 * pmax(1, abs(na)))
  }, bb$prior[bm], ff$prior[fm])
  st[st == "differs" & same_num] <- "same (as numbers)"
  # brms's dirichlet(1) on a mixture's weights is applied outside the
  # table: no frmtmb class names it (mix_theta_nlp())
  st[st %in% c("differs", "frmtmb lacks row") & grepl("^theta[|]", keys) &
       bb$prior[bm] %in% "dirichlet(1)"] <-
    if (HAS_MIX) "applied, not tabled (theta)" else "differs"
  # the reference weight of a predicted mixture: brms lists a row that
  # applies nothing (core's brms_prior_class_refusal())
  st[st == "frmtmb lacks row" & grepl("^theta[0-9]+[|]", keys)] <-
    "brms row applies nothing (reference weight)"
  # a smooth's sd: brms's class sds, frmtmb's class sd with the smooth's
  # label as group, the same density
  st[grepl("^sds[|]", keys) & st == "frmtmb lacks row"] <-
    "frmtmb spells it sd (smooth)"
  st[st == "brms lacks row" & grepl("^sd[|]", keys) &
       any(grepl("^sds[|]", keys))] <- "frmtmb spells it sd (smooth)"
  # brms estimates a Student-t latent's df; frmtmb fixes it (dist_nu)
  st[st == "frmtmb lacks row" & grepl("^df[|]", keys)] <-
    "no parameter (frmtmb fixes the df)"
  # frmtmb's raw internal class, which brms has no counterpart for
  st[st == "brms lacks row" & grepl("^theta[|]", keys)] <-
    "frmtmb-only class"
  # brms's per-coefficient sd rows: lane prilabel's
  st[st == "frmtmb lacks row" & grepl("^sd[|][^|]+[|]", keys) &
       bb$prior[bm] %in% "(flat)"] <- "prilabel (per-coef sd row)"
  data.frame(family = nm, variant = variant, key = keys,
             brms = bb$prior[bm], brms_lb = bb$lb[bm], brms_ub = bb$ub[bm],
             frmtmb = ff$prior[fm], frmtmb_lb = ff$lb[fm],
             frmtmb_ub = ff$ub[fm], status = st)
}

fams <- c(sort(names(gen)), "multinomial", "dirichlet", "dirichlet2",
          "logistic_normal", "wiener", "mixture")
all <- list()
report <- function(r, nm, v) {
  if (is.null(r$applied)) r$applied <- ""
  all[[length(all) + 1L]] <<- r
  cat(sprintf("%-28s %-22s %s
", nm, v, r$applied[1L]))
  cat(sprintf("%-28s %-22s rows %3d  same %3d  differs %3d  lacks %3d",
              nm, v, nrow(r), sum(startsWith(r$status, "same")),
              sum(r$status == "differs"),
              sum(r$status == "frmtmb lacks row")),
      if (any(r$status == "error")) " ERROR", "\n")
}
for (nm in fams) {
  dps <- tryCatch(setdiff(brms:::valid_dpars(brms_fam(nm)), "mu"),
                  error = function(e) character(0))
  # a categorical / simplex family's per-category mu<k> is a location
  dps <- dps[!grepl("^mu", dps)]
  # a mixture predicts its weights only all but one at a time
  if (identical(nm, "mixture")) dps <- setdiff(dps, "theta1")
  vs <- if (length(dps)) c(paste(dps, "~ x"), paste(dps, "~ 1"))
  for (v in c("base", vs)) {
    report(one(nm, v), nm, v)
  }
}

# Model terms with classes of their own, on one family each
dx <- D$gaussian
dx$t <- rep(1:10, each = 6)
dx$xo <- factor(sample(1:4, n, TRUE), ordered = TRUE)
dpos <- D$gamma
dmv <- data.frame(x = x, g = g, y1 = D$gaussian$y, y2 = D$student$y)
dthr <- D$cumulative
# drawn after every family's data, so no earlier response moves
dbig <- data.frame(x = x, g = g, y = 50 * D$gamma$y)
dwide <- data.frame(x = x, g = g, y = exp(rnorm(n, 0, 5)))
extras <- list(
  smooth = list("y ~ s(x)", "gaussian", dx),
  student_re = list("y ~ x + (1 | gr(g, dist = 'student'))", "gaussian", dx),
  cor = list("y ~ x + (x | g)", "gaussian", dx),
  ar = list("y ~ x + ar(time = t, gr = g)", "gaussian", dx),
  unstr = list("y ~ x + unstr(time = t, gr = g)", "gaussian", dx),
  mo = list("y ~ mo(xo)", "gaussian", dx),
  gaussian_log_link = list("y ~ x", c("gaussian", "log"), dpos),
  student_pos_sigma = list(c("y ~ x", "sigma ~ x"), "student", dpos),
  lognormal_sigma = list(c("y ~ x", "sigma ~ x"), "lognormal", dpos),
  equidistant = list("y ~ x", c("cumulative", "logit", "equidistant"), dthr),
  mixture3 = list("y ~ 1", "mixture3", D$mixture),
  mv_student_rescor = list("mvbind(y1, y2) ~ x", "mv_student", dmv),
  # sigma's scale read on the identity scale with a log mu link, and on
  # the log scale of a wide positive response under a sigma formula
  gaussian_log_link_wide = list("y ~ x", c("gaussian", "log"), dbig),
  student_wide_sigma = list(c("y ~ x", "sigma ~ x"), "student", dwide),
  # a mixture component's location on its own link, and no family hook
  mixture_poisson = list("y ~ x", "mixpois", D$poisson),
  mixture_negbinomial = list("y ~ x", "mixnb", D$negbinomial),
  mixture_cumulative = list("y ~ x", "mixcum", dthr))
efam <- function(spec, side) {
  if (identical(spec[1], "mixture3")) {
    return(if (side == "b") brms::mixture("gaussian", "gaussian", "gaussian")
           else frmtmb::mixture(stats::gaussian(), stats::gaussian(),
                                stats::gaussian()))
  }
  if (identical(spec[1], "mixpois")) {
    return(if (side == "b") brms::mixture("poisson", "poisson")
           else frmtmb::mixture(stats::poisson(), stats::poisson()))
  }
  if (identical(spec[1], "mixnb")) {
    return(if (side == "b") brms::mixture("negbinomial", "negbinomial")
           else frmtmb::mixture(frmtmb::brmsfamily("negbinomial"),
                                frmtmb::brmsfamily("negbinomial")))
  }
  if (identical(spec[1], "mixcum")) {
    return(if (side == "b") brms::mixture("cumulative", "cumulative")
           else frmtmb::mixture(frmtmb::cumulative(), frmtmb::cumulative()))
  }
  if (identical(spec[1], "mv_student")) {
    return(if (side == "b") brms::student() else frmtmb::student())
  }
  if (length(spec) == 3L) {
    return(if (side == "b") brms::cumulative(threshold = spec[3])
           else frmtmb::cumulative(threshold = spec[3]))
  }
  if (length(spec) == 2L) {
    return(if (side == "b") brms::brmsfamily(spec[1], link = spec[2])
           else frmtmb::brmsfamily(spec[1], link = spec[2]))
  }
  if (side == "b") brms::brmsfamily(spec) else frmtmb::brmsfamily(spec)
}
for (nm in names(extras)) {
  e <- extras[[nm]]
  fb <- e[[1]]
  r <- if (identical(nm, "mv_student_rescor")) {
    # set_rescor() is a formula component on both sides
    bb <- tryCatch(rows_of(suppressMessages(brms::default_prior(
      brms::bf(stats::as.formula(fb)) + brms::set_rescor(TRUE),
      data = e[[3]], family = brms::student()))),
      error = function(err) conditionMessage(err))
    ff <- tryCatch(rows_of(suppressMessages(frmtmb::default_prior(
      frmtmb::bf(stats::as.formula(fb)) + frmtmb::set_rescor(TRUE),
      data = e[[3]], family = frmtmb::student(), route = "sample"))),
      error = function(err) conditionMessage(err))
    rr <- compare_rows(nm, "extra", bb, ff)
    rr$applied <- apply_check(frmtmb::bf(stats::as.formula(fb)) +
                                frmtmb::set_rescor(TRUE), e[[3]],
                              frmtmb::student())
    rr
  } else {
    compare(nm, "extra", fb, fb, efam(e[[2]], "b"), efam(e[[2]], "f"),
            e[[3]])
  }
  report(r, nm, "extra")
}
res <- do.call(rbind, all)
utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(res), "rows to", out, "\n")
