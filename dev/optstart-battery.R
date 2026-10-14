# Lane optstart: a battery over families, designs and seeds. Per fit: the
# convergence code, the warnings and the log likelihood against a
# reference: lme4 (gaussian, and lognormal through log y), glmmTMB
# (student, Gamma, beta, negbinomial) or survival::survreg (weibull,
# no group term); for families without one (skew_normal, exgaussian,
# huber, asym_laplace) the reference is the best of this build's fit and
# a fit started at the true values. Run once per library; compare with
# dev/optstart-battery-sum.R.
#   Rscript dev/optstart-battery.R <lib or ""> <label> [families] [seeds]
# designs: ordinary (slope 0.5 in units of the residual scale, group sd
# half of it); hir2 (a slope 1e3 times the residual scale on the
# identity families, 3 on the log or logit scale); offset (an intercept
# 1e5 times it, or 10 on the log scale, 3 on the logit scale); rebig (a
# group sd 30 times the residual scale, or 2 on the log scale); offterm
# (an offset() term that varies across rows: up to 1e3 times the
# residual scale, a log exposure from 0 to log 1000, or sd 2 on the
# logit scale).
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb); library(lme4); library(glmmTMB); library(survival)
})
lab <- a[2]
fams <- if (length(a) > 2 && nzchar(a[3])) strsplit(a[3], ",")[[1]] else
  c("gaussian", "student", "lognormal", "skew_normal", "exgaussian",
    "huber", "asym_laplace", "Gamma", "weibull", "beta", "negbinomial")
seeds <- if (length(a) > 3) eval(parse(text = a[4])) else 1:10
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
reg <- frmtmb:::family_registry
spec_of <- function(fam) {
  switch(fam,
    gaussian = list(f = stats::gaussian(), link = "id", sc = c(sigma = 1)),
    student = list(f = frmtmb::student(), link = "id", sc = c(sigma = 1),
                   extra = list(nu = 5)),
    lognormal = list(f = frmtmb::lognormal(), link = "id",
                     sc = c(sigma = 0.4)),
    skew_normal = list(f = frmtmb::skew_normal(), link = "id",
                       sc = c(sigma = 1),
                       extra = list(alpha = 3)),
    exgaussian = list(f = frmtmb::exgaussian(), link = "id",
                      sc = c(sigma = 1.5, beta = 1.5)),
    huber = list(f = frmtmb::huber(), link = "id", sc = c(sigma = 1)),
    # quantile held at 0.5, the median regression the start is for; an
    # estimated quantile runs to a degenerate end of its range
    asym_laplace = list(f = frmtmb::asym_laplace(), link = "id",
                        sc = c(sigma = 1),
                        extra = list(quantile = 0.5)),
    Gamma = list(f = Gamma(link = "log"), link = "log", sc = c(shape = 4)),
    weibull = list(f = frmtmb::weibull(), link = "log", sc = c(shape = 2.5)),
    beta = list(f = frmtmb::Beta(), link = "logit", sc = c(phi = 15)),
    negbinomial = list(f = frmtmb::negbinomial(), link = "log",
                       sc = c(shape = 3)))
}
inv <- list(id = identity, log = exp, logit = stats::plogis)
design <- function(fam, des, s) {
  sp <- spec_of(fam)
  set.seed(s)
  n <- 400
  G <- 20
  g <- factor(rep(seq_len(G), length.out = n))
  x <- rnorm(n)
  # the residual scale on the scale of mu: sigma for the identity
  # families (sd of log y for lognormal), and 1 on the link scale
  rs <- if (sp$link == "id") sp$sc[[1]] else 1
  b0 <- 0.3
  b1 <- 0.5 * rs
  sdg <- 0.5 * rs
  if (sp$link != "id") sdg <- 0.3
  if (des == "hir2") b1 <- if (sp$link == "id") 1e3 * rs else 3
  if (des == "offset") b0 <- if (sp$link == "id") 1e5 * rs else
    if (sp$link == "log") 10 else 3
  if (des == "rebig") sdg <- if (sp$link == "id") 30 * rs else 2
  # weibull's reference (survreg) has no group term
  if (fam == "weibull") sdg <- 0
  # lognormal's mu is log y, so the multiples above overflow y; there
  # the slope is 50 and the intercept 300 times the residual scale
  if (fam == "lognormal" && des == "hir2") b1 <- 50 * rs
  if (fam == "lognormal" && des == "offset") b0 <- 300 * rs
  bg <- rnorm(G, 0, sdg)
  # offterm: a known offset that varies across rows and moves mu far
  # from the residual spread unless the offset is read
  o <- if (des != "offterm") rep(0, n) else if (sp$link == "id") {
    1e3 * rs * runif(n)
  } else if (sp$link == "log") log(runif(n, 1, 1000)) else rnorm(n, 0, 2)
  eta <- o + b0 + b1 * x + bg[g]
  mu <- inv[[sp$link]](eta)
  sim <- reg[[fam]]()$sim
  dpl <- c(list(mu = mu), lapply(c(as.list(sp$sc), sp$extra), rep, n))
  y <- if (fam == "huber") mu + rnorm(n, 0, sp$sc[[1]]) else
    sim(dpl, list(), n)
  if (fam == "beta") y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  d <- data.frame(y = y, x = x, g = g, o = o)
  fo <- if (fam == "weibull") y ~ x else y ~ x + (1 | g)
  if (des == "offterm") fo <- update(fo, . ~ . + offset(o))
  truth <- list(beta = c(`(Intercept)` = b0, x = b1),
                betad = stats::setNames(log(sp$sc),
                                        paste0(names(sp$sc), "_(Intercept)")),
                theta = if (fam != "weibull") log(sdg))
  # the truth start in the family's own convention: lognormal's mu is on
  # the log scale already; weibull's and Gamma's mu is the mean
  list(d = d, fo = fo, sp = sp, truth = truth)
}
ref_ll <- function(fam, m) {
  d <- m$d
  fo <- m$fo
  r <- tryCatch(suppressWarnings(suppressMessages(switch(fam,
    gaussian = as.numeric(logLik(lmer(fo, data = d, REML = FALSE,
      control = lmerControl(check.conv.singular = "ignore")))),
    lognormal = {
      d$ly <- log(d$y)
      as.numeric(logLik(lmer(update(fo, ly ~ .), data = d, REML = FALSE,
        control = lmerControl(check.conv.singular = "ignore")))) -
        sum(log(d$y))
    },
    student = as.numeric(logLik(glmmTMB(fo, data = d,
                                        family = t_family()))),
    Gamma = as.numeric(logLik(glmmTMB(fo, data = d,
                                      family = Gamma(link = "log")))),
    beta = as.numeric(logLik(glmmTMB(fo, data = d,
                                     family = beta_family()))),
    negbinomial = as.numeric(logLik(glmmTMB(fo, data = d,
                                            family = nbinom2()))),
    weibull = as.numeric(logLik(survreg(update(fo, Surv(y) ~ .), data = d,
                                        dist = "weibull"))),
    NA_real_))), error = function(e) NA_real_)
  r
}
one <- function(m, start = NULL) {
  w <- character()
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(suppressMessages(
    frm(if (identical(m$sp$f$family, "asym_laplace")) {
      bf(m$fo, quantile = 0.5)
    } else m$fo, data = m$d, family = m$sp$f, start = start)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  el <- proc.time()[[3]] - t0
  if (inherits(f, "error")) {
    return(data.frame(code = NA, ll = NA, nwarn = length(w), secs = el,
                      warn = substr(paste("ERROR", conditionMessage(f)), 1,
                                    80)))
  }
  data.frame(code = f$opt$convergence, ll = as.numeric(logLik(f)),
             nwarn = length(w), secs = el,
             warn = substr(gsub("[\t\n]", " ",
                                paste(unique(substr(w, 1, 40)),
                                      collapse = " | ")), 1, 120))
}
out <- list()
for (fam in fams) for (des in c("ordinary", "hir2", "offset", "rebig",
                                "offterm")) {
  for (s in seeds) {
    m <- design(fam, des, s)
    r <- one(m)
    tr <- tryCatch(one(m, start = Filter(Negate(is.null), m$truth)),
                   error = function(e) NULL)
    out[[length(out) + 1L]] <- data.frame(
      label = lab, family = fam, design = des, seed = s,
      code = r$code, ll = r$ll, nwarn = r$nwarn, secs = r$secs,
      ll_ref = ref_ll(fam, m), ll_truth = if (is.null(tr)) NA else tr$ll,
      warn = r$warn)
  }
  cat(fam, des, "done\n")
}
X <- do.call(rbind, out)
dir.create("dev/optstart-log/battery", showWarnings = FALSE)
write.table(X, sprintf("dev/optstart-log/battery/%s-%s.tsv", lab,
                       paste(fams, collapse = "_")),
            sep = "\t", quote = FALSE, row.names = FALSE)
