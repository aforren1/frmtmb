## Short check of the lane's fix round 1: F2 (prior row order), F3
## (one slot per coefficient), source equality of the installed build.
## Seed 21.   Rscript dev/review-brmsnames-fr1.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
if (arm == "lane") {
  src <- new.env()
  for (f in c("R/parse.R", "R/priors.R", "R/brms-shapes.R",
              "R/brms-names.R", "R/frame.R")) {
    sys.source(f, envir = src, keep.source = FALSE)
  }
  for (fn in c("unsub_specials", "parse_linpred", "parse_one_response",
               "prior_table", "prior_rows_sd_order", "prior_coef_brms",
               "prior_coef_match_brms", "prior_spec_slot", "sd_spec_reach",
               "sd_coef_shared", "resolve_priorlist", "brms_lp_prefix",
               "brms_fixef_rows", "check_cs_in_bar")) {
    cat("installed == source:", fn,
        identical(deparse(get(fn, ns)), deparse(get(fn, src))), "\n")
  }
  s2 <- new.env(parent = asNamespace("frmtmb.sample"))
  sys.source("extensions/frmtmb.sample/R/methods-draws.R", envir = s2,
             keep.source = FALSE)
  cat("installed == source: fixef.frmtmb_draws",
      identical(deparse(get("fixef.frmtmb_draws",
                            asNamespace("frmtmb.sample"))),
                deparse(get("fixef.frmtmb_draws", s2))), "\n")
}
set.seed(21)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), xp = runif(n, 0.5, 2),
                k = factor(sample(1:3, n, TRUE)), k0 = sample(0:3, n, TRUE),
                g = factor(rep(1:12, 20)), t = rep(1:20, each = 12),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$y <- 1 + d$x + 0.3 * d$z^2 + 0.5 * log(d$xp) + rnorm(12)[d$g] + rnorm(n)
d$y2 <- d$x + rnorm(n)
d$yo <- cut(d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)

## F3: two spellings of one coefficient, and different coefficients
cat("\n== F3\n")
tv <- function(lab, p, f, ...) {
  v <- tryCatch(q(validate_prior(p, f, data = d, ...)),
                error = function(e) e)
  cat(sprintf("%-44s ", lab))
  if (inherits(v, "error")) {
    cat("REFUSED:", substr(conditionMessage(v), 1, 150), "\n")
  } else {
    v <- as.data.frame(v)
    u <- v[v$source == "user", ]
    cat("accepted; user rows:", paste0("[", u$class, "|", u$coef, "|",
                                       u$prior, "]", collapse = " "), "\n")
  }
  invisible(v)
}
f1 <- y ~ I(x^2) + log(xp) + poly(z, 2) + factor(k) + (1 + I(x^2) | g)
tv("b: I(x^2) and IxE2", c(set_prior("normal(0, 1)", coef = "I(x^2)"),
                           set_prior("normal(0, 2)", coef = "IxE2")), f1)
tv("b: poly(z, 2)1 and polyz21",
   c(set_prior("normal(0, 1)", coef = "poly(z, 2)1"),
     set_prior("normal(0, 2)", coef = "polyz21")), f1)
tv("b: log(xp) and logxp", c(set_prior("normal(0, 1)", coef = "log(xp)"),
                             set_prior("normal(0, 2)", coef = "logxp")), f1)
tv("sd: I(x^2) and IxE2",
   c(set_prior("normal(0, 1)", class = "sd", coef = "I(x^2)", group = "g"),
     set_prior("normal(0, 2)", class = "sd", coef = "IxE2", group = "g")), f1)
tv("Intercept: (Intercept) and Intercept, 0+",
   c(set_prior("normal(0, 1)", coef = "(Intercept)"),
     set_prior("normal(0, 2)", coef = "Intercept")),
   bf(y ~ 0 + Intercept + x))
fs <- y ~ s(x) + z
tv("smooth: s(x).fx1 and sx_1", c(set_prior("normal(3, 0.001)",
                                             coef = "s(x).fx1"),
                                   set_prior("normal(-3, 0.001)",
                                             coef = "sx_1")), fs)
p5 <- c(set_prior("normal(1, 0.001)", coef = "logxp"),
        set_prior("normal(2, 0.001)", coef = "IxE2"),
        set_prior("normal(3, 0.001)", coef = "poly(z, 2)1"),
        set_prior("normal(4, 0.001)", coef = "polyz22"),
        set_prior("normal(5, 0.001)", coef = "factor(k)2"),
        set_prior("normal(6, 0.001)", coef = "factork3"))
tv("six distinct coefs, mixed spellings", p5, f1)
fit <- tryCatch(q(frm(f1, data = d, prior = p5)), error = function(e) e)
if (inherits(fit, "error")) {
  cat("fit with six priors:", conditionMessage(fit), "\n")
} else {
  fe <- fixef(fit)[, 1]
  cat("fit with six priors, estimates:", paste(names(fe), round(fe, 3),
                                               collapse = "; "), "\n")
}
fsm <- tryCatch(q(frm(fs, data = d, prior = c(
  set_prior("normal(3, 0.001)", coef = "s(x).fx1"),
  set_prior("normal(-3, 0.001)", coef = "sx_1")))), error = function(e) e)
cat("smooth fit, two spellings:", if (inherits(fsm, "error"))
  substr(conditionMessage(fsm), 1, 150) else
    paste(rownames(fixef(fsm)), round(fixef(fsm)[, 1], 3), collapse = "; "),
  "\n")
# a random-effect design whose columns rename to one name
fz <- tryCatch(q(frm(y ~ x + (1 + log(xp) + logxp | g),
                     data = transform(d, logxp = rnorm(n)))),
               error = function(e) e)
cat("Z columns log(xp) and logxp:", if (inherits(fz, "error"))
  substr(conditionMessage(fz), 1, 150) else
    paste(rownames(VarCorr(fz)$g$sd), collapse = " "), "\n")
bz <- tryCatch(q(brms::standata(y ~ x + (1 + log(xp) + logxp | g),
                                data = transform(d, logxp = rnorm(n)))),
               error = function(e) e)
cat("brms, same:", if (inherits(bz, "error"))
  substr(conditionMessage(bz), 1, 150) else "accepted", "\n")

## F2: full rows, saved for the cross-arm comparison, and the order
## against brms
cat("\n== F2\n")
models <- list(
  list(f = "y ~ exp(x) + z + (1 + z | g)", fam = "gaussian"),
  list(f = "bf(y ~ x + f, sigma ~ z)", fam = "gaussian"),
  list(f = paste("mvbf(bf(y ~ x + (1 | g)), bf(y2 ~ log(xp) + (1 | g)),",
                 "rescor = FALSE)"), fam = "gaussian"),
  list(f = "bf(y ~ a + b * x, a ~ exp(z) + (1 | g), b ~ 1 + z, nl = TRUE)",
       fam = "gaussian"),
  list(f = "y ~ s(x) + mo(k0) + (1 | g)", fam = "gaussian"),
  list(f = "bf(yo ~ x + cs(z))", fam = "acat"),
  list(f = "y ~ x + ar(time = t, gr = g)", fam = "gaussian"),
  list(f = "y ~ x + I(z^2) + (1 | g)", fam = "student")
)
up <- c(set_prior("normal(0, 2)", class = "b"),
        set_prior("exponential(1)", class = "sd"))
key <- function(p) with(p, paste(class, coef, group, resp, dpar, nlpar,
                                 sep = "/"))
out <- list()
for (m in models) {
  fr <- eval(parse(text = m$f))
  fb <- eval(parse(text = gsub("(mvbf|bf)\\(", "brms::\\1(", m$f)))
  fam <- get(m$fam, asNamespace("frmtmb"))()
  dp <- as.data.frame(q(default_prior(fr, data = d, family = fam)))
  has_sd <- any(dp$class == "sd")
  pr <- if (has_sd) up else up[1]
  vp <- tryCatch(as.data.frame(q(validate_prior(pr, fr, data = d,
                                                family = fam))),
                 error = function(e) conditionMessage(e))
  out[[m$f]] <- list(dp = dp, vp = vp)
  bd <- as.data.frame(q(brms::default_prior(fb, data = d, family = m$fam)))
  bpr <- q(brms::set_prior("normal(0, 2)", class = "b"))
  if (has_sd) bpr <- bpr + q(brms::set_prior("exponential(1)",
                                             class = "sd"))
  bv <- tryCatch(as.data.frame(q(brms::validate_prior(bpr, fb, data = d,
                                                      family = m$fam))),
                 error = function(e) conditionMessage(e))
  same <- function(a, b) {
    if (!is.data.frame(a) || !is.data.frame(b)) return(NA)
    ka <- key(a); kb <- key(b)
    identical(ka[ka %in% kb], kb[kb %in% ka])
  }
  cat(sprintf("%-70s default_prior order = brms: %s; validate_prior: %s\n",
              m$f, same(dp, bd), same(vp, bv)))
}
saveRDS(out, sprintf("dev/review-brmsnames-out/fr1-tables-%s.rds", arm))

## the lncor rows of two logistic_normal() responses
src <- readLines("tests/testthat/test-simplex-mvbf.R")
i <- grep("^smv_data <- function", src)
j <- i + which(src[(i + 1):length(src)] == "}")[1L]
eval(parse(text = src[i:j]))
dd <- smv_data(2)
bl <- as.data.frame(q(brms::default_prior(
  brms::mvbf(brms::bf(Y ~ x, family = brms::logistic_normal()),
             brms::bf(V ~ x, family = brms::logistic_normal())), data = dd)))
fl <- as.data.frame(q(default_prior(
  mvbf(bf(Y ~ x, family = logistic_normal()),
       bf(V ~ x, family = logistic_normal())), data = dd)))
cat("lncor resp order, brms:", bl$resp[bl$class == "lncor"],
    "; frmtmb:", fl$resp[fl$class == "lncor"], "\n")
cat("whole table order = brms on shared rows:",
    identical(key(fl)[key(fl) %in% key(bl)], key(bl)[key(bl) %in% key(fl)]),
    "\n")
