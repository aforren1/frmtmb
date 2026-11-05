# Blocker 1, measured. Bodies with a saddle at the default start, under ordinal,
# gaussian, poisson and hurdle families, frm() from the default start against
# frm() from a start, the model without the body, and brms's optimizing from 6
# starts under flat priors. Seeds 9202 (the review's data) and as named.
# Rscript dev/nlfix-saddle.R <one spec id or "all"> writes one row per case to
# dev/review-nlfix-saddle-out/<id>.txt.
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(Sys.getenv("NLFIX_LIB",
  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16"),
  ";", fixed = TRUE)[[1]], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                   c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
set.seed(1)
d$yh <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
set.seed(9205)
d$yg <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
d$yc <- rpois(n, 1.5 * exp(0.6 * d$x))
d$yzi <- ifelse(runif(n) < 0.25, 0L, rpois(n, 1.5 * exp(0.6 * d$x)))
d$yhu <- ifelse(runif(n) < 0.25, 0L, 1L + rpois(n, 1.5 * exp(0.6 * d$x)))
d$ypo <- rgamma(n, 3, 3 / (1.5 * exp(0.6 * d$x)))
d$ybi <- rbinom(n, 1, plogis(1.5 * exp(0.6 * d$x) - 2))
bodies <- list(
  exp = "b1 * exp(b2 * x)",
  sat = "b1 * (1 - exp(-b2 * x))",
  logi = "b1 / (1 + exp(b2 * x))",
  prod = "a + b1 * b2 * x")
fams <- list(
  cumulative = c("y", "cumulative()"), probit = c("y", "cumulative('probit')"),
  sratio = c("y", "sratio()"), cratio = c("y", "cratio()"), acat = c("y", "acat()"),
  hurdle_cumulative = c("yh", "hurdle_cumulative()"),
  gaussian = c("yg", "gaussian()"), poisson = c("yc", "poisson()"),
  poisson_id = c("yc", "poisson('identity')"),
  zip = c("yzi", "zero_inflated_poisson()"), hurdle_poisson = c("yhu", "hurdle_poisson()"),
  gamma_log = c("ypo", "Gamma('log')"), bernoulli = c("ybi", "bernoulli()"))
cases <- expand.grid(body = names(bodies), fam = names(fams),
                     stringsAsFactors = FALSE)
cases$id <- paste(cases$fam, cases$body, sep = "_")
mkf <- function(pkg, Y, b, nlp) {
  txt <- paste0("~ ", b)
  fs <- stats::as.formula(paste(Y, txt), env = globalenv())
  lst <- paste(nlp, collapse = " + ")
  if (identical(pkg, "brms")) {
    do.call(brms::bf, list(fs, stats::as.formula(paste(lst, "~ 1")), nl = TRUE))
  } else {
    do.call(frmtmb::bf, list(fs, stats::as.formula(paste(lst, "~ 1")), nl = TRUE))
  }
}
quiet <- function(expr) {
  w <- character(0); m <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") }),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "tterr"))
  attr(v, "w") <- w; attr(v, "m") <- m
  v
}
one <- function(id) {
  cs <- cases[cases$id == id, ]
  Y <- fams[[cs$fam]][1]
  famtxt <- fams[[cs$fam]][2]
  b <- bodies[[cs$body]]
  nlp <- intersect(c("a", "b1", "b2"), all.vars(parse(text = b)))
  f <- mkf("frm", Y, b, nlp)
  fam <- eval(parse(text = famtxt), envir = asNamespace("frmtmb"))
  ll <- function(r) if (inherits(r, "tterr")) NA_real_ else as.numeric(logLik(r))
  t0 <- proc.time()[["elapsed"]]
  r <- quiet(frm(f, family = fam, data = d))
  t_fit <- proc.time()[["elapsed"]] - t0
  st <- rep(0.3, length(nlp))
  names(st) <- paste0(nlp, "_(Intercept)")
  rs <- quiet(frm(f, family = fam, data = d, start = list(beta = st)))
  r0 <- quiet(frm(stats::as.formula(paste(Y, "~ 1")), family = fam, data = d))
  # brms from several starts, flat priors
  bm <- tryCatch({
    fb <- mkf("brms", Y, b, nlp)
    bfam <- eval(parse(text = famtxt), envir = asNamespace("brms"))
    pr <- E$brms_flat_prior(fb, data = d, family = bfam)
    code <- brms::make_stancode(fb, data = d, family = bfam, prior = pr)
    sdat <- brms::standata(fb, data = d, family = bfam, prior = pr)
    m <- E$brms_stan_model(code)
    vals <- numeric(0)
    for (i in 1:6) {
      o <- tryCatch(rstan::optimizing(m, data = sdat, seed = 100 + i, iter = 5000,
                                      hessian = FALSE, as_vector = FALSE),
                    error = function(e) NULL)
      if (!is.null(o)) vals <- c(vals, o$value)
    }
    if (length(vals)) max(vals) else NA_real_
  }, error = function(e) NA_real_)
  restarted <- any(grepl("restarted from moved", attr(r, "m"), fixed = TRUE))
  data.frame(
    id = id, family = cs$fam, body = b,
    default = ll(r), started = ll(rs), no_body = ll(r0), brms_best = bm,
    gap_brms = bm - ll(r), restarted = restarted,
    seconds = round(t_fit, 2), conv = if (inherits(r, "tterr")) NA else r$opt$convergence,
    err = if (inherits(r, "tterr")) substr(r$msg, 1, 150) else "",
    warn = substr(paste(attr(r, "w"), collapse = " | "), 1, 90),
    stringsAsFactors = FALSE)
}
id <- args[1]
if (identical(id, "list")) {
  cat(cases$id, sep = "\n")
} else {
  dir.create("dev/review-nlfix-saddle-out", showWarnings = FALSE)
  res <- one(id)
  utils::write.table(res, file.path("dev/review-nlfix-saddle-out", paste0(id, ".txt")),
                     sep = "\t", row.names = FALSE)
}
