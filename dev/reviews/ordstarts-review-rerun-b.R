# Review of lane ordstarts: an independent re-run of chosen seeds.
# Data written from the findings' description (mk() of
# dev/optima-csmix.R with mix = TRUE, then u0 ~ N(0, 0.5) per group and
# for the slope design u1 ~ N(0, 0.3) per group times x, drawn after the
# rest; groups rep_len(1:G, n)), compared with the lane's mkg() for
# identity. Per seed: default (A) and starts = 10 (B) fits, then at each
# point: the parameters the standard-error check loses, the smallest
# eigenvalue of the outer Hessian of the objective and its eigenvector's
# largest entries, the random-effect sds, and 8 refits from each point
# with every non-b coordinate moved by N(0, 0.5) (five times the lane's
# jitter). For the silent and moved seeds also starts = 30 (C).
#   Rscript dev/reviews/ordstarts-review-rerun.R <spec;spec;...> <out.txt>
#     spec: family,re,n,G,seed
args <- commandArgs(trailingOnly = TRUE)
specs <- strsplit(strsplit(args[1], ";")[[1]], ",")
out <- args[2]
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
sink(out, split = TRUE)
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
my_data <- function(seed, n, G, re) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  g <- rep_len(seq_len(G), n)
  lat <- lat + rnorm(G, 0, 0.5)[g]
  if (re == "slope") lat <- lat + rnorm(G, 0, 0.3)[g] * x
  y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  data.frame(x = x, z = z, y = y, g = factor(g))
}
# the lane's generator, for the identity check only
src <- readLines("dev/ordstarts-run.R")
mk_src <- readLines("dev/optima-csmix.R")
eval(parse(text = mk_src[grep("^mk <- ", mk_src):(grep("^one <- ",
                                                         mk_src) - 1L)]))
lane_body <- src[grep("^mkg <- ", src):(grep("^form <- ", src) - 1L)]
forms <- list(int = bf(y ~ x + z + (1 | g)),
              slope = bf(y ~ x + z + (1 + x | g)))
fams <- list(cum_cum = function() mixture(cumulative(), cumulative()),
             cum_sratio = function() mixture(cumulative(), sratio()))
q <- function(e) suppressWarnings(suppressMessages(e))
grab <- function(expr) {
  w <- character()
  m <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  list(v = v, w = w, m = m, el = proc.time()[["elapsed"]] - t0)
}
describe <- function(tag, f) {
  if (inherits(f, "error")) {
    cat(tag, "ERROR", conditionMessage(f), "\n")
    return(invisible())
  }
  ms <- f$opt$mixture_starts
  resp <- names(f$spec$responses)[1]
  cl <- ns$mixture_ord_start_class(f, resp)
  an <- tryCatch(q(ns$se_fit_analysis(f, force = TRUE)),
                 error = function(e) NULL)
  tse <- q(ns$mixture_ord_thres_se(f, resp))
  p <- f$opt$par
  H <- tryCatch(stats::optimHess(p, f$obj$fn, f$obj$gr),
                error = function(e) NULL)
  ev <- if (is.matrix(H) && all(is.finite(H))) eigen((H + t(H)) / 2, symmetric = TRUE) else NULL
  nm <- tryCatch(ns$outer_par_map(f)$names, error = function(e) NULL)
  if (length(nm) != length(p)) nm <- names(p)
  cat(sprintf("%s ll %.6f code %s starts %s gain %s tier %s cls %s\n", tag,
              as.numeric(logLik(f)), f$opt$convergence,
              if (is.null(ms)) "none" else ms[["starts"]],
              if (is.null(ms)) "-" else format(signif(ms[["gain"]], 4)),
              if (is.null(ms)) "-" else ms[["tier"]],
              if (length(cl$step)) "step" else if (length(cl$boundary))
                "flat" else "proper"))
  cat("   thres SE (forced):", format(signif(tse, 4)), "; lost:",
      paste(names(an$lost), collapse = ","), "\n")
  if (!is.null(ev)) {
    k <- length(ev$values)
    v <- ev$vectors[, k]
    top <- utils::head(order(-abs(v)), 3)
    cat("   min eigen", format(signif(ev$values[k], 3)), "along",
        paste0(nm[top], "(", round(v[top], 2), ")", collapse = " "),
        "; max |grad|", format(signif(max(abs(f$obj$gr(p))), 3)), "\n")
  }
  th <- grep("^theta|sd|tau", nm)
  cat("   outer:", paste0(nm, "=", signif(p, 3), collapse = " "), "\n")
  invisible()
}
jitter_check <- function(tag, f, d, sd = 0.5, K = 8) {
  if (inherits(f, "error")) return(invisible())
  f0 <- f
  f0$control$starts <- 0L
  ll0 <- as.numeric(logLik(f))
  res <- vapply(seq_len(K), function(k) {
    set.seed(9100 + k)
    st <- f$estimates
    for (nm in setdiff(names(st), "b")) {
      st[[nm]] <- st[[nm]] + stats::rnorm(length(st[[nm]]), 0, sd)
    }
    r <- tryCatch(q(refit(f0, d$y, start = st)), error = function(e) NULL)
    if (is.null(r) || !isTRUE(r$opt$convergence == 0)) return(NA_real_)
    as.numeric(logLik(r)) - ll0
  }, 0)
  cat(sprintf("   %s jitter sd %.1f: %d of %d converged; best - point %s;",
              tag, sd, sum(!is.na(res)), K,
              format(signif(suppressWarnings(max(res, na.rm = TRUE)), 4))),
      "ends above by > 1e-3:", sum(res > 1e-3, na.rm = TRUE), "\n")
}
for (s in specs) {
  fam <- s[1]
  re <- s[2]
  n <- as.integer(s[3])
  G <- as.integer(s[4])
  seed <- as.integer(s[5])
  cat("\n=====", fam, re, n, G, "seed", seed, "\n")
  d <- my_data(seed, n, G, re)
  e <- new.env()
  e$nobs <- n
  e$G <- G
  e$re <- re
  e$mk <- mk
  eval(parse(text = lane_body), envir = e)
  dl <- e$mkg(seed)
  cat("data identical to the lane's: y", identical(d$y, dl$y), "x",
      identical(d$x, dl$x), "g", identical(as.integer(d$g),
                                           as.integer(dl$g)), "\n")
  fit <- function(ctl) {
    frm(forms[[re]], family = fams[[fam]](), data = d, control = ctl)
  }
  a <- grab(fit(frmtmb_control()))
  b <- grab(fit(frmtmb_control(starts = 10)))
  cat(sprintf("A el %.1fs warnings %d messages %d | B el %.1fs warnings %d",
              a$el, length(a$w), length(a$m), b$el, length(b$w)),
      "messages", length(b$m), "\n")
  if (length(a$w)) cat("   A warns:", substr(a$w, 1, 140), sep = "\n    ")
  if (length(b$w)) cat("   B warns:", substr(b$w, 1, 140), sep = "\n    ")
  if (length(b$m)) cat("   B msgs:", substr(b$m, 1, 140), sep = "\n    ")
  cat("\n")
  describe("A", a$v)
  describe("B", b$v)
  if (!inherits(a$v, "error") && !inherits(b$v, "error")) {
    cat(sprintf("B - A = %.6f\n", as.numeric(logLik(b$v)) -
                  as.numeric(logLik(a$v))))
  }
  jitter_check("A", a$v, d)
  jitter_check("B", b$v, d)
  if (length(s) > 5 && s[6] == "c30") {
    c30 <- grab(fit(frmtmb_control(starts = 30)))
    describe("C30", c30$v)
    cat(sprintf("C30 el %.1fs; C30 - B = %.6f\n", c30$el,
                as.numeric(logLik(c30$v)) - as.numeric(logLik(b$v))))
  }
}
cat("DONE\n")
sink()
