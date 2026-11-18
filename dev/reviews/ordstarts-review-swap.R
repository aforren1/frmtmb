# Review of lane ordstarts: two questions on the (1 | g) designs, per
# seed, with the lane's data (dev/ordstarts-run.R).
# 1. Is the default fit A silent to a user who calls summary()? The lane
#    read warnings at the fit only, and the standard-error check of a
#    fit with random effects often waits for the first use. Here: the
#    warnings and messages of frm() and then of summary(A).
# 2. Does one start with the two components' roles swapped (mu1 <-> mu2,
#    tau_raw1 <-> tau_raw2, the group sds swapped, the weight's logit
#    negated) reach the point B that 10 jittered starts reach? Fitted by
#    refit(start = ) from A with starts = 0. Compared with the lane's
#    logLik of A and B (dev/ordstarts-out/main-*.tsv).
#   Rscript dev/reviews/ordstarts-review-swap.R <family> <n> <G> <seeds>
#     <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
nobs <- as.integer(args[2])
G <- as.integer(args[3])
seeds <- eval(parse(text = args[4]))
out <- args[5]
re <- "int"
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
lane <- readLines("dev/ordstarts-run.R")
eval(parse(text = lane[grep("^mkg <- ", lane):(grep("^form <- ", lane) -
                                                  1L)]))
fam <- switch(famname,
  cum_cum = function() mixture(cumulative(), cumulative()),
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
clean <- function(s) paste(substr(gsub("[\t\n\r]", " ", s), 1, 90),
                           collapse = " | ")
swap_start <- function(f) {
  p <- f$opt$par
  nm <- ns$outer_par_map(f)$names
  stopifnot(length(nm) == length(p))
  s <- p
  sw <- function(a, b) {
    ia <- match(a, nm)
    ib <- match(b, nm)
    ok <- !is.na(ia) & !is.na(ib)
    s[ia[ok]] <<- p[ib[ok]]
    s[ib[ok]] <<- p[ia[ok]]
  }
  m1 <- grep("^mu1_", nm, value = TRUE)
  sw(m1, sub("^mu1_", "mu2_", m1))
  t1 <- grep("^tau_raw1_", nm, value = TRUE)
  sw(t1, sub("^tau_raw1_", "tau_raw2_", t1))
  sw("theta_1", "theta_2")
  iw <- grep("^theta1_", nm)
  s[iw] <- -p[iw]
  st <- f$obj$env$parList(s)
  st[["b"]][] <- 0
  st
}
first <- TRUE
for (seed in seeds) {
  d <- mkg(seed)
  a <- grab(frm(bf(y ~ x + z + (1 | g)), family = fam(), data = d))
  row <- data.frame(family = famname, n = nobs, seed = seed,
                    elA = a$el, errA = inherits(a$v, "error"))
  if (!row$errA) {
    A <- a$v
    sm <- grab(summary(A))
    row$llA <- as.numeric(logLik(A))
    row$codeA <- A$opt$convergence
    row$fitwarnA <- length(a$w)
    row$summwarnA <- length(sm$w)
    row$summmsgA <- length(sm$m)
    row$summwarnsA <- clean(sm$w)
    row$summmsgsA <- clean(sm$m)
    st <- tryCatch(swap_start(A), error = function(e) e)
    if (inherits(st, "error")) {
      row$llS <- NA
      row$elS <- NA
      row$note <- clean(conditionMessage(st))
    } else {
      A0 <- A
      A0$control$starts <- 0L
      t0 <- proc.time()[["elapsed"]]
      S <- tryCatch(q(refit(A0, d$y, start = st)), error = function(e) e)
      row$elS <- proc.time()[["elapsed"]] - t0
      if (inherits(S, "error")) {
        row$llS <- NA
        row$note <- clean(conditionMessage(S))
      } else {
        row$llS <- as.numeric(logLik(S))
        row$codeS <- S$opt$convergence
        cl <- ns$mixture_ord_start_class(S, names(S$spec$responses)[1])
        row$clsS <- if (length(cl$step)) "step" else
          if (length(cl$boundary)) "flat" else "proper"
        q(ns$se_fit_analysis(S, force = TRUE))
        row$tseS <- max(q(ns$mixture_ord_thres_se(S,
          names(S$spec$responses)[1])), na.rm = TRUE)
        row$note <- ""
      }
    }
  }
  for (cn in c("llA", "codeA", "fitwarnA", "summwarnA", "summmsgA",
               "summwarnsA", "summmsgsA", "llS", "elS", "codeS", "clsS",
               "tseS", "note")) {
    if (is.null(row[[cn]])) row[[cn]] <- NA
  }
  utils::write.table(row, out, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = first, append = !first)
  first <- FALSE
  cat(famname, nobs, "seed", seed, "llA", row$llA, "llS", row$llS,
      "summwarn", row$summwarnA, "\n")
}
cat("DONE\n")
