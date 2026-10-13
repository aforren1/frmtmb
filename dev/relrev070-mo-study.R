# Reviewer copy of dev/optrap-mo-study.R (adds arm r8 and RELREV_INSTR).
# Lane optrap: dev/rel069-mo-study.R with arms r7 (rellib-r7, 0.69.0) and
# optrap (the lane's library before rellib-r7). Everything below the
# library lines is that script's.
#   Rscript dev/optrap-mo-study.R r7|optrap seeds out.tsv
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seeds <- eval(parse(text = args[2]))
out <- args[3]
libs <- switch(arm,
  r7 = "C:/Users/adf44/source/r/rellib-r7", r8 = "C:/Users/adf44/source/r/rellib-r8",
  optrap = c("C:/Users/adf44/source/r/wt-optrap-lib",
             "C:/Users/adf44/source/r/rellib-r7"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
INSTR <- Sys.getenv("RELREV_INSTR") == "1"
if (INSTR) source("dev/relrev070-instr.R")
# OPTRAP_NOKEEP=1 (only on a build that has it): optimize_obj() without
# lane optrap's rule, since removed, that kept
# a converged run over a troubled restart that gains nothing (the
# 0.69.0 rule, objective <= objective, in its place)
if (identical(Sys.getenv("OPTRAP_NOKEEP"), "1")) {
  ns <- asNamespace("frmtmb")
  src <- deparse(get("optimize_obj", ns))
  i <- grep("keep <- !inherits(opt2, \"error\") &&", src, fixed = TRUE)
  stopifnot(length(i) == 1L)
  src[i] <- "keep <- !inherits(opt2, \"error\") && isTRUE(opt2$objective <= opt$objective)"
  src <- src[-((i + 1L):(i + 3L))]
  f <- eval(parse(text = src))
  environment(f) <- ns
  utils::assignInNamespace("optimize_obj", f, "frmtmb")
  cat("keep rule removed\n")
}
# OPTIMA_NOSEARCH=1: the lane's simplex chart without mo_search(), to
# separate the two parts' effect and cost
if (identical(Sys.getenv("OPTIMA_NOSEARCH"), "1")) {
  utils::assignInNamespace("mo_search", function(obj, opt, ...) opt,
                           "frmtmb")
  cat("mo_search() disabled\n")
}
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")

mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}

# exact maximum over the four sign patterns
exact_max <- function(d) {
  code <- as.integer(d$income) - 1L
  S <- sapply(1:3, function(k) as.numeric(code >= k))
  X <- cbind(1, d$age, S, S * d$age)
  n <- nrow(X)
  best <- list(ll = -Inf)
  for (s1 in c(-1, 1)) for (s2 in c(-1, 1)) {
    A <- t(rbind(cbind(0, 0, diag(3) * s1, matrix(0, 3, 3)),
                 cbind(0, 0, matrix(0, 3, 3), diag(3) * s2)))
    qp <- quadprog::solve.QP(crossprod(X), crossprod(X, d$ls), A,
                             rep(0, 6))
    beta <- qp$solution
    r <- d$ls - X %*% beta
    ll <- -(n / 2) * (log(2 * pi * sum(r^2) / n) + 1)
    if (ll > best$ll) best <- list(ll = ll, beta = beta, s = c(s1, s2))
  }
  b <- best$beta
  d1 <- b[3:5]; d2 <- b[6:8]
  best$b_mo <- c(sum(d1), sum(d2)) / 3
  best$w1 <- if (sum(abs(d1)) > 0) d1 / sum(d1) else rep(NA, 3)
  best$w2 <- if (sum(abs(d2)) > 0) d2 / sum(d2) else rep(NA, 3)
  best
}

simplex_of <- function(z) {
  x <- exp(c(0, z) - max(c(0, z)))
  x / sum(x)
}

rows <- list()
for (s in seeds) {
  d <- mk(s)
  if (INSTR) relrev_log$seed <- s
  ex <- exact_max(d)
  t0 <- proc.time()[[3]]
  wl <- character()
  f <- withCallingHandlers(
    frm(ls ~ mo(income) * age, data = d),
    warning = function(w) {
      wl <<- c(wl, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  dt <- proc.time()[[3]] - t0
  ll <- as.numeric(logLik(f))
  fe <- fixef(f)
  bm <- fe[grepl("^moincome", rownames(fe)), 1]
  w <- if (!is.null(f$cache$mo_simplex)) f$cache$mo_simplex else {
    zt <- f$estimates[grepl("^zeta", names(f$estimates))]
    ms <- tryCatch(utils::getFromNamespace("mo_simplex", "frmtmb"),
                   error = function(e) simplex_of)
    lapply(zt, ms)
  }
  if (is.null(names(w))) names(w) <- paste0("zeta", seq_along(w))
  w <- unname(w)
  se <- fe[, 2]
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, ll_fit = ll, ll_exact = ex$ll, gap = ex$ll - ll,
    code = f$opt$convergence, sec = dt,
    nfev = f$opt$evals,
    b1_fit = bm[1], b2_fit = bm[2], b1_ex = ex$b_mo[1], b2_ex = ex$b_mo[2],
    w1_fit = paste(signif(w[[1]], 4), collapse = ","),
    w2_fit = paste(signif(w[[2]], 4), collapse = ","),
    w1_ex = paste(signif(ex$w1, 4), collapse = ","),
    w2_ex = paste(signif(ex$w2, 4), collapse = ","),
    min_w_fit = min(unlist(w)),
    nse_nan = sum(!is.finite(se)),
    nwarn = length(wl),
    warn = substr(paste(gsub("[\t\n]", " ", wl), collapse = " | "), 1, 300))
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("seeds", nrow(X), "\n")
if (INSTR) relrev_summary(paste0("mo-", arm, "-", sub("[.]tsv$", "", basename(out))))
