# Lane illse: dev/relrev070-segap.R with the arms base (rellib-r8),
# lane (the lane library) and dev (the worktree through load_all()).
#   Rscript dev/illse-segap.R base|lane|dev N
# Reviewer of the 0.70.0 consolidation: the backlog item closed as "9 of
# 9 named under both BLAS (dev/optima-segap.R)" (a mo() simplex with an
# unobserved interior category, so one direction is exactly flat), over
# more seeds. dev/optima-segap.R's generator (seeds 1..N, three
# formulas) and the reviewer's interplay generator (y ~ mo(x2), n 300).
# A fit is "named" when the SE check names a zeta, "silent" when every
# zeta SE is finite and nothing warns.
#   Rscript dev/relrev070-segap.R r8|seflags|optrap|r7 N
args <- commandArgs(TRUE)
arm <- args[1]
N <- as.integer(args[2])
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-illse-lib"),
            file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
cat("arm", arm, as.character(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
se_mo_gap <- function(seed) {
  set.seed(seed)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
  data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
}
src <- readLines("dev/relrev070-interplay.R")
i <- grep("^mk <- function", src)
j <- grep("^designs <- list", src)
eval(parse(text = src[i:(j - 1L)]))
cases <- list(
  list(tag = "segap mo", fo = ls ~ mo(inc), gen = se_mo_gap),
  list(tag = "segap mo+age", fo = ls ~ mo(inc) + age, gen = se_mo_gap),
  list(tag = "segap mo*age", fo = ls ~ mo(inc) * age, gen = se_mo_gap),
  list(tag = "inter mo(x2)", fo = y ~ mo(x2), gen = mk))
for (cs in cases) {
  named <- 0L; silent <- 0L; other <- 0L; codes <- integer(0)
  sil <- integer(0)
  oth <- character(0)
  for (s in seq_len(N)) {
    w <- character()
    f <- withCallingHandlers(frm(cs$fo, data = cs$gen(s)),
      warning = function(x) { w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning") },
      message = function(x) invokeRestart("muffleMessage"))
    sdr <- withCallingHandlers(frmtmb:::sdr_of(f),
      warning = function(x) { w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning") })
    ss <- summary(sdr, "fixed")
    zse <- ss[grepl("^zeta", rownames(ss)), 2]
    codes <- c(codes, f$opt$convergence)
    if (any(grepl("zeta", w))) named <- named + 1L else
      if (!length(w) && all(is.finite(zse))) {
        silent <- silent + 1L; sil <- c(sil, s)
      } else {
        other <- other + 1L
        oth <- c(oth, sprintf("%d: %s", s, if (length(w)) {
          substr(w[1], 1, 45)
        } else "no warning"))
      }
  }
  cat(sprintf("%-13s seeds %d | named %d | silent with finite zeta SEs %d (%s) | other %d | code != 0 %d\n",
              cs$tag, N, named, silent, paste(sil, collapse = ","), other,
              sum(codes != 0)))
  if (length(oth)) cat("   other:", oth, sep = "\n    ")
  if (length(oth)) cat("\n")
}
