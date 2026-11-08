# Fix round, m1: rows where one Gauss-Newton step (the `gn1` mutant) ends
# below two, with the second descent in place, to pin in a test. Rows: the
# lane96 set's C body (`exp(a) + b * z`, 24 rows) and the fresh set's E and F
# (24 rows), as dev/fitedge-r2-gpfresh.R builds them. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
orig <- get("gp_body_start", envir = ns)
txt <- paste(deparse(orig, width.cutoff = 500L), collapse = "\n")
stopifnot(lengths(regmatches(txt, gregexpr("for (step in 1:2)", txt,
                                           fixed = TRUE))) == 1L)
gn1 <- eval(parse(text = sub("for (step in 1:2)", "for (step in 1:1)", txt,
                             fixed = TRUE)))
environment(gn1) <- environment(orig)
use <- function(f) {
  unlockBinding("gp_body_start", ns)
  assign("gp_body_start", f, envir = ns)
  lockBinding("gp_body_start", ns)
}
src <- readLines(file.path(FE_WT, "dev/fitedge-r2-gpfresh.R"))
cut <- grep("^# theta starts in each design", src)
for (set in c("lane96", "fresh")) {
  env <- new.env()
  env$set <- set
  env$FE_WT <- FE_WT
  env$arm <- "x"
  eval(parse(text = src[grep("^if \\(set == \"fresh\"\\) \\{", src):(cut - 1L)]),
       env)
  keep <- env$rows$design %in% (if (set == "lane96") "C" else c("E", "F"))
  for (k in which(keep)) {
    r <- env$rows[k, ]
    d <- env$getd(r)
    sp <- env$getsp(r)
    use(orig)
    a <- fe_fit(sp$nl, data = d, family = sp$fam)
    use(gn1)
    b <- fe_fit(sp$nl, data = d, family = sp$fam)
    use(orig)
    la <- as.numeric(logLik(a$fit)); lb <- as.numeric(logLik(b$fit))
    cat(sprintf("G1 %s %s %s n %d seed %d two steps %.4f one step %.4f diff %.4f\n",
                set, r$design, env$lab(r), r$n, r$seed, la, lb, lb - la))
  }
}
cat("DONE\n")
