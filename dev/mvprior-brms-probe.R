# Lane wt-mvprior: what brms 2.23.0 does with each specification in
# dev/mvprior-cases.R. For each one: stancode() with the prior, reporting
# the lprior lines it writes or the error it raises, plus
# validate_prior()'s user rows. default_prior() per model at the end.
#
# Usage: Rscript dev/mvprior-brms-probe.R > dev/mvprior-log/brms-probe.txt
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages(library(brms))
source(file.path(MVPRIOR_ROOT, "dev/mvprior-cases.R"))
cat("brms", format(packageVersion("brms")), "\n")

brms_spec <- function(s) {
  brms::set_prior(s$prior, class = s$class, coef = s$coef, group = s$group,
                  resp = s$resp, dpar = s$dpar, nlpar = s$nlpar,
                  lb = if (is.na(s$lb)) NA else s$lb)
}
lab <- function(s) {
  f <- c(class = s$class, coef = s$coef, group = s$group, resp = s$resp,
         dpar = s$dpar, nlpar = s$nlpar)
  f <- f[nzchar(f)]
  paste0(s$prior, if (!is.na(s$lb)) paste0(" lb=", s$lb), " [",
         paste(names(f), f, sep = "=", collapse = ", "), "]")
}
flat <- function(x) gsub("[[:space:]]+", " ", x)

for (case in names(mvprior_cases)) {
  cc <- mvprior_cases[[case]]
  m <- mvprior_model(case, "brms")
  cat("\n== ", case, ": ", cc$what, " ==\n", sep = "")
  base_sc <- tryCatch(strsplit(stancode(m[[1]], data = m[[2]],
                                        family = m[[3]]), "\n")[[1]],
                      error = function(e) paste("ERROR", conditionMessage(e)))
  base_lp <- trimws(grep("lprior \\+=", base_sc, value = TRUE))
  for (s in cc$specs) {
    r <- tryCatch({
      pr <- brms_spec(s)
      sc <- strsplit(stancode(m[[1]], data = m[[2]], family = m[[3]],
                              prior = pr), "\n")[[1]]
      lp <- trimws(grep("lprior \\+=", sc, value = TRUE))
      bnd <- trimws(grep("<lower=", sc, value = TRUE))
      bnd0 <- trimws(grep("<lower=", base_sc, value = TRUE))
      newlp <- setdiff(lp, base_lp)
      newb <- setdiff(bnd, bnd0)
      paste0("ACCEPT: ",
             if (length(newlp)) paste(newlp, collapse = " | ") else
               "(no new lprior line)",
             if (length(newb)) paste0(" || bounds: ",
                                      paste(newb, collapse = " | ")))
    }, error = function(e) paste("REFUSE:", flat(conditionMessage(e))))
    cat(sprintf("%-58s %s\n", lab(s), r))
  }
  cat("-- default_prior() --\n")
  dp <- tryCatch(as.data.frame(default_prior(m[[1]], data = m[[2]],
                                             family = m[[3]])),
                 error = function(e) conditionMessage(e))
  if (is.data.frame(dp)) {
    print(dp[, c("prior", "class", "coef", "group", "resp", "dpar",
                 "nlpar")], row.names = FALSE)
  } else {
    cat(dp, "\n")
  }
}
