# brms 2.23.0's reading of every `||` case (lane dblbar): the sd and cor
# rows of default_prior(), the Z columns of standata() per block, and
# whether stancode() correlates a block. No compile.
# Output: dev/dblbar-rev-log/brms-probe.rds and brms-probe.txt.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
source(file.path(wt, "dev/dblbar-cases.R"))
suppressPackageStartupMessages(library(brms))
stopifnot(packageVersion("brms") == "2.23.0")
dd <- dblbar_data()
res <- list()
sink(file.path(wt, "dev/dblbar-rev-log/brms-probe.txt"), split = TRUE)
for (nm in names(dblbar_cases)) {
  cs <- dblbar_cases[[nm]]
  cat("\n==== ", nm, ": ", cs[1], "\n", sep = "")
  out <- list(case = cs)
  form <- tryCatch(eval(str2lang(cs[1])), error = function(e) e)
  fam <- eval(str2lang(cs[2]))
  if (inherits(form, "error")) {
    cat("FORMULA ERROR:", conditionMessage(form), "\n")
    out$error <- conditionMessage(form)
    res[[nm]] <- out
    next
  }
  pr <- tryCatch(default_prior(form, data = dd, family = fam),
                 error = function(e) e)
  if (inherits(pr, "error")) {
    cat("ERROR:", conditionMessage(pr), "\n")
    out$error <- conditionMessage(pr)
    res[[nm]] <- out
    next
  }
  pr <- as.data.frame(pr)
  prs <- pr[pr$class %in% c("sd", "cor"),
            c("class", "coef", "group", "resp", "dpar", "nlpar")]
  print(prs, row.names = FALSE)
  out$prior <- prs
  sdat <- standata(form, data = dd, family = fam)
  znm <- grep("^Z_", names(sdat), value = TRUE)
  out$Z <- sdat[znm]
  out$J <- sdat[grep("^J_", names(sdat), value = TRUE)]
  out$W <- sdat[grep("^W_", names(sdat), value = TRUE)]
  out$M <- unlist(sdat[grep("^M_", names(sdat), value = TRUE)])
  out$NC <- unlist(sdat[grep("^NC_", names(sdat), value = TRUE)])
  cat("M:", paste(names(out$M), out$M, collapse = " "), " NC:",
      paste(names(out$NC), out$NC, collapse = " "), "\n")
  cat("Z:", paste(znm, collapse = " "), "\n")
  scode <- stancode(form, data = dd, family = fam)
  cl <- grep("cor_|L_[0-9]|Cor_", strsplit(scode, "\n")[[1]],
             value = TRUE)
  out$corlines <- trimws(cl)
  cat("cor lines:", length(cl), "\n")
  if (length(cl)) cat(paste0("  ", trimws(head(cl, 4))), sep = "\n")
  # brms's own block table: group, coef, id, cor per row
  bt <- brms:::brmsterms(brms:::validate_formula(form, data = dd,
                                                 family = fam))
  rf <- tryCatch(brms:::frame_re(bt, dd), error = function(e) NULL)
  if (!is.null(rf)) {
    rf <- as.data.frame(rf)[, intersect(c("id", "group", "gn", "coef",
                                          "cn", "resp", "dpar", "nlpar",
                                          "cor", "type", "by"), names(rf))]
    print(rf, row.names = FALSE)
    out$frame_re <- rf
  }
  res[[nm]] <- out
}
sink()
saveRDS(res, file.path(wt, "dev/dblbar-rev-log/brms-probe.rds"))
