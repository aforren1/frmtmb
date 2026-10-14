# frmtmb's reading of every `||` case (lane dblbar): its group-level
# blocks, the Z matrix of each block per coefficient, its default_prior()
# sd and cor rows, and the sd and cor names of a fit.
# Usage: Rscript dev/dblbar-frm-probe.R base|lane
# Output: dev/dblbar-rev-log/x-frm-probe-<arm>.rds and .txt.
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-rev-data.R"))
source(file.path(wt, "dev/dblbar-rev-cases.R"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
dd <- dblbar_data()
res <- list()
sink(file.path(wt, sprintf("dev/dblbar-rev-log/x-frm-probe-%s.txt", arm)),
     split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
for (nm in names(dblbar_cases)) {
  cs <- dblbar_cases[[nm]]
  cat("\n==== ", nm, ": ", cs[1], "\n", sep = "")
  out <- list(case = cs)
  form <- eval(str2lang(cs[1]))
  fam <- eval(str2lang(cs[2]))
  fr <- tryCatch(frm(form, family = fam, data = dd, dry_run = "frame"),
                 error = function(e) e)
  if (inherits(fr, "error")) {
    cat("ERROR:", conditionMessage(fr), "\n")
    out$error <- conditionMessage(fr)
    res[[nm]] <- out
    next
  }
  blocks <- list()
  for (bk in fr$re_blocks) {
    if (bk$covstruct %in% c("smooth", "gp", "hsgp", "car", "spde")) next
    D <- bk$dim
    L <- bk$n_levels
    comps <- lapply(bk$components, function(cp) {
      Z <- as.matrix(fr$linpreds[[cp$lp_key]]$Z)
      mats <- lapply(seq_len(cp$dim), function(j) {
        idx <- bk$b_idx[(seq_len(L) - 1L) * D + cp$offset + j]
        M <- Z[, idx, drop = FALSE]
        colnames(M) <- bk$levels
        M
      })
      names(mats) <- cp$cnms
      list(lp_key = cp$lp_key, cnms = cp$cnms, mats = mats)
    })
    cat(sprintf("block %-8s %-5s dim %d group %s: %s\n", bk$covstruct,
                bk$dpar %||% "", D, bk$group_name,
                paste(vapply(comps, function(c) paste0(c$lp_key, "[",
                                                       paste(c$cnms,
                                                             collapse = ","),
                                                       "]"), ""),
                      collapse = " ")))
    blocks[[length(blocks) + 1L]] <- list(covstruct = bk$covstruct,
                                          dim = D, group = bk$group_name,
                                          levels = bk$levels,
                                          comps = comps)
  }
  out$blocks <- blocks
  pr <- tryCatch(as.data.frame(default_prior(form, family = fam, data = dd)),
                 error = function(e) NULL)
  if (!is.null(pr)) {
    prs <- pr[pr$class %in% c("sd", "cor"),
              c("class", "coef", "group", "resp", "dpar", "nlpar")]
    print(prs, row.names = FALSE)
    out$prior <- prs
  }
  fit <- tryCatch(suppressWarnings(frm(form, family = fam, data = dd)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat("FIT ERROR:", conditionMessage(fit), "\n")
    out$fit_error <- conditionMessage(fit)
  } else {
    v <- grep("^(sd|cor)_", variables(fit), value = TRUE)
    cat("names:", paste(v, collapse = " "), "\n")
    cat("logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")
    out$names <- v
    out$logLik <- as.numeric(logLik(fit))
    out$vc <- tryCatch(utils::capture.output(print(VarCorr(fit))),
                       error = function(e) conditionMessage(e))
  }
  res[[nm]] <- out
}
sink()
saveRDS(res, file.path(wt, sprintf("dev/dblbar-rev-log/x-frm-probe-%s.rds", arm)))
