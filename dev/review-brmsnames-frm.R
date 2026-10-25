## frmtmb's names for dev/review-brmsnames-cases.R on one arm.
##   Rscript dev/review-brmsnames-frm.R base|lane <case> ...
arm <- commandArgs(trailingOnly = TRUE)[1L]
ids <- as.integer(commandArgs(trailingOnly = TRUE)[-1L])
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")),
    "unsub:", exists("unsub_specials", asNamespace("frmtmb")), "\n")
source("dev/review-brmsnames-cases.R")
d <- rv_data()
cases <- rv_cases(frmtmb::bf, frmtmb::mvbf)
tn <- function(e) tryCatch(e, error = function(err) paste("ERROR:",
                                                         conditionMessage(err)))
fam_of <- function(s) get(s, envir = asNamespace("frmtmb"))()
for (i in ids) {
  cs <- cases[[i]]
  cat("== case", i, cs$id, "\n")
  ref <- readRDS(sprintf("dev/review-brmsnames-out/brms-%02d.rds", i))
  fam <- fam_of(cs$family)
  fit <- tryCatch(q(frm(cs$form, data = d, family = fam)),
                  error = function(e) e)
  out <- list(id = cs$id)
  if (inherits(fit, "error")) {
    out$error <- conditionMessage(fit)
  } else {
    out$fixef <- tn(rownames(fixef(fit)))
    out$vcov <- tn(rownames(stats::vcov(fit)))
    out$variables <- tn(variables(fit))
    out$confint <- tn(rownames(q(stats::confint(fit))))
    pr <- as.data.frame(q(default_prior(cs$form, data = d, family = fam)))
    out$prior <- pr[, c("class", "coef", "group", "resp", "dpar", "nlpar")]
    if (!is.null(ref$prior)) {
      rows <- which(nzchar(ref$prior$coef))
      out$prior_ok <- vapply(rows, function(r) {
        row <- ref$prior[r, ]
        a <- list(prior = if (row$class == "sd") "student_t(3, 0, 5)" else
          "normal(0, 5)", class = row$class, coef = row$coef)
        for (nm in c("group", "dpar", "nlpar", "resp")) {
          if (nzchar(row[[nm]])) a[[nm]] <- row[[nm]]
        }
        p <- tryCatch(do.call(frmtmb::set_prior, a), error = function(e) e)
        if (inherits(p, "error")) return(conditionMessage(p))
        v <- tryCatch(q(validate_prior(p, cs$form, data = d, family = fam)),
                      error = function(e) e)
        if (inherits(v, "error")) return(conditionMessage(v))
        v <- as.data.frame(v)
        if (!any(v$source == "user" & v$coef == a$coef)) return("not applied")
        "ok"
      }, "")
      names(out$prior_ok) <- paste(ref$prior$class[rows],
                                   ref$prior$dpar[rows],
                                   ref$prior$resp[rows],
                                   ref$prior$coef[rows])
    }
    hv <- grep("^(b|sd)_", if (is.character(ref$variables))
      ref$variables else character(0), value = TRUE)
    out$hyp <- vapply(hv, function(v) {
      cls <- sub("_.*$", "", v)
      rest <- sub("^[a-z]+_", "", v)
      h <- tryCatch(q(hypothesis(fit, paste0(rest, " = 0"), class = cls)),
                    error = function(e) e)
      if (inherits(h, "error")) conditionMessage(h) else "ok"
    }, "")
  }
  saveRDS(out, sprintf("dev/review-brmsnames-out/frm-%s-%02d.rds", arm, i))
  # side by side with brms
  cmp <- function(a, b) {
    if (is.character(a) && is.character(b) && length(a) == length(b) &&
        all(a == b)) "SAME" else if (is.character(a) && is.character(b) &&
          setequal(a, b)) "ORDER" else "DIFF"
  }
  cat("brms error:", ref$error %||% "-", "\n")
  cat("frm  error:", out$error %||% "-", "\n")
  if (is.null(out$error) && is.null(ref$error)) {
    vb <- grep("^(b|sd|cor|sigma|Intercept)", ref$variables, value = TRUE)
    vf <- grep("^(b|sd|cor|sigma|Intercept)", out$variables, value = TRUE)
    cat("fixef ", cmp(out$fixef, ref$fixef), "| brms:",
        paste(ref$fixef, collapse = " "), "| frm:",
        paste(out$fixef, collapse = " "), "\n")
    cat("vcov  ", cmp(out$vcov, ref$vcov), "\n")
    cat("vars  ", cmp(vf, vb), "| only brms:",
        paste(setdiff(vb, vf), collapse = " "), "| only frm:",
        paste(setdiff(vf, vb), collapse = " "), "\n")
    kb <- with(ref$prior, paste(class, coef, group, resp, dpar, nlpar,
                                sep = "|"))
    kf <- with(out$prior, paste(class, coef, group, resp, dpar, nlpar,
                                sep = "|"))
    cat("prior ", cmp(kf, kb), "| only brms:",
        paste(setdiff(kb, kf), collapse = " "), "| only frm:",
        paste(setdiff(kf, kb), collapse = " "), "\n")
    bad <- out$prior_ok[out$prior_ok != "ok"]
    cat("set_prior not ok:", length(bad), "\n")
    if (length(bad)) print(bad)
    badh <- out$hyp[out$hyp != "ok"]
    cat("hyp not ok:", length(badh), "\n")
    if (length(badh)) print(badh)
    cat("confint:", paste(out$confint, collapse = " "), "\n")
  }
}
