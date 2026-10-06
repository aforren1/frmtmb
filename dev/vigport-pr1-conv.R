# Punch round 1, B2: are the two vignette models whose frmtmb fit
# reports non-convergence (brms_distreg fit_smooth1, singular
# convergence; brms_multilevel fit_loss2, false convergence) at a wrong
# optimum, or only flagged?
#
#   Rscript dev/vigport-pr1-conv.R [lib]
#
# Rebuilds each fit exactly as the spell pass does (same extractor,
# transform, patches and per-expression seeds), then refits it under
# every optimizer frm_allfit() has and prints the log-likelihoods.
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
HERE <- file.path(here, "brms-port")
source(file.path(HERE, "port-lib.R"))
source(file.path(HERE, "patches.R"))
source(file.path(HERE, "shim.R"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")

rebuild <- function(vig, target) {
  env <- new.env(parent = make_shim(suppress_brms = TRUE))
  for (k in extract_vignette(vig)) {
    tr <- transform_code(k$code)
    for (j in seq_along(tr)) {
      id <- sprintf("%s.%d.%d", vig, k$idx, j)
      if (tr[[j]]$status %in% c("BRMS-ONLY", "PARSE-ERROR", "SETUP-SKIP") ||
          kind_of(tr[[j]]$src) == "post") next
      src <- if (!is.null(PATCH[[id]])) PATCH[[id]] else tr[[j]]$src
      set.seed(port_seed(id))
      v <- suppressWarnings(tryCatch(eval(parse(text = src), envir = env),
                                     error = function(e) NULL))
      if (id == target) return(v)
    }
  }
  NULL
}
grDevices::pdf(NULL)
for (t in list(c("brms_distreg", "brms_distreg.12.1"),
               c("brms_multilevel", "brms_multilevel.19.2"))) {
  f <- rebuild(t[1], t[2])
  cat("\n==", t[2], ": convergence", f$opt$convergence, "-", f$opt$message,
      "; logLik", format(as.numeric(logLik(f)), digits = 10), "\n")
  af <- tryCatch(frm_allfit(f), error = function(e) conditionMessage(e))
  if (is.character(af)) cat("frm_allfit:", af, "\n") else print(af)
}

# frm_allfit() refits with start = NULL (R/allfit.R), and a nonlinear
# model needs a start, so its table above says nothing about fit_loss2.
# Refit fit_loss2's model under each optimizer from the start the
# vignette port gives fit_loss1, and from that start moved by 10%.
loss <- read.csv(paste0("https://raw.githubusercontent.com/mages/",
                        "diesunddas/master/Data/ClarkTriangle.csv"))
nlform2 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
              ult ~ 1 + (1 | ID1 | AY), omega ~ 1 + (1 | ID1 | AY),
              theta ~ 1 + (1 | ID1 | AY), nl = TRUE)
cat("\n== fit_loss2's model, refit per optimizer with start =\n")
for (st in list(c(5000, 1, 45), c(5500, 1.1, 49.5))) {
  opts <- frmtmb:::allfit_optimizers()
  for (o in names(opts)) {
    r <- tryCatch({
      ctl <- frmtmb_control(optimizer = opts[[o]])
      g <- suppressWarnings(frm(nlform2, data = loss, family = gaussian(),
                                start = list(beta = st), control = ctl))
      sprintf("logLik %.4f  convergence %d", as.numeric(logLik(g)),
              g$opt$convergence)
    }, error = function(e) {
      paste("ERROR:", substr(conditionMessage(e), 1, 80))
    })
    cat(sprintf("start %-18s %-13s %s\n", paste(st, collapse = ","), o, r))
  }
}
