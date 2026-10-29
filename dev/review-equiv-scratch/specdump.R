# Reviewer: build every sampled spec and print what it spells, so a pair
# whose model omits one of its two features shows.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); requireNamespace("brms")})
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
r <- read.delim("dev/equiv-data/results.tsv")
out <- list()
for (i in seq_len(nrow(smp))) {
  key <- paste(smp$feature_a[i], smp$feature_b[i], sep = " x ")
  mk <- eq_env$eq_specs[[key]]
  if (is.null(mk)) { cat(smp$id[i], key, ": NO SPEC\n"); next }
  sp <- tryCatch(mk(), error = function(e) list(err = conditionMessage(e)))
  dl <- function(x) if (is.null(x)) "" else
    paste(deparse(x, width.cutoff = 500L), collapse = " ")
  txt <- c(bf = sp$bf_frm %||% sp$bf %||% "", bfb = sp$bf_brms %||% "",
           fam = dl(sp$family), famf = sp$family_frm %||% "",
           prior = sp$prior %||% "",
           frm_args = dl(sp$frm_args), ctrl = dl(sp$ctrl),
           method = paste(sp$method, collapse = ","),
           margs = substr(dl(lapply(sp$margs, function(m)
             if (is.data.frame(m)) "<df>" else m)), 1, 150),
           flags = paste(c(if (!is.null(sp$nobrms)) "nobrms",
                           if (!is.null(sp$decision)) "decision",
                           if (!is.null(sp$l1_without)) paste0("l1_without=",
                             paste(sp$l1_without, collapse = "/")),
                           if (isTRUE(sp$l4)) "l4", if (isTRUE(sp$no_l1)) "no_l1",
                           if (isTRUE(sp$no_points)) "no_points",
                           if (!is.null(sp$const)) paste0("const=",
                             signif(sp$const, 3)),
                           if (!is.null(sp$skip)) "SKIP"), collapse = " "))
  v <- r$verdict[r$id == smp$id[i]]
  cat(sprintf("%s [%s] %s => %s\n   %s\n", smp$id[i], smp$status[i], key, v,
              paste(names(txt)[nzchar(txt)], txt[nzchar(txt)], sep = ": ",
                    collapse = " | ")))
}
