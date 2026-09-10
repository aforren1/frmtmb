# The OTHER control: this lane edits frmtmb.eam, and item 1.0a's claim
# is that a model without ndt_group() is bitwise the model 0.6.0
# fitted. Nothing here may move that.
#
#   Rscript dev/rlddm-scripts/rlddm-eam-control.R <lib> <out.rds>
#
# Seed 4242. Five families, with and without a grouping, plus the
# fitted family's own bound record, compared under identical() by
# dev/rlddm-scripts/rlddm-compare.R.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
out <- args[[2L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("lib:", lib, " eam", format(packageVersion("frmtmb.eam")), "at",
    dirname(system.file(package = "frmtmb.eam")), "\n")

set.seed(4242)
ns <- 6L
nt <- 60L
u <- stats::rnorm(ns, 0, 0.12)
s <- rep(seq_len(ns), each = nt)
d <- ddm_simulate(ns * nt, mu = 1.2 + stats::rnorm(ns, 0, 0.3)[s],
                  bs = 1.5, ndt = 0.25 * exp(u[s]), bias = 0.5)
d$s <- factor(s)
d$cond <- factor(rep(rep(0:1, each = nt / 2L), times = ns),
                 labels = c("a", "b"))

got <- list(min_rt = min(d$rt), rt = d$rt, upper = d$upper)
grab <- function(key, fit) {
  got[[paste0("ll_", key)]] <<- as.numeric(stats::logLik(fit))
  got[[paste0("cf_", key)]] <<- unlist(fixef(fit))
  got[[paste0("par_", key)]] <<- as.numeric(fit$opt$par)
  # rdm() and lba() have no mean response time to report, and say so
  got[[paste0("fit_", key)]] <<- tryCatch(
    as.numeric(suppressWarnings(stats::fitted(fit))),
    error = function(e) conditionMessage(e))
  got[[paste0("ndt_", key)]] <<- as.numeric(suppressWarnings(
    stats::predict(fit, dpar = "ndt", type = "response")))
  invisible(NULL)
}

grab("wiener", frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1,
                      bias = 0.5), family = wiener(), data = d))
grab("wiener_sv", frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1,
                         bias = 0.5, sv ~ 1),
                      family = wiener(variability = "sv"), data = d))
grab("wiener_mx", frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1,
                         bias = 0.5),
                      family = wiener(max_ndt = 0.2), data = d))
dr <- data.frame(rt = d$rt, choice = d$upper + 1L, s = d$s)
grab("rdm", frm(bf(rt | vint(choice) ~ 1), family = rdm(2), data = dr))
grab("lba", frm(bf(rt | vint(choice) ~ 1), family = lba(2), data = dr))

# the grouped arm, which 1.0a added and this lane must not move either
fg <- frm(bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s), bs ~ 1,
             ndt ~ 1 + (1 | s), bias = 0.5), family = wiener(), data = d)
grab("wiener_grp", fg)
bd <- frmtmb::single_response(fg)[["family"]][["ndt_bound"]]
got$floors <- unname(bd[["floors"]])
got$floor_names <- names(bd[["floors"]])
got$sizes <- unname(bd[["sizes"]])
got$ub <- bd[["ub"]]
got$ndt_time_grp <- as.numeric(suppressWarnings(ndt_time(fg)))

# a fixed-parameter probe of each density, which is the density rather
# than an optimizer path
pin <- function(fam, form) {
  o <- frm(form, family = fam, data = d, dry_run = "objective")
  as.numeric(o$obj$fn(o$obj$par))
}
pin2 <- function(fam, dat, form) {
  o <- frm(form, family = fam, data = dat, dry_run = "objective")
  as.numeric(o$obj$fn(o$obj$par))
}
got$pin_wiener <- pin(wiener(max_ndt = 0.2),
                      bf(rt | dec(upper) ~ 1, bs = 1.5, ndt = 0.15,
                         bias = 0.5))
got$pin_rdm <- pin2(rdm(2, max_ndt = 0.2), dr,
                    bf(rt | vint(choice) ~ 1, ndt = 0.15))
got$pin_lba <- pin2(lba(2, max_ndt = 0.2), dr,
                    bf(rt | vint(choice) ~ 1, ndt = 0.15))

saveRDS(got, out)
cat("wrote", out, "  entries:", length(got), "\n")
for (k in grep("^ll_", names(got), value = TRUE)) {
  cat(sprintf("%-16s %.12f\n", k, got[[k]]))
}
