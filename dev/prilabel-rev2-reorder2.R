# Re-check of prilabel fix round 1, finding 4: the sorted term order
# keeps the model (base against lane), and start =, map =, priors and
# frmtmb.sample's theta columns follow the new positions.
#   arg 1: "base" or "lane". Data dev/dblbar-data.R seed 20261013, ga
#   seed 7; sampling seed 3.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
lg <- file.path(wt, "dev/prilabel-rev-log")
sink(file.path(lg, paste0("rev2-reorder2-", arm, ".txt")), split = TRUE)
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
dd <- dblbar_data()
set.seed(7)
dd$ga <- factor(sample(sprintf("a%d", 1:7), nrow(dd), TRUE))
dd$y <- dd$y + rnorm(7, 0, 0.5)[dd$ga]
form <- bf(y ~ x + (1 | h) + (1 + x | g) + (1 | ga) + (1 | mm(m1, m2)),
           sigma ~ (1 | h) + (1 | g))
fit <- frm(form, data = dd)
grp <- vapply(fit$frame$re_blocks, `[[`, "", "group_name")
cat("block order:", grp, "\n")
cat("theta_idx per block:",
    vapply(fit$frame$re_blocks, function(bk) paste(bk$theta_idx,
                                                    collapse = ","), ""), "\n")
cat(sprintf("logLik %.10f conv %d\n", as.numeric(logLik(fit)),
            fit$opt$convergence))
vc <- VarCorr(fit)
sds <- unlist(lapply(setdiff(names(vc), "residual__"), function(g) {
  v <- vc[[g]]$sd[, 1]
  stats::setNames(v, paste(g, names(v)))
}))
print(round(sds, 6))
if (FALSE) saveRDS(list(ll = as.numeric(logLik(fit)), sds = sds,
             fe = fixef(fit)[, 1]),
        file.path(lg, paste0("rev2-reorder-", arm, ".rds")))
if (arm == "lane") {
  base <- readRDS(file.path(lg, "rev2-reorder-base.rds"))
  cat(sprintf("\nlane - base logLik %.3e\n", as.numeric(logLik(fit)) - base$ll))
  common <- intersect(names(sds), names(base$sds))
  cat(sprintf("sds matched by name: %d of %d, max rel diff %.3e\n",
              length(common), length(sds),
              max(abs(sds[common] / base$sds[common] - 1))))
  cat(sprintf("fixef max rel diff %.3e\n",
              max(abs(fixef(fit)[, 1] / base$fe - 1))))

  cat("
## start = by position and by name
")
  pt <- par_template(fit)
  cat("par_template theta names:", names(pt$theta), "
")
  i_h <- fit$frame$re_blocks[[which(grp == "h")[1]]]$theta_idx
  st <- fit$estimates$theta
  st[i_h] <- log(0.5)
  ob <- frm(form, data = dd, start = list(theta = unname(st)), dry_run = "objective")
  cat(sprintf("positional start: theta %d of the objective = %.6f (block h)
", i_h, ob$obj$par[names(ob$obj$par) == "theta"][i_h]))
  if (!is.null(names(pt$theta))) { s2 <- pt$theta; s2[] <- 0; s2[i_h] <- log(0.5); ob2 <- frm(form, data = dd, start = list(theta = s2), dry_run = "objective"); cat("named start, theta par:", round(ob2$obj$par[names(ob2$obj$par) == "theta"], 4), "
") }
  cat("
## priors by group reach the block's theta
")
  des <- frmtmb:::prior_design(form, dd, NULL, list())
  for (gname in c("g", "ga", "h", "mmm1m2")) {
    e <- frmtmb:::resolve_priorlist(des, set_prior("normal(0, 0.5)",
                                                   class = "sd",
                                                   group = gname))$entries
    want <- unlist(lapply(Filter(function(bk) {
      identical(bk$group_name, gname) && identical(bk$dpar %||% "mu", "mu")
    }, des$frame$re_blocks), function(bk) bk$theta_idx[
      frmtmb:::block_sd_idx(bk)]))
    cat(sprintf("  group %-7s entries %s  block sds %s\n", gname,
                paste(vapply(e, `[[`, 0, "idx"), collapse = ","),
                paste(want, collapse = ",")))
  }

  cat("\n## frmtmb.sample: theta columns against the fit\n")
  fit2 <- frm(bf(y ~ x + (1 | h) + (1 | ga) + (1 | g)), data = dd)
  ds <- frm_sample(fit2, chains = 1, iter = 600, seed = 3, refresh = 0)
  m <- as.matrix(ds)
  tc <- grep("^theta", colnames(m), value = TRUE)
  cat("theta columns:", tc, "\n")
  cat("blocks:", vapply(fit2$frame$re_blocks, `[[`, "", "group_name"), "\n")
  cat("fit theta:   ", round(fit2$estimates$theta, 3), "\n")
  cat("draw medians:", round(apply(m[, tc, drop = FALSE], 2, median), 3), "\n")
  print(utils::head(as.data.frame(prior_summary(ds))[, c("prior", "class",
                                                        "group")], 8))
}
sink()
