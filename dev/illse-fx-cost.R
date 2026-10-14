# Lane illse fix round 1, item 7: where the design-coordinate Hessians
# of a wide fit are taken. The reviewer's design (dev/illse-rev-cost.R)
# at a smaller n: per call of reparam_hessian_refined_at() (the
# uncached work), its time, the caller and a hash of the point.
#   Rscript dev/illse-fx-cost.R base|lane|dev <re: 0|1> <n>
args <- commandArgs(TRUE)
arm <- args[1]
re <- as.integer(args[2])
n <- as.numeric(args[3])
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-illse-lib"),
            file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(1)
nc <- if (re) 91L else 92L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(sample(1:100, n, TRUE)))
for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
Z <- matrix(stats::rnorm(n * nc), n, nc)
colnames(Z) <- paste0("z", seq_len(nc))
d <- cbind(d, Z)
d$y <- sin(3 * d$x) + 0.05 * rowSums(Z[, 1:10]) +
  (if (re) stats::rnorm(100, 0, 0.3)[d$g] else 0) + stats::rnorm(n, 0, 0.5)
fo <- stats::as.formula(paste("y ~", paste(c(paste0("x", 1:6), colnames(Z)),
                                          collapse = " + "),
                              if (re) "+ (1 | g)"))
log <- list()
fn <- if (exists("reparam_hessian_refined_at", ns)) {
  "reparam_hessian_refined_at"
} else if (exists("reparam_hessian_refined", ns)) "reparam_hessian_refined"
if (!is.null(fn)) {
  trace(fn, where = ns, print = FALSE,
        tracer = quote({
          .fx_t0 <- proc.time()[["elapsed"]]
          .fx_who <- paste(vapply(sys.calls(), function(cl) {
            f <- cl[[1L]]
            if (is.name(f)) as.character(f) else ""
          }, ""), collapse = ">")
        }),
        exit = quote({
          assign("log", c(get("log", envir = globalenv()), list(list(
            t = proc.time()[["elapsed"]] - .fx_t0,
            who = .fx_who,
            p = sum(p * seq_along(p))))), envir = globalenv())
        }))
}
t0 <- proc.time()[["elapsed"]]
w <- character()
f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
  w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
})
t1 <- proc.time()[["elapsed"]] - t0
sdr <- ns$sdr_of(f)
cat(sprintf("arm %s | re %d | n %g | outer %d | frm %.1f s | code %d | lost %d | root %s | warnings %d\n",
            arm, re, n, length(f$opt$par), t1, f$opt$convergence,
            length(sdr$se_lost), !is.null(sdr$cov_fixed_root), length(w)))
cat(paste(" W:", substr(w, 1, 160)), sep = "\n")
for (e in log) {
  who <- strsplit(e$who, ">")[[1]]
  who <- who[who %in% c("newton_polish", "conv_headroom", "grad_headroom",
                        "se_reparam_cov", "check_convergence",
                        "fit_assembled", "sdr_of", "se_check",
                        "sdr_rescue_hessian", "cov_from_hessian")]
  cat(sprintf("  hessian %.1f s | point %.10g | %s\n", e$t, e$p,
              paste(who, collapse = ">")))
}
